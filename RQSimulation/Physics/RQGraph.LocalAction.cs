using System;
using System.Collections.Generic;
using System.Linq;

namespace RQSimulation
{
    /// <summary>
    /// Local Action Computation for RQ-Hypothesis Compliance.
    /// 
    /// PHYSICS PRINCIPLE: Locality of Action
    /// =====================================
    /// In General Relativity, the Lagrangian density L is local - it depends
    /// only on fields and their derivatives at a single spacetime point.
    /// On a graph, this means the action change from modifying edge (i,j)
    /// should depend only on:
    /// - The edge (i,j) itself
    /// - Triangles containing (i,j) (for curvature)
    /// - Fields at nodes i and j
    /// 
    /// This gives O(degree?) computation instead of O(N?) for global Hamiltonian.
    /// 
    /// FORMULA:
    /// ?S_local = ?S_geometry + ?S_matter + ?S_volume + ?S_field
    /// 
    /// where:
    /// - ?S_geometry = -G?? ? ?R_ij (Forman-Ricci change)
    /// - ?S_matter = T_ij ? ?w_ij (stress-energy coupling)
    /// - ?S_volume = ?(V - V_target) ? ?w_ij (QUADRATIC volume stabilization from CDT)
    /// - ?S_field = gradient energy contribution
    /// 
    /// RQ-HYPOTHESIS CHECKLIST FIX (Volume Stabilization):
    /// Changed from LINEAR: ?S_vol = ? ? ?w  (causes collapse/explosion)
    /// To QUADRATIC: ?S_vol = ?(V - V_target) ? ?w (restoring force, stable 4D)
    /// 
    /// This follows Causal Dynamical Triangulations (CDT) approach where
    /// a quadratic potential fixes the total volume around a target value,
    /// creating a "restoring force" that prevents both collapse and explosion.
    /// 
    /// NOTE: Volume constraint fields (_volumeLambda, _targetTotalWeight, etc.) are
    /// defined in RQGraph.VolumeStabilization.cs. This file uses those existing fields.
    /// </summary>
    public partial class RQGraph
    {
        // NOTE: Volume constraint fields are defined in RQGraph.VolumeStabilization.cs:
        // - _volumeLambda
        // - _targetTotalWeight  
        // - _targetEdgeCount
        // - _volumeConstraintInitialized
        // Methods CountEdges() and TotalEdgeWeight() are also defined there.
        /// <summary>
        /// Compute the change in local action due to modifying edge (i,j) weight.
        /// This is the core method for RQ-compliant Metropolis-Hastings.
        /// 
        /// Implements RQ-Hypothesis Checklist: Locality of Action (Step 1).
        /// 
        /// RQ-HYPOTHESIS CHECKLIST FIX (Volume Stabilization):
        /// Uses QUADRATIC volume potential S_vol = ?(V - V_target)? instead of
        /// linear ??V term. This creates a restoring force toward target volume,
        /// following CDT (Causal Dynamical Triangulations) approach for stable 4D emergence.
        /// </summary>
        /// <param name="i">First node of edge</param>
        /// <param name="j">Second node of edge</param>
        /// <param name="newWeight">Proposed new weight</param>
        /// <returns>Change in action ?S (positive = action increases = unfavorable)</returns>
        public double ComputeLocalActionChange(int i, int j, double newWeight)
        {
            if (i < 0 || i >= N || j < 0 || j >= N || i == j)
                return double.MaxValue; // Invalid edge

            double oldWeight = Edges[i, j] ? Weights[i, j] : 0.0;
            double deltaW = newWeight - oldWeight;
            
            if (Math.Abs(deltaW) < 1e-12)
                return 0.0; // No change

            // === 1. Geometry contribution: ?S_geom = -G?? ? ?R_ij ===
            // Curvature change for all triangles containing edge (i,j)
            double deltaS_geometry = ComputeLocalRicciChange(i, j, newWeight);

            // === 2. Matter contribution: ?S_matter = T_ij ? ?w_ij ===
            // RQ-FIX: GetStressEnergyTensor now uses unified NodeMasses.TotalMass
            double stressEnergy = GetStressEnergyTensor(i, j);
            double deltaS_matter = stressEnergy * deltaW * PhysicsConstants.CurvatureTermScale;

            // === 3. Volume contribution: QUADRATIC STABILIZATION (CDT approach) ===
            // RQ-HYPOTHESIS CHECKLIST FIX:
            // CHANGED FROM: deltaS_volume = ? ? ?w (linear - causes collapse/explosion)
            // TO: deltaS_volume = ?(V - V_target) ? ?w (quadratic - restoring force)
            //
            // Physics: S_vol = ?(V - V_target)? => dS/dw = 2?(V - V_target) ? dV/dw
            // For edge weight change: dV/dw = 1 (each edge contributes its weight to volume)
            // So: ?S_vol ? 2?(V_current - V_target) ? ?w
            //
            // Volume constraint fields (_volumeLambda, _targetTotalWeight, etc.) are
            // defined in RQGraph.VolumeStabilization.cs
            double deltaS_volume;
            if (_volumeConstraintInitialized && _volumeLambda > 0)
            {
                // Quadratic volume stabilization (CDT approach for stable 4D)
                // Use ComputeVolumePenaltyChange from VolumeStabilization.cs for efficient delta
                int edgeCreated = 0;
                if (!Edges[i, j] && newWeight > 0) edgeCreated = 1;  // Creating edge
                else if (Edges[i, j] && newWeight <= 0) edgeCreated = -1;  // Removing edge
                
                deltaS_volume = ComputeVolumePenaltyChange(edgeCreated, deltaW);
                
                // Also add small cosmological constant term for baseline
                deltaS_volume += PhysicsConstants.CosmologicalConstant * deltaW * 0.1;
            }
            else
            {
                // Fallback to linear cosmological constant (less stable)
                deltaS_volume = PhysicsConstants.CosmologicalConstant * deltaW;
            }

            // === 4. Field gradient contribution ===
            double deltaS_field = ComputeLocalFieldActionChange(i, j, newWeight);

            // Total action change
            return deltaS_geometry + deltaS_matter + deltaS_volume + deltaS_field;
        }

        /// <summary>
        /// Compute change in Forman-Ricci curvature contribution to action.
        /// Only considers triangles containing edge (i,j).
        /// 
        /// Einstein-Hilbert on graph: S_EH = ? w_ij ? R_ij
        /// Change: ?S_EH = w_new ? R_new - w_old ? R_old
        /// </summary>
        private double ComputeLocalRicciChange(int i, int j, double newWeight)
        {
            double oldWeight = Edges[i, j] ? Weights[i, j] : 0.0;
            
            // Compute curvature with old weight
            double R_old = ComputeFormanRicciCurvatureLocal(i, j, oldWeight);
            
            // Compute curvature with new weight (temporary calculation)
            double R_new = ComputeFormanRicciCurvatureLocal(i, j, newWeight);
            
            // Action change: ?S = -(1/G) ? ?(w ? R)
            // Negative sign because we MINIMIZE action (gravity is attractive)
            double S_old = oldWeight * R_old;
            double S_new = newWeight * R_new;
            
            double G_inv = 1.0 / PhysicsConstants.GravitationalCoupling;
            return -G_inv * (S_new - S_old);
        }

        /// <summary>
        /// Compute Forman-Ricci curvature for edge (i,j) with specified weight.
        /// Uses LOCAL computation - only considers triangles containing (i,j).
        /// 
        /// Formula: R_F(e) = w_e ? [?_? (w_ik ? w_jk)^(1/3) - ?(W_i + W_j - 2w_e)]
        /// </summary>
        private double ComputeFormanRicciCurvatureLocal(int i, int j, double w_e)
        {
            if (w_e <= 0)
                return 0.0;

            // Compute weighted degrees (excluding edge i-j for correct formula)
            double W_i = 0.0, W_j = 0.0;
            foreach (int k in Neighbors(i))
            {
                if (k != j)
                    W_i += Weights[i, k];
            }
            foreach (int k in Neighbors(j))
            {
                if (k != i)
                    W_j += Weights[j, k];
            }

            // Triangle contribution: for each common neighbor k
            double triangleSum = 0.0;
            foreach (int k in Neighbors(i))
            {
                if (k == j) continue;
                if (!Edges[j, k]) continue;

                // Triangle i-k-j exists
                double w_ik = Weights[i, k];
                double w_jk = Weights[j, k];
                triangleSum += Math.Pow(w_ik * w_jk * w_e, 1.0 / 3.0);
            }

            // Forman-Ricci formula
            double alpha = PhysicsConstants.DegreePenaltyFactor;
            double curvature = w_e * (triangleSum - alpha * (W_i + W_j));

            return curvature;
        }

        /// <summary>
        /// Compute change in field action from modifying edge weight.
        /// Gradient energy: E_grad = (1/2) ? w_ij |?_i - ?_j|?
        /// </summary>
        private double ComputeLocalFieldActionChange(int i, int j, double newWeight)
        {
            double oldWeight = Edges[i, j] ? Weights[i, j] : 0.0;
            double deltaW = newWeight - oldWeight;

            if (Math.Abs(deltaW) < 1e-12 || ScalarField == null)
                return 0.0;

            // Field gradient at edge (i,j)
            double phi_i = ScalarField[i];
            double phi_j = ScalarField[j];
            double gradSq = (phi_i - phi_j) * (phi_i - phi_j);

            // Gradient energy contribution: (1/2) ? ?w ? |??|?
            return 0.5 * deltaW * gradSq;
        }

        /// <summary>
        /// Get stress-energy tensor component T_ij for edge (i,j).
        /// This is the source term for gravity in Einstein equations.
        /// 
        /// T_ij = T_matter + T_field where:
        /// - T_matter = (m_i + m_j) / 2 (average mass at endpoints)
        /// - T_field = (1/2)|D_i ?||D_j ?| (field gradient product)
        /// 
        /// Implements RQ-Hypothesis Checklist: Replace _correlationMass with T_ij
        /// </summary>
        public double GetStressEnergyTensor(int i, int j)
        {
            double T_ij = 0.0;

            // === Matter contribution ===
            if (NodeMasses != null && i < NodeMasses.Length && j < NodeMasses.Length)
            {
                T_ij += 0.5 * (NodeMasses[i].TotalMass + NodeMasses[j].TotalMass);
            }
            else if (_correlationMass != null && i < _correlationMass.Length && j < _correlationMass.Length)
            {
                // Fallback to correlation mass
                T_ij += 0.5 * (_correlationMass[i] + _correlationMass[j]);
            }

            // === Scalar field contribution ===
            if (ScalarField != null && i < ScalarField.Length && j < ScalarField.Length)
            {
                double phi_i = ScalarField[i];
                double phi_j = ScalarField[j];

                // Kinetic energy density from field gradient
                // T_?? ~ (?_? ?)(?_? ?) - (1/2)g_?? L
                double gradPhiSq = (phi_i - phi_j) * (phi_i - phi_j);
                T_ij += 0.5 * PhysicsConstants.ScalarFieldEnergyWeight * gradPhiSq;

                // Potential energy contribution (if Mexican Hat)
                if (UseMexicanHatPotential)
                {
                    double V_i = -HiggsMuSquared * phi_i * phi_i + HiggsLambda * Math.Pow(phi_i, 4);
                    double V_j = -HiggsMuSquared * phi_j * phi_j + HiggsLambda * Math.Pow(phi_j, 4);
                    T_ij += 0.5 * (V_i + V_j);
                }
            }

            // === Spinor field contribution ===
            if (_spinorA != null && i < _spinorA.Length && j < _spinorA.Length)
            {
                double psi_i_sq = _spinorA[i].Magnitude * _spinorA[i].Magnitude;
                double psi_j_sq = _spinorA[j].Magnitude * _spinorA[j].Magnitude;
                T_ij += 0.5 * PhysicsConstants.FermionFieldEnergyWeight * (psi_i_sq + psi_j_sq);
            }

            // === Gauge field contribution ===
            if (_edgePhaseU1 != null)
            {
                // Electric field energy E? ~ (phase gradient)?
                double E_ij = _edgePhaseU1[i, j];
                T_ij += PhysicsConstants.GaugeFieldEnergyWeight * E_ij * E_ij;
            }

            return T_ij;
        }

        /// <summary>
        /// Storage for unified node mass models
        /// </summary>
        private Physics.NodeMassModel[]? _nodeMasses;

        /// <summary>
        /// Get or create node mass models array
        /// </summary>
        public Physics.NodeMassModel[] NodeMasses
        {
            get
            {
                if (_nodeMasses == null || _nodeMasses.Length != N)
                {
                    _nodeMasses = new Physics.NodeMassModel[N];
                    for (int i = 0; i < N; i++)
                    {
                        _nodeMasses[i] = new Physics.NodeMassModel();
                    }
                }
                return _nodeMasses;
            }
        }

        /// <summary>
        /// Update all node mass models from current field values.
        /// Should be called periodically to keep masses synchronized.
        /// </summary>
        public void UpdateNodeMassModels()
        {
            var masses = NodeMasses;

            for (int i = 0; i < N; i++)
            {
                masses[i].Reset();

                // Correlation mass (topological)
                if (_correlationMass != null && i < _correlationMass.Length)
                {
                    masses[i].CorrelationMass = _correlationMass[i];
                }

                // Fermion mass from spinor field
                if (_spinorA != null && i < _spinorA.Length)
                {
                    double spinorNorm = _spinorA[i].Magnitude;
                    if (_spinorB != null && i < _spinorB.Length)
                        spinorNorm += _spinorB[i].Magnitude;
                    masses[i].FermionMass = spinorNorm;
                }

                // Scalar field energy
                if (ScalarField != null && i < ScalarField.Length)
                {
                    double phi = ScalarField[i];
                    if (UseMexicanHatPotential)
                    {
                        masses[i].ScalarFieldEnergy = -HiggsMuSquared * phi * phi + HiggsLambda * Math.Pow(phi, 4);
                    }
                    else
                    {
                        masses[i].ScalarFieldEnergy = 0.5 * ScalarMass * ScalarMass * phi * phi;
                    }
                }

                // Gauge field energy at node
                double gaugeE = 0.0;
                foreach (int j in Neighbors(i))
                {
                    if (_edgePhaseU1 != null)
                    {
                        gaugeE += _edgePhaseU1[i, j] * _edgePhaseU1[i, j];
                    }
                }
                masses[i].GaugeFieldEnergy = 0.5 * gaugeE;

                // Vacuum energy (cosmological)
                masses[i].VacuumEnergy = PhysicsConstants.CosmologicalConstant;
            }
        }

        /// <summary>
        /// Optimized Metropolis step using local action calculation.
        /// Replaces O(N?) global Hamiltonian with O(degree?) local action.
        /// 
        /// Implements RQ-Hypothesis Checklist: Locality of Action (Step 1).
        /// Uses quadratic volume stabilization for stable 4D emergence.
        /// </summary>
        /// <returns>True if change was accepted</returns>
        public bool MetropolisEdgeStepLocalAction()
        {
            // Select random edge
            int i = _rng.Next(N);
            int j = _rng.Next(N);
            if (i == j) return false;

            // Store old state
            double w_old = Weights[i, j];
            bool edge_existed = Edges[i, j];

            // Propose change
            double proposal = _rng.NextDouble();
            double w_new;
            bool create_edge = false;
            bool remove_edge = false;

            if (proposal < 0.2 && !edge_existed)
            {
                // Causality check for new edge
                if (!AreCausallyConnected(i, j))
                    return false;

                // RQ-FIX: Bipartite check for Staggered Fermions (Issue 4)
                // Prevent edges within same sublattice (even-even or odd-odd)
                // This enforces a bipartite structure required for chiral fermions
                if ((i % 2) == (j % 2))
                    return false;

                create_edge = true;
                w_new = 0.1 + 0.3 * _rng.NextDouble();
            }
            else if (proposal < 0.4 && edge_existed && w_old < 0.15)
            {
                remove_edge = true;
                w_new = 0.0;
            }
            else if (edge_existed)
            {
                double delta = (_rng.NextDouble() * 2.0 - 1.0) * 0.1;
                w_new = Math.Clamp(w_old + delta, 0.05, 1.0);
            }
            else
            {
                return false;
            }

            // === KEY OPTIMIZATION: Compute LOCAL action change ===
            double deltaS = ComputeLocalActionChange(i, j, w_new);

            // Metropolis acceptance criterion
            bool accept;
            if (deltaS <= 0)
            {
                accept = true;
            }
            else
            {
                double prob = Math.Exp(-deltaS / _networkTemperature);
                accept = _rng.NextDouble() < prob;
            }

            if (accept)
            {
                // Apply topology change
                if (create_edge)
                {
                    Edges[i, j] = true;
                    Edges[j, i] = true;
                    _degree[i]++;
                    _degree[j]++;
                    InvalidateTopologyCache();
                }
                else if (remove_edge)
                {
                    Edges[i, j] = false;
                    Edges[j, i] = false;
                    _degree[i]--;
                    _degree[j]--;
                    InvalidateTopologyCache();
                }

                Weights[i, j] = w_new;
                Weights[j, i] = w_new;
            }

            return accept;
        }
    }
}
