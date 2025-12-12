using System;
using System.Linq;

namespace RQSimulation
{
    public partial class RQGraph
    {
        // Energy computation configuration - now using PhysicsConstants (Checklist B.2)
        // These getters allow runtime configuration while defaulting to PhysicsConstants
        private double _graphEnergyWeight = PhysicsConstants.GraphLinkEnergyWeight;
        private double _scalarFieldWeight = PhysicsConstants.ScalarFieldEnergyWeight;
        private double _fermionFieldWeight = PhysicsConstants.FermionFieldEnergyWeight;
        private double _gaugeFieldWeight = PhysicsConstants.GaugeFieldEnergyWeight;
        private double _yangMillsFieldWeight = PhysicsConstants.YangMillsFieldEnergyWeight;
        private double _gravityCurvatureWeight = PhysicsConstants.GravityCurvatureEnergyWeight;
        private double _clusterBindingWeight = PhysicsConstants.ClusterBindingEnergyWeight;
        
        /// <summary>Weight for graph link energy in unified Hamiltonian. Default from PhysicsConstants.</summary>
        public double GraphEnergyWeight { get => _graphEnergyWeight; set => _graphEnergyWeight = value; }
        
        /// <summary>Weight for scalar field energy in unified Hamiltonian. Default from PhysicsConstants.</summary>
        public double ScalarFieldWeight { get => _scalarFieldWeight; set => _scalarFieldWeight = value; }
        
        /// <summary>Weight for fermion field energy in unified Hamiltonian. Default from PhysicsConstants.</summary>
        public double FermionFieldWeight { get => _fermionFieldWeight; set => _fermionFieldWeight = value; }
        
        /// <summary>Weight for gauge field energy in unified Hamiltonian. Default from PhysicsConstants.</summary>
        public double GaugeFieldWeight { get => _gaugeFieldWeight; set => _gaugeFieldWeight = value; }
        
        /// <summary>Weight for Yang-Mills field energy in unified Hamiltonian. Default from PhysicsConstants.</summary>
        public double YangMillsFieldWeight { get => _yangMillsFieldWeight; set => _yangMillsFieldWeight = value; }
        
        /// <summary>Weight for gravity/curvature energy in unified Hamiltonian. Default from PhysicsConstants.</summary>
        public double GravityCurvatureWeight { get => _gravityCurvatureWeight; set => _gravityCurvatureWeight = value; }
        
        /// <summary>Weight for cluster binding energy in unified Hamiltonian. Default from PhysicsConstants.</summary>
        public double ClusterBindingWeight { get => _clusterBindingWeight; set => _clusterBindingWeight = value; }
        
        // Energy ledger for tracking energy conservation (checklist item 4)
        private EnergyLedger? _energyLedger;
        
        /// <summary>
        /// Get or create the energy ledger for tracking conservation.
        /// Implements checklist item 4: Unified energy ledger.
        /// </summary>
        public EnergyLedger GetEnergyLedger()
        {
            if (_energyLedger == null)
            {
                _energyLedger = new EnergyLedger();
                _energyLedger.Initialize(ComputeTotalEnergyUnified());
            }
            return _energyLedger;
        }
        
        /// <summary>
        /// Compute the unified Hamiltonian combining geometric and matter terms.
        /// H = H_geom + H_matter where:
        /// - H_geom = -Σ w_ij * R_ij (Einstein-Hilbert analog, curvature-weight product)
        /// - H_matter = ⟨ψ|D|ψ⟩ (Dirac operator expectation value)
        /// Implements checklist item 3.1: Unified action/Hamiltonian.
        /// </summary>
        /// <returns>Total Hamiltonian energy</returns>
        public double ComputeHamiltonian()
        {
            double H_geom = 0.0;
            double H_matter = 0.0;
            
            // H_geom = -Σ w_ij * R_ij (analog of Einstein-Hilbert action ∫R√g d⁴x)
            // Positive curvature regions contribute negative energy (attractive gravity)
            for (int i = 0; i < N; i++)
            {
                foreach (int j in Neighbors(i))
                {
                    if (j <= i) continue; // Count each edge once
                    
                    double weight = Weights[i, j];
                    double curvature = ComputeFormanRicciCurvature(i, j);
                    H_geom -= weight * curvature;
                }
            }
            
            // H_matter = ⟨ψ|D|ψ⟩ (Dirac operator expectation value with spinors)
            H_matter = ComputeDiracExpectationValue();
            
            return H_geom + H_matter;
        }
        
        /// <summary>
        /// Compute Forman-Ricci curvature for an edge.
        /// R(e) = w(e) * [#triangles - degree_penalty]
        /// Positive for clustered (sphere-like), negative for tree-like structures.
        /// Implements checklist item 3.1: Forman-Ricci curvature computation.
        /// </summary>
        /// <param name="i">First node of edge</param>
        /// <param name="j">Second node of edge</param>
        /// <returns>Forman-Ricci curvature for edge (i,j)</returns>
        public double ComputeFormanRicciCurvature(int i, int j)
        {
            if (!Edges[i, j])
                return 0.0;
            
            double w_edge = Weights[i, j];
            
            // Weighted degrees (excluding the edge itself)
            double w_i = 0.0, w_j = 0.0;
            foreach (int n in Neighbors(i))
                w_i += Weights[i, n];
            foreach (int n in Neighbors(j))
                w_j += Weights[j, n];
            w_i -= w_edge;
            w_j -= w_edge;
            
            // Triangle contribution (common neighbors)
            double triangleContribution = 0.0;
            foreach (int k in Neighbors(i))
            {
                if (k == j) continue;
                if (Edges[j, k]) // Triangle i-j-k
                {
                    double w_ik = Weights[i, k];
                    double w_jk = Weights[j, k];
                    // Geometric mean of triangle weights
                    triangleContribution += Math.Pow(w_ik * w_jk * w_edge, 1.0 / 3.0);
                }
            }
            
            // Forman-Ricci: R(e) = w(e) * (triangles - degree_penalty)
            // Using consistent degree penalty factor from PhysicsConstants
            double curvature = w_edge * (triangleContribution - PhysicsConstants.DegreePenaltyFactor * (w_i + w_j));
            
            return curvature;
        }
        
        /// <summary>
        /// Compute Dirac operator expectation value ⟨ψ|D|ψ⟩.
        /// The Dirac operator on the graph couples spinor components through links.
        /// D_ij = γ^μ * U_ij where U_ij is the gauge link.
        /// </summary>
        /// <returns>Dirac expectation value (matter energy)</returns>
        private double ComputeDiracExpectationValue()
        {
            if (_spinorA == null || _spinorA.Length != N)
                return 0.0;
            
            double expectation = 0.0;
            
            // For each node, compute D|ψ⟩ and inner product ⟨ψ|D|ψ⟩
            for (int i = 0; i < N; i++)
            {
                // Sum over neighbors (hopping term)
                System.Numerics.Complex Dpsi_A = System.Numerics.Complex.Zero;
                System.Numerics.Complex Dpsi_B = System.Numerics.Complex.Zero;
                
                foreach (int j in Neighbors(i))
                {
                    // Gauge link U_ij (U(1) phase)
                    System.Numerics.Complex U_ij = GetLinkVariable(i, j);
                    
                    // Dirac hopping: γ·∇ ~ σ^μ * (ψ_j - ψ_i)
                    // Simplified to just hopping with gauge link
                    if (_spinorA != null && j < _spinorA.Length)
                    {
                        Dpsi_A += Weights[i, j] * U_ij * _spinorA[j];
                        if (_spinorB != null && j < _spinorB.Length)
                            Dpsi_B += Weights[i, j] * U_ij * _spinorB[j];
                    }
                }
                
                // Mass term: m * ψ†γ⁰ψ
                double mass = (_correlationMass != null && i < _correlationMass.Length) ? _correlationMass[i] : 1.0;
                System.Numerics.Complex psi_A = _spinorA[i];
                System.Numerics.Complex psi_B = (_spinorB != null && i < _spinorB.Length) ? _spinorB[i] : System.Numerics.Complex.Zero;
                
                // ⟨ψ|D|ψ⟩ = ψ†Dψ
                double kinetic = (System.Numerics.Complex.Conjugate(psi_A) * Dpsi_A).Real 
                               + (System.Numerics.Complex.Conjugate(psi_B) * Dpsi_B).Real;
                double massTerm = mass * (psi_A.Magnitude * psi_A.Magnitude + psi_B.Magnitude * psi_B.Magnitude);
                
                expectation += kinetic + massTerm;
            }
            
            return expectation;
        }
        
        /// <summary>
        /// Compute total unified energy functional combining all contributions.
        /// Implements checklist item 4.2: H_total = H_matter + H_field + H_vacuum + K_geometry.
        /// 
        /// Includes:
        /// - Graph link energy (potential energy in edge weights)
        /// - Scalar field energy (kinetic + gradient + potential)
        /// - Fermion field energy
        /// - Gauge field energy (U(1))
        /// - Yang-Mills field energy (SU(3))
        /// - Gravity curvature energy
        /// - Cluster binding energy
        /// - Geometric kinetic energy K = Σ π_ij²/(2M) (Checklist F.3)
        /// </summary>
        public double ComputeTotalEnergyUnified()
        {
            double E_links = ComputeGraphLinkEnergy();
            double E_scalar = ComputeScalarFieldEnergy();
            double E_fermion = ComputeFermionFieldEnergy();
            double E_gauge = ComputeGaugeFieldEnergy();
            double E_yangMills = ComputeYangMillsFieldEnergy();
            double E_grav = ComputeGravityCurvatureEnergy();
            double E_bind = ComputeClusterBindingEnergy();
            
            // Checklist F.3: Add geometric kinetic energy K = Σ π²/(2M)
            // This represents gravitational wave energy and geometry momentum
            double E_geomKinetic = ComputeGeometryKineticEnergy();
            
            return GraphEnergyWeight * E_links
                 + ScalarFieldWeight * E_scalar
                 + FermionFieldWeight * E_fermion
                 + GaugeFieldWeight * E_gauge
                 + YangMillsFieldWeight * E_yangMills
                 + GravityCurvatureWeight * E_grav
                 + ClusterBindingWeight * E_bind
                 + E_geomKinetic;  // Kinetic energy of geometry (gravitational waves)
        }
        
        /// <summary>
        /// Compute Yang-Mills field energy (gluon, weak, hypercharge).
        /// Implements checklist item 4.2: Include Yang-Mills energy in total.
        /// This returns the action S = ∫ F² which represents field energy.
        /// </summary>
        private double ComputeYangMillsFieldEnergy()
        {
            // Delegate to existing Yang-Mills action computation
            return ComputeYangMillsAction();
        }
        
        /// <summary>
        /// Propose a weight change and check if it's affordable energy-wise.
        /// Implements checklist item 4.3: Ledger.CanAfford(deltaE) check.
        /// </summary>
        /// <param name="i">First node</param>
        /// <param name="j">Second node</param>
        /// <param name="newWeight">Proposed new weight</param>
        /// <returns>True if the change is accepted, false otherwise</returns>
        public bool ProposeWeightChangeWithEnergyCheck(int i, int j, double newWeight)
        {
            if (!Edges[i, j])
                return false;
            
            // Compute energy before
            double E_before = ComputeTotalEnergyUnified();
            double oldWeight = Weights[i, j];
            
            // Tentatively apply change
            Weights[i, j] = newWeight;
            Weights[j, i] = newWeight;
            
            // Compute energy after
            double E_after = ComputeTotalEnergyUnified();
            double deltaE = E_after - E_before;
            
            // Check with energy ledger
            var ledger = GetEnergyLedger();
            bool canAfford = ledger.CanAfford(deltaE);
            
            if (!canAfford)
            {
                // Revert the change
                Weights[i, j] = oldWeight;
                Weights[j, i] = oldWeight;
                return false;
            }
            
            // If energy increase, use Metropolis criterion
            if (deltaE > 0)
            {
                // Try to spend vacuum energy
                if (!ledger.TrySpendVacuumEnergy(deltaE))
                {
                    // Not enough vacuum energy, revert
                    Weights[i, j] = oldWeight;
                    Weights[j, i] = oldWeight;
                    return false;
                }
            }
            else
            {
                // Energy released, return to vacuum pool
                ledger.RegisterRadiation(-deltaE);
            }
            
            return true;
        }
        
        /// <summary>
        /// Graph link energy from edge weights
        /// </summary>
        private double ComputeGraphLinkEnergy()
        {
            double energy = 0;
            
            for (int i = 0; i < N; i++)
            {
                foreach (int j in Neighbors(i))
                {
                    if (j > i) // Count each edge once
                    {
                        // Energy proportional to weight squared (quadratic potential)
                        energy += Weights[i, j] * Weights[i, j];
                    }
                }
            }
            
            return energy;
        }
        
        // Note: ComputeScalarFieldEnergy() already exists in RQGraph.FieldTheory.cs
        
        /// <summary>
        /// Fermion field energy
        /// </summary>
        private double ComputeFermionFieldEnergy()
        {
            // Check if spinor fields are initialized (they're named _spinorA, _spinorB, _spinorC, _spinorD)
            if (_spinorA == null || _spinorA.Length != N)
                return 0.0;
            
            double energy = 0;
            
            // Energy from spinor field magnitude (4 components: A, B, C, D)
            for (int i = 0; i < N; i++)
            {
                var psiA = _spinorA[i];
                var psiB = _spinorB?[i] ?? System.Numerics.Complex.Zero;
                var psiC = _spinorC?[i] ?? System.Numerics.Complex.Zero;
                var psiD = _spinorD?[i] ?? System.Numerics.Complex.Zero;
                
                energy += psiA.Real * psiA.Real + psiA.Imaginary * psiA.Imaginary;
                energy += psiB.Real * psiB.Real + psiB.Imaginary * psiB.Imaginary;
                energy += psiC.Real * psiC.Real + psiC.Imaginary * psiC.Imaginary;
                energy += psiD.Real * psiD.Real + psiD.Imaginary * psiD.Imaginary;
            }
            
            return energy;
        }
        
        /// <summary>
        /// Gravity contribution from graph curvature
        /// </summary>
        private double ComputeGravityCurvatureEnergy()
        {
            double energy = 0;
            int edgeCount = 0;
            
            // Sum curvature contributions
            for (int i = 0; i < N; i++)
            {
                foreach (int j in Neighbors(i))
                {
                    if (j > i)
                    {
                        double curvature = CalculateGraphCurvature(i, j);
                        // Negative curvature increases energy (penalizes tree-like structures)
                        // Positive curvature decreases energy (favors clustering)
                        energy -= curvature;
                        edgeCount++;
                    }
                }
            }
            
            return edgeCount > 0 ? energy / edgeCount : 0;
        }
        
        /// <summary>
        /// Cluster binding energy - stabilizes clusters through triangle terms
        /// </summary>
        private double ComputeClusterBindingEnergy()
        {
            double energy = 0;
            
            // Energy reduction for closed triangles with strong edges
            for (int i = 0; i < N; i++)
            {
                foreach (int j in Neighbors(i))
                {
                    if (j <= i) continue;
                    
                    foreach (int k in Neighbors(i))
                    {
                        if (k <= j || !Edges[j, k]) continue;
                        
                        // Triangle i-j-k exists
                        double w1 = Weights[i, j];
                        double w2 = Weights[j, k];
                        double w3 = Weights[k, i];
                        
                        // Triangle strength (geometric mean)
                        double triStrength = Math.Pow(w1 * w2 * w3, 1.0 / 3.0);
                        
                        // Strong triangles reduce energy (bind clusters)
                        energy -= triStrength;
                    }
                }
            }
            
            return energy;
        }
        
        /// <summary>
        /// Check if topology change is acceptable based on energy criterion
        /// </summary>
        public bool AcceptTopologyChange(double energyBefore, double energyAfter, double temperature)
        {
            double dE = energyAfter - energyBefore;
            
            // Always accept if energy decreases
            if (dE < 0)
                return true;
            
            // Metropolis criterion for increases
            if (temperature > 0)
            {
                double probability = Math.Exp(-dE / temperature);
                return _rng.NextDouble() < probability;
            }
            
            return false;
        }
    }
}
