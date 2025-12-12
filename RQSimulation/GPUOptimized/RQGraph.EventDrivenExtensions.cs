using System;
using System.Linq;
using System.Numerics;
using RQSimulation.Gauge;

namespace RQSimulation
{
    /// <summary>
    /// Extensions to RQGraph for GPU-optimized event-driven simulation
    /// These methods support the new event-driven architecture
    /// </summary>
    public partial class RQGraph
    {
        /// <summary>
        /// Compute local proper time increment for a specific node
        /// Takes into account local field energy and curvature
        /// </summary>
        public double ComputeLocalProperTime(int nodeId)
        {
            if (nodeId < 0 || nodeId >= N)
                return 0.01;

            // Base time step
            double dt_base = 0.01;

            // Modify by local energy density (higher energy → slower time)
            if (LocalPotential != null && nodeId < LocalPotential.Length)
            {
                double energy = LocalPotential[nodeId];
                // Time dilation: dt_proper = dt_coordinate * sqrt(1 - 2GM/r)
                // Approximate: dt ∝ 1 / (1 + energy)
                double energyFactor = 1.0 / (1.0 + energy * 0.1);
                dt_base *= energyFactor;
            }

            return Math.Max(dt_base, 0.001); // Minimum time step
        }

        /// <summary>
        /// Update physics for a specific node and its local neighborhood
        /// This is the core of the event-driven evolution
        /// </summary>
        public void UpdateNodePhysics(int nodeId, double dt)
        {
            if (nodeId < 0 || nodeId >= N)
                return;

            // Update local scalar field (if using Mexican Hat potential)
            if (ScalarField != null && _scalarMomentum != null && nodeId < ScalarField.Length)
            {
                UpdateScalarFieldAtNode(nodeId, dt);
            }

            // Update local gauge fields
            if (_gluonField != null)
            {
                UpdateGaugeFieldsAtNode(nodeId, dt);
            }

            // Update spinor field (Dirac evolution)
            if (_spinorA != null && nodeId < _spinorA.Length)
            {
                UpdateSpinorFieldAtNode(nodeId, dt);
            }

            // Update state based on quantum dynamics
            UpdateNodeState(nodeId);
        }

        /// <summary>
        /// Get time dilation factor for a node based on local curvature and energy.
        /// Returns factor by which proper time differs from coordinate time.
        /// 
        /// FIX: Now uses unified mass model (including scalar field contribution)
        /// instead of hardcoded _correlationMass.
        /// 
        /// Physics: In GR, time dilation near a massive object is:
        ///   dτ/dt = sqrt(1 - 2GM/rc²)
        /// 
        /// In the RQ framework, mass emerges from multiple sources:
        ///   - Scalar field (Higgs mechanism): m_scalar ∝ |φ|
        ///   - Fermion density: m_fermion ∝ ψ†ψ
        ///   - Gauge field energy: E_gauge ∝ F²
        ///   - Correlation structure (geometric mass)
        /// 
        /// The total mass determines gravitational time dilation.
        /// </summary>
        public double GetTimeDilation(int nodeId)
        {
            if (nodeId < 0 || nodeId >= N)
                return 1.0;

            double dilationFactor = 1.0;

            // Factor 1: Energy density (higher energy → slower time)
            if (LocalPotential != null && nodeId < LocalPotential.Length)
            {
                double energy = LocalPotential[nodeId];
                // GR-like time dilation: sqrt(1 - 2GM/rc²)
                dilationFactor *= Math.Sqrt(1.0 / (1.0 + energy * 0.1));
            }

            // Factor 2: Mass from unified mass model (scalar field + fermion + gauge)
            // FIX: Use GetNodeTotalMass which includes scalar field (Higgs) contribution
            double totalMass = GetNodeTotalMass(nodeId);
            if (totalMass > 0)
            {
                // Time dilation from total mass/energy
                // This includes the dynamical mass from scalar field coupling
                dilationFactor *= Math.Sqrt(1.0 / (1.0 + totalMass * PhysicsConstants.TimeDilationMassCoupling));
            }

            // Factor 3: Local curvature contribution
            double curvature = GetLocalCurvatureNorm(nodeId);
            if (curvature > 0)
            {
                dilationFactor *= Math.Sqrt(1.0 / (1.0 + curvature * PhysicsConstants.TimeDilationCurvatureCoupling));
            }

            return Math.Max(dilationFactor, PhysicsConstants.MinTimeDilation); // Prevent extreme dilation
        }

        /// <summary>
        /// Update scalar field at a single node using Klein-Gordon equation
        /// </summary>
        private void UpdateScalarFieldAtNode(int nodeId, double dt)
        {
            if (ScalarField == null || _scalarMomentum == null)
                return;

            if (nodeId < 0 || nodeId >= ScalarField.Length)
                return;

            double phi = ScalarField[nodeId];
            double pi = _scalarMomentum[nodeId];

            // Compute Laplacian term (discrete)
            double laplacian = 0.0;
            int degree = 0;

            foreach (int j in Neighbors(nodeId))
            {
                if (j < ScalarField.Length)
                {
                    laplacian += (ScalarField[j] - phi) * Weights[nodeId, j];
                    degree++;
                }
            }

            if (degree > 0)
                laplacian /= degree;

            // Mexican Hat potential: V(φ) = -μ²φ² + λφ⁴
            double mu2 = PhysicsConstants.HiggsMuSquared;
            double lambda = PhysicsConstants.HiggsLambda;
            double dV_dphi = -2.0 * mu2 * phi + 4.0 * lambda * phi * phi * phi;

            // Klein-Gordon equation: d²φ/dt² = ∇²φ - dV/dφ
            double d2phi_dt2 = laplacian - dV_dphi;

            // Symplectic Euler integration
            _scalarMomentum[nodeId] += d2phi_dt2 * dt;
            ScalarField[nodeId] += _scalarMomentum[nodeId] * dt;
        }

        /// <summary>
        /// Update gauge fields at a single node for event-driven evolution.
        /// 
        /// CHECKLIST ITEM 8: Integrate gauge updates into local event loop.
        /// 
        /// This method evolves all gauge links emanating from the specified node
        /// using local neighbor information, enabling proper local time updates.
        /// </summary>
        private void UpdateGaugeFieldsAtNode(int nodeId, double dt)
        {
            if (nodeId < 0 || nodeId >= N) return;
            
            // Update SU(3) gluon field links
            if (_gluonField != null && _gluonFieldStrength != null)
            {
                UpdateGluonFieldAtNode(nodeId, dt);
            }
            
            // Update SU(2) weak field links
            if (_weakField != null && _weakFieldStrength != null)
            {
                UpdateWeakFieldAtNode(nodeId, dt);
            }
            
            // Update U(1) hypercharge field links
            if (_hyperchargeField != null && _hyperchargeFieldStrength != null)
            {
                UpdateHyperchargeFieldAtNode(nodeId, dt);
            }
        }
        
        /// <summary>
        /// Update SU(3) gluon field links emanating from a node.
        /// </summary>
        private void UpdateGluonFieldAtNode(int nodeId, double dt)
        {
            int[] scratch = new int[N];
            var neighbors = GetNeighborSpan(nodeId, ref scratch);
            
            foreach (int j in neighbors)
            {
                for (int a = 0; a < 8; a++)
                {
                    // Compute divergence of field strength at this edge
                    double divF = 0.0;
                    foreach (int k in neighbors)
                    {
                        if (Edges[k, j])
                        {
                            divF += _gluonFieldStrength![k, j, a] - _gluonFieldStrength[nodeId, j, a];
                        }
                    }
                    
                    // Self-interaction term: g f^{abc} A^b E^c
                    double selfInt = 0.0;
                    for (int b = 0; b < 8; b++)
                    {
                        double Ab = _gluonField![nodeId, j, b];
                        if (Math.Abs(Ab) < 1e-12) continue;
                        
                        for (int c = 0; c < 8; c++)
                        {
                            if (TryGetStructureConstant(a, b, c, out double fabc))
                            {
                                selfInt += StrongCoupling * fabc * Ab * _gluonFieldStrength![nodeId, j, c];
                            }
                        }
                    }
                    
                    // Color current from matter field
                    double J = ComputeColorCurrentCached(nodeId, j, a);
                    
                    // Update field: dA/dt = ∇×E + self-interaction - J
                    _gluonField![nodeId, j, a] += dt * (divF + selfInt - J);
                }
            }
        }
        
        /// <summary>
        /// Update SU(2) weak field links emanating from a node.
        /// Uses explicit Levi-Civita structure constants per checklist item 2.
        /// </summary>
        private void UpdateWeakFieldAtNode(int nodeId, double dt)
        {
            int[] scratch = new int[N];
            var neighbors = GetNeighborSpan(nodeId, ref scratch);
            
            foreach (int j in neighbors)
            {
                for (int a = 0; a < 3; a++)
                {
                    // Compute divergence
                    double divW = 0.0;
                    foreach (int k in neighbors)
                    {
                        if (Edges[k, j])
                        {
                            divW += _weakFieldStrength![k, j, a] - _weakFieldStrength[nodeId, j, a];
                        }
                    }
                    
                    // Self-interaction using explicit SU(2) structure constants
                    double selfInt = 0.0;
                    for (int b = 0; b < 3; b++)
                    {
                        for (int c = 0; c < 3; c++)
                        {
                            double fabc = Gauge.PauliMatrices.GetStructureConstant(a, b, c);
                            if (fabc != 0.0)
                            {
                                selfInt += WeakCoupling * fabc * _weakField![nodeId, j, b] * _weakFieldStrength![nodeId, j, c];
                            }
                        }
                    }
                    
                    // Weak current
                    double Jw = ComputeWeakCurrent(nodeId, j, a);
                    
                    // Update field
                    _weakField![nodeId, j, a] += dt * (divW + selfInt - Jw);
                }
            }
        }
        
        /// <summary>
        /// Update U(1) hypercharge field links emanating from a node.
        /// </summary>
        private void UpdateHyperchargeFieldAtNode(int nodeId, double dt)
        {
            int[] scratch = new int[N];
            var neighbors = GetNeighborSpan(nodeId, ref scratch);
            
            foreach (int j in neighbors)
            {
                // Compute divergence
                double divB = 0.0;
                foreach (int k in neighbors)
                {
                    if (Edges[k, j])
                    {
                        divB += _hyperchargeFieldStrength![k, j] - _hyperchargeFieldStrength[nodeId, j];
                    }
                }
                
                // Hypercharge current
                double Jh = ComputeHyperchargeCurrent(nodeId, j);
                
                // Update field
                _hyperchargeField![nodeId, j] += dt * (divB - Jh);
            }
        }

        /// <summary>
        /// Update spinor field at a single node using local Dirac equation.
        /// 
        /// FIX: Replaces the empty placeholder with proper local spinor evolution.
        /// 
        /// Physics: The Dirac equation on a graph for node i reads:
        ///   dψ_i/dt = -i/ℏ * H * ψ_i
        /// where H includes:
        ///   1. Kinetic term: sum over neighbors j of (γ·∇)_ij with gauge-covariant derivative
        ///   2. Mass term: m(φ_i) * γ⁰ψ_i where mass comes from scalar field (Yukawa coupling)
        /// 
        /// The mass is NOT hardcoded but computed dynamically from the local scalar field value
        /// using the Yukawa coupling: m_i = g_Y * φ_i (Higgs mechanism).
        /// </summary>
        private void UpdateSpinorFieldAtNode(int nodeId, double dt)
        {
            if (_spinorA == null || nodeId >= _spinorA.Length)
                return;

            double hbar = VectorMath.HBar;
            double c = VectorMath.SpeedOfLight;

            // Get dynamical fermion mass from scalar field (Yukawa coupling)
            // FIX: Mass is NOT manually set but emerges from scalar field interaction
            double mass = ComputeDynamicalFermionMass(nodeId);

            Complex deltaA = Complex.Zero, deltaB = Complex.Zero;
            Complex deltaC = Complex.Zero, deltaD = Complex.Zero;

            bool isEvenSite = (nodeId % 2 == 0);

            // Kinetic term: sum over neighbors with gauge-covariant parallel transport
            foreach (int j in Neighbors(nodeId))
            {
                double weight = Weights[nodeId, j];
                if (weight < 1e-12) continue;

                // Gauge-covariant parallel transport for U(1)
                Complex parallelTransport = Complex.One;
                if (_edgePhaseU1 != null)
                {
                    double phase = _edgePhaseU1[nodeId, j];
                    parallelTransport = Complex.FromPolarCoordinates(1.0, -phase); // U† = e^{-iθ}
                }

                // Staggered fermion sign
                bool isNeighborEven = (j % 2 == 0);
                double sign = (isEvenSite != isNeighborEven) ? 1.0 : -1.0;

                // Direction flavor based on edge parity
                int edgeDirection = Math.Abs(nodeId - j) % 2;

                // Apply parallel transport to neighbor spinor
                Complex gaugedA_j = _spinorA[j] * parallelTransport;
                Complex gaugedB_j = _spinorB![j] * parallelTransport;
                Complex gaugedC_j = _spinorC![j] * parallelTransport;
                Complex gaugedD_j = _spinorD![j] * parallelTransport;

                // Hopping terms with staggered signs
                if (edgeDirection == 0)
                {
                    // X-like direction: couples A↔B and C↔D
                    deltaB += sign * weight * (gaugedA_j - _spinorA[nodeId]);
                    deltaA += sign * weight * (gaugedB_j - _spinorB[nodeId]);
                    deltaD += sign * weight * (gaugedC_j - _spinorC[nodeId]);
                    deltaC += sign * weight * (gaugedD_j - _spinorD[nodeId]);
                }
                else
                {
                    // Y-like direction: couples A↔B with i and C↔D with -i
                    Complex iSign = Complex.ImaginaryOne * sign;
                    deltaB += iSign * weight * (gaugedA_j - _spinorA[nodeId]);
                    deltaA -= iSign * weight * (gaugedB_j - _spinorB[nodeId]);
                    deltaD -= iSign * weight * (gaugedC_j - _spinorC[nodeId]);
                    deltaC += iSign * weight * (gaugedD_j - _spinorD[nodeId]);
                }
            }

            // Mass term: couples left and right handed components
            // m_dynamic comes from Yukawa coupling to scalar field
            double mc = mass * c;
            Complex massTermA = -Complex.ImaginaryOne * mc / hbar * _spinorC![nodeId];
            Complex massTermB = -Complex.ImaginaryOne * mc / hbar * _spinorD![nodeId];
            Complex massTermC = -Complex.ImaginaryOne * mc / hbar * _spinorA[nodeId];
            Complex massTermD = -Complex.ImaginaryOne * mc / hbar * _spinorB![nodeId];

            // Dirac evolution: dψ/dt = -i/ℏ * H * ψ
            double factor = -1.0 / hbar;

            // Apply local update (simple Euler for single node - full step uses leapfrog)
            _spinorA[nodeId] += dt * factor * (c * deltaA + massTermA);
            _spinorB[nodeId] += dt * factor * (c * deltaB + massTermB);
            _spinorC[nodeId] += dt * factor * (c * deltaC + massTermC);
            _spinorD[nodeId] += dt * factor * (c * deltaD + massTermD);
        }

        /// <summary>
        /// Compute dynamical fermion mass from scalar field via Yukawa coupling.
        /// 
        /// FIX: Replaces hardcoded _correlationMass with proper Higgs mechanism.
        /// 
        /// Physics: In the Standard Model, fermion masses arise from Yukawa coupling:
        ///   m_f = g_Y * ⟨φ⟩
        /// where g_Y is the Yukawa coupling constant and ⟨φ⟩ is the Higgs VEV.
        /// 
        /// Here we use the local scalar field value φ_i as the dynamical "Higgs" field.
        /// The mass smoothly varies as the scalar field evolves.
        /// </summary>
        private double ComputeDynamicalFermionMass(int nodeId)
        {
            // Yukawa coupling constant
            const double g_Yukawa = 0.1;

            // Get scalar field value (Higgs-like field)
            double phi = 0.0;
            if (ScalarField != null && nodeId < ScalarField.Length)
            {
                phi = ScalarField[nodeId];
            }

            // Higgs VEV for comparison
            double v = PhysicsConstants.HiggsVEV;

            // Mass from Yukawa coupling: m = g_Y * |φ|
            // Near the vacuum (φ ≈ v), this gives the "natural" fermion mass
            double mass = g_Yukawa * Math.Abs(phi);

            // Add small bare mass for numerical stability
            const double bareMass = 0.001;
            return bareMass + mass;
        }

        /// <summary>
        /// Update node state based on local physics AND neighbor excitation.
        /// 
        /// FIX: Original code only checked LocalPotential which stays near 0
        /// in Event-Based mode, causing all nodes to remain in Rest state.
        /// 
        /// Now uses same logic as StepExcitableMedium:
        /// - Rest -> Excited: based on excited neighbors (weighted sum)
        /// - Excited -> Refractory: immediate transition
        /// - Refractory -> Rest: after counter expires
        /// </summary>
        private void UpdateNodeState(int nodeId)
        {
            if (nodeId < 0 || nodeId >= N)
                return;

            // === REFRACTORY STATE: count down and transition to Rest ===
            if (State[nodeId] == NodeState.Refractory)
            {
                _refractoryCounter[nodeId]--;
                if (_refractoryCounter[nodeId] <= 0)
                {
                    State[nodeId] = NodeState.Rest;
                    _refractoryCounter[nodeId] = 0;
                }
                return;
            }

            // === EXCITED STATE: transition to Refractory ===
            if (State[nodeId] == NodeState.Excited)
            {
                State[nodeId] = NodeState.Refractory;
                // Refractory period scales with correlation mass (heavy nodes stay refractory longer)
                double massFactor = (_correlationMass != null && _correlationMass.Length == N && _avgCorrelationMass > 0)
                    ? _correlationMass[nodeId] / _avgCorrelationMass
                    : 0.0;
                int extraSteps = (int)Math.Round(massFactor);
                _refractoryCounter[nodeId] = DynamicBaseRefractorySteps + Math.Max(0, extraSteps);
                return;
            }

            // === REST STATE: probabilistic excitation based on neighbors ===
            // Count excited neighbors and compute weighted excitation sum
            int excitedNeighbors = 0;
            double weightedExcitation = 0.0;
            int degree = 0;

            foreach (int j in Neighbors(nodeId))
            {
                degree++;
                if (State[j] == NodeState.Excited)
                {
                    excitedNeighbors++;
                    weightedExcitation += Weights[nodeId, j];
                }
            }

            // Spontaneous excitation probability based on local curvature
            // Higher curvature (gravity well) = higher spontaneous excitation
            double spontProb = 0.0;
            {
                double kNorm = GetLocalCurvatureNorm(nodeId);
                spontProb = 1.0 - Math.Exp(-kNorm);
                spontProb = Math.Clamp(spontProb, 0.0, 0.5);
            }

            // Neighbor-driven excitation probability
            double neighborProb = 0.0;
            if (excitedNeighbors > 0 && degree > 0)
            {
                double meanW = weightedExcitation / excitedNeighbors;
                double density = (double)excitedNeighbors / degree;
                neighborProb = 1.0 - Math.Exp(-density * meanW);
            }

            // Combined probability with LocalPotential boost if available
            double totalProb = spontProb + neighborProb;

            // LocalPotential boost (if initialized)
            if (LocalPotential != null && nodeId < LocalPotential.Length && LocalPotential[nodeId] > 0)
            {
                totalProb *= (1.0 + LocalPotential[nodeId]);
            }

            // Clamp and apply
            totalProb = Math.Clamp(totalProb, 0.0, 0.99);

            // Probabilistic transition
            if (_rng.NextDouble() < totalProb)
            {
                State[nodeId] = NodeState.Excited;
            }
        }

        /// <summary>
        /// Get gauge field component along edge (i,j)
        /// Used for electric field divergence calculation
        /// </summary>
        public double GetGaugeFieldComponent(int i, int j)
        {
            if (!Edges[i, j])
                return 0.0;

            // Return U(1) phase as field component
            if (_edgePhaseU1 != null && i < _edgePhaseU1.GetLength(0) && j < _edgePhaseU1.GetLength(1))
            {
                return _edgePhaseU1[i, j];
            }

            // Fallback: use gluon field if available
            if (_gluonField != null && i < _gluonField.GetLength(0) && j < _gluonField.GetLength(1))
            {
                // Return first color component
                return _gluonField[i, j, 0];
            }

            return 0.0;
        }

        /// <summary>
        /// Apply gauge transformation to edge (i,j) consistently.
        /// 
        /// FIX: Original code only updated one direction of the edge phase,
        /// violating gauge invariance (Wilson loop unitarity).
        /// 
        /// Physics: Under a U(1) gauge transformation with parameter χ:
        ///   θ_ij → θ_ij - (χ_i - χ_j)
        ///   θ_ji → θ_ji - (χ_j - χ_i) = θ_ji + (χ_i - χ_j)
        /// 
        /// Since θ_ji = -θ_ij for antisymmetric gauge links, both directions
        /// must be updated consistently to preserve:
        ///   θ_ij + θ_ji = 0  (antisymmetry)
        ///   Π_{plaquette} exp(iθ) = const  (Wilson loop invariance)
        /// 
        /// chi_diff = χ_i - χ_j (difference of gauge parameters at nodes)
        /// </summary>
        public void ApplyGaugeTransformation(int i, int j, double chi_diff)
        {
            if (!Edges[i, j])
                return;

            // Apply transformation to U(1) phase - both directions
            if (_edgePhaseU1 != null && i < _edgePhaseU1.GetLength(0) && j < _edgePhaseU1.GetLength(1))
            {
                // θ_ij → θ_ij - χ_diff
                _edgePhaseU1[i, j] -= chi_diff;

                // θ_ji → θ_ji + χ_diff (antisymmetric: θ_ji = -θ_ij)
                // FIX: This was missing! Both directions must be updated.
                _edgePhaseU1[j, i] += chi_diff;

                // Normalize both to [-π, π] for numerical stability
                _edgePhaseU1[i, j] = NormalizeAngle(_edgePhaseU1[i, j]);
                _edgePhaseU1[j, i] = NormalizeAngle(_edgePhaseU1[j, i]);

                // Verify antisymmetry (diagnostic)
                // After normalization, θ_ij + θ_ji should be ≈ 0 or ±2π
            }

            // Also apply to SU(3) gluon field if present
            if (_gluonField != null && i < _gluonField.GetLength(0) && j < _gluonField.GetLength(1))
            {
                // For SU(3), gauge transformation is more complex: U → g_i U g_j†
                // Here we apply a simplified Abelian subgroup transformation
                double scale = 0.1; // Coupling to Abelian component
                for (int a = 0; a < 8; a++)
                {
                    // Transform both directions consistently
                    _gluonField[i, j, a] -= chi_diff * scale;
                    _gluonField[j, i, a] += chi_diff * scale;
                }
            }
        }

        /// <summary>
        /// Normalize angle to [-π, π] range.
        /// </summary>
        private static double NormalizeAngle(double angle)
        {
            while (angle > Math.PI) angle -= 2 * Math.PI;
            while (angle < -Math.PI) angle += 2 * Math.PI;
            return angle;
        }

        /// <summary>
        /// Update local geometry (gravity) for a node and its edges.
        /// Uses Ollivier-Ricci curvature and volume constraint.
        /// </summary>
        public void UpdateLocalGeometry(int node, double dt)
        {
            if (node < 0 || node >= N) return;

            // Update edges connected to this node
            // We iterate neighbors and update the edge (node, neighbor)
            // Note: EvolveGeometry updates symmetric weights
            double learningRate = PhysicsConstants.GravitationalCoupling * dt;
            
            foreach (int neighbor in Neighbors(node))
            {
                // Only update if node < neighbor to avoid double counting in a full sweep,
                // but here we are updating a specific node's local geometry.
                // To be safe and local, we update all connected edges.
                // EvolveGeometry handles symmetry.
                EvolveGeometry(node, neighbor, learningRate);
            }
        }

        /// <summary>
        /// Wrapper for UpdateNodePhysics to match the requested API.
        /// Updates scalar, gauge, spinor fields and node state.
        /// </summary>
        public void UpdateLocalFields(int node, double dt)
        {
            UpdateNodePhysics(node, dt);
        }

        /// <summary>
        /// Perform topological updates using Metropolis-Hastings with local action.
        /// Includes bipartite check for staggered fermions.
        /// </summary>
        public void MetropolisTopologyUpdate_Checked()
        {
            // Perform a batch of local updates
            // Number of attempts proportional to N to ensure coverage
            int attempts = Math.Max(10, N / 10);
            
            for (int i = 0; i < attempts; i++)
            {
                MetropolisEdgeStepLocalAction();
            }
        }

        /// <summary>
        /// Calculate total energy of the system from all components.
        /// Used for conservation validation.
        /// </summary>
        public double CalculateTotalEnergy()
        {
            double total = 0.0;
            
            // 1. Matter Energy (from NodeMasses)
            if (_nodeMasses != null)
            {
                for (int i = 0; i < N; i++)
                {
                    total += _nodeMasses[i].TotalMass;
                }
            }
            else
            {
                // Fallback if NodeMasses not initialized
                total += ComputeTotalEnergy();
            }
            
            // 2. Geometry Energy (Kinetic)
            total += ComputeGeometryKineticEnergy();
            
            // 3. Vacuum Pool (from Ledger)
            total += Ledger.VacuumPool;
            
            return total;
        }
    }
}
