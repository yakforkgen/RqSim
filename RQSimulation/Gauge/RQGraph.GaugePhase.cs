using System;
using System.Numerics;

namespace RQSimulation
{
    public partial class RQGraph
    {
        // Edge gauge phases for U(1) gauge field
        private double[,]? _edgePhaseU1;
        
        // Configuration constants
        private const double GaugeCouplingConstant = 0.1;
        private const double PlaquetteWeight = 0.05;
        
        /// <summary>
        /// Initialize gauge phases on edges
        /// </summary>
        public void InitEdgeGaugePhases()
        {
            _edgePhaseU1 = new double[N, N];
            
            // Initialize with small random phases
            var rng = new Random();
            for (int i = 0; i < N; i++)
            {
                for (int j = i + 1; j < N; j++)
                {
                    if (Edges[i, j])
                    {
                        double phase = (rng.NextDouble() - 0.5) * 0.1; // Small random phase
                        _edgePhaseU1[i, j] = phase;
                        _edgePhaseU1[j, i] = -phase; // Antisymmetric
                    }
                }
            }
        }
        
        /// <summary>
        /// Get edge gauge data for edge (i,j)
        /// </summary>
        public EdgeGaugeData GetEdgeGaugeData(int i, int j)
        {
            if (_edgePhaseU1 == null)
                InitEdgeGaugePhases();
            
            double weight = Edges[i, j] ? Weights[i, j] : 0.0;
            double phase = _edgePhaseU1?[i, j] ?? 0.0;
            
            return new EdgeGaugeData(weight, phase);
        }
        
        /// <summary>
        /// Update gauge phases based on currents and field strength
        /// </summary>
        public void UpdateGaugePhases(double dt)
        {
            if (_edgePhaseU1 == null)
                InitEdgeGaugePhases();
            
            double[,] phaseUpdates = new double[N, N];
            
            // Compute plaquette contributions (field curvature)
            for (int i = 0; i < N; i++)
            {
                foreach (int j in Neighbors(i))
                {
                    if (i >= j) continue; // Only process each edge once
                    
                    // Find minimal cycles (plaquettes) through this edge
                    double curvatureSum = 0;
                    int plaquetteCount = 0;
                    
                    foreach (int k in Neighbors(i))
                    {
                        if (k == j) continue;
                        if (Edges[j, k]) // Triangle i-j-k
                        {
                            // Sum phases around plaquette
                            double phaseSum = _edgePhaseU1[i, j] 
                                            + _edgePhaseU1[j, k] 
                                            + _edgePhaseU1[k, i];
                            
                            // Field strength (curvature) proportional to sin(phase sum)
                            curvatureSum += Math.Sin(phaseSum);
                            plaquetteCount++;
                        }
                    }
                    
                    // Average curvature contribution
                    double avgCurvature = plaquetteCount > 0 ? curvatureSum / plaquetteCount : 0;
                    
                    // Compute current through edge (from fermion density if available)
                    double current = ComputeEdgeCurrent(i, j);
                    
                    // Checklist E.2: Add scalar field back-reaction current
                    // J_ij = g * φ_i * φ_j * sin(θ_ij) from scalar field coupling to gauge field
                    // This ensures the scalar field affects gauge phase evolution (Higgs mechanism)
                    double scalarCurrent = ComputeScalarFieldCurrent(i, j);
                    current += scalarCurrent;
                    
                    // Phase evolution: dφ/dt = -g * (J_total + curl)
                    // J_total = J_fermion + J_scalar (back-reaction from matter fields)
                    double dPhase = -GaugeCouplingConstant * (current + PlaquetteWeight * avgCurvature) * dt;
                    
                    phaseUpdates[i, j] = dPhase;
                    phaseUpdates[j, i] = -dPhase; // Antisymmetric
                }
            }
            
            // Apply updates and wrap to [-π, π]
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    if (Edges[i, j])
                    {
                        _edgePhaseU1[i, j] += phaseUpdates[i, j];
                        
                        // Wrap to [-π, π]
                        while (_edgePhaseU1[i, j] > Math.PI) _edgePhaseU1[i, j] -= 2 * Math.PI;
                        while (_edgePhaseU1[i, j] < -Math.PI) _edgePhaseU1[i, j] += 2 * Math.PI;
                    }
                }
            }
        }
        
        /// <summary>
        /// Compute gauge field energy (sum over plaquettes)
        /// </summary>
        public double ComputeGaugeFieldEnergy()
        {
            if (_edgePhaseU1 == null)
                return 0.0;
            
            double energy = 0;
            int plaquetteCount = 0;
            
            // Sum over all triangular plaquettes
            for (int i = 0; i < N; i++)
            {
                foreach (int j in Neighbors(i))
                {
                    if (j <= i) continue;
                    
                    foreach (int k in Neighbors(i))
                    {
                        if (k <= j || !Edges[j, k]) continue;
                        
                        // Triangle i-j-k forms a plaquette
                        double phaseSum = _edgePhaseU1[i, j] 
                                        + _edgePhaseU1[j, k] 
                                        + _edgePhaseU1[k, i];
                        
                        // Energy: E_gauge = sum(1 - cos(phase around plaquette))
                        energy += 1.0 - Math.Cos(phaseSum);
                        plaquetteCount++;
                    }
                }
            }
            
            return energy;
        }
        
        /// <summary>
        /// Compute current through edge based on fermion/excitation flow
        /// </summary>
        private double ComputeEdgeCurrent(int i, int j)
        {
            // Current based on density difference
            double densityI = State[i] == NodeState.Excited ? 1.0 : 0.0;
            double densityJ = State[j] == NodeState.Excited ? 1.0 : 0.0;
            
            // Add contribution from local potential if available
            if (LocalPotential != null && i < LocalPotential.Length && j < LocalPotential.Length)
            {
                densityI += LocalPotential[i] * 0.1;
                densityJ += LocalPotential[j] * 0.1;
            }
            
            // Current flows from high to low density
            return (densityI - densityJ) * Weights[i, j];
        }

        // NOTE: ComputeScalarFieldCurrent is implemented in RQGraph.FieldTheory.cs (line 229)
        // It computes: J_ij = g * φ_i * φ_j * sin(θ_ij)
        // Do not duplicate here - use the method from FieldTheory which is gauge-covariant
    }
}
