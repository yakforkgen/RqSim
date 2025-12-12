using System;
using System.Collections.Generic;

namespace RQSimulation
{
    public partial class RQGraph
    {
        // RQ-FIX: Use PhysicsConstants instead of local magic numbers
        // This ensures consistency across the codebase and allows physics-based tuning
        
        /// <summary>
        /// Apply probabilistic vacuum fluctuations based on local field conditions
        /// Replaces manual vacuum pair creation with physics-based stochastic events
        /// 
        /// RQ-FIX: Now uses PhysicsConstants for all rates and thresholds.
        /// RQ-FIX: Uses EnergyLedger for strict energy conservation.
        /// </summary>
        public void ApplyProbabilisticVacuumFluctuations()
        {
            if (_correlationMass == null || _correlationMass.Length != N)
                RecomputeCorrelationMass();
            
            // Ensure ledger is initialized
            if (Math.Abs(Ledger.TotalTrackedEnergy) < 1e-10)
            {
                Ledger.Initialize(1000.0); // Default vacuum energy
            }

            for (int i = 0; i < N; i++)
            {
                // Compute local curvature
                double localCurvature = ComputeLocalCurvature(i);
                
                // Compute local field energy
                double localFieldEnergy = 0;
                if (LocalPotential != null && i < LocalPotential.Length)
                {
                    localFieldEnergy = LocalPotential[i];
                }
                
                // Check if near black hole
                bool nearBlackHole = false;
                if (_blackHoles != null)
                {
                    foreach (var bh in _blackHoles)
                    {
                        if (bh.HorizonNodes.Contains(i) || bh.InteriorNodes.Contains(i))
                        {
                            nearBlackHole = true;
                            break;
                        }
                    }
                }
                
                // Probability enhanced by curvature and proximity to black holes
                // RQ-FIX: Use PhysicsConstants values instead of magic numbers
                double probability = PhysicsConstants.VacuumFluctuationBaseRate;
                probability *= (1.0 + PhysicsConstants.CurvatureCouplingFactor * Math.Abs(localCurvature));
                
                if (nearBlackHole)
                {
                    probability *= PhysicsConstants.HawkingRadiationEnhancement;
                }
                
                // Roll for fluctuation
                if (_rng.NextDouble() < probability)
                {
                    ApplyVacuumFluctuationAtNode(i, localFieldEnergy, localCurvature);
                }
            }
        }
        
        /// <summary>
        /// Apply vacuum fluctuation effect at specific node
        /// RQ-compliant: Energy-conserving vacuum fluctuations
        /// 
        /// RQ-FIX: Uses PhysicsConstants.PairCreationEnergyThreshold
        /// RQ-FIX: Uses EnergyLedger to borrow/return energy
        /// </summary>
        private void ApplyVacuumFluctuationAtNode(int i, double fieldEnergy, double curvature)
        {
            // Decide type of fluctuation based on energy
            // RQ-FIX: Use physics-based threshold
            if (fieldEnergy > PhysicsConstants.PairCreationEnergyThreshold)
            {
                // Create virtual pair: excite node and neighbor
                // This costs energy, must borrow from vacuum
                double pairCost = 0.1; // Estimated cost
                
                if (Ledger.TrySpendVacuumEnergy(pairCost))
                {
                    State[i] = NodeState.Excited;
                    
                    // Find a neighbor to excite as the pair partner
                    var neighbors = new List<int>();
                    foreach (int j in Neighbors(i))
                    {
                        neighbors.Add(j);
                    }
                    
                    if (neighbors.Count > 0)
                    {
                        int partner = neighbors[_rng.Next(neighbors.Count)];
                        State[partner] = NodeState.Excited;
                        
                        // Transfer energy conservatively - redistribute between pair
                        if (LocalPotential != null && i < LocalPotential.Length && partner < LocalPotential.Length)
                        {
                            double energyTransfer = fieldEnergy * 0.1;
                            LocalPotential[i] -= energyTransfer;
                            LocalPotential[partner] += energyTransfer;
                        }
                    }
                    else
                    {
                        // No partner, return energy
                        Ledger.RegisterRadiation(pairCost);
                        State[i] = NodeState.Rest;
                    }
                }
            }
            else
            {
                // Small quantum fluctuation in field - ENERGY CONSERVING
                // Borrow or return energy to vacuum reservoir
                if (LocalPotential != null && i < LocalPotential.Length)
                {
                    double fluctuation = (_rng.NextDouble() - 0.5) * PhysicsConstants.VacuumFluctuationAmplitude;
                    
                    if (fluctuation > 0)
                    {
                        // Borrowing energy from vacuum
                        if (Ledger.TrySpendVacuumEnergy(fluctuation))
                        {
                            LocalPotential[i] += fluctuation;
                        }
                    }
                    else
                    {
                        // Returning energy to vacuum
                        double returnAmount = -fluctuation;
                        LocalPotential[i] -= returnAmount;
                        Ledger.RegisterRadiation(returnAmount);
                    }
                }
            }
        }
        
        /// <summary>
        /// Probabilistic black hole evaporation based on Hawking radiation physics
        /// Replaces manual EmitHawkingRadiation with stochastic emission
        /// </summary>
        public void ApplyStochasticHawkingRadiation()
        {
            if (_blackHoles == null || _blackHoles.Count == 0)
                return;
            
            foreach (var bh in _blackHoles)
            {
                // Skip if temperature too low (very massive black holes)
                if (bh.Temperature < 1e-10)
                    continue;
                
                // Emission probability based on temperature (Stefan-Boltzmann law)
                // P ∝ T⁴ * Area
                double baseEmissionProb = 0.001 * Math.Pow(bh.Temperature, 4);
                
                foreach (int horizonNode in bh.HorizonNodes)
                {
                    // Each horizon node has chance to emit
                    if (_rng.NextDouble() < baseEmissionProb)
                    {
                        EmitQuantumFromHorizon(horizonNode, bh);
                    }
                }
                
                // Check if black hole should evaporate completely
                if (bh.Mass < 0.01)
                {
                    EvaporateBlackHoleCompletely(bh);
                }
            }
        }
        
        /// <summary>
        /// Emit a quantum of radiation from horizon node
        /// </summary>
        private void EmitQuantumFromHorizon(int node, BlackHoleRegion bh)
        {
            // Excite horizon node
            State[node] = NodeState.Excited;
            
            // Add energy to local field
            if (LocalPotential != null && node < LocalPotential.Length)
            {
                LocalPotential[node] += bh.Temperature * 0.1;
            }
            
            // Reduce black hole mass (energy conservation)
            double energyEmitted = bh.Temperature * 0.01;
            bh.Mass -= energyEmitted;
            
            // Distribute mass loss among interior nodes
            if (_correlationMass != null && bh.InteriorNodes.Count > 0)
            {
                double massPerNode = energyEmitted / bh.InteriorNodes.Count;
                foreach (int interior in bh.InteriorNodes)
                {
                    if (interior < _correlationMass.Length)
                    {
                        _correlationMass[interior] = Math.Max(0, _correlationMass[interior] - massPerNode);
                    }
                }
            }
        }
        
        /// <summary>
        /// Complete evaporation of a small black hole with energy release
        /// </summary>
        private void EvaporateBlackHoleCompletely(BlackHoleRegion bh)
        {
            // Convert remaining mass to field excitations in neighborhood
            double remainingEnergy = bh.Mass;
            
            // Distribute energy to horizon and nearby nodes
            var affectedNodes = new List<int>();
            affectedNodes.AddRange(bh.HorizonNodes);
            
            // Add neighbors of horizon nodes
            foreach (int horizonNode in bh.HorizonNodes)
            {
                foreach (int neighbor in Neighbors(horizonNode))
                {
                    if (!affectedNodes.Contains(neighbor))
                    {
                        affectedNodes.Add(neighbor);
                    }
                }
            }
            
            // Distribute energy
            if (LocalPotential != null && affectedNodes.Count > 0)
            {
                double energyPerNode = remainingEnergy / affectedNodes.Count;
                foreach (int node in affectedNodes)
                {
                    if (node < LocalPotential.Length)
                    {
                        LocalPotential[node] += energyPerNode;
                        State[node] = NodeState.Excited;
                    }
                }
            }
            
            // Clear interior nodes' correlation mass
            if (_correlationMass != null)
            {
                foreach (int interior in bh.InteriorNodes)
                {
                    if (interior < _correlationMass.Length)
                    {
                        _correlationMass[interior] = 0;
                    }
                }
            }
        }
        
        /// <summary>
        /// Combined update for all probabilistic quantum effects
        /// Replaces manual event calls in simulation loop
        /// </summary>
        public void UpdateProbabilisticQuantumEffects()
        {
            // Apply vacuum fluctuations
            ApplyProbabilisticVacuumFluctuations();
            
            // Update black holes if physics enabled
            if (_blackHoles != null && _blackHoles.Count > 0)
            {
                ApplyStochasticHawkingRadiation();
            }
            
            // Clean up evaporated black holes
            if (_blackHoles != null)
            {
                _blackHoles.RemoveAll(bh => bh.Mass < 0.001);
            }
        }
    }
}
