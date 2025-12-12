// RQSimulation/Experiments/Definitions/BioFoldingExperiment.cs
using System;

namespace RQSimulation.Experiments.Definitions
{
    /// <summary>
    /// Experiment C: "DNA Folding" (ДНК-Фолдинг / Bio-Graphity)
    /// 
    /// A special mode modeling topological self-organization of a linear chain
    /// under "gravity" (analogous to Van der Waals / hydrogen bond forces).
    /// 
    /// Sequence: DNA hairpin
    /// Code: 5'-GCGCAAAAGCGC-3' (12 bases)
    /// - GCGC (stem, strong bond)
    /// - AAAA (loop, flexible)
    /// - GCGC (stem, complementary to start)
    /// 
    /// What to observe:
    /// - Spectral Dimension: Start at dS ≈ 1.0 (line), end at dS ≈ 2.0-3.0 (folded)
    /// - Average Path Length: Start at 11, end at ~4-5 (if hairpin folds properly)
    /// - New edges: Check for 0-11, 1-10 connections (complementary pairing)
    /// </summary>
    public class BioFoldingExperiment : IExperiment
    {
        public string Name => "Bio-Folding (DNA Hairpin)";
        
        public string Description => 
            "Models topological self-organization of a 12-nucleotide DNA hairpin " +
            "(5'-GCGCAAAAGCGC-3') under graph 'gravity'. Starts as linear chain (dS ≈ 1), " +
            "should fold into hairpin structure (dS ≈ 2-3). Watch for complementary base pairing.";
        
        public StartupConfig GetConfig()
        {
            // Physical intent: keep the DNA backbone intact while gentle gravity nudges bases into a hairpin without crushing it.
            return new StartupConfig
            {
                // Exact number of nucleotides
                NodeCount = 12,
                
                // Allow more time for folding without increasing heat
                TotalSteps = 2000,
                
                // IMPORTANT: No random edges - use CustomInitializer
                InitialEdgeProb = 0.0,
                
                // Moderate coupling so hydrogen bonds form but do not collapse backbone
                GravitationalCoupling = 0.5,
                
                // Cold start - don't break the chain, only bend it
                HotStartTemperature = 0.1,
                
                // No cooling schedule – structure evolves at fixed near-zero T
                AnnealingCoolingRate = 1.0,
                
                // No spontaneous excitation
                InitialExcitedProb = 0.0,
                
                // Linear chain has degree 2 (except endpoints)
                TargetDegree = 2,
                Temperature = 0.5,
                LambdaState = 0.3,
                EdgeTrialProbability = 0.05, // Allow new edge formation for folding
                DecoherenceRate = 0.001,
                WarmupDuration = 50,
                GravityTransitionDuration = 100,
                
                // Enable spectral geometry for dS measurement
                UseSpectralGeometry = true,
                UseNetworkGravity = true,
                UseHotStartAnnealing = false, // Keep cold
                UseQuantumDrivenStates = true,
                UseVacuumFluctuations = false, // No random fluctuations
                UseSpinorField = false,
                UseTopologicalProtection = true, // Prevents backbone breakage during folding
                
                // No fractal topology
                FractalLevels = 0,
                FractalBranchFactor = 0
            };
        }
        
        public void ApplyPhysicsOverrides()
        {
            // Physics override: External pressure to help folding
            // CosmologicalConstant = 0.05 provides inward pressure
            // Note: PhysicsConstants.CosmologicalConstant is readonly = 0.001
            // The effect is simulated through high GravitationalCoupling
        }
        
        /// <summary>
        /// Custom initializer creates a linear chain topology:
        /// 0-1-2-3-4-5-6-7-8-9-10-11
        /// with weight 1.0 (unbreakable covalent backbone)
        /// </summary>
        public Action<RQGraph>? CustomInitializer => InitializeLinearChain;
        
        private static void InitializeLinearChain(RQGraph graph)
        {
            int n = graph.N;
            if (n < 2) return;
            
            // Clear any existing edges (graph may have random initialization)
            // Note: For this experiment n=12, so O(n²) = O(144) is negligible.
            // For larger graphs, consider using FlatEdgesFrom/FlatEdgesTo arrays.
            for (int i = 0; i < n; i++)
            {
                for (int j = i + 1; j < n; j++) // Only check upper triangle (undirected graph)
                {
                    if (graph.Edges[i, j])
                    {
                        graph.RemoveEdge(i, j);
                    }
                }
            }
            
            // Create linear chain: 0-1-2-3-...-n-1
            for (int i = 0; i < n - 1; i++)
            {
                graph.AddEdge(i, i + 1);
                // Set high weight (covalent bond - should not break)
                graph.Weights[i, i + 1] = 1.0;
                graph.Weights[i + 1, i] = 1.0;
            }
            
            // Initialize coordinates in a line for visualization.
            // IMPORTANT: Coordinates property is marked obsolete because it should
            // NOT be used for physics calculations (use graph distances instead).
            // However, it IS the correct API for visualization setup, which is
            // exactly what we're doing here. The obsolete warning is intentionally
            // suppressed for this legitimate rendering use case.
            double spacing = 0.8 / Math.Max(1, n - 1);
            for (int i = 0; i < n; i++)
            {
#pragma warning disable CS0618 // Coordinates are used correctly here for visualization setup
                if (graph.Coordinates != null && graph.Coordinates.Length == n)
                {
                    graph.Coordinates[i] = (-0.4 + i * spacing, 0.0);
                }
#pragma warning restore CS0618
            }
        }
    }
}
