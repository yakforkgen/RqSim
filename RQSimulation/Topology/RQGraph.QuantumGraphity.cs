using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;

namespace RQSimulation
{
    /// <summary>
    /// Quantum Graphity dynamics: topology evolution through action minimization.
    /// Uses Metropolis-Hastings algorithm instead of deterministic state rules.
    /// </summary>
    public partial class RQGraph
    {
        // Network Hamiltonian parameters
        private double _linkCostCoeff = 0.1;    // Penalty for having links
        private double _lengthCostCoeff = 0.05; // Penalty for long links
        private double _matterCouplingCoeff = 0.2; // Matter-geometry coupling

        /// <summary>
        /// Cosmological constant (Lambda) - prevents graph collapse/explosion.
        /// Positive Lambda favors expansion (de Sitter), negative favors contraction.
        /// The term Lambda * V is added to the Hamiltonian where V is effective volume.
        /// REDUCED from 0.01 to 0.001 to prevent excessive edge deletion.
        /// </summary>
        private double _cosmologicalConstant = 0.001;

        /// <summary>
        /// Gets or sets the cosmological constant (Lambda).
        /// Positive values resist collapse, negative values resist expansion.
        /// </summary>
        public double CosmologicalConstant
        {
            get => _cosmologicalConstant;
            set => _cosmologicalConstant = value;
        }

        // Effective temperature for Metropolis-Hastings
        private double _networkTemperature = 1.0;

        // Track energy for optimization
        private double _lastNetworkEnergy = double.MaxValue;

        /// <summary>
        /// Network temperature for Metropolis-Hastings (controls acceptance of uphill moves)
        /// </summary>
        public double NetworkTemperature
        {
            get => _networkTemperature;
            set => _networkTemperature = Math.Max(0.001, value);
        }

        /// <summary>
        /// Compute the network Hamiltonian H = H_links + H_nodes
        /// H_links: Cost of having edges (prefers sparse graphs)
        /// H_nodes: Matter contribution (correlation mass)
        /// </summary>
        public double ComputeNetworkHamiltonian()
        {
            double H_links = 0.0;
            double H_nodes = 0.0;

            // H_links: Sum over all edges of (1 - w_ij) weighted by link cost
            // This penalizes weak links more than strong ones
            int edgeCount = 0;
            for (int i = 0; i < N; i++)
            {
                for (int j = i + 1; j < N; j++)
                {
                    if (!Edges[i, j]) continue;
                    edgeCount++;

                    double w = Weights[i, j];

                    // Link existence cost (sparse graph preference)
                    H_links += _linkCostCoeff * (1.0 - w * w);

                    // Length cost (strong links = short effective distance)
                    double effectiveLength = w > 1e-10 ? 1.0 / w : 10.0;
                    H_links += _lengthCostCoeff * effectiveLength;
                }
            }

            // H_nodes: Matter contribution
            if (_correlationMass != null && _correlationMass.Length == N)
            {
                for (int i = 0; i < N; i++)
                {
                    // Mass-curvature coupling (Einstein-like)
                    double mass = _correlationMass[i];
                    double curvature = GetLocalCurvature(i);

                    // Matter wants to curve space (positive contribution when uncurved)
                    H_nodes += _matterCouplingCoeff * mass * (1.0 - Math.Abs(curvature));
                }
            }

            // Add string energy contribution if present
            if (_stringEnergy != null && _stringEnergy.GetLength(0) == N)
            {
                for (int i = 0; i < N; i++)
                {
                    for (int j = i + 1; j < N; j++)
                    {
                        H_links += _stringEnergy[i, j];
                    }
                }
            }

            return H_links + H_nodes;
        }

        /// <summary>
        /// Compute local energy contribution for a single edge.
        /// Used for efficient delta-energy calculation in Metropolis step.
        /// </summary>
        private double ComputeEdgeEnergy(int i, int j)
        {
            if (!Edges[i, j]) return 0.0;

            double w = Weights[i, j];

            double energy = _linkCostCoeff * (1.0 - w * w);
            double effectiveLength = w > 1e-10 ? 1.0 / w : 10.0;
            energy += _lengthCostCoeff * effectiveLength;

            if (_stringEnergy != null && i < N && j < N)
                energy += _stringEnergy[i, j];

            return energy;
        }

        /// <summary>
        /// Compute local Hamiltonian contribution for a single edge (i,j) and its neighbors.
        /// This includes: gravity (Ricci curvature), matter coupling, and cosmological constant.
        /// Used for efficient O(degree) delta-energy calculation instead of O(N^2) full recalculation.
        /// </summary>
        /// <param name="i">First node index</param>
        /// <param name="j">Second node index</param>
        /// <returns>Local energy contribution from edge (i,j) and affected neighbors</returns>
        public double ComputeLocalHamiltonian(int i, int j)
        {
            double H_local = 0.0;

            // 1. Edge energy contribution (link cost, length cost)
            if (Edges[i, j])
            {
                double w = Weights[i, j];

                // Link existence cost
                H_local += _linkCostCoeff * (1.0 - w * w);

                // Length cost (strong links = short effective distance)
                double effectiveLength = w > 1e-10 ? 1.0 / w : 10.0;
                H_local += _lengthCostCoeff * effectiveLength;

                // String energy contribution
                if (_stringEnergy != null && i < N && j < N)
                    H_local += _stringEnergy[i, j];

                // 2. Gravity contribution: Ricci curvature for this edge
                double ricci = CalculateApproximateRicci(i, j);
                double G_inv = 1.0 / 0.1;
                H_local -= G_inv * ricci;
            }

            // 3. Matter coupling for affected nodes (i and j)
            if (_correlationMass != null && _correlationMass.Length == N)
            {
                // Node i contribution
                double mass_i = _correlationMass[i];
                double curvature_i = GetLocalCurvature(i);
                H_local += _matterCouplingCoeff * mass_i * (1.0 - Math.Abs(curvature_i));

                // Node j contribution
                double mass_j = _correlationMass[j];
                double curvature_j = GetLocalCurvature(j);
                H_local += _matterCouplingCoeff * mass_j * (1.0 - Math.Abs(curvature_j));
            }

            // 4. Cosmological constant term: Lambda * local_volume
            // Local volume is approximated as sum of weights around nodes i and j
            double localVolume = 0.0;
            foreach (int k in Neighbors(i))
                localVolume += Weights[i, k];
            foreach (int k in Neighbors(j))
                localVolume += Weights[j, k];
            // Avoid double counting edge (i,j)
            if (Edges[i, j])
                localVolume -= Weights[i, j];

            H_local += _cosmologicalConstant * localVolume;

            return H_local;
        }

        /// <summary>
        /// OBSOLETE: Uses global Hamiltonian which violates RQ-Hypothesis Locality.
        /// 
        /// Propose a random change to edge weight and accept/reject via Metropolis criterion.
        /// Returns true if the change was accepted.
        /// 
        /// USE INSTEAD: MetropolisEdgeStepLocalAction() in RQGraph.LocalAction.cs
        /// which computes LOCAL action change in O(degree²) time instead of O(N²).
        /// 
        /// RQ-Hypothesis Compliant (Item 2): Enforces strict light cone constraints.
        /// New edges can only be created between causally connected nodes
        /// within the light cone defined by c*dt distance.
        /// </summary>
        [Obsolete("Uses global Hamiltonian. Use MetropolisEdgeStepLocalAction() for RQ-compliant local action.", error: true)]
        public bool MetropolisEdgeStep()
        {
            // LEGACY CODE REMOVED - Use MetropolisEdgeStepLocalAction() instead
            throw new NotSupportedException(
                "MetropolisEdgeStep() is obsolete and violates RQ-Hypothesis locality. " +
                "Use MetropolisEdgeStepLocalAction() from RQGraph.LocalAction.cs instead.");
        }

        /// <summary>
        /// OBSOLETE: Uses ComputeLocalHamiltonian which is a partial fix but not fully RQ-compliant.
        /// 
        /// Optimized Metropolis step using local delta-H calculation instead of full Hamiltonian.
        /// This reduces complexity from O(N^2) to O(degree) per step, solving the freeze hazard
        /// for graphs with more than 100 nodes.
        /// 
        /// USE INSTEAD: MetropolisEdgeStepLocalAction() in RQGraph.LocalAction.cs
        /// which implements proper local action computation following RQ-Hypothesis.
        /// 
        /// RQ-Hypothesis Compliant (Item 2): Enforces strict light cone constraints.
        /// New edges can only be created between causally connected nodes
        /// within the light cone defined by c*dt distance.
        /// </summary>
        /// <returns>True if the proposed change was accepted</returns>
        [Obsolete("Use MetropolisEdgeStepLocalAction() for proper RQ-compliant local action.", error: true)]
        public bool MetropolisEdgeStepOptimized()
        {
            // LEGACY CODE REMOVED - Use MetropolisEdgeStepLocalAction() instead
            throw new NotSupportedException(
                "MetropolisEdgeStepOptimized() is obsolete. " +
                "Use MetropolisEdgeStepLocalAction() from RQGraph.LocalAction.cs instead.");
        }

        /// <summary>
        /// OBSOLETE: Uses MetropolisEdgeStep which violates RQ-Hypothesis Locality.
        /// 
        /// Optimize network topology using Metropolis-Hastings.
        /// Replaces deterministic ComputeNextState logic.
        /// 
        /// USE INSTEAD: A loop calling MetropolisEdgeStepLocalAction() directly.
        /// Example:
        ///   for (int s = 0; s < steps; s++)
        ///       graph.MetropolisEdgeStepLocalAction();
        /// </summary>
        [Obsolete("Uses global Hamiltonian. Use MetropolisEdgeStepLocalAction() in a loop instead.", error: true)]
        public void OptimizeNetworkTopology(int steps = 10)
        {
            // LEGACY CODE REMOVED - Use loop with MetropolisEdgeStepLocalAction() instead
            throw new NotSupportedException(
                "OptimizeNetworkTopology() is obsolete and violates RQ-Hypothesis locality. " +
                "Use a loop calling MetropolisEdgeStepLocalAction() instead:\n" +
                "  for (int s = 0; s < steps; s++)\n" +
                "      graph.MetropolisEdgeStepLocalAction();");
        }

        /// <summary>
        /// Simulated annealing: gradually reduce temperature to find ground state.
        /// Uses local action Metropolis step for RQ-Hypothesis compliance.
        /// </summary>
        public void SimulatedAnnealingStep(double coolingRate = 0.99)
        {
            // RQ-FIX: Use local action Metropolis step instead of obsolete global method
            int steps = Math.Max(1, N / 10);
            int accepted = 0;

            for (int s = 0; s < steps; s++)
            {
                if (MetropolisEdgeStepLocalAction())
                    accepted++;
            }

            // Update edge delays after topology changes
            if (accepted > 0)
            {
                UpdateEdgeDelaysFromDistances();
            }

            _networkTemperature *= coolingRate;

            // Don't let temperature go too low (numerical issues)
            if (_networkTemperature < 0.001)
                _networkTemperature = 0.001;
        }

        /// <summary>
        /// Quantum evolution via exp(-iHt) where H is the graph Laplacian.
        /// This replaces ad-hoc QuantumDiffusion with principled evolution.
        /// </summary>
        public void EvolveQuantumUnitary(double dt)
        {
            if (_waveMulti == null) return;

            int d = GaugeDimension;
            int len = N * d;

            // Build Hamiltonian matrix (graph Laplacian + potential)
            var H = ComputeGraphLaplacian();

            // Add potential diagonal terms from correlation mass
            if (_correlationMass != null && _correlationMass.Length == N)
            {
                for (int i = 0; i < N; i++)
                {
                    H[i, i] += _correlationMass[i]; // Mass as potential
                }
            }

            // For each gauge component, apply evolution operator
            // Using first-order approximation: psi(t+dt) ≈ (I - i*H*dt) * psi(t)
            // Then normalize for unitarity

            var psiNew = new Complex[len];

            for (int a = 0; a < d; a++)
            {
                for (int i = 0; i < N; i++)
                {
                    int idx = i * d + a;
                    Complex psi_i = _waveMulti[idx];

                    // Apply -iH operator
                    Complex Hpsi = Complex.Zero;
                    for (int j = 0; j < N; j++)
                    {
                        int jdx = j * d + a;
                        Hpsi += H[i, j] * _waveMulti[jdx];
                    }

                    // Evolution: psi_new = psi - i*dt*H*psi
                    psiNew[idx] = psi_i - Complex.ImaginaryOne * dt * Hpsi;
                }
            }

            // Normalize to preserve probability
            double norm = 0.0;
            for (int i = 0; i < len; i++)
            {
                double mag = psiNew[i].Magnitude;
                norm += mag * mag;
            }

            if (norm > 1e-10)
            {
                double invNorm = 1.0 / Math.Sqrt(norm);
                for (int i = 0; i < len; i++)
                {
                    psiNew[i] *= invNorm;
                }
            }

            _waveMulti = psiNew;
        }

        /// <summary>
        /// Combined physics step following Quantum Graphity principles:
        /// 1. Quantum state evolves unitarily with Hamiltonian
        /// 2. Topology optimizes to minimize action
        /// 3. Classical states emerge from quantum measurement
        /// 
        /// RQ-FIX: Now uses MetropolisEdgeStepLocalAction() for topology optimization
        /// instead of obsolete global Hamiltonian method.
        /// </summary>
        public void QuantumGraphityStep()
        {
            double dt = ComputeRelationalDt();
            EvolveQuantumUnitary(dt);

            // РЕКОМЕНДАЦИЯ: Обновить массы, так как поля только что изменились
            UpdateNodeMassModels();


            // 2. Optimize network topology using LOCAL action Metropolis-Hastings
            // RQ-FIX: Use local action step instead of obsolete OptimizeNetworkTopology
            int steps = Math.Max(1, N / 10);
            int accepted = 0;
            for (int s = 0; s < steps; s++)
            {
                if (MetropolisEdgeStepLocalAction())
                    accepted++;
            }

            // Update edge delays after topology changes
            if (accepted > 0)
            {
                UpdateEdgeDelaysFromDistances();
            }

            // 3. Update correlation mass from new topology
            RecomputeCorrelationMass();

            // 4. Update spectral geometry (emergent coordinates)
            if (_rng.NextDouble() < 0.1) // Expensive, do occasionally
            {
                UpdateSpectralCoordinates();
                SyncCoordinatesFromSpectral();
            }

            // 5. Update clock correlations (internal time)
            UpdateClockCorrelations();
            AdvanceInternalClock();

            // 6. Classical state emerges from quantum + clock
            UpdateStatesFromClockCondProb();

            // 7. Update correlation weights (Hebbian-like learning)
            UpdateCorrelationWeights();
        }

        /// <summary>
        /// Get current network energy (computed from local action sum)
        /// </summary>
        public double GetNetworkEnergy() => _lastNetworkEnergy;

        /// <summary>
        /// Get acceptance rate estimate for monitoring convergence.
        /// Uses local action Metropolis step for RQ-Hypothesis compliance.
        /// </summary>
        public double GetAcceptanceRate()
        {
            // Do a test batch using local action step
            int accepted = 0;
            int trials = 100;

            for (int t = 0; t < trials; t++)
            {
                if (MetropolisEdgeStepLocalAction())
                    accepted++;
            }

            return (double)accepted / trials;
        }


























        // Внутри RQGraph.QuantumGraphity.cs или нового файла RQGraph.Ricci.cs

        public double CalculateApproximateRicci(int i, int j)
        {
            // Упрощенная Форман-Риччи кривизна для взвешенных графов
            // Ric(e) ~ 4 - deg(i) - deg(j) + 3 * (triangles containing e)
            // Для взвешенных графов используем сумму весов вместо степени.

            if (!Edges[i, j]) return 0.0;

            double w_e = Weights[i, j];
            double w_i = 0.0; // Сумма весов соседей i
            double w_j = 0.0; // Сумма весов соседей j

            // Считаем "взвешенные степени", исключая само ребро (i,j)
            foreach (var n in Neighbors(i)) w_i += Weights[i, n];
            foreach (var n in Neighbors(j)) w_j += Weights[j, n];

            w_i -= w_e;
            w_j -= w_e;

            // Считаем треугольники (корреляции соседей)
            double triangles = 0.0;
            // Находим общих соседей
            foreach (var n_i in Neighbors(i))
            {
                if (n_i == j) continue;
                if (Edges[j, n_i]) // Треугольник i-j-n_i
                {
                    // Вклад треугольника зависит от силы связей
                    double w_in = Weights[i, n_i];
                    double w_jn = Weights[j, n_i];
                    triangles += Math.Sqrt(w_in * w_jn); // Геометрическое среднее
                }
            }

            // Эвристическая формула кривизны (базируется на Forman's Ricci curvature)
            // Положительная кривизна = много треугольников (кластер).
            // Отрицательная = древовидная структура.
            return w_e * (triangles - (w_i + w_j) * 0.1);
        }

    }
}