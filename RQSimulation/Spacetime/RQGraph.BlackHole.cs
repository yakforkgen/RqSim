using System;
using System.Collections.Generic;
using System.Numerics;
using System.Runtime.CompilerServices;

namespace RQSimulation
{
    /// <summary>
    /// Implements black hole physics based on correlation density singularities.
    /// When correlation density exceeds a critical threshold, an event horizon forms
    /// and Hawking-like radiation emerges from the boundary.
    /// </summary>
    public partial class RQGraph
    {
        /// <summary>
        /// Represents a detected black hole region in the graph.
        /// </summary>
        public class BlackHoleRegion
        {
            public int CenterNode { get; init; }
            public List<int> InteriorNodes { get; init; } = new();
            public List<int> HorizonNodes { get; init; } = new();
            public double Mass { get; set; }
            public double SchwarzschildRadius { get; set; }
            public double Temperature { get; set; }
            public double Entropy { get; set; }
            public double AngularMomentum { get; set; }
            public double Charge { get; set; }
        }

        // Critical correlation density for horizon formation
        public double HorizonDensityThreshold { get; set; } = 2.0;

        // List of detected black hole regions
        private List<BlackHoleRegion>? _blackHoles;

        // Hawking radiation accumulator
        private double[]? _hawkingRadiation;

        // Information paradox tracking: entanglement with interior
        private double[,]? _horizonEntanglement;

        /// <summary>
        /// Initialize black hole detection and Hawking radiation tracking.
        /// </summary>
        public void InitBlackHolePhysics()
        {
            _blackHoles = new List<BlackHoleRegion>();
            _hawkingRadiation = new double[N];
            _horizonEntanglement = new double[N, N];
        }

        /// <summary>
        /// Detects black hole regions based on correlation density exceeding threshold.
        /// </summary>
        public List<BlackHoleRegion> DetectBlackHoles()
        {
            if (_blackHoles == null) InitBlackHolePhysics();
            _blackHoles!.Clear();

            // Find nodes with supercritical correlation density
            var supercritical = new HashSet<int>();
            for (int i = 0; i < N; i++)
            {
                double density = GetLocalCorrelationDensity(i);
                if (density >= HorizonDensityThreshold)
                {
                    supercritical.Add(i);
                }
            }

            if (supercritical.Count == 0) return _blackHoles;

            // Cluster supercritical nodes into black hole regions
            var visited = new HashSet<int>();
            foreach (int seed in supercritical)
            {
                if (visited.Contains(seed)) continue;

                var bh = new BlackHoleRegion { CenterNode = seed };
                var queue = new Queue<int>();
                queue.Enqueue(seed);
                visited.Add(seed);

                while (queue.Count > 0)
                {
                    int current = queue.Dequeue();
                    bh.InteriorNodes.Add(current);

                    foreach (int nb in Neighbors(current))
                    {
                        if (visited.Contains(nb)) continue;

                        if (supercritical.Contains(nb))
                        {
                            visited.Add(nb);
                            queue.Enqueue(nb);
                        }
                        else
                        {
                            // Boundary node = horizon
                            bh.HorizonNodes.Add(nb);
                        }
                    }
                }

                // Compute black hole properties
                ComputeBlackHoleProperties(bh);
                _blackHoles.Add(bh);
            }

            return _blackHoles;
        }

        /// <summary>
        /// Computes physical properties of a black hole region.
        /// Uses spectral (graph-based) coordinates instead of external coordinates
        /// to comply with RQ-hypothesis principles.
        /// </summary>
        private void ComputeBlackHoleProperties(BlackHoleRegion bh)
        {
            // Mass from correlation mass of interior
            double totalMass = 0;
            foreach (int node in bh.InteriorNodes)
            {
                double m = _correlationMass != null ? _correlationMass[node] : 1.0;
                totalMass += m;
            }
            bh.Mass = totalMass;

            // Center of mass using spectral coordinates (RQ-compliant)
            var (cx, cy, cz) = ComputeSpectralCenterOfMassWeighted(bh.InteriorNodes);

            // Schwarzschild radius: r_s = 2GM/c²
            double G = 1.0;
            double c2 = VectorMath.SpeedOfLight * VectorMath.SpeedOfLight;
            bh.SchwarzschildRadius = 2.0 * G * totalMass / c2;

            // Hawking temperature: T_H = ℏc³/(8πGMk_B)
            // In natural units with k_B = 1:
            double c3 = VectorMath.SpeedOfLight * c2;
            bh.Temperature = totalMass > 0.1
                ? VectorMath.HBar * c3 / (8.0 * Math.PI * G * totalMass)
                : 100.0;  // Small black holes are very hot

            // Bekenstein-Hawking entropy: S = A/(4l_P²) where A = 4πr_s²
            // In natural units: S = π r_s² / (G ℏ)
            double rs = bh.SchwarzschildRadius;
            bh.Entropy = Math.PI * rs * rs / (G * VectorMath.HBar);

            // Angular momentum from cluster rotation (using spectral coordinates)
            bh.AngularMomentum = ComputeClusterAngularMomentumSpectral(bh.InteriorNodes, cx, cy);

            // Electric charge from node charges
            bh.Charge = 0;
            if (_charges != null)
            {
                foreach (int node in bh.InteriorNodes)
                {
                    if (node < _charges.Length)
                        bh.Charge += _charges[node];
                }
            }
        }

        /// <summary>
        /// Computes angular momentum of a cluster using spectral coordinates.
        /// RQ-compliant: uses only graph-derived coordinates.
        /// </summary>
        private double ComputeClusterAngularMomentumSpectral(List<int> nodes, double cx, double cy)
        {
            if (_correlationMass == null) return 0;
            
            // Ensure spectral coordinates are computed
            if (_spectralX == null || _spectralX.Length != N)
            {
                UpdateSpectralCoordinates();
            }

            double L = 0;
            foreach (int node in nodes)
            {
                if (node < 0 || node >= N) continue;
                
                double rx = _spectralX![node] - cx;
                double ry = (_spectralY != null && node < _spectralY.Length ? _spectralY[node] : 0) - cy;

                // Approximate velocity from proper time gradient (relational)
                double vx = 0.1 * ry;  // Rough rotational velocity
                double vy = -0.1 * rx;

                double m = node < _correlationMass.Length ? _correlationMass[node] : 1.0;
                L += m * (rx * vy - ry * vx);
            }
            return L;
        }

        /// <summary>
        /// Emits Hawking radiation from black hole horizons.
        /// </summary>
        public void EmitHawkingRadiation()
        {
            if (_blackHoles == null || _hawkingRadiation == null) return;

            foreach (var bh in _blackHoles)
            {
                double T = bh.Temperature;
                if (T < 1e-10) continue;

                // Stefan-Boltzmann emission rate: P ∝ σT⁴A
                double emissionRate = 0.001 * Math.Pow(T, 4) * bh.HorizonNodes.Count;

                foreach (int horizonNode in bh.HorizonNodes)
                {
                    // Probability to emit a quantum
                    double emissionProb = emissionRate / bh.HorizonNodes.Count;

                    if (_rng.NextDouble() < emissionProb)
                    {
                        // Emit energy as excitation
                        State[horizonNode] = NodeState.Excited;
                        _hawkingRadiation[horizonNode] += T;

                        // Reduce black hole mass (very small amount)
                        if (bh.Mass > 0.01)
                        {
                            bh.Mass -= 0.001;
                            // Distribute mass loss to interior nodes
                            foreach (int interior in bh.InteriorNodes)
                            {
                                if (_correlationMass != null && interior < _correlationMass.Length)
                                {
                                    _correlationMass[interior] *= 0.9999;
                                }
                            }
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Checks if a node is inside a black hole (beyond the horizon).
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public bool IsInsideBlackHole(int node)
        {
            if (_blackHoles == null) return false;

            foreach (var bh in _blackHoles)
            {
                if (bh.InteriorNodes.Contains(node)) return true;
            }
            return false;
        }

        /// <summary>
        /// Checks if a node is on a black hole horizon.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public bool IsOnHorizon(int node)
        {
            if (_blackHoles == null) return false;

            foreach (var bh in _blackHoles)
            {
                if (bh.HorizonNodes.Contains(node)) return true;
            }
            return false;
        }

        /// <summary>
        /// Gets the total Hawking radiation received by a node.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public double GetHawkingRadiation(int node)
        {
            return _hawkingRadiation != null ? _hawkingRadiation[node] : 0;
        }

        /// <summary>
        /// Computes the Penrose process energy extraction potential at a node.
        /// </summary>
        public double PenroseEnergyPotential(int node)
        {
            if (_blackHoles == null) return 0;

            foreach (var bh in _blackHoles)
            {
                if (!bh.HorizonNodes.Contains(node)) continue;

                // In Kerr black hole, energy can be extracted from ergosphere
                // E_extract ≤ (1 - √((1+√(1-a²))/2)) * M
                // where a = J/(Mc)
                double a = Math.Abs(bh.AngularMomentum) / (bh.Mass * VectorMath.SpeedOfLight + 0.01);
                a = Math.Min(a, 0.99);  // Limit to sub-extremal

                double sqrtTerm = Math.Sqrt((1 + Math.Sqrt(1 - a * a)) / 2);
                double efficiency = 1 - sqrtTerm;

                return efficiency * bh.Mass;
            }
            return 0;
        }

        /// <summary>
        /// Tracks information scrambling across the black hole horizon.
        /// Implements "fast scrambling" where information spreads in time O(log S).
        /// </summary>
        public void ScrambleHorizonInformation()
        {
            if (_blackHoles == null || _horizonEntanglement == null) return;

            foreach (var bh in _blackHoles)
            {
                // Scrambling time: t_s ∝ log(S)/T
                double scramblingRate = bh.Temperature / (Math.Log(bh.Entropy + 1) + 1);

                // Mix entanglement between horizon nodes
                for (int i = 0; i < bh.HorizonNodes.Count; i++)
                {
                    int node1 = bh.HorizonNodes[i];
                    for (int j = i + 1; j < bh.HorizonNodes.Count; j++)
                    {
                        int node2 = bh.HorizonNodes[j];

                        // Increase entanglement
                        _horizonEntanglement[node1, node2] +=
                            scramblingRate * (1 - _horizonEntanglement[node1, node2]);
                        _horizonEntanglement[node2, node1] = _horizonEntanglement[node1, node2];
                    }
                }
            }
        }

        /// <summary>
        /// Checks for black hole evaporation (complete mass loss).
        /// </summary>
        public List<BlackHoleRegion> CheckEvaporation()
        {
            var evaporated = new List<BlackHoleRegion>();

            if (_blackHoles == null) return evaporated;

            // Iterate in reverse to safely remove items
            for (int idx = _blackHoles.Count - 1; idx >= 0; idx--)
            {
                var bh = _blackHoles[idx];
                if (bh.Mass < 0.1)
                {
                    // Black hole has evaporated!
                    evaporated.Add(bh);

                    // Release remaining energy as burst
                    foreach (int node in bh.InteriorNodes)
                    {
                        State[node] = NodeState.Excited;
                        if (_nodeEnergy != null && node < _nodeEnergy.Length)
                        {
                            _nodeEnergy[node] += bh.Mass / bh.InteriorNodes.Count;
                        }
                    }

                    _blackHoles.RemoveAt(idx);
                }
            }

            return evaporated;
        }

        /// <summary>
        /// Computes the surface gravity at a black hole horizon.
        /// κ = c⁴/(4GM) for Schwarzschild
        /// </summary>
        public double SurfaceGravity(BlackHoleRegion bh)
        {
            double G = 1.0;
            double c4 = Math.Pow(VectorMath.SpeedOfLight, 4);
            return c4 / (4.0 * G * Math.Max(0.01, bh.Mass));
        }

        /// <summary>
        /// Implements black hole mergers when two horizons overlap.
        /// </summary>
        public void CheckBlackHoleMergers()
        {
            if (_blackHoles == null || _blackHoles.Count < 2) return;

            for (int i = 0; i < _blackHoles.Count; i++)
            {
                for (int j = i + 1; j < _blackHoles.Count; j++)
                {
                    var bh1 = _blackHoles[i];
                    var bh2 = _blackHoles[j];

                    // Check for horizon overlap
                    bool overlap = false;
                    foreach (int node in bh1.HorizonNodes)
                    {
                        if (bh2.InteriorNodes.Contains(node) || bh2.HorizonNodes.Contains(node))
                        {
                            overlap = true;
                            break;
                        }
                    }

                    if (overlap)
                    {
                        // Merge: combine into bh1
                        bh1.Mass += bh2.Mass;
                        bh1.AngularMomentum += bh2.AngularMomentum;
                        bh1.Charge += bh2.Charge;

                        foreach (int node in bh2.InteriorNodes)
                        {
                            if (!bh1.InteriorNodes.Contains(node))
                                bh1.InteriorNodes.Add(node);
                        }
                        foreach (int node in bh2.HorizonNodes)
                        {
                            if (!bh1.HorizonNodes.Contains(node))
                                bh1.HorizonNodes.Add(node);
                        }

                        // Recompute properties
                        ComputeBlackHoleProperties(bh1);

                        // Remove merged black hole
                        _blackHoles.RemoveAt(j);
                        j--;

                        // Emit gravitational wave energy (10% of mass)
                        double gwEnergy = 0.1 * bh2.Mass;
                        EmitGravitationalWave(bh1.CenterNode, gwEnergy);
                    }
                }
            }
        }

        /// <summary>
        /// Emits gravitational wave energy from a point.
        /// </summary>
        private void EmitGravitationalWave(int source, double energy)
        {
            if (_nodeEnergy == null) return;

            // Distribute energy radially
            double remaining = energy;
            var visited = new HashSet<int> { source };
            var frontier = new Queue<int>();
            frontier.Enqueue(source);

            int ring = 0;
            while (frontier.Count > 0 && remaining > 0.01)
            {
                int ringSize = frontier.Count;
                double energyPerNode = remaining / (ringSize * 2);

                for (int k = 0; k < ringSize; k++)
                {
                    int node = frontier.Dequeue();
                    if (node < _nodeEnergy.Length)
                    {
                        _nodeEnergy[node] += energyPerNode;
                        remaining -= energyPerNode;
                    }

                    foreach (int nb in Neighbors(node))
                    {
                        if (!visited.Contains(nb))
                        {
                            visited.Add(nb);
                            frontier.Enqueue(nb);
                        }
                    }
                }
                ring++;
            }
        }

        /// <summary>
        /// Gets all detected black holes.
        /// </summary>
        public IReadOnlyList<BlackHoleRegion> GetBlackHoles()
        {
            return _blackHoles ?? new List<BlackHoleRegion>();
        }
        
        // =====================================================================
        // Black Hole Saturation Mechanism (Checklist item 5.2)
        // =====================================================================
        
        /// <summary>
        /// Maximum entropy constant (bits per horizon unit area).
        /// In Bekenstein-Hawking formula: S = A/(4*l_p²)
        /// For our discrete model, we use a scaled version.
        /// </summary>
        private const double BekensteinConstant = 1.0;
        
        /// <summary>
        /// Compute horizon area from node's connection density.
        /// Area ∝ number of links (discrete analog of r²).
        /// </summary>
        /// <param name="nodeId">Node to analyze</param>
        /// <returns>Effective horizon area</returns>
        public double GetHorizonArea(int nodeId)
        {
            if (nodeId < 0 || nodeId >= N)
                return 0.0;
            
            // Area is proportional to link count and their total weight
            double area = 0.0;
            foreach (int j in Neighbors(nodeId))
            {
                area += Weights[nodeId, j];
            }
            
            return area;
        }
        
        /// <summary>
        /// Calculate maximum entropy (maximum number of links) for given horizon area.
        /// Implements Bekenstein bound: S_max = A/(4*l_p²)
        /// </summary>
        /// <param name="horizonArea">Horizon area of the node</param>
        /// <returns>Maximum allowed link count</returns>
        public int MaxEntropy(double horizonArea)
        {
            // S_max = (BekensteinConstant * A)
            // Each link carries roughly 1 bit of entropy
            int maxLinks = (int)(BekensteinConstant * horizonArea * N / 4.0);
            return Math.Max(1, maxLinks);
        }
        
        /// <summary>
        /// Get current link count for a node
        /// </summary>
        /// <param name="nodeId">Node index</param>
        /// <returns>List of neighboring nodes (current links)</returns>
        public List<int> CurrentLinks(int nodeId)
        {
            return Neighbors(nodeId).ToList();
        }
        
        /// <summary>
        /// Trigger Hawking radiation from an oversaturated node.
        /// Called when a node exceeds its entropy bound.
        /// </summary>
        /// <param name="nodeId">Node that is radiating</param>
        private void TriggerHawkingRadiation(int nodeId)
        {
            if (_hawkingRadiation == null)
                _hawkingRadiation = new double[N];
            
            // Find the weakest link to radiate
            int weakestNeighbor = -1;
            double minWeight = double.MaxValue;
            
            foreach (int j in Neighbors(nodeId))
            {
                if (Weights[nodeId, j] < minWeight)
                {
                    minWeight = Weights[nodeId, j];
                    weakestNeighbor = j;
                }
            }
            
            if (weakestNeighbor >= 0)
            {
                // Emit energy as Hawking radiation
                double radiatedEnergy = minWeight * 0.1;
                _hawkingRadiation[weakestNeighbor] += radiatedEnergy;
                
                // Weaken the oversaturated link
                Weights[nodeId, weakestNeighbor] *= 0.9;
                Weights[weakestNeighbor, nodeId] *= 0.9;
                
                // Very weak links may be removed entirely
                if (Weights[nodeId, weakestNeighbor] < 0.01)
                {
                    RemoveEdge(nodeId, weakestNeighbor);
                }
                
                // Also excite the radiating neighbor
                State[weakestNeighbor] = NodeState.Excited;
            }
        }
        
        /// <summary>
        /// Check if a new link can be added to a node without violating entropy bounds.
        /// Implements RQ-hypothesis checklist item 5.2 (black hole saturation).
        /// If the node is saturated, triggers Hawking radiation instead.
        /// </summary>
        /// <param name="nodeId">Node to check</param>
        /// <returns>True if link can be added, false if saturated (will radiate instead)</returns>
        public bool CanAddLink(int nodeId)
        {
            if (nodeId < 0 || nodeId >= N)
                return false;
            
            double horizonArea = GetHorizonArea(nodeId);
            int currentLinkCount = CurrentLinks(nodeId).Count;
            int maxLinks = MaxEntropy(horizonArea);
            
            if (currentLinkCount > maxLinks)
            {
                // Entropy bound exceeded: trigger Hawking radiation (evaporation instead of absorption)
                TriggerHawkingRadiation(nodeId);
                return false;
            }
            
            return true;
        }
        
        /// <summary>
        /// Safe edge addition that respects entropy bounds.
        /// Will not add edge if it would violate Bekenstein bound.
        /// </summary>
        /// <param name="nodeA">First node</param>
        /// <param name="nodeB">Second node</param>
        /// <returns>True if edge was added, false if blocked by entropy limit</returns>
        public bool TryAddEdgeSafe(int nodeA, int nodeB)
        {
            if (!CanAddLink(nodeA) || !CanAddLink(nodeB))
                return false;
            
            AddEdge(nodeA, nodeB);
            return true;
        }
    }
}
