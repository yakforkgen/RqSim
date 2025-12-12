using System;
using ComputeSharp;

namespace RQSimulation.GPUOptimized
{
    /// <summary>
    /// GPU-accelerated random walk engine for spectral dimension computation.
    /// Launches thousands of independent random walkers in parallel.
    /// 
    /// Spectral dimension d_s is computed from return probability P(t):
    /// P(t) ~ t^(-d_s/2) for large t
    /// 
    /// Each GPU thread manages one walker, enabling 10,000+ walkers
    /// to be simulated simultaneously.
    /// 
    /// FIX: Added topology version tracking to prevent stale topology issues.
    /// FIX: Added UpdateTopologyFromGraph method for easier synchronization.
    /// </summary>
    public class SpectralWalkEngine : IDisposable
    {
        private readonly GraphicsDevice _device;

        // Walker state
        private ReadWriteBuffer<int>? _walkerPositions;
        private ReadOnlyBuffer<int>? _startPositions;
        private ReadWriteBuffer<int>? _returnCounts;

        // Topology in CSR format
        private ReadOnlyBuffer<int>? _adjOffsets;
        private ReadOnlyBuffer<int>? _adjNeighbors;
        private ReadOnlyBuffer<float>? _cumulativeWeights;

        private int _walkerCount;
        private int _nodeCount;
        private uint _currentSeed;
        private bool _initialized;

        // Topology version tracking to detect stale data
        private int _topologyVersion = -1;

        /// <summary>
        /// Current cached topology version. Compare with graph.TopologyVersion to detect staleness.
        /// </summary>
        public int TopologyVersion => _topologyVersion;

        /// <summary>
        /// Whether the engine is initialized and ready for computation.
        /// </summary>
        public bool IsInitialized => _initialized;

        /// <summary>
        /// Number of walker threads.
        /// </summary>
        public int WalkerCount => _walkerCount;

        /// <summary>
        /// Number of nodes in the cached topology.
        /// </summary>
        public int NodeCount => _nodeCount;

        public SpectralWalkEngine()
        {
            _device = GraphicsDevice.GetDefault();
            _currentSeed = (uint)Environment.TickCount;
        }

        /// <summary>
        /// Initialize the random walk engine.
        /// </summary>
        /// <param name="walkerCount">Number of parallel random walkers</param>
        /// <param name="nodeCount">Number of nodes in the graph</param>
        /// <param name="totalEdges">Total number of directed edges (2 * undirected)</param>
        public void Initialize(int walkerCount, int nodeCount, int totalEdges)
        {
            _walkerCount = walkerCount;
            _nodeCount = nodeCount;

            // Dispose old buffers
            _walkerPositions?.Dispose();
            _startPositions?.Dispose();
            _returnCounts?.Dispose();
            _adjOffsets?.Dispose();
            _adjNeighbors?.Dispose();
            _cumulativeWeights?.Dispose();

            _walkerPositions = _device.AllocateReadWriteBuffer<int>(walkerCount);
            _startPositions = _device.AllocateReadOnlyBuffer<int>(walkerCount);
            _returnCounts = _device.AllocateReadWriteBuffer<int>(1);

            _adjOffsets = _device.AllocateReadOnlyBuffer<int>(nodeCount + 1);
            _adjNeighbors = _device.AllocateReadOnlyBuffer<int>(totalEdges);
            _cumulativeWeights = _device.AllocateReadOnlyBuffer<float>(totalEdges);
        }

        /// <summary>
        /// Update topology buffers from RQGraph.
        /// Call this when graph.TopologyVersion changes or before first computation.
        /// 
        /// This method builds CSR topology with weighted edges, matching the
        /// CPU random walk implementation.
        /// </summary>
        /// <param name="graph">The RQGraph to read topology from</param>
        /// <param name="walkerCount">Number of parallel walkers (default 10000)</param>
        public void UpdateTopologyFromGraph(RQGraph graph, int walkerCount = 10000)
        {
            ArgumentNullException.ThrowIfNull(graph);

            // Build CSR format
            graph.BuildSoAViews();

            int nodeCount = graph.N;
            int[] offsets = graph.CsrOffsets;
            int[] neighbors = graph.CsrIndices;
            int totalEdges = neighbors.Length;

            // Build weights array matching CSR order
            float[] weights = new float[totalEdges];
            for (int i = 0; i < nodeCount; i++)
            {
                int start = offsets[i];
                int end = offsets[i + 1];
                for (int k = start; k < end; k++)
                {
                    int j = neighbors[k];
                    weights[k] = (float)graph.Weights[i, j];
                }
            }

            // Initialize buffers if needed
            if (_walkerCount != walkerCount || _nodeCount != nodeCount || _adjNeighbors == null ||
                _adjNeighbors.Length != totalEdges)
            {
                Initialize(walkerCount, nodeCount, totalEdges);
            }

            // Update topology
            UpdateTopology(offsets, neighbors, weights);

            // Store topology version
            _topologyVersion = graph.TopologyVersion;

            Console.WriteLine($"[SpectralWalkEngine] Updated topology: N={nodeCount}, E={totalEdges/2}, version={_topologyVersion}");
        }

        /// <summary>
        /// Update topology buffers.
        /// Also precomputes cumulative weights for efficient weighted sampling.
        /// </summary>
        public void UpdateTopology(int[] offsets, int[] neighbors, float[] weights)
        {
            if (_adjOffsets == null || _adjNeighbors == null || _cumulativeWeights == null)
            {
                throw new InvalidOperationException("Buffers not initialized.");
            }

            // Precompute cumulative weights for each node's neighbors
            float[] cumulative = new float[weights.Length];
            for (int i = 0; i < _nodeCount; i++)
            {
                int start = offsets[i];
                int end = offsets[i + 1];
                float sum = 0;

                for (int k = start; k < end; k++)
                {
                    sum += weights[k];
                    cumulative[k] = sum;
                }
            }

            _adjOffsets.CopyFrom(offsets);
            _adjNeighbors.CopyFrom(neighbors);
            _cumulativeWeights.CopyFrom(cumulative);
            _initialized = true;
        }

        /// <summary>
        /// Initialize walkers at random starting positions.
        /// </summary>
        /// <param name="random">Random number generator</param>
        public void InitializeWalkersRandom(Random random)
        {
            if (_walkerPositions == null || _startPositions == null)
            {
                throw new InvalidOperationException("Buffers not initialized.");
            }

            int[] positions = new int[_walkerCount];
            for (int i = 0; i < _walkerCount; i++)
            {
                positions[i] = random.Next(_nodeCount);
            }

            _walkerPositions.CopyFrom(positions);
            
            // Need to recreate start positions buffer with new data
            _startPositions.Dispose();
            _startPositions = _device.AllocateReadOnlyBuffer(positions);
        }

        /// <summary>
        /// Initialize all walkers at the same starting node.
        /// Useful for measuring return probability from a specific node.
        /// </summary>
        public void InitializeWalkersAt(int startNode)
        {
            if (_walkerPositions == null || _startPositions == null)
            {
                throw new InvalidOperationException("Buffers not initialized.");
            }

            int[] positions = new int[_walkerCount];
            Array.Fill(positions, startNode);

            _walkerPositions.CopyFrom(positions);
            
            _startPositions.Dispose();
            _startPositions = _device.AllocateReadOnlyBuffer(positions);
        }

        /// <summary>
        /// Initialize walkers uniformly distributed across all nodes.
        /// </summary>
        public void InitializeWalkersUniform()
        {
            if (_walkerPositions == null || _startPositions == null)
            {
                throw new InvalidOperationException("Buffers not initialized.");
            }

            int[] positions = new int[_walkerCount];
            for (int i = 0; i < _walkerCount; i++)
            {
                positions[i] = i % _nodeCount;
            }

            _walkerPositions.CopyFrom(positions);
            
            _startPositions.Dispose();
            _startPositions = _device.AllocateReadOnlyBuffer(positions);
        }

        /// <summary>
        /// Perform one step of random walk for all walkers.
        /// Returns the number of walkers that returned to their starting position.
        /// </summary>
        public int Step()
        {
            if (!_initialized || _walkerPositions == null || _startPositions == null ||
                _returnCounts == null || _adjOffsets == null || _adjNeighbors == null ||
                _cumulativeWeights == null)
            {
                throw new InvalidOperationException("Engine not properly initialized.");
            }

            // Reset return counter
            int[] zero = [0];
            _returnCounts.CopyFrom(zero);

            // Advance seed for this step
            _currentSeed = _currentSeed * 1103515245u + 12345u;

            var shader = new RandomWalkShader(
                _walkerPositions,
                _startPositions,
                _returnCounts,
                _adjOffsets,
                _adjNeighbors,
                _cumulativeWeights,
                _currentSeed,
                _walkerCount,
                _nodeCount);

            _device.For(_walkerCount, shader);

            // Read return count
            int[] result = new int[1];
            _returnCounts.CopyTo(result);

            return result[0];
        }

        /// <summary>
        /// Perform multiple steps and collect return statistics.
        /// Returns array of return counts at each step.
        /// </summary>
        public int[] RunSteps(int numSteps)
        {
            int[] returns = new int[numSteps];

            for (int t = 0; t < numSteps; t++)
            {
                returns[t] = Step();
            }

            return returns;
        }

        /// <summary>
        /// Compute spectral dimension from return probability data.
        /// Uses log-log linear fit of P(t) ~ t^(-d_s/2)
        /// </summary>
        /// <param name="returns">Return counts at each time step</param>
        /// <param name="skipInitial">Number of initial steps to skip (thermalization)</param>
        public double ComputeSpectralDimension(int[] returns, int skipInitial = 10)
        {
            // P(t) = returns[t] / walkerCount
            // log(P(t)) = -d_s/2 * log(t) + const
            // Slope of log-log plot gives -d_s/2

            double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0;
            int count = 0;

            // Skip late steps too (noisy)
            int skipLate = Math.Min(10, returns.Length / 5);

            for (int t = skipInitial; t < returns.Length - skipLate; t++)
            {
                if (returns[t] <= 0) continue;

                double logT = Math.Log(t);
                double logP = Math.Log((double)returns[t] / _walkerCount);

                sumX += logT;
                sumY += logP;
                sumXY += logT * logP;
                sumX2 += logT * logT;
                count++;
            }

            if (count < 2)
            {
                Console.WriteLine("[SpectralWalkEngine] Insufficient data points for regression");
                return double.NaN; // Signal invalid computation
            }

            // Linear regression: slope = (n*?xy - ?x*?y) / (n*?x? - (?x)?)
            double denominator = count * sumX2 - sumX * sumX;
            if (Math.Abs(denominator) < 1e-10)
            {
                Console.WriteLine("[SpectralWalkEngine] Degenerate regression (zero denominator)");
                return double.NaN;
            }

            double slope = (count * sumXY - sumX * sumY) / denominator;

            Console.WriteLine($"[SpectralWalkEngine] Fit: count={count}, slope={slope:F6}");

            // If slope is non-negative, walkers are not diffusing (trapped or disconnected)
            if (slope >= -1e-4)
            {
                Console.WriteLine("[SpectralWalkEngine] Non-negative slope detected. Graph may be disconnected.");
                return double.NaN;
            }

            // d_s = -2 * slope
            double spectralDim = -2.0 * slope;

            // Clamp to reasonable range [1, 10]
            return Math.Clamp(spectralDim, 1.0, 10.0);
        }

        /// <summary>
        /// Compute spectral dimension with automatic topology synchronization.
        /// This is the recommended entry point for computing d_s on GPU.
        /// 
        /// FIX: Automatically checks topology version and updates if stale.
        /// </summary>
        /// <param name="graph">The RQGraph to compute d_s for</param>
        /// <param name="numSteps">Number of random walk steps</param>
        /// <param name="walkerCount">Number of parallel walkers</param>
        /// <param name="skipInitial">Number of initial steps to skip</param>
        /// <returns>Spectral dimension estimate, or NaN if computation failed</returns>
        public double ComputeSpectralDimensionWithSyncCheck(
            RQGraph graph,
            int numSteps = 100,
            int walkerCount = 10000,
            int skipInitial = 10)
        {
            ArgumentNullException.ThrowIfNull(graph);

            // Check for stale topology
            if (!_initialized || _topologyVersion != graph.TopologyVersion)
            {
                Console.WriteLine($"[SpectralWalkEngine] Topology stale (cached={_topologyVersion}, graph={graph.TopologyVersion}). Updating...");
                UpdateTopologyFromGraph(graph, walkerCount);
            }

            // Initialize walkers at random positions
            InitializeWalkersRandom(new Random());

            // Run random walks
            int[] returns = RunSteps(numSteps);

            // Compute spectral dimension
            double ds = ComputeSpectralDimension(returns, skipInitial);

            Console.WriteLine($"[SpectralWalkEngine] d_S = {ds:F4}");

            return ds;
        }

        /// <summary>
        /// Get current positions of all walkers.
        /// </summary>
        public int[] GetWalkerPositions()
        {
            if (_walkerPositions == null)
            {
                throw new InvalidOperationException("Buffers not initialized.");
            }

            int[] positions = new int[_walkerCount];
            _walkerPositions.CopyTo(positions);

            return positions;
        }

        public void Dispose()
        {
            _walkerPositions?.Dispose();
            _startPositions?.Dispose();
            _returnCounts?.Dispose();
            _adjOffsets?.Dispose();
            _adjNeighbors?.Dispose();
            _cumulativeWeights?.Dispose();
        }
    }

    /// <summary>
    /// GPU shader for weighted random walk.
    /// Each thread manages one walker.
    /// 
    /// Uses PCG hash for pseudo-random number generation.
    /// Weighted neighbor selection via cumulative weights.
    /// </summary>
    [ThreadGroupSize(64, 1, 1)]
    [GeneratedComputeShaderDescriptor]
    public readonly partial struct RandomWalkShader : IComputeShader
    {
        public readonly ReadWriteBuffer<int> walkerPositions;
        public readonly ReadOnlyBuffer<int> startPositions;
        public readonly ReadWriteBuffer<int> returnCounts;
        public readonly ReadOnlyBuffer<int> offsets;
        public readonly ReadOnlyBuffer<int> neighbors;
        public readonly ReadOnlyBuffer<float> cumulativeWeights;
        public readonly uint seed;
        public readonly int walkerCount;
        public readonly int nodeCount;

        public RandomWalkShader(
            ReadWriteBuffer<int> walkerPositions,
            ReadOnlyBuffer<int> startPositions,
            ReadWriteBuffer<int> returnCounts,
            ReadOnlyBuffer<int> offsets,
            ReadOnlyBuffer<int> neighbors,
            ReadOnlyBuffer<float> cumulativeWeights,
            uint seed,
            int walkerCount,
            int nodeCount)
        {
            this.walkerPositions = walkerPositions;
            this.startPositions = startPositions;
            this.returnCounts = returnCounts;
            this.offsets = offsets;
            this.neighbors = neighbors;
            this.cumulativeWeights = cumulativeWeights;
            this.seed = seed;
            this.walkerCount = walkerCount;
            this.nodeCount = nodeCount;
        }

        public void Execute()
        {
            int walkerIdx = ThreadIds.X;
            if (walkerIdx >= walkerCount) return;

            int currentNode = walkerPositions[walkerIdx];

            // 1. Generate pseudo-random number using PCG hash
            // Use both walkerIdx and seed for better entropy
            uint state = (uint)(walkerIdx * 73856093) ^ seed ^ (uint)(currentNode * 19349663);
            state = state * 747796405u + 2891336453u;
            uint word = ((state >> ((int)(state >> 28) + 4)) ^ state) * 277803737u;
            uint result = (word >> 22) ^ word;
            float rnd = (float)result / 4294967295.0f; // Normalize to [0, 1]

            // 2. Get neighbor range for current node
            int start = offsets[currentNode];
            int end = offsets[currentNode + 1];

            // BUG FIX: Track whether walker actually moved
            // Isolated nodes (degree=0) should NOT count as returns
            bool didMove = false;
            int nextNode = currentNode;

            if (end > start)
            {
                // 3. Select neighbor using cumulative weights
                float totalWeight = cumulativeWeights[end - 1];

                // BUG FIX: Only move if there's actual weight
                if (totalWeight > 0.0f)
                {
                    float targetWeight = rnd * totalWeight;

                    // Linear search for neighbor selection
                    // For typical graphs (degree ~8), this is fine
                    for (int k = start; k < end; k++)
                    {
                        if (cumulativeWeights[k] >= targetWeight)
                        {
                            nextNode = neighbors[k];
                            didMove = true;
                            break;
                        }
                    }

                    // Fallback to last neighbor if loop didn't select (floating point edge case)
                    if (!didMove && end > start)
                    {
                        nextNode = neighbors[end - 1];
                        didMove = true;
                    }
                }
            }

            // 4. Update position
            walkerPositions[walkerIdx] = nextNode;

            // 5. Check if returned to start and atomically increment counter
            // BUG FIX: Only count as return if walker actually moved this step
            // Otherwise isolated nodes would inflate return probability to 100%
            if (didMove && nextNode == startPositions[walkerIdx])
            {
                Hlsl.InterlockedAdd(ref returnCounts[0], 1);
            }
        }
    }
}
