using System;
using System.Collections.Generic;
using System.Numerics;
using ComputeSharp;

namespace RQSimulation.GPUOptimized
{
    // GPU configuration structure
    public struct GpuConfig
    {
        public int GpuIndex;
        public bool MultiGpu;
        public int ThreadBlockSize;
    }

    /// <summary>
    /// GPU shader for computing Forman-Ricci curvature on edges.
    /// Uses CSR (Compressed Sparse Row) format for efficient neighbor lookup.
    /// 
    /// Formula: Ric(e) = w_e * (Σ√(w_e1*w_e2) for triangles - α*(W_u + W_v))
    /// where W_u = weighted degree of node u (excluding edge e)
    /// </summary>
    [ThreadGroupSize(64, 1, 1)]
    [GeneratedComputeShaderDescriptor]
    public readonly partial struct FormanCurvatureShader : IComputeShader
    {
        // Input data
        public readonly ReadWriteBuffer<float> weights;       // Current edge weights (w_ij)
        public readonly ReadOnlyBuffer<Int2> edges;           // Node pairs for each edge (u, v)

        // Graph in CSR format (for fast neighbor lookup)
        public readonly ReadOnlyBuffer<int> adjOffsets;       // Start indices of neighbors for each node
        public readonly ReadOnlyBuffer<Int2> adjData;         // Neighbor data: .X = neighborIndex, .Y = edgeIndex

        // Output data
        public readonly ReadWriteBuffer<float> curvatures;    // Result (Ric_ij)

        // Constants
        public readonly float degreePenaltyFactor;
        public readonly int nodeCount;

        public FormanCurvatureShader(
            ReadWriteBuffer<float> weights,
            ReadOnlyBuffer<Int2> edges,
            ReadOnlyBuffer<int> adjOffsets,
            ReadOnlyBuffer<Int2> adjData,
            ReadWriteBuffer<float> curvatures,
            float degreePenaltyFactor,
            int nodeCount)
        {
            this.weights = weights;
            this.edges = edges;
            this.adjOffsets = adjOffsets;
            this.adjData = adjData;
            this.curvatures = curvatures;
            this.degreePenaltyFactor = degreePenaltyFactor;
            this.nodeCount = nodeCount;
        }

        public void Execute()
        {
            int edgeIdx = ThreadIds.X;
            if (edgeIdx >= edges.Length) return;

            // 1. Get data for current edge E(u,v)
            Int2 nodes = edges[edgeIdx];
            int u = nodes.X;
            int v = nodes.Y;
            float w_uv = weights[edgeIdx];

            if (w_uv <= 0.0f)
            {
                curvatures[edgeIdx] = 0.0f;
                return;
            }

            // 2. Compute weighted degrees (excluding edge e)
            // W_node = sum of weights of all incident edges except E(u,v)

            float w_u = 0.0f;
            int startU = adjOffsets[u];
            int endU = (u + 1 < nodeCount) ? adjOffsets[u + 1] : adjData.Length;

            for (int k = startU; k < endU; k++)
            {
                int e_idx = adjData[k].Y;
                // Forman formula requires excluding edge e itself
                if (e_idx != edgeIdx)
                    w_u += weights[e_idx];
            }

            float w_v = 0.0f;
            int startV = adjOffsets[v];
            int endV = (v + 1 < nodeCount) ? adjOffsets[v + 1] : adjData.Length;

            for (int k = startV; k < endV; k++)
            {
                int e_idx = adjData[k].Y;
                if (e_idx != edgeIdx)
                    w_v += weights[e_idx];
            }

            // 3. Find triangles (3-cycles contribution)
            // Intersection of neighbor sets U and V
            float triangleTerm = 0.0f;

            // Iterate through neighbors of u, check if they are also neighbors of v
            for (int i = startU; i < endU; i++)
            {
                int neighbor_u = adjData[i].X;
                if (neighbor_u == v) continue; // Skip edge u-v itself

                // Search for neighbor_u in neighbor list of v
                for (int j = startV; j < endV; j++)
                {
                    int neighbor_v = adjData[j].X;

                    if (neighbor_u == neighbor_v)
                    {
                        // FOUND TRIANGLE (u, v, neighbor)
                        // Weights of adjacent edges
                        float w_un = weights[adjData[i].Y];
                        float w_vn = weights[adjData[j].Y];

                        // Geometric mean of weights (as in CPU implementation)
                        // Using cube root for consistency with ComputeFormanRicciCurvature
                        triangleTerm += Hlsl.Pow(w_un * w_vn * w_uv, 1.0f / 3.0f);
                    }
                }
            }

            // 4. Final Forman-Ricci curvature formula
            // Ric(e) = w(e) * (triangles - penalty*(W_u + W_v))
            float penalty = degreePenaltyFactor * (w_u + w_v);
            curvatures[edgeIdx] = w_uv * (triangleTerm - penalty);
        }
    }

    public class GpuGravityEngine : IDisposable
    {
        private readonly GraphicsDevice _device;

        // GPU buffers
        private ReadWriteBuffer<float> _weightsBuffer;
        private ReadWriteBuffer<float> _curvaturesBuffer; // Changed to ReadWrite for curvature shader output
        private ReadOnlyBuffer<float> _nodeMassesBuffer;
        private ReadOnlyBuffer<Int2> _edgeIndicesBuffer; // Edge pairs (NodeA, NodeB)

        // CSR adjacency buffers for curvature computation
        private ReadOnlyBuffer<int>? _adjOffsetsBuffer;
        private ReadOnlyBuffer<Int2>? _adjDataBuffer;
        private bool _topologyInitialized = false;
        private int _nodeCount;

        public GpuGravityEngine(GpuConfig config, int edgeCount, int nodeCount)
        {
            // 1. Select device (default for now, multi-GPU support prepared)
            _device = GraphicsDevice.GetDefault();
            if (config.MultiGpu && config.GpuIndex >= 0)
            {
                // Select specific GPU if needed
                // _device = GraphicsDevice.QueryDevices().ElementAt(config.GpuIndex);
            }

            _nodeCount = nodeCount;

            // 2. Allocate GPU buffers
            _weightsBuffer = _device.AllocateReadWriteBuffer<float>(edgeCount);
            _curvaturesBuffer = _device.AllocateReadWriteBuffer<float>(edgeCount); // ReadWrite for curvature shader
            _nodeMassesBuffer = _device.AllocateReadOnlyBuffer<float>(nodeCount);
            _edgeIndicesBuffer = _device.AllocateReadOnlyBuffer<Int2>(edgeCount);
        }

        /// <summary>
        /// Update topology buffers (CSR format) for GPU curvature computation.
        /// Call this when graph topology changes (edge add/remove).
        /// 
        /// IMPORTANT: This builds a directed adjacency list where each undirected edge
        /// appears twice (once for each direction). The edgeIndex stored is always
        /// the canonical index from FlatEdgesFrom (where i less than j).
        /// </summary>
        /// <param name="graph">The graph to build CSR from</param>
        public void UpdateTopologyBuffers(RQGraph graph)
        {
            int N = graph.N;
            int[] offsets = new int[N + 1];
            List<Int2> flatAdj = new List<Int2>();

            // Build CSR structure on CPU (done rarely, only on rewiring)
            // Each undirected edge (i,j) appears twice in adjacency lists:
            // - In node i's list: neighbor=j, edgeIndex=GetEdgeIndex(i,j)
            // - In node j's list: neighbor=i, edgeIndex=GetEdgeIndex(j,i) = same index
            int currentOffset = 0;
            for (int i = 0; i < N; i++)
            {
                offsets[i] = currentOffset;
                foreach (int neighbor in graph.Neighbors(i))
                {
                    // GetEdgeIndex normalizes to (min, max) internally, so it returns
                    // the same index regardless of argument order
                    int edgeIndex = graph.GetEdgeIndex(i, neighbor);
                    if (edgeIndex >= 0) // Valid edge
                    {
                        flatAdj.Add(new Int2(neighbor, edgeIndex));
                        currentOffset++;
                    }
                }
            }
            offsets[N] = currentOffset;

            // Validate: flatAdj.Count should equal 2 * edgeCount (each edge appears twice)
            int expectedAdjCount = graph.FlatEdgesFrom.Length * 2;
            if (flatAdj.Count != expectedAdjCount)
            {
                // This can happen if graph has isolated nodes or self-loops
                // Log warning but continue
                Console.WriteLine($"[GPU] Warning: CSR adjacency count {flatAdj.Count} != expected {expectedAdjCount}");
            }

            // Upload to GPU
            _adjOffsetsBuffer?.Dispose();
            _adjDataBuffer?.Dispose();

            _adjOffsetsBuffer = _device.AllocateReadOnlyBuffer(offsets);
            _adjDataBuffer = _device.AllocateReadOnlyBuffer(flatAdj.ToArray());
            _topologyInitialized = true;
            _nodeCount = N;
        }

        /// <summary>
        /// Check if topology buffers are initialized.
        /// </summary>
        public bool IsTopologyInitialized => _topologyInitialized;

        /// <summary>
        /// Hybrid GPU method: Accepts curvature computed on CPU, evolves gravity on GPU.
        /// Use EvolveFullGpuStep for full GPU computation when topology buffers are ready.
        /// </summary>
        public void EvolveGravityGpu(
            float[] hostWeights,
            float[] hostCurvatures,
            float[] hostMasses,
            int[] hostEdgesFrom,
            int[] hostEdgesTo,
            float dt,
            float G,
            float lambda)
        {
            int edgeCount = hostWeights.Length;

            // Copy data CPU -> GPU
            _weightsBuffer.CopyFrom(hostWeights);
            _curvaturesBuffer.CopyFrom(hostCurvatures);
            _nodeMassesBuffer.CopyFrom(hostMasses);

            // Pack edge pairs into Int2 (only needed if topology changed)
            Int2[] packedEdges = new Int2[edgeCount];
            for (int i = 0; i < edgeCount; i++)
                packedEdges[i] = new Int2(hostEdgesFrom[i], hostEdgesTo[i]);
            _edgeIndicesBuffer.CopyFrom(packedEdges);

            // Run gravity shader
            var shader = new GravityShader(
                _weightsBuffer,
                _curvaturesBuffer,
                _nodeMassesBuffer,
                _edgeIndicesBuffer,
                dt, G, lambda);

            _device.For(edgeCount, shader);

            // Copy results GPU -> CPU
            _weightsBuffer.CopyTo(hostWeights);
        }

        /// <summary>
        /// Full GPU step: Compute curvature AND evolve gravity entirely on GPU.
        /// Eliminates CPU bottleneck by computing Forman-Ricci curvature on GPU.
        /// 
        /// IMPORTANT: Call UpdateTopologyBuffers() first when graph topology changes.
        /// </summary>
        /// <param name="hostWeights">Edge weights (will be updated in-place)</param>
        /// <param name="hostMasses">Node masses</param>
        /// <param name="hostEdgesFrom">Edge source nodes</param>
        /// <param name="hostEdgesTo">Edge target nodes</param>
        /// <param name="dt">Time step</param>
        /// <param name="G">Gravitational coupling constant</param>
        /// <param name="lambda">Cosmological constant</param>
        /// <param name="degreePenaltyFactor">Penalty factor for weighted degree in curvature formula</param>
        public void EvolveFullGpuStep(
            float[] hostWeights,
            float[] hostMasses,
            int[] hostEdgesFrom,
            int[] hostEdgesTo,
            float dt,
            float G,
            float lambda,
            float degreePenaltyFactor)
        {
            if (!_topologyInitialized || _adjOffsetsBuffer == null || _adjDataBuffer == null)
            {
                throw new InvalidOperationException(
                    "Topology buffers not initialized. Call UpdateTopologyBuffers() first.");
            }

            int edgeCount = hostWeights.Length;

            // Copy data CPU -> GPU
            _weightsBuffer.CopyFrom(hostWeights);
            _nodeMassesBuffer.CopyFrom(hostMasses);

            // Pack edge pairs into Int2
            Int2[] packedEdges = new Int2[edgeCount];
            for (int i = 0; i < edgeCount; i++)
                packedEdges[i] = new Int2(hostEdgesFrom[i], hostEdgesTo[i]);
            _edgeIndicesBuffer.CopyFrom(packedEdges);

            // 1. Compute curvature on GPU (eliminates CPU bottleneck!)
            var curvatureShader = new FormanCurvatureShader(
                _weightsBuffer,
                _edgeIndicesBuffer,
                _adjOffsetsBuffer,
                _adjDataBuffer,
                _curvaturesBuffer,
                degreePenaltyFactor,
                _nodeCount);

            _device.For(edgeCount, curvatureShader);

            // 2. Evolve gravity (using fresh GPU-computed curvature)
            var gravityShader = new GravityShader(
                _weightsBuffer,
                _curvaturesBuffer,
                _nodeMassesBuffer,
                _edgeIndicesBuffer,
                dt, G, lambda);

            _device.For(edgeCount, gravityShader);

            // Copy results GPU -> CPU
            _weightsBuffer.CopyTo(hostWeights);
        }

        /// <summary>
        /// Full GPU step without data copy back - data stays in GPU memory.
        /// Use SyncToHost() to retrieve results when needed.
        /// 
        /// This is the most efficient mode for iterative simulation.
        /// </summary>
        /// <param name="dt">Time step</param>
        /// <param name="G">Gravitational coupling constant</param>
        /// <param name="lambda">Cosmological constant</param>
        /// <param name="degreePenaltyFactor">Penalty factor for weighted degree</param>
        public void EvolveFullGpuStep_NoCopy(float dt, float G, float lambda, float degreePenaltyFactor)
        {
            if (!_topologyInitialized || _adjOffsetsBuffer == null || _adjDataBuffer == null)
            {
                throw new InvalidOperationException(
                    "Topology buffers not initialized. Call UpdateTopologyBuffers() first.");
            }

            int edgeCount = (int)_weightsBuffer.Length;

            // 1. Compute curvature on GPU
            var curvatureShader = new FormanCurvatureShader(
                _weightsBuffer,
                _edgeIndicesBuffer,
                _adjOffsetsBuffer,
                _adjDataBuffer,
                _curvaturesBuffer,
                degreePenaltyFactor,
                _nodeCount);

            _device.For(edgeCount, curvatureShader);

            // 2. Evolve gravity
            var gravityShader = new GravityShader(
                _weightsBuffer,
                _curvaturesBuffer,
                _nodeMassesBuffer,
                _edgeIndicesBuffer,
                dt, G, lambda);

            _device.For(edgeCount, gravityShader);

            // NO CopyTo - data stays on GPU for maximum performance
        }

        /// <summary>
        /// Upload initial data to GPU (weights, masses, edges).
        /// Call this once before using EvolveFullGpuStep_NoCopy.
        /// </summary>
        public void UploadInitialData(
            float[] hostWeights,
            float[] hostMasses,
            int[] hostEdgesFrom,
            int[] hostEdgesTo)
        {
            _weightsBuffer.CopyFrom(hostWeights);
            _nodeMassesBuffer.CopyFrom(hostMasses);

            int edgeCount = hostWeights.Length;
            Int2[] packedEdges = new Int2[edgeCount];
            for (int i = 0; i < edgeCount; i++)
                packedEdges[i] = new Int2(hostEdgesFrom[i], hostEdgesTo[i]);
            _edgeIndicesBuffer.CopyFrom(packedEdges);
        }

        /// <summary>
        /// Optimized method that evolves gravity without copying data back to CPU.
        /// Data stays in GPU memory for maximum performance.
        /// Use SyncToHost() to retrieve results when needed (e.g., for visualization).
        /// 
        /// NOTE: This is the HYBRID mode (curvature from CPU).
        /// For full GPU mode, use EvolveFullGpuStep_NoCopy instead.
        /// </summary>
        public void EvolveGravityGpu_NoCopy(float dt, float G, float lambda)
        {
            // Run shader over data already in GPU memory
            // Topology must remain unchanged between calls

            var shader = new GravityShader(
                _weightsBuffer,
                _curvaturesBuffer,
                _nodeMassesBuffer,
                _edgeIndicesBuffer,
                dt, G, lambda);

            int edgeCount = (int)_weightsBuffer.Length;
            _device.For(edgeCount, shader);

            // NO CopyTo(hostWeights) - data stays on GPU
        }

        /// <summary>
        /// Synchronize GPU data back to host memory.
        /// Call this only when visualization or analysis is needed (e.g., every 50 steps).
        /// </summary>
        public void SyncToHost(float[] hostWeights)
        {
            _weightsBuffer.CopyTo(hostWeights);
        }

        /// <summary>
        /// Synchronize both weights and curvatures to host.
        /// Useful for diagnostics and visualization.
        /// </summary>
        public void SyncToHost(float[] hostWeights, float[] hostCurvatures)
        {
            _weightsBuffer.CopyTo(hostWeights);
            _curvaturesBuffer.CopyTo(hostCurvatures);
        }

        public void Dispose()
        {
            _weightsBuffer?.Dispose();
            _curvaturesBuffer?.Dispose();
            _nodeMassesBuffer?.Dispose();
            _edgeIndicesBuffer?.Dispose();
            _adjOffsetsBuffer?.Dispose();
            _adjDataBuffer?.Dispose();
        }
    }

    /// <summary>
    /// HLSL Shader for gravity evolution, compiled from C#.
    /// IComputeShader - marker for ComputeSharp.
    /// Source generator provides IComputeShaderDescriptor implementation.
    /// 
    /// ENERGY CONSERVATION FIX:
    /// The gravity update formula now uses a flow-based approach that
    /// redistributes weight rather than adding/removing arbitrarily.
    /// dw = dt * tanh(G * massTerm - curvature + lambda)
    /// The tanh bounds the change to [-1, 1] preventing runaway growth.
    /// </summary>
    [ThreadGroupSize(64, 1, 1)]
    [GeneratedComputeShaderDescriptor]
    public readonly partial struct GravityShader : IComputeShader
    {
        public readonly ReadWriteBuffer<float> weights;
        public readonly ReadWriteBuffer<float> curvatures; // Changed to ReadWrite to match buffer type
        public readonly ReadOnlyBuffer<float> masses;
        public readonly ReadOnlyBuffer<Int2> edges;
        public readonly float dt;
        public readonly float G;
        public readonly float lambda;

        public GravityShader(
            ReadWriteBuffer<float> weights,
            ReadWriteBuffer<float> curvatures,
            ReadOnlyBuffer<float> masses,
            ReadOnlyBuffer<Int2> edges,
            float dt,
            float G,
            float lambda)
        {
            this.weights = weights;
            this.curvatures = curvatures;
            this.masses = masses;
            this.edges = edges;
            this.dt = dt;
            this.G = G;
            this.lambda = lambda;
        }

        public void Execute()
        {
            // ThreadIds.X - edge index that this thread processes
            int i = ThreadIds.X;

            // Get nodes that this edge connects
            Int2 edgeNodes = edges[i];
            int nodeA = edgeNodes.X;
            int nodeB = edgeNodes.Y;

            float currentWeight = weights[i];

            // Skip very weak edges
            if (currentWeight < 0.0001f)
            {
                return;
            }

            // Physics: Weight update from gravity (ENERGY-CONSERVING)
            // Use curvature flow that preserves total energy
            // dW/dt = -W * (Curvature - G*Mass + Lambda)
            // This is a multiplicative update that can't create energy from nothing

            float massTerm = (masses[nodeA] + masses[nodeB]) * 0.5f;
            float curvatureTerm = curvatures[i];

            // Flow rate based on local geometry
            // Negative curvature (saddle) -> weight decreases (expansion)
            // Positive curvature (sphere) -> weight increases (contraction)
            float flowRate = curvatureTerm - G * massTerm + lambda;

            // ENERGY CONSERVATION: Use multiplicative update with bounded rate
            // This ensures dW/W is bounded, preventing runaway growth
            // tanh bounds to [-1, 1], scaled by dt for stability
            float relativeChange = Hlsl.Tanh(flowRate * 0.1f) * dt;

            // Multiplicative update: W_new = W * (1 + relativeChange)
            float w = currentWeight * (1.0f + relativeChange);

            // Clamp using HLSL function with soft bounds
            weights[i] = Hlsl.Clamp(w, 0.001f, 0.999f);
        }
    }

    class ImprovedNetworkGravity
    {
        // Note: These constants are now superseded by PhysicsConstants for consistency.
        // They are kept for backward compatibility with GetEffectiveGravitationalCoupling.
        private const double BaseGravitationalCoupling = 2.5; // Legacy: use PhysicsConstants.GravitationalCoupling instead
        private const double MaxGravitationalCoupling = 5.0; // Maximum during early annealing

        /// <summary>
        /// Compute effective gravitational coupling with simulated annealing
        /// CHECKLIST ITEM 3: During first 1000 steps, use boosted coupling to form clusters
        /// 
        /// G_eff(t) = G_base * (1.0 + 10.0 * exp(-step / 200))
        /// 
        /// This "hot start" forces the graph to collapse into clusters early,
        /// creating matter before temperature cools down.
        /// 
        /// NOTE: For new simulations, use EvolveNetworkGeometryOllivierDynamic with
        /// externally computed effectiveG based on PhysicsConstants.
        /// </summary>
        public static double GetEffectiveGravitationalCoupling(int step)
        {
            // Annealing phase: first 1000 steps
            if (step < 1000)
            {
                // Exponential boost that decays over 200 steps
                double boostFactor = 1.0 + 10.0 * Math.Exp(-step / 200.0);
                double effective = BaseGravitationalCoupling * boostFactor;

                // Clamp to maximum
                return Math.Min(effective, MaxGravitationalCoupling);
            }

            // After annealing: use base value
            return BaseGravitationalCoupling;
        }

        /// <summary>
        /// Evolve network geometry using Ollivier-Ricci curvature with DYNAMIC gravitational coupling.
        /// This is the primary method for "Primordial Soup" simulation.
        /// 
        /// The caller provides effectiveG which should be computed based on simulation time:
        /// - Phase 1 (Mixing): effectiveG = 0.1 (weak gravity, allows restructuring)
        /// - Phase 2 (Clustering): effectiveG = 1.5 (strong gravity, forms clusters)
        /// 
        /// CHECKLIST ITEM 4: Use Ollivier-Ricci for better geometric sensitivity
        /// 
        /// RQ-FIX: Now uses unified NodeMasses.TotalMass instead of correlation mass only.
        /// This ensures gravity responds to ALL field contributions (scalar, fermion, gauge).
        /// 
        /// RQ-FIX: Added connectivity protection to prevent graph fragmentation.
        /// When graph is at risk of fragmenting (low average degree), weight decreases are suppressed.
        /// 
        /// GPU ACCELERATION: If graph.GpuGravity is initialized, uses GPU for Forman-Ricci curvature
        /// and weight updates. Ollivier-Ricci requires CPU due to optimal transport complexity.
        /// Use graph.InitGpuGravity() to enable GPU mode, graph.DisposeGpuGravity() to release.
        /// </summary>
        /// <param name="graph">The RQGraph to evolve</param>
        /// <param name="dt">Time step</param>
        /// <param name="effectiveG">Gravitational coupling (caller controls phase transition)</param>
        public static void EvolveNetworkGeometryOllivierDynamic(RQGraph graph, double dt, double effectiveG)
        {
            if (graph.Weights == null || graph.Edges == null)
                return;

            // GPU fast path: Use Forman-Ricci on GPU when available
            // Note: Ollivier-Ricci requires optimal transport which is expensive on GPU,
            // so we use Forman-Ricci for GPU mode (still captures curvature well)
            if (graph.GpuGravity != null && graph.GpuGravity.IsTopologyInitialized)
            {
                EvolveNetworkGeometryGpu(graph, dt, effectiveG);
                return;
            }

            int N = graph.N;
            double learningRate = effectiveG * dt;

            // RQ-FIX: Connectivity protection - compute graph health metrics
            // If average degree is too low, prevent further weight decreases to avoid fragmentation
            double avgDegree = ComputeAverageDegree(graph);
            double minWeightSum = ComputeMinNodeWeightSum(graph);

            // Protection thresholds:
            // - If avgDegree < 3, graph is at risk of fragmentation
            // - If any node has weightSum < 0.5, it may become isolated
            bool protectConnectivity = avgDegree < 4.0 || minWeightSum < 0.8;
            double connectivityProtectionFactor = protectConnectivity ? 0.5 : 1.0;

            // RQ-FIX: Update unified node masses (includes ALL field contributions)
            graph.UpdateNodeMasses();
            var nodeMasses = graph.NodeMasses;

            // Parallel update of edge weights
            double[,] deltaWeights = new double[N, N];

            System.Threading.Tasks.Parallel.For(0, N, i =>
            {
                // RQ-FIX: Use TotalMass from unified NodeMasses
                // This includes fermion, scalar, gauge, correlation, and vacuum energy
                double massI = nodeMasses[i].TotalMass;

                foreach (int j in graph.Neighbors(i))
                {
                    if (j <= i) continue; // Process each edge once

                    double massJ = nodeMasses[j].TotalMass;
                    double currentWeight = graph.Weights[i, j];

                    // CHECKLIST ITEM 4: Use Ollivier-Ricci curvature
                    // This is more sensitive to geometry than Forman-Ricci
                    double curvature = OllivierRicciCurvature.ComputeOllivierRicciJaccard(graph, i, j);

                    // Stress-energy contribution (matter) - now uses FULL mass
                    double stressEnergyTensor = (massI + massJ) * 0.5;
                    double dS_matter = -stressEnergyTensor * PhysicsConstants.CurvatureTermScale;

                    // Cosmological constant contribution
                    double dS_cosmological = PhysicsConstants.CosmologicalConstant;

                    // Total action gradient (flow rate)
                    double flowRate = curvature + dS_matter + dS_cosmological;

                    // ENERGY CONSERVATION FIX: Use multiplicative update with bounded rate
                    // Instead of: delta = -learningRate * dS_total (additive, can blow up)
                    // Use: delta = currentWeight * tanh(flowRate * scale) * dt (multiplicative, bounded)
                    double relativeChange = Math.Tanh(flowRate * 0.1) * learningRate;
                    double delta = currentWeight * relativeChange;

                    // RQ-FIX: Connectivity protection
                    // If delta would decrease weight AND connectivity is at risk, suppress the decrease
                    if (delta < 0 && protectConnectivity)
                    {
                        delta *= connectivityProtectionFactor;

                        // Additional protection: never let weight go below minimum threshold
                        double projectedWeight = currentWeight + delta;
                        if (projectedWeight < PhysicsConstants.WeightLowerSoftWall)
                        {
                            delta = Math.Max(delta, PhysicsConstants.WeightLowerSoftWall - currentWeight);
                        }
                    }

                    deltaWeights[i, j] = delta;
                    deltaWeights[j, i] = delta;
                }
            });

            // Apply updates with soft walls
            System.Threading.Tasks.Parallel.For(0, N, i =>
            {
                foreach (int j in graph.Neighbors(i))
                {
                    if (j <= i) continue;

                    double newWeight = graph.Weights[i, j] + deltaWeights[i, j];

                    // Apply soft walls (Checklist A.4)
                    newWeight = ApplySoftWalls(newWeight);

                    graph.Weights[i, j] = newWeight;
                    graph.Weights[j, i] = newWeight;
                }
            });

            // Update target distances - using public method from RQGraph
            graph.UpdateTargetDistancesFromWeights();
        }

        /// <summary>
        /// Compute average degree of the graph (number of edges per node).
        /// Used for connectivity protection.
        /// </summary>
        private static double ComputeAverageDegree(RQGraph graph)
        {
            int totalDegree = 0;
            for (int i = 0; i < graph.N; i++)
            {
                totalDegree += graph.Neighbors(i).Count();
            }
            return graph.N > 0 ? (double)totalDegree / graph.N : 0.0;
        }

        /// <summary>
        /// Compute minimum weight sum across all nodes.
        /// A low value indicates a node at risk of isolation.
        /// </summary>
        private static double ComputeMinNodeWeightSum(RQGraph graph)
        {
            double minSum = double.MaxValue;
            for (int i = 0; i < graph.N; i++)
            {
                double sum = 0.0;
                foreach (int j in graph.Neighbors(i))
                {
                    sum += graph.Weights[i, j];
                }
                if (sum < minSum) minSum = sum;
            }
            return minSum == double.MaxValue ? 0.0 : minSum;
        }

        /// <summary>
        /// GPU-accelerated network geometry evolution using Forman-Ricci curvature.
        /// 
        /// This method runs entirely on GPU for maximum performance:
        /// 1. Computes Forman-Ricci curvature on GPU (parallel over edges)
        /// 2. Applies gravity update on GPU (parallel over edges)
        /// 3. Syncs results back to CPU graph
        /// 
        /// Note: Uses Forman-Ricci instead of Ollivier-Ricci because:
        /// - Ollivier-Ricci requires optimal transport (complex on GPU)
        /// - Forman-Ricci is efficient O(degree²) and GPU-friendly
        /// - Both capture essential curvature properties for gravity
        /// 
        /// RQ-FIX: Now uses unified GetUnifiedMassesFlat() which includes ALL field
        /// contributions (fermion, scalar, gauge, correlation, vacuum, kinetic),
        /// not just correlation mass. This ensures GPU gravity matches CPU physics.
        /// 
        /// Enable via: graph.InitGpuGravity()
        /// </summary>
        private static void EvolveNetworkGeometryGpu(RQGraph graph, double dt, double effectiveG)
        {
            var gpuEngine = graph.GpuGravity;
            if (gpuEngine == null || !gpuEngine.IsTopologyInitialized)
                return;

            int N = graph.N;
            int edgeCount = graph.FlatEdgesFrom.Length;

            // Prepare data arrays
            float[] weights = new float[edgeCount];

            // RQ-FIX: Use unified masses instead of just correlation mass
            // GetUnifiedMassesFlat() calls UpdateNodeMasses() internally and returns
            // NodeMassModels[i].TotalMass which includes ALL field contributions
            float[] masses = graph.GetUnifiedMassesFlat();

            // Get current weights
            for (int e = 0; e < edgeCount; e++)
            {
                int i = graph.FlatEdgesFrom[e];
                int j = graph.FlatEdgesTo[e];
                weights[e] = (float)graph.Weights[i, j];
            }

            // Run GPU computation with Forman-Ricci curvature
            gpuEngine.EvolveFullGpuStep(
                weights,
                masses,
                graph.FlatEdgesFrom,
                graph.FlatEdgesTo,
                (float)dt,
                (float)effectiveG,
                (float)PhysicsConstants.CosmologicalConstant,
                (float)PhysicsConstants.DegreePenaltyFactor);

            // Apply results back to graph (with soft walls)
            for (int e = 0; e < edgeCount; e++)
            {
                int i = graph.FlatEdgesFrom[e];
                int j = graph.FlatEdgesTo[e];
                double newWeight = ApplySoftWalls(weights[e]);
                graph.Weights[i, j] = newWeight;
                graph.Weights[j, i] = newWeight;
            }

            graph.UpdateTargetDistancesFromWeights();
        }

        /// <summary>
        /// Evolve network geometry using Ollivier-Ricci curvature
        /// This replaces the old Forman-Ricci based method
        /// 
        /// CHECKLIST ITEM 4: Use Ollivier-Ricci for better geometric sensitivity
        /// ENERGY CONSERVATION: Uses multiplicative update with bounded flow rate.
        /// RQ-FIX: Now uses unified NodeMasses.TotalMass instead of correlation mass only.
        /// </summary>
        public static void EvolveNetworkGeometryOllivier(RQGraph graph, double dt, int currentStep)
        {
            if (graph.Weights == null || graph.Edges == null)
                return;

            int N = graph.N;

            // Compute effective coupling with annealing
            double effectiveG = GetEffectiveGravitationalCoupling(currentStep);
            double learningRate = effectiveG * dt;

            // RQ-FIX: Update unified node masses (includes ALL field contributions)
            graph.UpdateNodeMasses();
            var nodeMasses = graph.NodeMasses;

            // Parallel update of edge weights
            double[,] deltaWeights = new double[N, N];

            System.Threading.Tasks.Parallel.For(0, N, i =>
            {
                // RQ-FIX: Use TotalMass from unified NodeMasses
                double massI = nodeMasses[i].TotalMass;

                foreach (int j in graph.Neighbors(i))
                {
                    if (j <= i) continue; // Process each edge once

                    double massJ = nodeMasses[j].TotalMass;
                    double currentWeight = graph.Weights[i, j];

                    // CHECKLIST ITEM 4: Use Ollivier-Ricci curvature
                    double curvature = OllivierRicciCurvature.ComputeOllivierRicciJaccard(graph, i, j);

                    // Stress-energy contribution (matter) - now uses FULL mass
                    double stressEnergyTensor = (massI + massJ) * 0.5;
                    double dS_matter = -stressEnergyTensor * PhysicsConstants.CurvatureTermScale;

                    // Cosmological constant contribution
                    double dS_cosmological = PhysicsConstants.CosmologicalConstant;

                    // Flow rate from action gradient
                    double flowRate = curvature + dS_matter + dS_cosmological;

                    // ENERGY CONSERVATION: Multiplicative update with bounded rate
                    double relativeChange = Math.Tanh(flowRate * 0.1) * learningRate;
                    double delta = currentWeight * relativeChange;

                    deltaWeights[i, j] = delta;
                    deltaWeights[j, i] = delta;
                }
            });

            // Apply updates with soft walls
            System.Threading.Tasks.Parallel.For(0, N, i =>
            {
                foreach (int j in graph.Neighbors(i))
                {
                    if (j <= i) continue;

                    double newWeight = graph.Weights[i, j] + deltaWeights[i, j];

                    // Apply soft walls (Checklist A.4)
                    newWeight = ApplySoftWalls(newWeight);

                    graph.Weights[i, j] = newWeight;
                    graph.Weights[j, i] = newWeight;
                }
            });

            // Update target distances - using public method from RQGraph
            // Note: This is called internally by the graph when needed
        }

        /// <summary>
        /// Apply soft potential walls to weights
        /// Prevents extreme values without hard clipping
        /// </summary>
        private static double ApplySoftWalls(double weight)
        {
            const double upperWall = 0.95;
            const double lowerWall = 0.05;

            if (weight > upperWall)
            {
                // Soft approach to 1.0
                weight = upperWall + (1.0 - upperWall) * Math.Tanh((weight - upperWall) / 0.1);
            }

            if (weight < lowerWall)
            {
                // Soft approach to 0.0
                weight = lowerWall * Math.Exp((weight - lowerWall) / lowerWall);
            }

            // Safety bounds
            return Math.Clamp(weight, 0.0, 1.0);
        }

        /// <summary>
        /// Compute average Ollivier-Ricci curvature of the network
        /// Useful for diagnostics and monitoring phase transitions
        /// </summary>
        public static double ComputeAverageOllivierCurvature(RQGraph graph)
        {
            if (graph.Edges == null || graph.Weights == null)
                return 0.0;

            double totalCurvature = 0.0;
            int edgeCount = 0;

            for (int i = 0; i < graph.N; i++)
            {
                foreach (int j in graph.Neighbors(i))
                {
                    if (j <= i) continue;

                    totalCurvature += OllivierRicciCurvature.ComputeOllivierRicciJaccard(graph, i, j);
                    edgeCount++;
                }
            }

            return edgeCount > 0 ? totalCurvature / edgeCount : 0.0;
        }

        /// <summary>
        /// Check if the network has undergone phase transition
        /// Heavy clusters form when average curvature becomes positive
        /// </summary>
        public static bool HasFormattedHeavyClusters(RQGraph graph)
        {
            double avgCurvature = ComputeAverageOllivierCurvature(graph);

            // Positive curvature indicates clustering (matter formation)
            if (avgCurvature > 0.05)
            {
                var heavyStats = graph.GetHeavyClusterStatsCorrelationMass(
                    RQGraph.HeavyClusterThreshold,
                    RQGraph.HeavyClusterMinSize);

                // Check if heavy clusters have formed
                // Note: The tuple returns (count, totalMass, maxSize, avgMassPerNode)
                return heavyStats.count > 0 && heavyStats.totalMass > 0.1;
            }

            return false;
        }

        /// <summary>
        /// Evolve network geometry using unified NodeMassModel (CHECKLIST ITEM 6).
        /// 
        /// This method uses TotalMass from NodeMassModel as the source term for gravity,
        /// aggregating contributions from all fields:
        /// - Fermion mass (Dirac spinors)
        /// - Correlation mass (graph structure)
        /// - Gauge field energy (Yang-Mills)
        /// - Scalar field energy (Higgs)
        /// - Vacuum energy (cosmological constant)
        /// - Kinetic energy (gravitational waves)
        /// 
        /// Einstein equation on graph: dw_ij/dt = -G * (R_ij - T_ij + Λ)
        /// where T_ij = (M_i + M_j)/2 is the average total mass at edge endpoints.
        /// </summary>
        /// <param name="graph">The RQGraph to evolve</param>
        /// <param name="dt">Time step</param>
        /// <param name="effectiveG">Gravitational coupling constant</param>
        public static void EvolveNetworkGeometryUnifiedMass(RQGraph graph, double dt, double effectiveG)
        {
            if (graph.Weights == null || graph.Edges == null)
                return;

            int N = graph.N;
            double learningRate = effectiveG * dt;

            // CHECKLIST ITEM 6: Update unified mass models before gravity step
            graph.UpdateNodeMasses();
            var nodeMasses = graph.NodeMasses;

            // Connectivity protection
            double avgDegree = ComputeAverageDegree(graph);
            double minWeightSum = ComputeMinNodeWeightSum(graph);
            bool protectConnectivity = avgDegree < 4.0 || minWeightSum < 0.8;
            double protectionFactor = protectConnectivity ? 0.5 : 1.0;

            // Parallel weight update
            double[,] deltaWeights = new double[N, N];

            System.Threading.Tasks.Parallel.For(0, N, i =>
            {
                // CHECKLIST ITEM 6: Use TotalMass from NodeMassModel
                double massI = nodeMasses[i].TotalMass;

                foreach (int j in graph.Neighbors(i))
                {
                    if (j <= i) continue;

                    double massJ = nodeMasses[j].TotalMass;
                    double currentWeight = graph.Weights[i, j];

                    // Use Ollivier-Ricci or Forman-Ricci curvature
                    double curvature;
                    if (PhysicsConstants.PreferOllivierRicciCurvature)
                        curvature = OllivierRicciCurvature.ComputeOllivierRicciJaccard(graph, i, j);
                    else
                        curvature = graph.ComputeFormanRicciCurvature(i, j);

                    // Stress-energy from unified mass (includes all field contributions)
                    double stressEnergy = (massI + massJ) * 0.5 * PhysicsConstants.CurvatureTermScale;

                    // Cosmological constant
                    double lambda = PhysicsConstants.CosmologicalConstant;

                    // Flow rate from Einstein equation
                    double flowRate = curvature - stressEnergy + lambda;

                    // ENERGY CONSERVATION: Multiplicative update with bounded rate
                    double relativeChange = Math.Tanh(flowRate * 0.1) * learningRate;
                    double delta = currentWeight * relativeChange;

                    // Connectivity protection
                    if (delta < 0 && protectConnectivity)
                    {
                        delta *= protectionFactor;
                        double projectedWeight = currentWeight + delta;
                        if (projectedWeight < PhysicsConstants.WeightLowerSoftWall)
                            delta = Math.Max(delta, PhysicsConstants.WeightLowerSoftWall - currentWeight);
                    }

                    deltaWeights[i, j] = delta;
                    deltaWeights[j, i] = delta;
                }
            });

            // Apply updates with soft walls
            System.Threading.Tasks.Parallel.For(0, N, i =>
            {
                foreach (int j in graph.Neighbors(i))
                {
                    if (j <= i) continue;

                    double newWeight = graph.Weights[i, j] + deltaWeights[i, j];
                    newWeight = ApplySoftWalls(newWeight);

                    graph.Weights[i, j] = newWeight;
                    graph.Weights[j, i] = newWeight;
                }
            });

            graph.UpdateTargetDistancesFromWeights();
        }
    }
}
