using System;
using System.Collections.Generic;
using System.Linq;
using ComputeSharp;

namespace RQSimulation.GPUOptimized
{
    /// <summary>
    /// Curvature computations for graph geometry: Ollivier-Ricci and Forman-Ricci.
    /// 
    /// CHECKLIST ITEM 4: Both Ollivier and Forman curvatures are provided.
    /// 
    /// Ollivier-Ricci (more accurate but slower):
    /// - Based on optimal transport between probability measures
    /// - Captures geodesic deviation like Einstein's gravity
    /// 
    /// Forman-Ricci (faster approximation):
    /// - Based on local combinatorics (degrees and triangles)
    /// - O(degree²) per edge vs O(N²) for full Wasserstein
    /// - Good approximation for Ricci flow on graphs
    /// </summary>
    public static class OllivierRicciCurvature
    {
        /// <summary>
        /// Compute Ollivier-Ricci curvature for edge (i,j)
        /// 
        /// κ(i,j) = 1 - W₁(μᵢ, μⱼ) / d(i,j)
        /// 
        /// where:
        /// - W₁ is the Wasserstein-1 distance between probability measures
        /// - μᵢ is the probability distribution on i's neighborhood
        /// - d(i,j) is the edge weight (distance)
        /// 
        /// Positive curvature → space contracts (like sphere)
        /// Negative curvature → space expands (like hyperbolic space)
        /// </summary>
        public static double ComputeOllivierRicci(RQGraph graph, int i, int j)
        {
            if (!graph.Edges[i, j])
                return 0.0;

            double edgeWeight = graph.Weights[i, j];
            if (edgeWeight <= 0)
                return 0.0;

            // Get weighted probability distributions on neighborhoods
            var distI = GetWeightedNeighborhood(graph, i);
            var distJ = GetWeightedNeighborhood(graph, j);

            // Compute Wasserstein-1 distance
            double w1Distance = ComputeWassersteinDistance(graph, distI, distJ);

            // Ollivier-Ricci curvature formula
            double curvature = 1.0 - (w1Distance / edgeWeight);

            return curvature;
        }

        /// <summary>
        /// Compute simplified Ollivier-Ricci using Jaccard index approximation
        /// This is much faster than full Wasserstein distance computation
        /// 
        /// Uses: W₁ ≈ d * (1 - Jaccard_weighted)
        /// where Jaccard_weighted measures neighborhood overlap
        /// </summary>
        public static double ComputeOllivierRicciJaccard(RQGraph graph, int i, int j)
        {
            if (!graph.Edges[i, j])
                return 0.0;

            double edgeWeight = graph.Weights[i, j];
            if (edgeWeight <= 0)
                return 0.0;

            // Get neighborhoods
            var neighborsI = new HashSet<int>(graph.Neighbors(i));
            var neighborsJ = new HashSet<int>(graph.Neighbors(j));

            // Compute weighted Jaccard index
            double intersection = 0.0;
            double union_sum = 0.0;

            // Add weights from i's neighborhood
            foreach (int k in neighborsI)
            {
                double w_ik = graph.Weights[i, k];
                union_sum += w_ik;

                if (neighborsJ.Contains(k))
                {
                    // k is in both neighborhoods
                    double w_jk = graph.Weights[j, k];
                    intersection += Math.Min(w_ik, w_jk);
                }
            }

            // Add weights from j's neighborhood (not already counted)
            foreach (int k in neighborsJ)
            {
                if (!neighborsI.Contains(k))
                {
                    double w_jk = graph.Weights[j, k];
                    union_sum += w_jk;
                }
            }

            double jaccard = union_sum > 0 ? intersection / union_sum : 0.0;

            // Approximate Wasserstein distance
            double w1Approx = edgeWeight * (1.0 - jaccard);

            // Ollivier-Ricci curvature
            double curvature = 1.0 - (w1Approx / edgeWeight);

            return curvature;
        }

        /// <summary>
        /// Get weighted probability distribution on node's neighborhood
        /// Returns dictionary: neighbor -> probability
        /// </summary>
        private static Dictionary<int, double> GetWeightedNeighborhood(RQGraph graph, int node)
        {
            var distribution = new Dictionary<int, double>();
            double totalWeight = 0.0;

            // Self-loop with small weight (lazy random walk)
            const double alpha = 0.1; // Probability of staying at current node
            distribution[node] = alpha;
            totalWeight += alpha;

            // Neighbors with weights proportional to edge weights
            foreach (int neighbor in graph.Neighbors(node))
            {
                double weight = graph.Weights[node, neighbor];
                distribution[neighbor] = (1.0 - alpha) * weight;
                totalWeight += (1.0 - alpha) * weight;
            }

            // Normalize to probability distribution
            if (totalWeight > 0)
            {
                var normalized = new Dictionary<int, double>();
                foreach (var kvp in distribution)
                {
                    normalized[kvp.Key] = kvp.Value / totalWeight;
                }
                return normalized;
            }

            return distribution;
        }

        /// <summary>
        /// Compute Wasserstein-1 (Earth Mover's Distance) between two discrete distributions
        /// 
        /// Uses linear programming formulation:
        /// W₁(μ, ν) = min Σᵢⱼ d(i,j) * T[i,j]
        /// subject to: row sums = μ, column sums = ν
        /// 
        /// For graph distributions, we use a greedy approximation:
        /// - Match closest pairs first (by graph distance)
        /// - Subtract matched mass
        /// - Repeat until all mass transported
        /// </summary>
        private static double ComputeWassersteinDistance(
            RQGraph graph,
            Dictionary<int, double> distI,
            Dictionary<int, double> distJ)
        {
            // Create mass arrays
            var massI = new Dictionary<int, double>(distI);
            var massJ = new Dictionary<int, double>(distJ);

            double totalCost = 0.0;

            // Greedy matching: transport mass along shortest paths
            while (massI.Count > 0 && massJ.Count > 0)
            {
                // Find closest pair (minimum transport distance)
                double minDist = double.MaxValue;
                int bestSrc = -1;
                int bestDst = -1;

                foreach (var src in massI.Keys)
                {
                    foreach (var dst in massJ.Keys)
                    {
                        double dist = ComputeGraphDistance(graph, src, dst);
                        if (dist < minDist)
                        {
                            minDist = dist;
                            bestSrc = src;
                            bestDst = dst;
                        }
                    }
                }

                if (bestSrc == -1 || bestDst == -1)
                    break;

                // Transport mass along this edge
                double transportMass = Math.Min(massI[bestSrc], massJ[bestDst]);
                totalCost += transportMass * minDist;

                // Update remaining mass
                massI[bestSrc] -= transportMass;
                massJ[bestDst] -= transportMass;

                if (massI[bestSrc] < 1e-10)
                    massI.Remove(bestSrc);
                if (massJ[bestDst] < 1e-10)
                    massJ.Remove(bestDst);
            }

            return totalCost;
        }

        /// <summary>
        /// Compute graph distance between two nodes (shortest path)
        /// Uses Dijkstra's algorithm with edge weights as distances
        /// </summary>
        private static double ComputeGraphDistance(RQGraph graph, int source, int target)
        {
            if (source == target)
                return 0.0;

            int N = graph.N;
            double[] dist = new double[N];
            bool[] visited = new bool[N];

            for (int i = 0; i < N; i++)
                dist[i] = double.MaxValue;

            dist[source] = 0.0;

            // Dijkstra's algorithm
            for (int iter = 0; iter < N; iter++)
            {
                // Find unvisited node with minimum distance
                int u = -1;
                double minDist = double.MaxValue;

                for (int i = 0; i < N; i++)
                {
                    if (!visited[i] && dist[i] < minDist)
                    {
                        minDist = dist[i];
                        u = i;
                    }
                }

                if (u == -1 || u == target)
                    break;

                visited[u] = true;

                // Update neighbors
                foreach (int v in graph.Neighbors(u))
                {
                    if (!visited[v])
                    {
                        double alt = dist[u] + graph.Weights[u, v];
                        if (alt < dist[v])
                        {
                            dist[v] = alt;
                        }
                    }
                }
            }

            return dist[target];
        }
    }

    /// <summary>
    /// GPU-accelerated Ollivier-Ricci curvature computation engine
    /// Uses ComputeSharp for parallel edge curvature calculations
    /// </summary>
    public class GpuCurvatureEngine : IDisposable
    {
        private readonly GraphicsDevice _device;
        private ReadWriteBuffer<float>? _curvaturesBuffer;
        private ReadOnlyBuffer<float>? _weightsBuffer;
        private ReadOnlyBuffer<Int2>? _edgesBuffer;
        private ReadOnlyBuffer<int>? _neighborOffsetsBuffer;
        private ReadOnlyBuffer<int>? _neighborIndicesBuffer;
        private int _edgeCount;
        private int _nodeCount;

        public GpuCurvatureEngine()
        {
            _device = GraphicsDevice.GetDefault();
        }

        /// <summary>
        /// Compute curvatures for all edges on GPU
        /// Returns array of curvatures indexed by edge
        /// </summary>
        public float[] ComputeAllCurvaturesGpu(
            float[] weights,
            int[] edgesFrom,
            int[] edgesTo,
            int[] neighborOffsets,
            int[] neighborIndices,
            int nodeCount)
        {
            _edgeCount = edgesFrom.Length;
            _nodeCount = nodeCount;

            // Allocate GPU buffers
            _curvaturesBuffer = _device.AllocateReadWriteBuffer<float>(_edgeCount);
            _weightsBuffer = _device.AllocateReadOnlyBuffer(weights);

            // Pack edges into Int2
            Int2[] packedEdges = new Int2[_edgeCount];
            for (int i = 0; i < _edgeCount; i++)
            {
                packedEdges[i] = new Int2(edgesFrom[i], edgesTo[i]);
            }
            _edgesBuffer = _device.AllocateReadOnlyBuffer(packedEdges);

            // CSR format for neighbors
            _neighborOffsetsBuffer = _device.AllocateReadOnlyBuffer(neighborOffsets);
            _neighborIndicesBuffer = _device.AllocateReadOnlyBuffer(neighborIndices);

            // Create and run shader
            var shader = new OllivierRicciJaccardShader(
                _curvaturesBuffer,
                _weightsBuffer,
                _edgesBuffer,
                _neighborOffsetsBuffer,
                _neighborIndicesBuffer,
                nodeCount);

            _device.For(_edgeCount, shader);

            // Copy results back
            float[] results = new float[_edgeCount];
            _curvaturesBuffer.CopyTo(results);

            return results;
        }

        public void Dispose()
        {
            _curvaturesBuffer?.Dispose();
            _weightsBuffer?.Dispose();
            _edgesBuffer?.Dispose();
            _neighborOffsetsBuffer?.Dispose();
            _neighborIndicesBuffer?.Dispose();
        }
    }

    /// <summary>
    /// GPU shader for computing Ollivier-Ricci curvature using Jaccard approximation
    /// Each thread processes one edge
    /// </summary>
    [ThreadGroupSize(64, 1, 1)]
    [GeneratedComputeShaderDescriptor]
    public readonly partial struct OllivierRicciJaccardShader : IComputeShader
    {
        public readonly ReadWriteBuffer<float> curvatures;
        public readonly ReadOnlyBuffer<float> weights;
        public readonly ReadOnlyBuffer<Int2> edges;
        public readonly ReadOnlyBuffer<int> neighborOffsets;
        public readonly ReadOnlyBuffer<int> neighborIndices;
        public readonly int nodeCount;

        public OllivierRicciJaccardShader(
            ReadWriteBuffer<float> curvatures,
            ReadOnlyBuffer<float> weights,
            ReadOnlyBuffer<Int2> edges,
            ReadOnlyBuffer<int> neighborOffsets,
            ReadOnlyBuffer<int> neighborIndices,
            int nodeCount)
        {
            this.curvatures = curvatures;
            this.weights = weights;
            this.edges = edges;
            this.neighborOffsets = neighborOffsets;
            this.neighborIndices = neighborIndices;
            this.nodeCount = nodeCount;
        }

        public void Execute()
        {
            int edgeIdx = ThreadIds.X;

            Int2 edge = edges[edgeIdx];
            int i = edge.X;
            int j = edge.Y;

            // Get edge weight (used as distance)
            float edgeWeight = GetWeight(i, j);
            if (edgeWeight <= 0.0f)
            {
                curvatures[edgeIdx] = 0.0f;
                return;
            }

            // Compute weighted Jaccard index between neighborhoods
            float intersection = 0.0f;
            float unionSum = 0.0f;

            // Get neighbor ranges for both nodes
            int startI = neighborOffsets[i];
            int endI = neighborOffsets[i + 1];
            int startJ = neighborOffsets[j];
            int endJ = neighborOffsets[j + 1];

            // Process neighbors of i
            for (int idx = startI; idx < endI; idx++)
            {
                int k = neighborIndices[idx];
                float w_ik = GetWeight(i, k);
                unionSum += w_ik;

                // Check if k is also neighbor of j
                bool inJ = false;
                for (int idx2 = startJ; idx2 < endJ; idx2++)
                {
                    if (neighborIndices[idx2] == k)
                    {
                        inJ = true;
                        break;
                    }
                }

                if (inJ)
                {
                    float w_jk = GetWeight(j, k);
                    intersection += Hlsl.Min(w_ik, w_jk);
                }
            }

            // Process neighbors of j not in i
            for (int idx = startJ; idx < endJ; idx++)
            {
                int k = neighborIndices[idx];

                // Check if k is NOT neighbor of i
                bool inI = false;
                for (int idx2 = startI; idx2 < endI; idx2++)
                {
                    if (neighborIndices[idx2] == k)
                    {
                        inI = true;
                        break;
                    }
                }

                if (!inI)
                {
                    float w_jk = GetWeight(j, k);
                    unionSum += w_jk;
                }
            }

            // Compute Jaccard index
            float jaccard = unionSum > 0.0f ? intersection / unionSum : 0.0f;

            // Ollivier-Ricci curvature: κ = 1 - W₁/d ≈ 1 - (1 - Jaccard) = Jaccard
            curvatures[edgeIdx] = jaccard;
        }

        /// <summary>
        /// Get weight between nodes i and j from flat weights array
        /// Assumes weights stored as nodeCount x nodeCount matrix flattened row-major
        /// </summary>
        private float GetWeight(int i, int j)
        {
            // Index into flattened NxN matrix
            int idx = i * nodeCount + j;
            return weights[idx];
        }
    }
    
    // ================================================================
    // CHECKLIST ITEM 4: FORMAN-RICCI CURVATURE
    // ================================================================
    
    /// <summary>
    /// Forman-Ricci curvature computation for graph geometry.
    /// 
    /// CHECKLIST ITEM 4: Fast local curvature approximation.
    /// 
    /// Forman-Ricci curvature for an edge e = (i,j):
    /// 
    ///   R_F(e) = 4 - deg(i) - deg(j) + 3 * Σ_△ (w_ik * w_jk)^(1/3)
    /// 
    /// where the sum is over triangles △ containing edge e, with nodes k.
    /// 
    /// Interpretation:
    /// - Base curvature = 4 (maximal for isolated edge)
    /// - Degree penalty: high-degree nodes reduce curvature (more "spreading")
    /// - Triangle bonus: shared neighbors increase curvature (more "clustering")
    /// 
    /// Positive curvature → sphere-like (clustering)
    /// Negative curvature → hyperbolic (tree-like, spreading)
    /// 
    /// Advantage over Ollivier-Ricci:
    /// - O(degree²) per edge instead of O(N²) for Wasserstein
    /// - Local computation, no global optimization
    /// - Sufficient for Ricci flow gravity simulations
    /// </summary>
    public static class FormanRicciCurvature
    {
        /// <summary>
        /// Compute Forman-Ricci curvature for edge (i,j).
        /// 
        /// Formula: R_F(e) = 4 - deg(i) - deg(j) + 3 * Σ_△ (w_ik * w_jk)^(1/3)
        /// </summary>
        public static double ComputeFormanRicci(RQGraph graph, int i, int j)
        {
            if (!graph.Edges[i, j])
                return 0.0;
                
            // Get degrees
            int degI = 0, degJ = 0;
            foreach (int _ in graph.Neighbors(i)) degI++;
            foreach (int _ in graph.Neighbors(j)) degJ++;
            
            // Base curvature minus degree penalty
            double curvature = 4.0 - degI - degJ;
            
            // Triangle contribution: for each common neighbor k
            // Add 3 * (w_ik * w_jk)^(1/3)
            foreach (int k in graph.Neighbors(i))
            {
                if (k == j) continue;
                if (!graph.Edges[k, j]) continue;
                
                // k is a common neighbor: triangle i-k-j exists
                double w_ik = graph.Weights[i, k];
                double w_jk = graph.Weights[j, k];
                
                // Weight contribution from triangle
                double triangleWeight = Math.Pow(w_ik * w_jk, 1.0 / 3.0);
                curvature += 3.0 * triangleWeight;
            }
            
            return curvature;
        }
        
        /// <summary>
        /// Compute Forman-Ricci curvature with weighted edge contribution.
        /// 
        /// Extended formula including edge weight w_ij:
        ///   R_F(e) = w_e * [4 - α*(W_i + W_j - 2*w_e) + 3 * Σ_△ (w_ik * w_jk)^(1/3)]
        /// 
        /// where:
        /// - W_i = weighted degree of node i = Σ_k w_ik
        /// - α = 1/⟨k⟩ = degree penalty factor (from PhysicsConstants.DegreePenaltyFactor)
        /// 
        /// This version is from the RQ-hypothesis documentation and includes
        /// the edge weight prefactor for consistency with Regge calculus.
        /// </summary>
        public static double ComputeFormanRicciWeighted(RQGraph graph, int i, int j)
        {
            if (!graph.Edges[i, j])
                return 0.0;
                
            double w_ij = graph.Weights[i, j];
            if (w_ij <= 0)
                return 0.0;
                
            // Compute weighted degrees
            double W_i = 0.0, W_j = 0.0;
            foreach (int k in graph.Neighbors(i)) W_i += graph.Weights[i, k];
            foreach (int k in graph.Neighbors(j)) W_j += graph.Weights[j, k];
            
            // Degree penalty factor: α ≈ 1/⟨k⟩
            double alpha = PhysicsConstants.DegreePenaltyFactor;
            
            // Triangle contribution
            double triangleSum = 0.0;
            foreach (int k in graph.Neighbors(i))
            {
                if (k == j) continue;
                if (!graph.Edges[k, j]) continue;
                
                double w_ik = graph.Weights[i, k];
                double w_jk = graph.Weights[j, k];
                triangleSum += Math.Pow(w_ik * w_jk, 1.0 / 3.0);
            }
            
            // Full Forman curvature with edge weight
            double curvature = w_ij * (4.0 - alpha * (W_i + W_j - 2.0 * w_ij) + 3.0 * triangleSum);
            
            return curvature;
        }
        
        /// <summary>
        /// Compute Forman-Ricci curvature for all edges in the graph.
        /// Returns array indexed by flattened edge coordinates.
        /// </summary>
        public static double[,] ComputeAllFormanRicci(RQGraph graph)
        {
            int N = graph.N;
            var curvatures = new double[N, N];
            
            for (int i = 0; i < N; i++)
            {
                foreach (int j in graph.Neighbors(i))
                {
                    if (j <= i) continue; // Only compute once per edge
                    
                    double R = ComputeFormanRicciWeighted(graph, i, j);
                    curvatures[i, j] = R;
                    curvatures[j, i] = R; // Symmetric
                }
            }
            
            return curvatures;
        }
        
        /// <summary>
        /// Compute scalar curvature at a node (average of incident edge curvatures).
        /// 
        /// R(v) = (1/deg(v)) * Σ_e R_F(e) for edges e incident to v
        /// </summary>
        public static double ComputeScalarCurvature(RQGraph graph, int node)
        {
            double sum = 0.0;
            int count = 0;
            
            foreach (int neighbor in graph.Neighbors(node))
            {
                sum += ComputeFormanRicciWeighted(graph, node, neighbor);
                count++;
            }
            
            return count > 0 ? sum / count : 0.0;
        }
        
        /// <summary>
        /// Compute average scalar curvature over the entire graph.
        /// This is analogous to the Einstein-Hilbert action integrand.
        /// </summary>
        public static double ComputeAverageScalarCurvature(RQGraph graph)
        {
            double totalCurvature = 0.0;
            int edgeCount = 0;
            
            for (int i = 0; i < graph.N; i++)
            {
                foreach (int j in graph.Neighbors(i))
                {
                    if (j <= i) continue;
                    totalCurvature += ComputeFormanRicciWeighted(graph, i, j);
                    edgeCount++;
                }
            }
            
            return edgeCount > 0 ? totalCurvature / edgeCount : 0.0;
        }
    }
}
