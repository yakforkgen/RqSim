using System;
using ComputeSharp;

namespace RQSimulation.GPUOptimized
{
    /// <summary>
    /// GPU-accelerated scalar field diffusion engine.
    /// Implements Klein-Gordon equation: d?/dt = D??? - m??
    /// 
    /// This is a classic SpMV (sparse matrix-vector multiplication) operation
    /// where each GPU thread processes one node and its neighbors.
    /// </summary>
    public class ScalarFieldEngine : IDisposable
    {
        private readonly GraphicsDevice _device;
        private ReadWriteBuffer<float>? _fieldBuffer;
        private ReadWriteBuffer<float>? _deltaBuffer;

        // Topology in CSR format (for fast neighbor access)
        private ReadOnlyBuffer<int>? _adjOffsets;
        private ReadOnlyBuffer<int>? _adjNeighbors;
        private ReadOnlyBuffer<float>? _adjWeights;

        private int _nodeCount;
        private bool _topologyInitialized;

        public ScalarFieldEngine()
        {
            _device = GraphicsDevice.GetDefault();
        }

        /// <summary>
        /// Initialize buffers for a graph of given size.
        /// Call this before UpdateTopology and UpdateField.
        /// </summary>
        /// <param name="nodeCount">Number of nodes in the graph</param>
        /// <param name="totalEdges">Total number of directed edges (2 * undirected edges)</param>
        public void Initialize(int nodeCount, int totalEdges)
        {
            _nodeCount = nodeCount;

            // Dispose old buffers if reinitializing
            _fieldBuffer?.Dispose();
            _deltaBuffer?.Dispose();
            _adjOffsets?.Dispose();
            _adjNeighbors?.Dispose();
            _adjWeights?.Dispose();

            _fieldBuffer = _device.AllocateReadWriteBuffer<float>(nodeCount);
            _deltaBuffer = _device.AllocateReadWriteBuffer<float>(nodeCount);
            _adjOffsets = _device.AllocateReadOnlyBuffer<int>(nodeCount + 1);
            _adjNeighbors = _device.AllocateReadOnlyBuffer<int>(totalEdges);
            _adjWeights = _device.AllocateReadOnlyBuffer<float>(totalEdges);
        }

        /// <summary>
        /// Update topology buffers (CSR format).
        /// Call this when graph structure changes.
        /// </summary>
        /// <param name="offsets">CSR row offsets array (length = nodeCount + 1)</param>
        /// <param name="neighbors">CSR column indices (neighbor node IDs)</param>
        /// <param name="weights">Edge weights corresponding to neighbors</param>
        public void UpdateTopology(int[] offsets, int[] neighbors, float[] weights)
        {
            if (_adjOffsets == null || _adjNeighbors == null || _adjWeights == null)
            {
                throw new InvalidOperationException(
                    "Buffers not initialized. Call Initialize() first.");
            }

            _adjOffsets.CopyFrom(offsets);
            _adjNeighbors.CopyFrom(neighbors);
            _adjWeights.CopyFrom(weights);
            _topologyInitialized = true;
        }

        /// <summary>
        /// Evolve scalar field by one timestep using GPU.
        /// Implements Klein-Gordon: d? = dt * (D??? - m??)
        /// </summary>
        /// <param name="hostField">Field values at each node (updated in-place)</param>
        /// <param name="dt">Time step</param>
        /// <param name="diffusionRate">Diffusion coefficient D</param>
        /// <param name="mass">Scalar field mass m</param>
        public void UpdateField(float[] hostField, float dt, float diffusionRate, float mass)
        {
            if (!_topologyInitialized || _fieldBuffer == null || _deltaBuffer == null ||
                _adjOffsets == null || _adjNeighbors == null || _adjWeights == null)
            {
                throw new InvalidOperationException(
                    "Engine not properly initialized. Call Initialize() and UpdateTopology() first.");
            }

            if (hostField.Length != _nodeCount)
            {
                throw new ArgumentException(
                    $"Field array length ({hostField.Length}) doesn't match node count ({_nodeCount})");
            }

            // Upload field to GPU
            _fieldBuffer.CopyFrom(hostField);

            // Step 1: Compute Laplacian and delta (change in field)
            var laplacianShader = new ScalarLaplacianShader(
                _fieldBuffer,
                _deltaBuffer,
                _adjOffsets,
                _adjNeighbors,
                _adjWeights,
                dt,
                diffusionRate,
                mass,
                _nodeCount);

            _device.For(_nodeCount, laplacianShader);

            // Step 2: Apply delta (Euler integration step)
            var applyShader = new ApplyScalarDeltaShader(_fieldBuffer, _deltaBuffer, _nodeCount);
            _device.For(_nodeCount, applyShader);

            // Download results
            _fieldBuffer.CopyTo(hostField);
        }

        /// <summary>
        /// Evolve field without copying back to CPU.
        /// Use SyncToHost() to retrieve results when needed.
        /// </summary>
        public void UpdateFieldNoCopy(float dt, float diffusionRate, float mass)
        {
            if (!_topologyInitialized || _fieldBuffer == null || _deltaBuffer == null ||
                _adjOffsets == null || _adjNeighbors == null || _adjWeights == null)
            {
                throw new InvalidOperationException(
                    "Engine not properly initialized. Call Initialize() and UpdateTopology() first.");
            }

            var laplacianShader = new ScalarLaplacianShader(
                _fieldBuffer,
                _deltaBuffer,
                _adjOffsets,
                _adjNeighbors,
                _adjWeights,
                dt,
                diffusionRate,
                mass,
                _nodeCount);

            _device.For(_nodeCount, laplacianShader);

            var applyShader = new ApplyScalarDeltaShader(_fieldBuffer, _deltaBuffer, _nodeCount);
            _device.For(_nodeCount, applyShader);
        }

        /// <summary>
        /// Upload initial field data to GPU.
        /// Use with UpdateFieldNoCopy for batch processing.
        /// </summary>
        public void UploadField(float[] hostField)
        {
            if (_fieldBuffer == null)
            {
                throw new InvalidOperationException("Buffers not initialized.");
            }

            _fieldBuffer.CopyFrom(hostField);
        }

        /// <summary>
        /// Download field data from GPU.
        /// </summary>
        public void SyncToHost(float[] hostField)
        {
            if (_fieldBuffer == null)
            {
                throw new InvalidOperationException("Buffers not initialized.");
            }

            _fieldBuffer.CopyTo(hostField);
        }

        public void Dispose()
        {
            _fieldBuffer?.Dispose();
            _deltaBuffer?.Dispose();
            _adjOffsets?.Dispose();
            _adjNeighbors?.Dispose();
            _adjWeights?.Dispose();
        }
    }

    /// <summary>
    /// GPU shader for computing discrete Laplacian and Klein-Gordon update.
    /// Each thread processes one node.
    /// 
    /// Discrete Laplacian: (???)_i = ?_j w_ij * (?_j - ?_i)
    /// Klein-Gordon: d? = dt * (D * ??? - m? * ?)
    /// </summary>
    [ThreadGroupSize(64, 1, 1)]
    [GeneratedComputeShaderDescriptor]
    public readonly partial struct ScalarLaplacianShader : IComputeShader
    {
        public readonly ReadWriteBuffer<float> field;
        public readonly ReadWriteBuffer<float> delta;
        public readonly ReadOnlyBuffer<int> offsets;
        public readonly ReadOnlyBuffer<int> neighbors;
        public readonly ReadOnlyBuffer<float> weights;
        public readonly float dt;
        public readonly float diffusion;
        public readonly float mass;
        public readonly int nodeCount;

        public ScalarLaplacianShader(
            ReadWriteBuffer<float> field,
            ReadWriteBuffer<float> delta,
            ReadOnlyBuffer<int> offsets,
            ReadOnlyBuffer<int> neighbors,
            ReadOnlyBuffer<float> weights,
            float dt,
            float diffusion,
            float mass,
            int nodeCount)
        {
            this.field = field;
            this.delta = delta;
            this.offsets = offsets;
            this.neighbors = neighbors;
            this.weights = weights;
            this.dt = dt;
            this.diffusion = diffusion;
            this.mass = mass;
            this.nodeCount = nodeCount;
        }

        public void Execute()
        {
            int i = ThreadIds.X;
            if (i >= nodeCount) return;

            float phi_i = field[i];
            float laplacian = 0.0f;

            // Iterate over neighbors via CSR format
            int start = offsets[i];
            int end = offsets[i + 1];

            for (int k = start; k < end; k++)
            {
                int neighborIdx = neighbors[k];
                float w = weights[k];
                float phi_j = field[neighborIdx];

                // Discrete Laplacian: ? w_ij * (?_j - ?_i)
                laplacian += w * (phi_j - phi_i);
            }

            // Klein-Gordon equation: d?/dt = D??? - m??
            float dPhi = diffusion * laplacian - (mass * mass * phi_i);
            delta[i] = dPhi * dt;
        }
    }

    /// <summary>
    /// GPU shader for applying delta to field: ?_new = ? + ??
    /// </summary>
    [ThreadGroupSize(64, 1, 1)]
    [GeneratedComputeShaderDescriptor]
    public readonly partial struct ApplyScalarDeltaShader : IComputeShader
    {
        public readonly ReadWriteBuffer<float> field;
        public readonly ReadWriteBuffer<float> delta;
        public readonly int nodeCount;

        public ApplyScalarDeltaShader(
            ReadWriteBuffer<float> field,
            ReadWriteBuffer<float> delta,
            int nodeCount)
        {
            this.field = field;
            this.delta = delta;
            this.nodeCount = nodeCount;
        }

        public void Execute()
        {
            int i = ThreadIds.X;
            if (i < nodeCount)
            {
                field[i] += delta[i];
            }
        }
    }
}
