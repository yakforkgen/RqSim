using System;
using ComputeSharp;

namespace RQSimulation.GPUOptimized
{
    /// <summary>
    /// GPU-accelerated statistics engine for fast aggregation operations.
    /// Implements parallel reduction algorithms for sum, min, max, and histogram.
    /// 
    /// Key operations:
    /// - Energy summation (total energy of all nodes/edges)
    /// - Mass aggregation (total mass, mass distribution)
    /// - Weight statistics (min, max, mean, variance)
    /// 
    /// Note: GPU atomics only support int, so we use scaled integers for accumulation
    /// and convert back to float on CPU.
    /// </summary>
    public class StatisticsEngine : IDisposable
    {
        private readonly GraphicsDevice _device;

        // Input/output buffers
        private ReadOnlyBuffer<float>? _inputBuffer;
        private ReadWriteBuffer<int>? _partialSumsBuffer;  // Int for atomics
        private ReadWriteBuffer<int>? _outputBuffer;       // Int for atomics

        // Histogram buffers
        private ReadWriteBuffer<int>? _histogramBuffer;

        private int _maxInputSize;
        private int _blockCount;
        private const int BlockSize = 256;
        private const float Scale = 1000000.0f; // Scale factor for float->int conversion

        public StatisticsEngine()
        {
            _device = GraphicsDevice.GetDefault();
        }

        /// <summary>
        /// Initialize buffers for a maximum input size.
        /// </summary>
        public void Initialize(int maxInputSize)
        {
            _maxInputSize = maxInputSize;
            _blockCount = (maxInputSize + BlockSize - 1) / BlockSize;

            _inputBuffer?.Dispose();
            _partialSumsBuffer?.Dispose();
            _outputBuffer?.Dispose();
            _histogramBuffer?.Dispose();

            _inputBuffer = _device.AllocateReadOnlyBuffer<float>(maxInputSize);
            _partialSumsBuffer = _device.AllocateReadWriteBuffer<int>(_blockCount);
            _outputBuffer = _device.AllocateReadWriteBuffer<int>(4); // For [sum, min, max, count]
            _histogramBuffer = _device.AllocateReadWriteBuffer<int>(256); // 256-bin histogram
        }

        /// <summary>
        /// Compute sum of all values using parallel reduction.
        /// </summary>
        public float Sum(float[] values)
        {
            if (_inputBuffer == null || _partialSumsBuffer == null)
            {
                throw new InvalidOperationException("Engine not initialized.");
            }

            int N = values.Length;
            if (N > _maxInputSize)
            {
                throw new ArgumentException($"Input size ({N}) exceeds maximum ({_maxInputSize})");
            }

            // Recreate input buffer with actual data
            _inputBuffer.Dispose();
            _inputBuffer = _device.AllocateReadOnlyBuffer(values);

            // First pass: block-level reduction
            int currentBlockCount = (N + BlockSize - 1) / BlockSize;
            
            // Clear partial sums
            int[] zeros = new int[currentBlockCount];
            _partialSumsBuffer.CopyFrom(zeros);
            
            var sumShader = new BlockSumShader(_inputBuffer, _partialSumsBuffer, N, Scale);
            _device.For(N, sumShader);

            // Second pass: final reduction on CPU
            int[] partialSums = new int[currentBlockCount];
            _partialSumsBuffer.CopyTo(partialSums);

            long total = 0;
            for (int i = 0; i < currentBlockCount; i++)
            {
                total += partialSums[i];
            }

            return (float)(total / Scale);
        }

        /// <summary>
        /// Compute sum, count, and approximate min/max.
        /// Note: Min/max use integer approximation.
        /// </summary>
        public (float sum, int count) ComputeSumAndCount(float[] values)
        {
            if (_inputBuffer == null || _outputBuffer == null)
            {
                throw new InvalidOperationException("Engine not initialized.");
            }

            int N = values.Length;
            if (N == 0) return (0, 0);

            _inputBuffer.Dispose();
            _inputBuffer = _device.AllocateReadOnlyBuffer(values);

            // Initialize output: [sum=0, unused, unused, count=0]
            int[] initOutput = [0, 0, 0, 0];
            _outputBuffer.CopyFrom(initOutput);

            // Run statistics shader
            var statsShader = new SumCountShader(_inputBuffer, _outputBuffer, N, Scale);
            _device.For(N, statsShader);

            // Read results
            int[] result = new int[4];
            _outputBuffer.CopyTo(result);

            return ((float)(result[0] / Scale), result[3]);
        }

        /// <summary>
        /// Compute weighted sum: ? values[i] * weights[i]
        /// </summary>
        public float WeightedSum(float[] values, float[] weights)
        {
            if (values.Length != weights.Length)
            {
                throw new ArgumentException("Values and weights must have same length");
            }

            // Compute products on CPU, sum on GPU
            float[] products = new float[values.Length];
            for (int i = 0; i < values.Length; i++)
            {
                products[i] = values[i] * weights[i];
            }

            return Sum(products);
        }

        /// <summary>
        /// Compute mean and variance (CPU fallback for variance due to float atomics limitation).
        /// </summary>
        public (float mean, float variance) ComputeMeanVariance(float[] values)
        {
            var (sum, count) = ComputeSumAndCount(values);
            if (count == 0) return (0, 0);

            float mean = sum / count;

            // Second pass for variance on CPU (GPU float atomics not available)
            float sumSquares = 0;
            for (int i = 0; i < values.Length; i++)
            {
                float d = values[i] - mean;
                sumSquares += d * d;
            }

            float variance = sumSquares / count;

            return (mean, variance);
        }

        /// <summary>
        /// Compute histogram of values in [0, 1] range.
        /// Returns 256-bin histogram.
        /// </summary>
        public int[] ComputeHistogram(float[] values)
        {
            if (_inputBuffer == null || _histogramBuffer == null)
            {
                throw new InvalidOperationException("Engine not initialized.");
            }

            int N = values.Length;
            
            _inputBuffer.Dispose();
            _inputBuffer = _device.AllocateReadOnlyBuffer(values);

            // Clear histogram
            int[] zeros = new int[256];
            _histogramBuffer.CopyFrom(zeros);

            // Compute histogram
            var histShader = new HistogramShader(_inputBuffer, _histogramBuffer, N);
            _device.For(N, histShader);

            // Read results
            int[] histogram = new int[256];
            _histogramBuffer.CopyTo(histogram);

            return histogram;
        }

        /// <summary>
        /// Compute total energy from edge weights and curvatures.
        /// E = ? w_ij * (1 + ?_ij)
        /// </summary>
        public float ComputeTotalEnergy(float[] weights, float[] curvatures)
        {
            if (weights.Length != curvatures.Length)
            {
                throw new ArgumentException("Weights and curvatures must have same length");
            }

            float[] energyTerms = new float[weights.Length];
            for (int i = 0; i < weights.Length; i++)
            {
                energyTerms[i] = weights[i] * (1.0f + curvatures[i]);
            }

            return Sum(energyTerms);
        }

        /// <summary>
        /// Count how many values are above/below thresholds.
        /// Returns (belowLow, inRange, aboveHigh)
        /// </summary>
        public (int belowLow, int inRange, int aboveHigh) CountThresholds(
            float[] values, float lowThreshold, float highThreshold)
        {
            if (_inputBuffer == null || _outputBuffer == null)
            {
                throw new InvalidOperationException("Engine not initialized.");
            }

            int N = values.Length;
            
            _inputBuffer.Dispose();
            _inputBuffer = _device.AllocateReadOnlyBuffer(values);

            // Output: [belowCount, inRangeCount, aboveCount, unused]
            int[] initOutput = [0, 0, 0, 0];
            _outputBuffer.CopyFrom(initOutput);

            var countShader = new ThresholdCountShader(
                _inputBuffer, _outputBuffer, lowThreshold, highThreshold, N);
            _device.For(N, countShader);

            int[] result = new int[4];
            _outputBuffer.CopyTo(result);

            return (result[0], result[1], result[2]);
        }

        public void Dispose()
        {
            _inputBuffer?.Dispose();
            _partialSumsBuffer?.Dispose();
            _outputBuffer?.Dispose();
            _histogramBuffer?.Dispose();
        }
    }

    /// <summary>
    /// GPU shader for block-level sum reduction using atomic int operations.
    /// Converts float to scaled int for accumulation.
    /// </summary>
    [ThreadGroupSize(64, 1, 1)]
    [GeneratedComputeShaderDescriptor]
    public readonly partial struct BlockSumShader : IComputeShader
    {
        public readonly ReadOnlyBuffer<float> input;
        public readonly ReadWriteBuffer<int> partialSums;
        public readonly int inputLength;
        public readonly float scale;

        public BlockSumShader(
            ReadOnlyBuffer<float> input,
            ReadWriteBuffer<int> partialSums,
            int inputLength,
            float scale)
        {
            this.input = input;
            this.partialSums = partialSums;
            this.inputLength = inputLength;
            this.scale = scale;
        }

        public void Execute()
        {
            int gId = ThreadIds.X;
            if (gId >= inputLength) return;

            int blockId = gId / 64; // Block size is 64 for this shader

            // Load value and scale to int
            float val = input[gId];
            int scaledVal = (int)(val * scale);

            // Atomic add to block result
            Hlsl.InterlockedAdd(ref partialSums[blockId], scaledVal);
        }
    }

    /// <summary>
    /// GPU shader for computing sum and count using atomic int operations.
    /// </summary>
    [ThreadGroupSize(64, 1, 1)]
    [GeneratedComputeShaderDescriptor]
    public readonly partial struct SumCountShader : IComputeShader
    {
        public readonly ReadOnlyBuffer<float> input;
        public readonly ReadWriteBuffer<int> output; // [scaledSum, unused, unused, count]
        public readonly int inputLength;
        public readonly float scale;

        public SumCountShader(
            ReadOnlyBuffer<float> input,
            ReadWriteBuffer<int> output,
            int inputLength,
            float scale)
        {
            this.input = input;
            this.output = output;
            this.inputLength = inputLength;
            this.scale = scale;
        }

        public void Execute()
        {
            int i = ThreadIds.X;
            if (i >= inputLength) return;

            float val = input[i];
            int scaledVal = (int)(val * scale);

            // Atomic add for sum
            Hlsl.InterlockedAdd(ref output[0], scaledVal);

            // Atomic increment for count
            Hlsl.InterlockedAdd(ref output[3], 1);
        }
    }

    /// <summary>
    /// GPU shader for computing histogram of values in [0, 1].
    /// </summary>
    [ThreadGroupSize(64, 1, 1)]
    [GeneratedComputeShaderDescriptor]
    public readonly partial struct HistogramShader : IComputeShader
    {
        public readonly ReadOnlyBuffer<float> input;
        public readonly ReadWriteBuffer<int> histogram;
        public readonly int inputLength;

        public HistogramShader(
            ReadOnlyBuffer<float> input,
            ReadWriteBuffer<int> histogram,
            int inputLength)
        {
            this.input = input;
            this.histogram = histogram;
            this.inputLength = inputLength;
        }

        public void Execute()
        {
            int i = ThreadIds.X;
            if (i >= inputLength) return;

            float val = input[i];

            // Clamp to [0, 1] and compute bin
            val = Hlsl.Clamp(val, 0.0f, 0.9999f);
            int bin = (int)(val * 256.0f);

            // Atomic increment
            Hlsl.InterlockedAdd(ref histogram[bin], 1);
        }
    }

    /// <summary>
    /// GPU shader for counting values in threshold ranges.
    /// </summary>
    [ThreadGroupSize(64, 1, 1)]
    [GeneratedComputeShaderDescriptor]
    public readonly partial struct ThresholdCountShader : IComputeShader
    {
        public readonly ReadOnlyBuffer<float> input;
        public readonly ReadWriteBuffer<int> output; // [below, inRange, above, unused]
        public readonly float lowThreshold;
        public readonly float highThreshold;
        public readonly int inputLength;

        public ThresholdCountShader(
            ReadOnlyBuffer<float> input,
            ReadWriteBuffer<int> output,
            float lowThreshold,
            float highThreshold,
            int inputLength)
        {
            this.input = input;
            this.output = output;
            this.lowThreshold = lowThreshold;
            this.highThreshold = highThreshold;
            this.inputLength = inputLength;
        }

        public void Execute()
        {
            int i = ThreadIds.X;
            if (i >= inputLength) return;

            float val = input[i];

            if (val < lowThreshold)
            {
                Hlsl.InterlockedAdd(ref output[0], 1);
            }
            else if (val > highThreshold)
            {
                Hlsl.InterlockedAdd(ref output[2], 1);
            }
            else
            {
                Hlsl.InterlockedAdd(ref output[1], 1);
            }
        }
    }
}
