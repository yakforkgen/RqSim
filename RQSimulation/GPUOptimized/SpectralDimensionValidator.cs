using System;
using System.Collections.Generic;
using System.Linq;
using ComputeSharp;

namespace RQSimulation.GPUOptimized
{
    /// <summary>
    /// Spectral Dimension Validation and Spinor Compatibility Checker
    /// 
    /// CHECKLIST ITEM 5: Verify that spinor fields are compatible with actual spectral dimension
    /// 
    /// Problem: The code uses 4-component spinors (assuming 3+1 dimensions), but if the graph
    /// is fractal with d_S ≈ 2.5, the Dirac operator is incorrect.
    /// 
    /// Solution: Compute spectral dimension periodically. If d_S deviates from 4 by >10%,
    /// issue warning or suppress spinor evolution until geometry stabilizes.
    /// 
    /// ADDITIONAL VALIDATION (Action 3):
    /// If d_S = 1.000 exactly for many steps, this indicates walkers are trapped in small
    /// clusters (1D chain) and the calculation may be unreliable.
    /// </summary>
    public static class SpectralDimensionValidator
    {
        private const double TargetDimension = 4.0;
        private const double TolerancePercent = 0.10; // 10% tolerance
        private const double MinAcceptableDimension = TargetDimension * (1.0 - TolerancePercent);
        private const double MaxAcceptableDimension = TargetDimension * (1.0 + TolerancePercent);

        // Suspicious dimension detection (Action 3)
        private const double SuspiciousExactDimension = 1.0;
        private const double SuspiciousDimensionTolerance = 0.001; // If d_S = 1.000 ± 0.001
        private const int SuspiciousStabilityWindow = 100; // Number of steps to track
        private const int SuspiciousStabilityThreshold = 80; // If 80% of recent values are exactly 1.0

        // History for detecting suspiciously stable d_S
        private static readonly Queue<double> _recentDimensions = new();
        private static int _consecutiveExactOne = 0;

        /// <summary>
        /// Check if spinor fields should be active given current spectral dimension
        /// </summary>
        public static bool ShouldEnableSpinorFields(double spectralDimension)
        {
            // Spinors are only valid for d_S ≈ 4 (within 10% tolerance)
            return spectralDimension >= MinAcceptableDimension &&
                   spectralDimension <= MaxAcceptableDimension;
        }

        /// <summary>
        /// Check if spectral dimension is suspiciously stable at exactly 1.0
        /// This indicates walkers are trapped in small clusters (1D chain structure)
        /// </summary>
        public static SuspiciousDimensionStatus CheckForSuspiciousStability(double spectralDimension)
        {
            var status = new SuspiciousDimensionStatus
            {
                CurrentDimension = spectralDimension,
                IsSuspicious = false,
                ConsecutiveExactOneCount = 0,
                RecentExactOnePercent = 0.0,
                Diagnosis = string.Empty
            };

            // Check if current dimension is exactly 1.0 (within tolerance)
            bool isExactlyOne = Math.Abs(spectralDimension - SuspiciousExactDimension) < SuspiciousDimensionTolerance;

            // Track consecutive exact-1.0 values
            if (isExactlyOne)
            {
                _consecutiveExactOne++;
            }
            else
            {
                _consecutiveExactOne = 0;
            }

            // Add to history
            _recentDimensions.Enqueue(spectralDimension);
            while (_recentDimensions.Count > SuspiciousStabilityWindow)
            {
                _recentDimensions.Dequeue();
            }

            // Calculate percentage of recent values that are exactly 1.0
            int exactOneCount = _recentDimensions.Count(d => Math.Abs(d - SuspiciousExactDimension) < SuspiciousDimensionTolerance);
            double exactOnePercent = _recentDimensions.Count > 0
                ? (double)exactOneCount / _recentDimensions.Count
                : 0.0;

            status.ConsecutiveExactOneCount = _consecutiveExactOne;
            status.RecentExactOnePercent = exactOnePercent;

            // Determine if suspicious
            if (_consecutiveExactOne >= 10 || exactOnePercent >= 0.8)
            {
                status.IsSuspicious = true;
                status.Diagnosis = DiagnoseSuspiciousDimension();
            }

            return status;
        }

        /// <summary>
        /// Diagnose why spectral dimension might be stuck at 1.0
        /// </summary>
        private static string DiagnoseSuspiciousDimension()
        {
            return "[SUSPICIOUS] d_S = 1.000 is too stable - possible causes:\n" +
                   "  1. Random walkers are trapped in isolated 1D chains\n" +
                   "  2. Graph has fragmented into disconnected linear components\n" +
                   "  3. EdgeCreationBarrier may still be too high for long-range links\n" +
                   "  4. GravitationalCoupling may be cementing 1D structure\n" +
                   "RECOMMENDED ACTIONS:\n" +
                   "  - Check graph connectivity (number of connected components)\n" +
                   "  - Verify LargestCluster size is reasonable\n" +
                   "  - Consider increasing TopologyTunnelingRate further\n" +
                   "  - Ensure initial graph is not too sparse (initialEdgeProb)";
        }

        /// <summary>
        /// Reset the history tracking (call when starting new simulation)
        /// </summary>
        public static void ResetHistory()
        {
            _recentDimensions.Clear();
            _consecutiveExactOne = 0;
        }

        /// <summary>
        /// Validate spinor compatibility and return status
        /// Also checks for suspicious d_S = 1.0 stability
        /// </summary>
        public static SpinorCompatibilityStatus ValidateSpinorCompatibility(RQGraph graph)
        {
            double spectralDim = graph.ComputeSpectralDimension();

            // Check for suspicious stability at d_S = 1.0
            var suspiciousStatus = CheckForSuspiciousStability(spectralDim);

            var status = new SpinorCompatibilityStatus
            {
                SpectralDimension = spectralDim,
                IsCompatible = ShouldEnableSpinorFields(spectralDim),
                DeviationPercent = Math.Abs(spectralDim - TargetDimension) / TargetDimension,
                Recommendation = GetRecommendation(spectralDim),
                IsSuspiciouslyStable = suspiciousStatus.IsSuspicious,
                SuspiciousDiagnosis = suspiciousStatus.Diagnosis
            };

            return status;
        }

        /// <summary>
        /// Extended validation with walker spread analysis
        /// Use this for detailed diagnostics
        /// </summary>
        public static ExtendedValidationResult ValidateWithWalkerAnalysis(
            RQGraph graph,
            int numWalkers = 100,
            int walkSteps = 50)
        {
            var result = new ExtendedValidationResult();

            // Run random walks and track spread
            int N = graph.N;
            var random = new Random();
            var visitedNodes = new HashSet<int>();
            var walkerFinalPositions = new int[numWalkers];
            var walkerDistances = new double[numWalkers];

            for (int w = 0; w < numWalkers; w++)
            {
                int startNode = random.Next(N);
                int currentNode = startNode;

                for (int step = 0; step < walkSteps; step++)
                {
                    var neighbors = graph.Neighbors(currentNode).ToList();
                    if (neighbors.Count == 0) break;

                    // Random walk step
                    currentNode = neighbors[random.Next(neighbors.Count)];
                    visitedNodes.Add(currentNode);
                }

                walkerFinalPositions[w] = currentNode;
                walkerDistances[w] = graph.ShortestPathDistance(startNode, currentNode);
            }

            // Analyze walker spread
            result.UniqueNodesVisited = visitedNodes.Count;
            result.UniqueNodesFraction = (double)visitedNodes.Count / N;
            result.AverageWalkerDistance = walkerDistances.Average();
            result.MaxWalkerDistance = walkerDistances.Max();

            // Check for trapped walkers (distance = 0 or very small)
            int trappedWalkers = walkerDistances.Count(d => d <= 2);
            result.TrappedWalkerFraction = (double)trappedWalkers / numWalkers;

            // Compute spectral dimension
            result.SpectralDimension = graph.ComputeSpectralDimension(numWalkers, walkSteps);

            // Determine if walkers are locked
            result.WalkersAreLocked = result.TrappedWalkerFraction > 0.5 ||
                                       result.UniqueNodesFraction < 0.1;

            if (result.WalkersAreLocked)
            {
                result.Diagnosis = $"[LOCKED WALKERS] {result.TrappedWalkerFraction * 100:F1}% of walkers " +
                                   $"traveled ≤2 hops. Only {result.UniqueNodesFraction * 100:F1}% of nodes reached.\n" +
                                   "This invalidates spectral dimension calculation.\n" +
                                   "Graph may be disconnected or have bottleneck topology.";
            }
            else
            {
                result.Diagnosis = $"Walker spread OK: {result.UniqueNodesFraction * 100:F1}% nodes visited, " +
                                   $"avg distance {result.AverageWalkerDistance:F2}";
            }

            return result;
        }

        /// <summary>
        /// Get recommendation based on spectral dimension
        /// </summary>
        private static string GetRecommendation(double spectralDim)
        {
            if (spectralDim < MinAcceptableDimension)
            {
                return $"Spectral dimension {spectralDim:F2} < {MinAcceptableDimension:F2}. " +
                       "Graph is too fractal for 4D spinors. Suppress spinor evolution " +
                       "until geometry crystallizes to 4D.";
            }
            else if (spectralDim > MaxAcceptableDimension)
            {
                return $"Spectral dimension {spectralDim:F2} > {MaxAcceptableDimension:F2}. " +
                       "Graph has too many dimensions. Check for over-connectivity.";
            }
            else
            {
                return $"Spectral dimension {spectralDim:F2} is compatible with 4D spinors.";
            }
        }

        /// <summary>
        /// Export diagnostics with spectral dimension and spinor status
        /// Now includes suspicious stability detection
        /// </summary>
        public static void ExportDiagnostics(
            RQGraph graph,
            int step,
            List<string> diagnosticsExport)
        {
            var status = ValidateSpinorCompatibility(graph);

            string statusFlag = status.IsCompatible ? "OK" : "INCOMPATIBLE";
            if (status.IsSuspiciouslyStable)
            {
                statusFlag = "SUSPICIOUS";
            }

            string line = $"{step},{status.SpectralDimension:F4},{statusFlag}," +
                         $"{status.DeviationPercent:F4},{status.Recommendation}";

            diagnosticsExport.Add(line);

            // Log warning if incompatible (to standard console to avoid UI dependency)
            if (!status.IsCompatible)
            {
                Console.WriteLine($"[WARNING] Step {step}: {status.Recommendation}");
            }

            // Log suspicious stability warning
            if (status.IsSuspiciouslyStable)
            {
                Console.WriteLine($"[SUSPICIOUS] Step {step}: d_S = {status.SpectralDimension:F3} is suspiciously stable");
                Console.WriteLine(status.SuspiciousDiagnosis);
            }
        }

        /// <summary>
        /// Decide whether to suppress spinor field evolution
        /// </summary>
        public static bool ShouldSuppressSpinorEvolution(RQGraph graph)
        {
            double spectralDim = graph.ComputeSpectralDimension();
            return !ShouldEnableSpinorFields(spectralDim);
        }

        /// <summary>
        /// Monitor spectral dimension evolution and detect transitions
        /// </summary>
        public static SpectralTransitionInfo MonitorTransition(
            double previousDimension,
            double currentDimension)
        {
            var info = new SpectralTransitionInfo
            {
                PreviousDimension = previousDimension,
                CurrentDimension = currentDimension,
                DimensionChange = currentDimension - previousDimension,
                HasCrossedThreshold = false,
                TransitionType = TransitionType.None
            };

            // Detect transition from fractal (d < 3) to spacetime (d ≥ 4)
            if (previousDimension < 3.0 && currentDimension >= 3.8)
            {
                info.HasCrossedThreshold = true;
                info.TransitionType = TransitionType.FractalToSpacetime;
                Console.WriteLine($"[SUCCESS] Spectral dimension crossed from {previousDimension:F2} to {currentDimension:F2}!");
                Console.WriteLine("[SUCCESS] Graph has crystallized into 4D spacetime!");
            }
            // Detect transition from spacetime to fractal (collapse)
            else if (previousDimension >= 3.8 && currentDimension < 3.0)
            {
                info.HasCrossedThreshold = true;
                info.TransitionType = TransitionType.SpacetimeToFractal;
                Console.WriteLine($"[WARNING] Spectral dimension collapsed from {previousDimension:F2} to {currentDimension:F2}");
                Console.WriteLine("[WARNING] Graph has become fractal - spacetime structure lost!");
            }
            // Detect gradual increase toward 4D
            else if (currentDimension > previousDimension && currentDimension < 4.5)
            {
                info.TransitionType = TransitionType.Crystallizing;
            }
            // Detect gradual decrease (fragmentation)
            else if (currentDimension < previousDimension && currentDimension > 1.5)
            {
                info.TransitionType = TransitionType.Fragmenting;
            }

            return info;
        }
    }

    /// <summary>
    /// Status of spinor field compatibility with spectral dimension
    /// </summary>
    public class SpinorCompatibilityStatus
    {
        public double SpectralDimension { get; set; }
        public bool IsCompatible { get; set; }
        public double DeviationPercent { get; set; }
        public string Recommendation { get; set; } = string.Empty;

        /// <summary>
        /// True if d_S has been suspiciously stable at exactly 1.0
        /// </summary>
        public bool IsSuspiciouslyStable { get; set; }

        /// <summary>
        /// Diagnosis message if suspicious stability detected
        /// </summary>
        public string SuspiciousDiagnosis { get; set; } = string.Empty;
    }

    /// <summary>
    /// Status of suspicious dimension detection
    /// </summary>
    public class SuspiciousDimensionStatus
    {
        public double CurrentDimension { get; set; }
        public bool IsSuspicious { get; set; }
        public int ConsecutiveExactOneCount { get; set; }
        public double RecentExactOnePercent { get; set; }
        public string Diagnosis { get; set; } = string.Empty;
    }

    /// <summary>
    /// Extended validation result with walker spread analysis
    /// </summary>
    public class ExtendedValidationResult
    {
        public double SpectralDimension { get; set; }
        public int UniqueNodesVisited { get; set; }
        public double UniqueNodesFraction { get; set; }
        public double AverageWalkerDistance { get; set; }
        public double MaxWalkerDistance { get; set; }
        public double TrappedWalkerFraction { get; set; }
        public bool WalkersAreLocked { get; set; }
        public string Diagnosis { get; set; } = string.Empty;
    }

    /// <summary>
    /// Information about spectral dimension transitions
    /// </summary>
    public class SpectralTransitionInfo
    {
        public double PreviousDimension { get; set; }
        public double CurrentDimension { get; set; }
        public double DimensionChange { get; set; }
        public bool HasCrossedThreshold { get; set; }
        public TransitionType TransitionType { get; set; }
    }

    /// <summary>
    /// Types of spectral dimension transitions
    /// </summary>
    public enum TransitionType
    {
        None,
        FractalToSpacetime,    // d: 2 → 4 (SUCCESS: crystallization)
        SpacetimeToFractal,    // d: 4 → 2 (FAILURE: collapse)
        Crystallizing,         // d increasing toward 4
        Fragmenting           // d decreasing
    }

    /// <summary>
    /// GPU-accelerated spectral dimension computation using heat kernel trace
    /// 
    /// Spectral dimension d_S is computed via:
    /// d_S = -2 * d(log P(t)) / d(log t)
    /// 
    /// where P(t) = Tr(e^{-tL}) is the heat kernel trace and L is the graph Laplacian
    /// 
    /// FIX: Added topology version tracking to prevent stale data issues.
    /// FIX: Added CPU comparison mode for debugging CPU vs GPU discrepancies.
    /// </summary>
    public class GpuSpectralEngine : IDisposable
    {
        private readonly GraphicsDevice _device;
        private ReadWriteBuffer<float>? _heatBuffer;      // Heat distribution
        private ReadWriteBuffer<float>? _newHeatBuffer;   // New heat after diffusion step
        private ReadOnlyBuffer<float>? _weightsBuffer;    // Graph weights (CSR order)
        private ReadOnlyBuffer<int>? _neighborOffsets;    // CSR format offsets
        private ReadOnlyBuffer<int>? _neighborIndices;
        private ReadOnlyBuffer<float>? _invSqrtDegrees;   // 1/sqrt(degree) per node (normalized Laplacian)

        // Topology version tracking to detect stale data
        private int _topologyVersion = -1;
        private int _nodeCount;
        private int _edgeCount;
        private bool _isInitialized;

        // Cached topology data for CPU comparison
        private float[]? _cachedWeights;
        private int[]? _cachedOffsets;
        private int[]? _cachedIndices;
        private float[]? _cachedInvSqrtDegrees;

        /// <summary>
        /// Current topology version. Compare with graph.TopologyVersion to detect staleness.
        /// </summary>
        public int TopologyVersion => _topologyVersion;

        /// <summary>
        /// Whether the engine has valid topology data loaded.
        /// </summary>
        public bool IsInitialized => _isInitialized;

        public GpuSpectralEngine()
        {
            _device = GraphicsDevice.GetDefault();
        }

        /// <summary>
        /// Update topology from graph. Call this when graph.TopologyVersion changes.
        /// </summary>
        /// <param name="graph">The RQGraph to read topology from</param>
        public void UpdateTopology(RQGraph graph)
        {
            ArgumentNullException.ThrowIfNull(graph);

            // Build CSR format from graph
            graph.BuildSoAViews();

            _nodeCount = graph.N;
            _cachedOffsets = graph.CsrOffsets;
            _cachedIndices = graph.CsrIndices;
            _edgeCount = _cachedIndices.Length;

            // Build weights array matching CSR order
            _cachedWeights = new float[_edgeCount];
            for (int i = 0; i < _nodeCount; i++)
            {
                int start = _cachedOffsets[i];
                int end = _cachedOffsets[i + 1];
                for (int k = start; k < end; k++)
                {
                    int j = _cachedIndices[k];
                    _cachedWeights[k] = (float)graph.Weights[i, j];
                }
            }

            // Precompute inverse square roots of weighted degrees
            _cachedInvSqrtDegrees = new float[_nodeCount];
            for (int node = 0; node < _nodeCount; node++)
            {
                int start = _cachedOffsets[node];
                int end = _cachedOffsets[node + 1];
                double degree = 0.0;
                for (int idx = start; idx < end; idx++)
                {
                    degree += _cachedWeights[idx];
                }

                _cachedInvSqrtDegrees[node] = degree > 1e-10
                    ? (float)(1.0 / Math.Sqrt(degree))
                    : 0.0f;
            }

            // Upload to GPU
            _weightsBuffer?.Dispose();
            _neighborOffsets?.Dispose();
            _neighborIndices?.Dispose();
            _invSqrtDegrees?.Dispose();
            _heatBuffer?.Dispose();
            _newHeatBuffer?.Dispose();

            _weightsBuffer = _device.AllocateReadOnlyBuffer(_cachedWeights);
            _neighborOffsets = _device.AllocateReadOnlyBuffer(_cachedOffsets);
            _neighborIndices = _device.AllocateReadOnlyBuffer(_cachedIndices);
            _invSqrtDegrees = _device.AllocateReadOnlyBuffer(_cachedInvSqrtDegrees);
            _heatBuffer = _device.AllocateReadWriteBuffer<float>(_nodeCount);
            _newHeatBuffer = _device.AllocateReadWriteBuffer<float>(_nodeCount);

            _topologyVersion = graph.TopologyVersion;
            _isInitialized = true;
        }

        /// <summary>
        /// Compute spectral dimension using GPU-accelerated heat kernel method.
        /// 
        /// IMPORTANT: Call UpdateTopology() first if graph has changed since last call.
        /// </summary>
        /// <param name="graph">The graph (for topology version check)</param>
        /// <param name="dt">Time step for heat diffusion (default 0.01)</param>
        /// <param name="numSteps">Number of diffusion steps (default 100)</param>
        /// <param name="numProbeVectors">Number of Rademacher vectors for trace estimation (default 8)</param>
        /// <param name="enableCpuComparison">If true, also compute CPU version for debugging</param>
        /// <returns>Spectral dimension estimate</returns>
        public float ComputeSpectralDimension(
            RQGraph graph,
            float dt = 0.01f,
            int numSteps = 100,
            int numProbeVectors = 8,
            bool enableCpuComparison = false)
        {
            ArgumentNullException.ThrowIfNull(graph);

            // Check for stale topology
            if (!_isInitialized || _topologyVersion != graph.TopologyVersion)
            {
                Console.WriteLine($"[GPU d_S] Topology version mismatch: cached={_topologyVersion}, current={graph.TopologyVersion}. Updating...");
                UpdateTopology(graph);
            }

            // Run GPU computation
            float gpuResult = ComputeSpectralDimensionGpuInternal(dt, numSteps, numProbeVectors);

            // Optional CPU comparison for debugging
            if (enableCpuComparison)
            {
                double cpuResult = ComputeSpectralDimensionCpuReference(dt, numSteps, numProbeVectors);
                double diff = Math.Abs(gpuResult - cpuResult);
                Console.WriteLine($"[d_S COMPARE] GPU={gpuResult:F4}, CPU={cpuResult:F4}, diff={diff:F4}");
                if (diff > 0.5)
                {
                    Console.WriteLine($"[d_S WARNING] Large CPU/GPU discrepancy detected! Check topology sync.");
                }
            }

            return gpuResult;
        }

        /// <summary>
        /// Compute spectral dimension using GPU-accelerated heat kernel method
        /// </summary>
        public float ComputeSpectralDimensionGpu(
            float[] flatWeights,
            int[] neighborOffsets,
            int[] neighborIndices,
            int nodeCount,
            float dt = 0.01f,
            int numSteps = 100,
            int numProbeVectors = 8)
        {
            int N = nodeCount;

            // Cache the input data
            _cachedWeights = flatWeights;
            _cachedOffsets = neighborOffsets;
            _cachedIndices = neighborIndices;
            _nodeCount = N;

            // Allocate buffers
            _heatBuffer?.Dispose();
            _newHeatBuffer?.Dispose();
            _weightsBuffer?.Dispose();
            _neighborOffsets?.Dispose();
            _neighborIndices?.Dispose();
            _invSqrtDegrees?.Dispose();

            _heatBuffer = _device.AllocateReadWriteBuffer<float>(N);
            _newHeatBuffer = _device.AllocateReadWriteBuffer<float>(N);
            _weightsBuffer = _device.AllocateReadOnlyBuffer(flatWeights);
            _neighborOffsets = _device.AllocateReadOnlyBuffer(neighborOffsets);
            _neighborIndices = _device.AllocateReadOnlyBuffer(neighborIndices);

            // Precompute inverse square roots of degrees to match CPU normalized Laplacian
            _cachedInvSqrtDegrees = new float[nodeCount];
            for (int node = 0; node < nodeCount; node++)
            {
                int start = neighborOffsets[node];
                int end = neighborOffsets[node + 1];
                double degree = 0.0;
                for (int idx = start; idx < end; idx++)
                {
                    degree += flatWeights[idx];
                }

                _cachedInvSqrtDegrees[node] = degree > 1e-10
                    ? (float)(1.0 / Math.Sqrt(degree))
                    : 0.0f;
            }

            _invSqrtDegrees = _device.AllocateReadOnlyBuffer(_cachedInvSqrtDegrees);
            _isInitialized = true;

            return ComputeSpectralDimensionGpuInternal(dt, numSteps, numProbeVectors);
        }

        /// <summary>
        /// Internal GPU computation using pre-uploaded buffers.
        /// </summary>
        private float ComputeSpectralDimensionGpuInternal(float dt, int numSteps, int numProbeVectors)
        {
            if (!_isInitialized || _heatBuffer == null || _newHeatBuffer == null ||
                _weightsBuffer == null || _neighborOffsets == null ||
                _neighborIndices == null || _invSqrtDegrees == null)
            {
                throw new InvalidOperationException("GPU buffers not initialized. Call UpdateTopology or ComputeSpectralDimensionGpu first.");
            }

            int N = _nodeCount;
            numProbeVectors = Math.Max(1, numProbeVectors);

            var random = new Random(42); // Fixed seed for reproducibility
            var heatScratch = new float[N];
            var probeVector = new float[N];
            var traceSums = new double[numSteps];

            for (int probe = 0; probe < numProbeVectors; probe++)
            {
                // Generate Rademacher vector (±1) and copy to GPU
                for (int i = 0; i < N; i++)
                {
                    probeVector[i] = random.NextDouble() < 0.5 ? -1f : 1f;
                }

                _heatBuffer.CopyFrom(probeVector);

                // Diffuse and accumulate z^T e^{-tL} z estimates
                for (int step = 0; step < numSteps; step++)
                {
                    var diffusionShader = new HeatDiffusionShader(
                        _newHeatBuffer,
                        _heatBuffer,
                        _weightsBuffer,
                        _neighborOffsets,
                        _neighborIndices,
                        _invSqrtDegrees,
                        N,
                        dt);

                    _device.For(N, diffusionShader);
                    (_heatBuffer, _newHeatBuffer) = (_newHeatBuffer, _heatBuffer);

                    _heatBuffer.CopyTo(heatScratch);
                    double dot = 0.0;
                    for (int i = 0; i < N; i++)
                    {
                        dot += probeVector[i] * heatScratch[i];
                    }

                    traceSums[step] += dot;
                }
            }

            // Fit log-log slope
            return FitSpectralDimension(traceSums, numProbeVectors, dt, numSteps, "GPU");
        }

        /// <summary>
        /// CPU reference implementation for debugging.
        /// Uses the same algorithm as GPU but runs on CPU.
        /// </summary>
        private double ComputeSpectralDimensionCpuReference(float dt, int numSteps, int numProbeVectors)
        {
            if (_cachedWeights == null || _cachedOffsets == null ||
                _cachedIndices == null || _cachedInvSqrtDegrees == null)
            {
                return 4.0; // Default if not initialized
            }

            int N = _nodeCount;
            var random = new Random(42); // Same seed as GPU
            var heat = new double[N];
            var newHeat = new double[N];
            var probeVector = new double[N];
            var traceSums = new double[numSteps];

            for (int probe = 0; probe < numProbeVectors; probe++)
            {
                // Generate same Rademacher vector as GPU
                for (int i = 0; i < N; i++)
                {
                    probeVector[i] = random.NextDouble() < 0.5 ? -1.0 : 1.0;
                    heat[i] = probeVector[i];
                }

                // Diffuse using CPU
                for (int step = 0; step < numSteps; step++)
                {
                    for (int i = 0; i < N; i++)
                    {
                        int start = _cachedOffsets[i];
                        int end = _cachedOffsets[i + 1];
                        double invSqrtDegI = _cachedInvSqrtDegrees[i];

                        if (invSqrtDegI <= 1e-10 || start == end)
                        {
                            newHeat[i] = heat[i];
                            continue;
                        }

                        // Compute off-diagonal sum (same formula as GPU shader)
                        double offDiagonalSum = 0.0;
                        for (int k = start; k < end; k++)
                        {
                            int j = _cachedIndices[k];
                            double invSqrtDegJ = _cachedInvSqrtDegrees[j];
                            if (invSqrtDegJ <= 1e-10) continue;

                            double weightTerm = _cachedWeights[k] * invSqrtDegI * invSqrtDegJ;
                            offDiagonalSum += weightTerm * heat[j];
                        }

                        // Same formula: h_new = (1-dt)*h + dt*sum
                        newHeat[i] = (1.0 - dt) * heat[i] + dt * offDiagonalSum;
                        newHeat[i] = Math.Max(newHeat[i], 0.0);
                    }

                    // Swap buffers
                    (heat, newHeat) = (newHeat, heat);

                    // Accumulate trace
                    double dot = 0.0;
                    for (int i = 0; i < N; i++)
                    {
                        dot += probeVector[i] * heat[i];
                    }
                    traceSums[step] += dot;
                }
            }

            return FitSpectralDimension(traceSums, numProbeVectors, dt, numSteps, "CPU");
        }

        /// <summary>
        /// Fit spectral dimension from trace sums using log-log linear regression.
        /// </summary>
        private float FitSpectralDimension(double[] traceSums, int numProbeVectors, float dt, int numSteps, string label)
        {
            // Collect averaged heat kernel traces at different times
            List<(float t, float trace)> tracePoints = new();
            float currentTime = dt;
            for (int step = 0; step < numSteps; step++)
            {
                double avgTrace = traceSums[step] / numProbeVectors;
                double traceMagnitude = Math.Abs(avgTrace);
                if (traceMagnitude > 1e-10)
                {
                    tracePoints.Add((currentTime, (float)traceMagnitude));
                }

                currentTime += dt;
            }

            // Compute spectral dimension from log-log slope
            // d_S = -2 * d(log P) / d(log t)
            if (tracePoints.Count < 10)
            {
                Console.WriteLine($"[{label} d_S] Insufficient data points: {tracePoints.Count}");
                return 4.0f; // Default to 4D if insufficient data
            }

            // Use linear regression on log-log scale
            // Skip early points (thermalization) and late points (noise)
            int skipEarly = Math.Min(5, tracePoints.Count / 10);
            int skipLate = Math.Min(10, tracePoints.Count / 5);

            float sumLogT = 0f, sumLogP = 0f, sumLogT2 = 0f, sumLogTLogP = 0f;
            int count = 0;

            for (int i = skipEarly; i < tracePoints.Count - skipLate; i++)
            {
                var (time, trace) = tracePoints[i];
                if (trace > 0 && time > 0)
                {
                    float logT = MathF.Log(time);
                    float logP = MathF.Log(trace);
                    sumLogT += logT;
                    sumLogP += logP;
                    sumLogT2 += logT * logT;
                    sumLogTLogP += logT * logP;
                    count++;
                }
            }

            if (count < 2)
            {
                Console.WriteLine($"[{label} d_S] Not enough valid points after filtering: {count}");
                return 4.0f;
            }

            // Linear regression: log P = a + b * log t
            // Slope b = d(log P) / d(log t)
            float denominator = count * sumLogT2 - sumLogT * sumLogT;
            if (Math.Abs(denominator) < 1e-10)
            {
                Console.WriteLine($"[{label} d_S] Degenerate regression (zero denominator)");
                return 4.0f;
            }

            float slope = (count * sumLogTLogP - sumLogT * sumLogP) / denominator;

            // Log raw fit parameters for CPU vs GPU debugging
            Console.WriteLine(
                $"[{label} d_S] HeatKernelFit: count={count}, dt={dt}, steps={numSteps}, " +
                $"sumLogT={sumLogT:F6}, sumLogP={sumLogP:F6}, slope={slope:F6}");

            // If slope is non-negative, heat trace is not decaying as expected.
            // This indicates either a disconnected graph or incorrect parameters.
            if (slope >= -1e-4f)
            {
                Console.WriteLine($"[{label} d_S] Non-decaying heat trace detected (slope={slope:F6} >= 0). Check graph connectivity.");
                // Instead of returning 1.0, return NaN to signal invalid computation
                return float.NaN;
            }

            // Spectral dimension: d_S = -2 * slope
            float spectralDim = -2.0f * slope;

            Console.WriteLine($"[{label} d_S] RawSpectralDim={spectralDim:F6}");

            // Clamp to reasonable range
            float clamped = Math.Clamp(spectralDim, 1.0f, 10.0f);
            if (Math.Abs(clamped - spectralDim) > 1e-4f)
            {
                Console.WriteLine($"[{label} d_S] ClampedSpectralDim={clamped:F6} (was {spectralDim:F6})");
            }

            return clamped;
        }

        public void Dispose()
        {
            _heatBuffer?.Dispose();
            _newHeatBuffer?.Dispose();
            _weightsBuffer?.Dispose();
            _neighborOffsets?.Dispose();
            _neighborIndices?.Dispose();
            _invSqrtDegrees?.Dispose();
        }
    }

    /// <summary>
    /// GPU shader for heat diffusion step using normalized Laplacian.
    /// 
    /// Heat kernel evolution: h(t+dt) = (I - dt * L_norm) * h(t)
    /// 
    /// Normalized Laplacian L_norm = I - D^{-1/2} A D^{-1/2} where:
    /// - L_norm_ii = 1 (diagonal)
    /// - L_norm_ij = -w_ij / sqrt(deg_i * deg_j) (off-diagonal)
    /// 
    /// Therefore: (L_norm * h)_i = h_i - sum_j [w_ij * invSqrtDeg_i * invSqrtDeg_j * h_j]
    /// 
    /// And: h_new_i = h_i - dt * (L_norm * h)_i
    ///              = h_i - dt * h_i + dt * sum_j [...]
    ///              = (1 - dt) * h_i + dt * sum_j [w_ij * D_i^{-1/2} * D_j^{-1/2} * h_j]
    /// 
    /// FIX: This matches the CPU ComputeSpectralDimensionHeatKernel formula exactly.
    /// weights[idx] is in CSR format and corresponds to edge neighborIndices[idx].
    /// </summary>
    [ThreadGroupSize(64, 1, 1)]
    [GeneratedComputeShaderDescriptor]
    public readonly partial struct HeatDiffusionShader : IComputeShader
    {
        public readonly ReadWriteBuffer<float> newHeat;
        public readonly ReadWriteBuffer<float> heat;
        public readonly ReadOnlyBuffer<float> weights;
        public readonly ReadOnlyBuffer<int> neighborOffsets;
        public readonly ReadOnlyBuffer<int> neighborIndices;
        public readonly ReadOnlyBuffer<float> invSqrtDegrees;
        public readonly int nodeCount;
        public readonly float dt;

        public HeatDiffusionShader(
            ReadWriteBuffer<float> newHeat,
            ReadWriteBuffer<float> heat,
            ReadOnlyBuffer<float> weights,
            ReadOnlyBuffer<int> neighborOffsets,
            ReadOnlyBuffer<int> neighborIndices,
            ReadOnlyBuffer<float> invSqrtDegrees,
            int nodeCount,
            float dt)
        {
            this.newHeat = newHeat;
            this.heat = heat;
            this.weights = weights;
            this.neighborOffsets = neighborOffsets;
            this.neighborIndices = neighborIndices;
            this.invSqrtDegrees = invSqrtDegrees;
            this.nodeCount = nodeCount;
            this.dt = dt;
        }

        public void Execute()
        {
            int i = ThreadIds.X;
            if (i >= nodeCount) return;

            int start = neighborOffsets[i];
            int end = neighborOffsets[i + 1];

            float invSqrtDegI = invSqrtDegrees[i];

            // Nodes without neighbors retain their heat (matching CPU fallback)
            if (invSqrtDegI <= 0.0f || start == end)
            {
                newHeat[i] = heat[i];
                return;
            }

            // Compute off-diagonal sum: sum_j [w_ij * D_i^{-1/2} * D_j^{-1/2} * h_j]
            // This is the weighted averaging term from normalized Laplacian
            float offDiagonalSum = 0.0f;
            for (int idx = start; idx < end; idx++)
            {
                int j = neighborIndices[idx];
                float invSqrtDegJ = invSqrtDegrees[j];
                if (invSqrtDegJ <= 0.0f)
                {
                    continue;
                }

                // L_norm_ij = -w_ij / sqrt(deg_i * deg_j)
                // Contribution to sum is +w_ij * invSqrtDegI * invSqrtDegJ * heat[j]
                // (positive because we subtract -L_ij)
                float weightTerm = weights[idx] * invSqrtDegI * invSqrtDegJ;
                offDiagonalSum += weightTerm * heat[j];
            }

            // Heat diffusion: h_new = h - dt * L_norm * h
            // L_norm * h = h_i (diagonal) - offDiagonalSum
            // Therefore: h_new = h_i - dt * (h_i - offDiagonalSum)
            //                  = (1 - dt) * h_i + dt * offDiagonalSum
            float h_current = heat[i];
            float h_new = (1.0f - dt) * h_current + dt * offDiagonalSum;

            // Ensure non-negative (heat cannot be negative)
            newHeat[i] = Hlsl.Max(h_new, 0.0f);
        }
    }
}
