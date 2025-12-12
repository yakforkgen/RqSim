using RQSimulation;
using RQSimulation.GPUOptimized;
using System.Text.Json.Serialization;

namespace RqSimForms.Forms.Interfaces;

/// <summary>
/// Simulation API - contains all simulation logic separated from UI.
/// Used by Form_Main to run simulations without coupling to WinForms controls.
/// 
/// Threading model:
/// - Calculation thread: AboveNormal priority, writes to MetricsDispatcher
/// - UI thread: Normal priority, reads from MetricsDispatcher (lock-free)
/// - MetricsDispatcher: double-buffered, auto-decimates large datasets
/// </summary>
public class FormSimAPI
{
    // === Simulation Engine ===
    public SimulationEngine? SimulationEngine { get; private set; }
    public SimulationConfig? LastConfig { get; private set; }
    public SimulationResult? LastResult { get; private set; }
    public RQSimulation.ExampleModernSimulation.ScenarioResult? ModernResult { get; private set; }

    // === GPU Engines ===
    public bool GpuAvailable { get; set; }
    public SpectralWalkEngine? GpuSpectralWalkEngine { get; private set; }
    public GpuSpectralEngine? GpuHeatKernelEngine { get; private set; }
    public StatisticsEngine? GpuStatisticsEngine { get; private set; }
    public OptimizedGpuSimulationEngine? OptimizedGpuEngine { get; private set; }

    // === GPU Performance Tracking ===
    public GpuDiagnostics GpuStats { get; } = new();

    /// <summary>
    /// Holds GPU performance and shader diagnostics
    /// </summary>
    public class GpuDiagnostics
    {
        public bool IsActive { get; set; }
        public string? DeviceName { get; set; }
        public int KernelLaunches { get; set; }
        public double TotalGpuTimeMs { get; set; }
        public double TotalCopyTimeMs { get; set; }
        public int TopologyRebuilds { get; set; }
        public int WeightSyncs { get; set; }

        /// <summary>
        /// GPU-accelerated operations (shaders used per step)
        /// </summary>
        public string[] AcceleratedOperations { get; } =
        [
            "FormanCurvatureShader - Ricci curvature on edges (O(E) parallel)",
            "GravityShader - Weight evolution via curvature flow (O(E) parallel)",
            "ScalarLaplacianShader - Klein-Gordon field diffusion (O(N) parallel)",
            "ApplyScalarDeltaShader - Field update accumulation (O(N) parallel)",
            "SpectralWalkShader - Random walks for d_S (50K walkers parallel)"
        ];

        /// <summary>
        /// CPU-bound operations (not GPU accelerated)
        /// </summary>
        public string[] CpuBoundOperations { get; } =
        [
            "Topology updates (Step, BuildSoAViews) - every 100 steps",
            "Cluster detection (Union-Find) - every 10+ steps",
            "Metric collection (excited count, cluster stats) - every 10+ steps",
            "Quantum state updates - every step",
            "Spectral dimension computation - every 200 steps"
        ];

        /// <summary>
        /// Optimized intervals for performance (steps)
        /// </summary>
        public int TopologyUpdateInterval { get; set; } = 100;
        public int WeightSyncInterval { get; set; } = 50;
        public int MetricsInterval { get; set; } = 10;
        public int SpectralDimInterval { get; set; } = 200;

        public void Reset()
        {
            IsActive = false;
            KernelLaunches = 0;
            TotalGpuTimeMs = 0;
            TotalCopyTimeMs = 0;
            TopologyRebuilds = 0;
            WeightSyncs = 0;
        }

        public void RecordKernelLaunch(int count = 4) => KernelLaunches += count;
        public void RecordTopologyRebuild() => TopologyRebuilds++;
        public void RecordWeightSync() => WeightSyncs++;
    }

    // === Simulation State ===
    public bool IsModernRunning { get; set; }
    public bool SpectrumLoggingEnabled { get; set; }
    public DateTime SimulationWallClockStart { get; private set; }

    // === Metrics Dispatcher (thread-safe, double-buffered) ===
    public MetricsDispatcher Dispatcher { get; } = new();

    // === Legacy Live Metrics (for compatibility, delegate to Dispatcher) ===
    public int LiveStep => Dispatcher.LiveStep;
    public int LiveExcited => Dispatcher.LiveExcited;
    public double LiveHeavyMass => Dispatcher.LiveHeavyMass;
    public int LiveLargestCluster => Dispatcher.LiveLargestCluster;
    public int LiveStrongEdges => Dispatcher.LiveStrongEdges;
    public double LiveQNorm => Dispatcher.LiveQNorm;
    public double LiveEntanglement => Dispatcher.LiveEntanglement;
    public double LiveCorrelation => Dispatcher.LiveCorrelation;
    public double LiveSpectralDim => Dispatcher.LiveSpectralDim;
    public double LiveTemp => Dispatcher.LiveTemp;
    public double LiveEffectiveG => Dispatcher.LiveEffectiveG;
    public double LiveAdaptiveThreshold => Dispatcher.LiveAdaptiveThreshold;
    public int LiveTotalSteps => Dispatcher.LiveTotalSteps;

    // Cluster statistics accessors
    public int LiveHeavyClusterCount => Dispatcher.LiveHeavyClusterCount;
    public int LiveTotalClusters => Dispatcher.LiveTotalClusters;
    public double LiveAvgClusterMass => Dispatcher.LiveAvgClusterMass;
    public double LiveMaxClusterMass => Dispatcher.LiveMaxClusterMass;
    public double LiveAvgDegree => Dispatcher.LiveAvgDegree;

    public volatile bool SimulationComplete;

    // === Final State Values ===
    public double FinalSpectralDimension { get; private set; }
    public double FinalNetworkTemperature { get; private set; }

    // === Time Series Buffers (full resolution, for export) ===
    public readonly List<int> SeriesSteps = [];
    public readonly List<int> SeriesExcited = [];
    public readonly List<double> SeriesHeavyMass = [];
    public readonly List<int> SeriesLargestCluster = [];
    public readonly List<double> SeriesEnergy = [];
    public readonly List<int> SeriesHeavyCount = [];
    public readonly List<double> SeriesAvgDist = [];
    public readonly List<double> SeriesDensity = [];
    public readonly List<double> SeriesCorr = [];
    public readonly List<int> SeriesStrongEdges = [];
    public readonly List<double> SeriesQNorm = [];
    public readonly List<double> SeriesQEnergy = [];
    public readonly List<double> SeriesEntanglement = [];
    public readonly List<double> SeriesSpectralDimension = [];
    public readonly List<double> SeriesNetworkTemperature = [];
    public readonly List<double> SeriesEffectiveG = [];
    public readonly List<double> SeriesAdaptiveThreshold = [];

    // === Synthesis Analysis ===
    public List<(int volume, double deltaMass)>? SynthesisData { get; set; }
    public int SynthesisCount { get; set; }
    public int FissionCount { get; set; }

    // === Session History ===
    public readonly List<SimulationSession> SessionHistory = [];
    public SimulationSession? CurrentSession { get; set; }

    // === Callbacks for UI updates ===
    public Action<string>? OnConsoleLog { get; set; }

    // === Live Config (thread-safe, for runtime parameter updates) ===
    /// <summary>
    /// Thread-safe container for runtime-adjustable simulation parameters.
    /// UI thread writes, calculation thread reads.
    /// Note: double reads/writes are atomic on x64, no volatile needed.
    /// </summary>
    public class LiveConfigData
    {
        // Physics Constants (grpPhysicsConstants)
        public double GravitationalCoupling = 0.05;
        public double VacuumEnergyScale = 0.00005;
        public double AnnealingCoolingRate = 0.995;
        public double DecoherenceRate = 0.001;
        public double HotStartTemperature = 3.0;
        public double InitialEdgeProb = 0.02;
        public double AdaptiveThresholdSigma = 1.5;
        public double WarmupDuration = 200;
        public double GravityTransitionDuration = 137;

        // Simulation Parameters (grpSimParams)
        public volatile int TargetDegree = 8;
        public double InitialExcitedProb = 0.02;
        public double LambdaState = 0.5;
        public double Temperature = 10.0;
        public double EdgeTrialProb = 0.02;
        public double MeasurementThreshold = 0.3;
        public volatile int FractalLevels = 0;
        public volatile int FractalBranchFactor = 0;

        /// <summary>
        /// Timestamp of last update (for change detection)
        /// </summary>
        private long _lastUpdateTicks;
        public long LastUpdateTicks
        {
            get => Interlocked.Read(ref _lastUpdateTicks);
            private set => Interlocked.Exchange(ref _lastUpdateTicks, value);
        }

        /// <summary>
        /// Marks config as updated
        /// </summary>
        public void MarkUpdated() => LastUpdateTicks = DateTime.UtcNow.Ticks;
    }

    /// <summary>
    /// Live configuration that can be modified during simulation run.
    /// Thread-safe: UI writes, calculation thread reads.
    /// </summary>
    public LiveConfigData LiveConfig { get; } = new();

    // === Auto-Tuning ===
    /// <summary>
    /// Whether auto-tuning is enabled (UI checkbox state)
    /// </summary>
    public volatile bool AutoTuningEnabled;

    /// <summary>
    /// Random generator for auto-tuning perturbations
    /// </summary>
    private readonly Random _autoTuneRng = new(42);

    /// <summary>
    /// Last step when auto-tuning was applied
    /// </summary>
    private int _lastAutoTuneStep;

    /// <summary>
    /// Counter for consecutive steps where giant cluster persists despite intervention
    /// Used to trigger more aggressive measures
    /// </summary>
    private int _persistentGiantClusterCount;

    /// <summary>
    /// Auto-tuning interval (steps between adjustments)
    /// Reduced from 500 to 100 for faster response to d_S collapse
    /// </summary>
    public int AutoTuneInterval { get; set; } = 100;

    /// <summary>
    /// Request to trigger topology tunneling (checked by main loop)
    /// </summary>
    public volatile bool RequestTopologyTunneling;

    // === Graph Health Live Config (from UI NumericUpDown controls) ===
    // These mirror PhysicsConstants but can be adjusted via UI at runtime
    // Note: double reads/writes are atomic on x64, no volatile needed
    public class GraphHealthConfig
    {
        public double GiantClusterThreshold = PhysicsConstants.GiantClusterThreshold;
        public double EmergencyGiantClusterThreshold = PhysicsConstants.EmergencyGiantClusterThreshold;
        public double GiantClusterDecoherenceRate = PhysicsConstants.GiantClusterDecoherenceRate;
        public double MaxDecoherenceEdgesFraction = PhysicsConstants.MaxDecoherenceEdgesFraction;
        public double CriticalSpectralDimension = PhysicsConstants.CriticalSpectralDimension;
        public double WarningSpectralDimension = PhysicsConstants.WarningSpectralDimension;
    }

    /// <summary>
    /// Live Graph Health configuration (updated from UI controls).
    /// Thread-safe: UI writes, calculation thread reads.
    /// </summary>
    public GraphHealthConfig GraphHealthLive { get; } = new();

    /// <summary>
    /// Performs auto-tuning based on current simulation metrics.
    /// Called from calculation thread every AutoTuneInterval steps.
    /// 
    /// RQ-HYPOTHESIS COMPLIANT AUTO-TUNING:
    /// Uses Graph Health parameters from UI (GraphHealthLive) instead of constants.
    /// 
    /// Priority order:
    /// 1. FRAGMENTATION (d_S less than Critical) - Most urgent, add edges + reduce G
    /// 2. GIANT CLUSTER (greater than threshold) - Break up via decoherence
    /// 3. LOW d_S (less than Warning) - Preventive gravity reduction
    /// 4. HIGH d_S (greater than 5) - Increase gravity to compact
    /// 5. CLUSTER TUNING - Adjust threshold for detection
    /// 6. ACTIVITY BALANCE - Decoherence for frozen/hyperactive states
    /// </summary>
    /// <returns>Description of adjustments made, or null if no tuning needed</returns>
    public string? PerformAutoTuning(int step, double spectralDim, int excitedCount, int clusterCount,
        int largestCluster, double heavyMass, int nodeCount)
    {
        if (!AutoTuningEnabled) return null;
        if (step - _lastAutoTuneStep < AutoTuneInterval) return null;

        _lastAutoTuneStep = step;
        var adjustments = new List<string>();

        // Read live Graph Health thresholds from UI
        double criticalDim = GraphHealthLive.CriticalSpectralDimension;
        double warningDim = GraphHealthLive.WarningSpectralDimension;
        double giantThreshold = GraphHealthLive.GiantClusterThreshold;
        double emergencyThreshold = GraphHealthLive.EmergencyGiantClusterThreshold;
        double decoherenceRate = GraphHealthLive.GiantClusterDecoherenceRate;

        double clusterRatio = (double)largestCluster / Math.Max(1, nodeCount);
        double excitedRatio = (double)excitedCount / Math.Max(1, nodeCount);

        // === PRIORITY 1: FRAGMENTATION DETECTION (d_S less than Critical) ===
        if (spectralDim > 0 && spectralDim <= criticalDim)
        {
            // CRITICAL: Graph has collapsed! Emergency recovery
            double oldG = LiveConfig.GravitationalCoupling;
            LiveConfig.GravitationalCoupling *= 0.2; // Aggressive reduction
            LiveConfig.GravitationalCoupling = Math.Max(0.0005, LiveConfig.GravitationalCoupling);

            double oldEdgeProb = LiveConfig.EdgeTrialProb;
            LiveConfig.EdgeTrialProb = Math.Min(0.5, LiveConfig.EdgeTrialProb * 3.0);

            double oldTemp = LiveConfig.HotStartTemperature;
            LiveConfig.HotStartTemperature = Math.Min(20.0, LiveConfig.HotStartTemperature * 2.0);

            adjustments.Add($"FRAGMENTED d_S={spectralDim:F2}<={criticalDim}: G {oldG:F4}->{LiveConfig.GravitationalCoupling:F4}, " +
                           $"EdgeProb {oldEdgeProb:F3}->{LiveConfig.EdgeTrialProb:F3}, Temp->{LiveConfig.HotStartTemperature:F1}");

            LiveConfig.MarkUpdated();
            return string.Join("; ", adjustments);
        }

        // === PRIORITY 2: GIANT CLUSTER DETECTION ===
        // Extreme threshold (>70%): Decoherence alone is insufficient
        const double extremeClusterThreshold = 0.70;

        if (clusterRatio >= extremeClusterThreshold)
        {
            _persistentGiantClusterCount++;

            // EXTREME: Cluster dominates despite previous interventions
            double oldDecoherence = LiveConfig.DecoherenceRate;
            LiveConfig.DecoherenceRate = 0.15; // Maximum decoherence

            double oldG = LiveConfig.GravitationalCoupling;
            LiveConfig.GravitationalCoupling = 0.0005; // Near-zero gravity

            // After 3 consecutive extreme clusters, trigger topology tunneling
            if (_persistentGiantClusterCount >= 3)
            {
                RequestTopologyTunneling = true;
                adjustments.Add($"EXTREME CLUSTER {clusterRatio:P0}>=70% (persist={_persistentGiantClusterCount}): " +
                               $"TOPOLOGY TUNNELING REQUESTED, Decoherence->0.15, G->0.0005");
            }
            else
            {
                adjustments.Add($"EXTREME CLUSTER {clusterRatio:P0}>=70% (persist={_persistentGiantClusterCount}): " +
                               $"Decoherence {oldDecoherence:F4}->0.15, G {oldG:F4}->0.0005");
            }

            // Also reduce edge trial probability to prevent reconnection
            LiveConfig.EdgeTrialProb = Math.Max(0.001, LiveConfig.EdgeTrialProb * 0.3);
            // Immediately apply decoherence to the graph if available
            try
            {
                var graph = SimulationEngine?.Graph;
                if (graph != null)
                {
                    int weakened = graph.ApplyGiantClusterDecoherence();
                    adjustments.Add($"AppliedGraphDecoherence: weakened {weakened} edges");

                    // If persistent extreme cluster, perform topology tunneling
                    if (_persistentGiantClusterCount >= 3)
                    {
                        int removed = graph.PerformTopologyTunneling(removalFraction: 0.30);
                        adjustments.Add($"TopologyTunneling: removed {removed} edges");
                    }
                }
            }
            catch (Exception ex)
            {
                OnConsoleLog?.Invoke($"[AUTO-TUNE] Error applying decoherence/tunneling: {ex.Message}\n");
            }
        }
        else if (clusterRatio >= emergencyThreshold)
        {
            // Reset persistence counter when cluster is breaking up
            if (clusterRatio < extremeClusterThreshold - 0.05)
                _persistentGiantClusterCount = Math.Max(0, _persistentGiantClusterCount - 1);

            // EMERGENCY: Single giant cluster dominates - aggressive decoherence
            double oldDecoherence = LiveConfig.DecoherenceRate;
            LiveConfig.DecoherenceRate = Math.Min(0.1, LiveConfig.DecoherenceRate * 3.0 + decoherenceRate);

            double oldG = LiveConfig.GravitationalCoupling;
            LiveConfig.GravitationalCoupling *= 0.5;
            LiveConfig.GravitationalCoupling = Math.Max(0.001, LiveConfig.GravitationalCoupling);

            adjustments.Add($"EMERGENCY CLUSTER {clusterRatio:P0}>={emergencyThreshold:P0}: " +
                           $"Decoherence {oldDecoherence:F4}->{LiveConfig.DecoherenceRate:F4}, G->{LiveConfig.GravitationalCoupling:F4}");

            // Apply decoherence immediately
            try
            {
                var graph = SimulationEngine?.Graph;
                if (graph != null)
                {
                    int weakened = graph.ApplyGiantClusterDecoherence();
                    adjustments.Add($"AppliedGraphDecoherence: weakened {weakened} edges");

                    // Inject additional decoherence noise into giant clusters when emergency
                    var clusters = graph.GetStrongCorrelationClusters(graph.GetAdaptiveHeavyThreshold());
                    int giantCount = (int)(nodeCount * emergencyThreshold);
                    foreach (var cluster in clusters.Where(c => c.Count >= giantCount))
                    {
                        // noise amplitude proportional to LiveConfig.DecoherenceRate (clamped)
                        double amp = Math.Clamp(LiveConfig.DecoherenceRate * 10.0, 0.05, 0.25);
                        graph.InjectDecoherenceIntoCluster(cluster, noiseAmplitude: amp);
                        adjustments.Add($"Injected noise into cluster size {cluster.Count} (amp={amp:F3})");
                    }
                }
            }
            catch (Exception ex)
            {
                OnConsoleLog?.Invoke($"[AUTO-TUNE] Error applying emergency decoherence: {ex.Message}\n");
            }
        }
        else if (clusterRatio >= giantThreshold)
        {
            // Reset persistence counter
            _persistentGiantClusterCount = 0;

            // Giant cluster forming - apply decoherence
            double oldDecoherence = LiveConfig.DecoherenceRate;
            LiveConfig.DecoherenceRate = Math.Min(0.05, LiveConfig.DecoherenceRate * 1.5 + decoherenceRate * 0.5);

            adjustments.Add($"GIANT CLUSTER {clusterRatio:P0}>={giantThreshold:P0}: Decoherence->{LiveConfig.DecoherenceRate:F4}");

            // Apply a gentle decoherence pass immediately
            try
            {
                var graph = SimulationEngine?.Graph;
                if (graph != null)
                {
                    int weakened = graph.ApplyGiantClusterDecoherence();
                    adjustments.Add($"AppliedGraphDecoherence: weakened {weakened} edges");
                }
            }
            catch (Exception ex)
            {
                OnConsoleLog?.Invoke($"[AUTO-TUNE] Error applying giant-cluster decoherence: {ex.Message}\n");
            }
        }
        else
        {
            // Cluster is under control - reset persistence
            _persistentGiantClusterCount = 0;
        }

        // === PRIORITY 3: LOW SPECTRAL DIMENSION (Warning Zone) ===
        if (spectralDim > criticalDim && spectralDim < warningDim)
        {
            // Approaching fragmentation - preventive measures
            double oldG = LiveConfig.GravitationalCoupling;
            LiveConfig.GravitationalCoupling *= 0.6;
            LiveConfig.GravitationalCoupling = Math.Max(0.001, LiveConfig.GravitationalCoupling);

            double oldEdgeProb = LiveConfig.EdgeTrialProb;
            LiveConfig.EdgeTrialProb = Math.Min(0.3, LiveConfig.EdgeTrialProb * 1.5);

            adjustments.Add($"d_S WARNING {spectralDim:F2}<{warningDim}: G {oldG:F4}->{LiveConfig.GravitationalCoupling:F4}, " +
                           $"EdgeProb->{LiveConfig.EdgeTrialProb:F3}");
        }
        // === HEALTHY ZONE: Restore parameters gradually ===
        else if (spectralDim >= 2.0 && spectralDim <= 4.0)
        {
            // Healthy dimension - restore suppressed parameters
            bool restored = false;

            if (LiveConfig.GravitationalCoupling < 0.03)
            {
                double oldG = LiveConfig.GravitationalCoupling;
                LiveConfig.GravitationalCoupling = Math.Min(0.1, LiveConfig.GravitationalCoupling * 1.3);
                adjustments.Add($"d_S HEALTHY: G restore {oldG:F4}->{LiveConfig.GravitationalCoupling:F4}");
                restored = true;
            }

            if (LiveConfig.DecoherenceRate > 0.02 && clusterRatio < giantThreshold * 0.5)
            {
                double oldDec = LiveConfig.DecoherenceRate;
                LiveConfig.DecoherenceRate = Math.Max(0.001, LiveConfig.DecoherenceRate * 0.8);
                if (!restored) adjustments.Add($"d_S HEALTHY: Decoherence {oldDec:F4}->{LiveConfig.DecoherenceRate:F4}");
            }
        }
        // === HIGH d_S (Hyperbolic) ===
        else if (spectralDim > 5.0)
        {
            double oldG = LiveConfig.GravitationalCoupling;
            LiveConfig.GravitationalCoupling = Math.Min(0.5, LiveConfig.GravitationalCoupling * 1.5);
            adjustments.Add($"d_S HIGH {spectralDim:F1}>5: G {oldG:F4}->{LiveConfig.GravitationalCoupling:F4}");
        }

        // === PRIORITY 4: CLUSTER FORMATION TUNING ===
        if (clusterCount < 3 && heavyMass < 50 && step > 200 && spectralDim >= warningDim)
        {
            double oldThreshold = LiveConfig.AdaptiveThresholdSigma;
            LiveConfig.AdaptiveThresholdSigma = Math.Max(0.3, LiveConfig.AdaptiveThresholdSigma * 0.85);
            adjustments.Add($"No clusters: Threshold {oldThreshold:F2}->{LiveConfig.AdaptiveThresholdSigma:F2}");
        }

        // === PRIORITY 5: ACTIVITY BALANCE ===
        if (excitedRatio > 0.6 && spectralDim >= warningDim)
        {
            double oldDecoherence = LiveConfig.DecoherenceRate;
            LiveConfig.DecoherenceRate = Math.Min(0.1, LiveConfig.DecoherenceRate * 1.3);
            adjustments.Add($"Hyperactive {excitedRatio:P0}: Decoherence->{LiveConfig.DecoherenceRate:F4}");
        }
        else if (excitedRatio < 0.03 && step > 300)
        {
            double oldDecoherence = LiveConfig.DecoherenceRate;
            LiveConfig.DecoherenceRate = Math.Max(0.0001, LiveConfig.DecoherenceRate * 0.7);
            adjustments.Add($"Frozen {excitedRatio:P0}: Decoherence->{LiveConfig.DecoherenceRate:F4}");
        }

        // === EXPLORATION (5% chance, only if stable) ===
        if (_autoTuneRng.NextDouble() < 0.05 && adjustments.Count == 0 &&
            spectralDim >= 2.0 && spectralDim <= 4.0 && clusterRatio < giantThreshold)
        {
            int param = _autoTuneRng.Next(3);
            double factor = 0.95 + _autoTuneRng.NextDouble() * 0.1;

            switch (param)
            {
                case 0:
                    double oldG = LiveConfig.GravitationalCoupling;
                    LiveConfig.GravitationalCoupling = Math.Clamp(LiveConfig.GravitationalCoupling * factor, 0.005, 0.3);
                    adjustments.Add($"Explore: G {oldG:F4}->{LiveConfig.GravitationalCoupling:F4}");
                    break;
                case 1:
                    double oldTemp = LiveConfig.HotStartTemperature;
                    LiveConfig.HotStartTemperature = Math.Clamp(LiveConfig.HotStartTemperature * factor, 2.0, 15.0);
                    adjustments.Add($"Explore: Temp {oldTemp:F2}->{LiveConfig.HotStartTemperature:F2}");
                    break;
                case 2:
                    double oldDec = LiveConfig.DecoherenceRate;
                    LiveConfig.DecoherenceRate = Math.Clamp(LiveConfig.DecoherenceRate * factor, 0.0005, 0.05);
                    adjustments.Add($"Explore: Decoherence {oldDec:F4}->{LiveConfig.DecoherenceRate:F4}");
                    break;
            }
        }

        if (adjustments.Count == 0) return null;

        LiveConfig.MarkUpdated();
        return string.Join("; ", adjustments);
    }

    /// <summary>
    /// Initializes LiveConfig from SimulationConfig at simulation start
    /// </summary>
    public void InitializeLiveConfig(SimulationConfig config)
    {
        LiveConfig.GravitationalCoupling = config.GravitationalCoupling;
        LiveConfig.VacuumEnergyScale = config.VacuumEnergyScale;
        LiveConfig.AnnealingCoolingRate = config.AnnealingCoolingRate;
        LiveConfig.DecoherenceRate = config.DecoherenceRate;
        LiveConfig.HotStartTemperature = config.HotStartTemperature;
        LiveConfig.InitialEdgeProb = config.InitialEdgeProb;
        LiveConfig.TargetDegree = config.TargetDegree;
        LiveConfig.InitialExcitedProb = config.InitialExcitedProb;
        LiveConfig.LambdaState = config.LambdaState;
        LiveConfig.Temperature = config.Temperature;
        LiveConfig.EdgeTrialProb = config.EdgeTrialProbability;
        LiveConfig.MeasurementThreshold = config.MeasurementThreshold;
        LiveConfig.FractalLevels = config.FractalLevels;
        LiveConfig.FractalBranchFactor = config.FractalBranchFactor;
        LiveConfig.MarkUpdated();
    }

    /// <summary>
    /// Clears all time series buffers and dispatcher
    /// </summary>
    public void ClearTimeSeries()
    {
        SeriesSteps.Clear();
        SeriesExcited.Clear();
        SeriesHeavyMass.Clear();
        SeriesLargestCluster.Clear();
        SeriesEnergy.Clear();
        SeriesHeavyCount.Clear();
        SeriesStrongEdges.Clear();
        SeriesCorr.Clear();
        SeriesQNorm.Clear();
        SeriesEntanglement.Clear();
        SeriesSpectralDimension.Clear();
        SeriesNetworkTemperature.Clear();
        SeriesEffectiveG.Clear();
        SeriesAdaptiveThreshold.Clear();
        SeriesAvgDist.Clear();
        SeriesDensity.Clear();
        SeriesQEnergy.Clear();

        // Clear dispatcher for new simulation
        Dispatcher.Clear();
    }

    // === Custom Initializer for Experiments ===
    private Action<RQGraph>? _pendingCustomInitializer;

    /// <summary>
    /// Sets a custom graph initializer to be used in the next simulation.
    /// Called by Form_Main when loading an experiment with custom topology.
    /// </summary>
    public void SetCustomInitializer(Action<RQGraph>? initializer)
    {
        _pendingCustomInitializer = initializer;
    }

    /// <summary>
    /// Initializes simulation engine with given config
    /// </summary>
    public void InitializeSimulation(SimulationConfig config)
    {
        LastConfig = config;
        SimulationEngine = new SimulationEngine(config, _pendingCustomInitializer);
        _pendingCustomInitializer = null; // Clear after use
        Dispatcher.LiveTotalSteps = config.TotalSteps;
        SimulationWallClockStart = DateTime.UtcNow;
    }

    /// <summary>
    /// Initializes GPU engines for the current graph
    /// </summary>
    public bool InitializeGpuEngines()
    {
        if (SimulationEngine?.Graph == null) return false;

        var graph = SimulationEngine.Graph;
        int edgeCount = graph.FlatEdgesFrom.Length;
        if (edgeCount == 0) return false;

        try
        {
            graph.BuildSoAViews();
            OptimizedGpuEngine = new OptimizedGpuSimulationEngine(graph);
            OptimizedGpuEngine.Initialize();
            OptimizedGpuEngine.UploadState();

            // Also initialize the standalone GpuGravityEngine for ImprovedNetworkGravity
            // This enables GPU-accelerated gravity in CPU simulation paths too
            bool gravityGpuOk = graph.InitGpuGravity();
            if (gravityGpuOk)
            {
                OnConsoleLog?.Invoke("[GPU] GpuGravityEngine initialized for Ollivier/Forman curvature\n");
            }

            int totalDirectedEdges = graph.CsrOffsets[graph.N];
            float[] csrWeights = new float[totalDirectedEdges];
            for (int n = 0; n < graph.N; n++)
            {
                int start = graph.CsrOffsets[n];
                int end = graph.CsrOffsets[n + 1];
                for (int k = start; k < end; k++)
                {
                    int to = graph.CsrIndices[k];
                    csrWeights[k] = (float)graph.Weights[n, to];
                }
            }

            GpuSpectralWalkEngine = new SpectralWalkEngine();
            int walkerCount = Math.Clamp(graph.N * 100, 50000, 200000);
            GpuSpectralWalkEngine.Initialize(walkerCount, graph.N, totalDirectedEdges);
            GpuSpectralWalkEngine.UpdateTopology(graph.CsrOffsets, graph.CsrIndices, csrWeights);
            GpuSpectralWalkEngine.InitializeWalkersUniform();

            // Initialize Heat Kernel GPU engine (alternative d_S computation method)
            // Heat Kernel is better for dense graphs, Random Walk for sparse
            GpuHeatKernelEngine = new GpuSpectralEngine();
            GpuHeatKernelEngine.UpdateTopology(graph);

            GpuStatisticsEngine = new StatisticsEngine();
            GpuStatisticsEngine.Initialize(Math.Max(graph.N, edgeCount));

            return true;
        }
        catch (Exception ex)
        {
            OnConsoleLog?.Invoke($"[GPU] Error: {ex.Message}, using CPU fallback\n");
            DisposeGpuEngines();
            return false;
        }
    }

    /// <summary>
    /// Disposes GPU engines
    /// </summary>
    public void DisposeGpuEngines()
    {
        GpuSpectralWalkEngine?.Dispose();
        GpuSpectralWalkEngine = null;
        GpuHeatKernelEngine?.Dispose();
        GpuHeatKernelEngine = null;
        GpuStatisticsEngine?.Dispose();
        GpuStatisticsEngine = null;
        OptimizedGpuEngine?.Dispose();
        OptimizedGpuEngine = null;

        // Also dispose standalone gravity engine
        SimulationEngine?.Graph?.DisposeGpuGravity();
    }

    /// <summary>
    /// Compute spectral dimension using the best available GPU method.
    /// 
    /// Method selection:
    /// - Heat Kernel (GpuHeatKernelEngine): Better for dense graphs (density > 10%)
    /// - Random Walk (GpuSpectralWalkEngine): Better for sparse graphs
    /// 
    /// Both methods automatically sync topology when graph.TopologyVersion changes.
    /// </summary>
    /// <param name="graph">The RQGraph to compute d_S for</param>
    /// <param name="enableCpuComparison">If true, also runs CPU version for debugging</param>
    /// <returns>Spectral dimension estimate, or NaN if computation failed</returns>
    public double ComputeSpectralDimensionGpu(RQGraph graph, bool enableCpuComparison = false)
    {
        ArgumentNullException.ThrowIfNull(graph);

        // Compute graph density to choose method
        int edgeCount = 0;
        for (int i = 0; i < graph.N; i++)
        {
            foreach (int j in graph.Neighbors(i))
            {
                if (j > i) edgeCount++;
            }
        }
        int maxEdges = graph.N * (graph.N - 1) / 2;
        double density = maxEdges > 0 ? (double)edgeCount / maxEdges : 0;

        // Choose method based on density
        if (density > 0.10 && GpuHeatKernelEngine != null)
        {
            // Dense graph: Heat Kernel is more stable
            float ds = GpuHeatKernelEngine.ComputeSpectralDimension(
                graph,
                dt: 0.01f,
                numSteps: 100,
                numProbeVectors: 8,
                enableCpuComparison: enableCpuComparison);

            if (!float.IsNaN(ds))
            {
                OnConsoleLog?.Invoke($"[GPU d_S] HeatKernel method: d_S={ds:F4} (density={density:P1})\n");
                return ds;
            }
        }

        // Sparse graph or HeatKernel fallback: Random Walk
        if (GpuSpectralWalkEngine != null)
        {
            double ds = GpuSpectralWalkEngine.ComputeSpectralDimensionWithSyncCheck(
                graph,
                numSteps: 100,
                walkerCount: 10000,
                skipInitial: 10);

            if (!double.IsNaN(ds))
            {
                OnConsoleLog?.Invoke($"[GPU d_S] RandomWalk method: d_S={ds:F4} (density={density:P1})\n");
                return ds;
            }
        }

        // CPU fallback
        OnConsoleLog?.Invoke($"[GPU d_S] GPU methods failed, using CPU fallback\n");
        return graph.ComputeSpectralDimension(t_max: 100, num_walkers: 50);
    }

    /// <summary>
    /// Runs one physics step with GPU or CPU
    /// </summary>
    public void RunPhysicsStep(int step, double dt, double effectiveG, bool useGpu)
    {
        if (SimulationEngine?.Graph == null) return;
        var graph = SimulationEngine.Graph;

        if (effectiveG > 0)
        {
            if (useGpu && OptimizedGpuEngine != null)
            {
                if (step == 0) OptimizedGpuEngine.UploadState();

                OptimizedGpuEngine.StepGpu(
                    dt: (float)dt,
                    G: (float)effectiveG,
                    lambda: (float)PhysicsConstants.CosmologicalConstant,
                    degreePenalty: (float)PhysicsConstants.DegreePenaltyFactor,
                    diffusionRate: (float)PhysicsConstants.FieldDiffusionRate,
                    scalarMass: (float)PhysicsConstants.KleinGordonMass);

                // Track GPU kernel launches (4 shaders per step)
                GpuStats.RecordKernelLaunch(4);

                // Sync weights less frequently (every 50 steps instead of 10)
                // This reduces GPU?CPU transfer overhead significantly
                if (step % 50 == 0)
                {
                    OptimizedGpuEngine.SyncWeightsToGraph();
                    OptimizedGpuEngine.SyncScalarFieldToGraph();
                    GpuStats.RecordWeightSync();
                }
            }
            else
            {
                ImprovedNetworkGravity.EvolveNetworkGeometryOllivierDynamic(
                    graph, dt: dt, effectiveG: effectiveG);
            }
        }
    }

    /// <summary>
    /// Updates topology and syncs GPU buffers if needed
    /// </summary>
    public void UpdateTopology(bool useGpu)
    {
        if (SimulationEngine?.Graph == null) return;
        var graph = SimulationEngine.Graph;

        graph.Step();

        if (useGpu && OptimizedGpuEngine != null)
        {
            graph.BuildSoAViews();
            OptimizedGpuEngine.UpdateTopologyBuffers();
            OptimizedGpuEngine.UploadState();
            GpuStats.RecordTopologyRebuild();
        }
    }

    /// <summary>
    /// Collects metrics from current graph state (full version - expensive)
    /// </summary>
    public (int excited, double heavyMass, int heavyCount, int largestCluster, double energy,
            int strongEdges, double correlation, double qNorm, double entanglement,
            int totalClusters, double avgClusterMass, double maxClusterMass, double avgDegree)
        CollectMetrics()
    {
        if (SimulationEngine?.Graph == null)
            return (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

        var graph = SimulationEngine.Graph;

        int excited = graph.State.Count(s => s == NodeState.Excited);
        double effectiveThreshold = Math.Min(graph.GetAdaptiveHeavyThreshold(), RQGraph.HeavyClusterThreshold);
        var heavyStats = graph.GetHeavyClusterStatsCorrelationMass(effectiveThreshold, RQGraph.HeavyClusterMinSize);
        double heavyMass = heavyStats.totalMass;
        int heavyCount = heavyStats.count;

        var clusters = graph.GetStrongCorrelationClusters(effectiveThreshold);
        int largestCluster = clusters.Count > 0 ? clusters.Max(c => c.Count) : 0;
        int totalClusters = clusters.Count;

        // Compute cluster mass statistics
        double avgClusterMass = 0.0;
        double maxClusterMass = 0.0;
        if (clusters.Count > 0)
        {
            double totalMassSum = 0.0;
            foreach (var cluster in clusters)
            {
                double clusterMass = 0.0;
                foreach (int u in cluster)
                {
                    foreach (int v in cluster)
                    {
                        if (u < v && graph.Edges[u, v])
                            clusterMass += graph.Weights[u, v];
                    }
                }
                totalMassSum += clusterMass;
                if (clusterMass > maxClusterMass)
                    maxClusterMass = clusterMass;
            }
            avgClusterMass = totalMassSum / clusters.Count;
        }

        double energy = graph.ComputeTotalEnergy();

        var weightStats = graph.GetWeightStats(effectiveThreshold);
        int strongEdges = weightStats.strongEdges;
        double correlation = weightStats.avgWeight;
        double qNorm = graph.GetQuantumNorm();

        // Compute average degree
        var edgeStats = graph.GetEdgeStats();
        double avgDegree = edgeStats.avgDegree;

        var largestClusterNodes = clusters.Count > 0 ? clusters.OrderByDescending(c => c.Count).First() : new List<int>();
        double entanglement = largestClusterNodes.Count > 0 ? graph.ComputeEntanglementEntropy(largestClusterNodes) : 0.0;

        return (excited, heavyMass, heavyCount, largestCluster, energy, strongEdges, correlation, qNorm, entanglement,
                totalClusters, avgClusterMass, maxClusterMass, avgDegree);
    }

    /// <summary>
    /// Collects only lightweight metrics (fast - O(N) instead of O(N?))
    /// Use this for steps where full metrics are not needed
    /// </summary>
    public int CollectExcitedCount()
    {
        if (SimulationEngine?.Graph == null) return 0;
        return SimulationEngine.Graph.State.Count(s => s == NodeState.Excited);
    }

    /// <summary>
    /// Stores metrics to time series buffers and dispatcher.
    /// Called from calculation thread - minimal lock time.
    /// </summary>
    public void StoreMetrics(int step, int excited, double heavyMass, int heavyCount, int largestCluster,
                             double energy, int strongEdges, double correlation, double qNorm, double entanglement,
                             double spectralDim, double temp, double effectiveG, double threshold,
                             int totalClusters = 0, double avgClusterMass = 0.0, double maxClusterMass = 0.0, double avgDegree = 0.0,
                             int edgeCount = 0, int componentCount = 1)
    {
        // Store to full-resolution lists (for export)
        SeriesSteps.Add(step);
        SeriesExcited.Add(excited);
        SeriesHeavyMass.Add(heavyMass);
        SeriesLargestCluster.Add(largestCluster);
        SeriesEnergy.Add(energy);
        SeriesHeavyCount.Add(heavyCount);
        SeriesStrongEdges.Add(strongEdges);
        SeriesCorr.Add(correlation);
        SeriesQNorm.Add(qNorm);
        SeriesEntanglement.Add(entanglement);
        SeriesSpectralDimension.Add(spectralDim);
        SeriesNetworkTemperature.Add(temp);
        SeriesEffectiveG.Add(effectiveG);
        SeriesAdaptiveThreshold.Add(threshold);

        // Get node count for ratios
        int nodeCount = SimulationEngine?.Graph?.N ?? 0;

        // Update dispatcher (thread-safe, will be decimated for UI)
        // NOTE: Pass Dispatcher.LiveTotalSteps instead of LastConfig.TotalSteps
        // because in event-based mode, LiveTotalSteps is set to sweepCount, not config steps
        int totalStepsForUi = Dispatcher.LiveTotalSteps > 0 ? Dispatcher.LiveTotalSteps : (LastConfig?.TotalSteps ?? 0);
        Dispatcher.UpdateLiveMetrics(step, excited, heavyMass, largestCluster, strongEdges,
            qNorm, entanglement, correlation, spectralDim, temp, effectiveG, totalStepsForUi, threshold,
            heavyCount, totalClusters, avgClusterMass, maxClusterMass, avgDegree,
            edgeCount, componentCount, nodeCount);

        Dispatcher.AppendTimeSeriesPoint(step, excited, heavyMass, heavyCount, largestCluster,
            energy, strongEdges, correlation, qNorm, entanglement, spectralDim, temp, effectiveG, threshold);
    }

    /// <summary>
    /// Computes annealing temperature for given step using dynamic ?.
    /// 
    /// RQ-FIX: Uses PhysicsConstants.ComputeAnnealingTimeConstant(totalSteps)
    /// instead of fixed ?=18779 which was too slow for typical simulations.
    /// 
    /// Formula: T(t) = T_f + (T_i - T_f) ? exp(-t/?)
    /// where ? = totalSteps / 5 ensures ~99% cooling by simulation end.
    /// </summary>
    public double ComputeAnnealingTemperature(int step, double startTemp, int totalSteps)
    {
        double finalTemp = PhysicsConstants.FinalAnnealingTemperature;

        // RQ-FIX: Use dynamic ? based on simulation length
        // This ensures proper cooling regardless of TotalSteps
        double annealingTau = PhysicsConstants.ComputeAnnealingTimeConstant(totalSteps);

        return finalTemp + (startTemp - finalTemp) * Math.Exp(-step / annealingTau);
    }

    /// <summary>
    /// Computes effective gravitational coupling with warmup
    /// </summary>
    public double ComputeEffectiveG(int step, double targetG)
    {
        double warmupG = PhysicsConstants.WarmupGravitationalCoupling;
        double warmupEnd = PhysicsConstants.WarmupDuration;
        double transitionDuration = PhysicsConstants.GravityTransitionDuration;

        if (step < warmupEnd)
            return warmupG;
        else if (step < warmupEnd + transitionDuration)
            return warmupG + (targetG - warmupG) * ((step - warmupEnd) / transitionDuration);
        else
            return targetG;
    }

    /// <summary>
    /// Finalizes simulation and stores final metrics
    /// </summary>
    public void FinalizeSimulation(List<int> excitedHistory)
    {
        if (SimulationEngine?.Graph == null) return;
        var graph = SimulationEngine.Graph;

        // Get final spectral dimension - compute fresh if series is empty or last value is 0
        double finalSpectralDim = SeriesSpectralDimension.Count > 0
            ? SeriesSpectralDimension.Where(d => d != 0).LastOrDefault()
            : 0;

        // If still 0, compute it now
        if (finalSpectralDim == 0)
        {
            finalSpectralDim = graph.ComputeSpectralDimension(t_max: 100, num_walkers: 100);
        }

        FinalSpectralDimension = finalSpectralDim;
        FinalNetworkTemperature = graph.NetworkTemperature;

        int finalExcited = graph.State.Count(s => s == NodeState.Excited);
        int maxExcited = excitedHistory.Count > 0 ? excitedHistory.Max() : 0;
        double avgExcitedFinal = excitedHistory.Count > 0 ? excitedHistory.Average() : 0.0;
        var finalHeavy = graph.ComputeHeavyClustersEnergy();
        var finalClusters = graph.GetStrongCorrelationClusters(graph.GetAdaptiveHeavyThreshold());
        int finalHeavyCount = finalClusters.Count(c => c.Count >= RQGraph.HeavyClusterMinSize);
        double scalarEnergy = graph.ComputeScalarFieldEnergy();
        double simulationTime = SeriesSteps.Count * 0.01;

        ModernResult = new RQSimulation.ExampleModernSimulation.ScenarioResult
        {
            FinalTime = simulationTime,
            ExcitedCount = finalExcited,
            HeavyClusterCount = finalHeavyCount,
            HeavyClusterTotalMass = finalHeavy.totalMass,
            HeavyClusterMaxMass = finalHeavy.maxMass,
            HeavyClusterMeanMass = finalHeavyCount > 0 ? finalHeavy.totalMass / finalHeavyCount : 0.0,
            ScalarFieldEnergy = scalarEnergy,
            HiggsFieldEnergy = scalarEnergy * 0.5
        };

        LastResult = new SimulationResult
        {
            AverageExcited = avgExcitedFinal,
            MaxExcited = maxExcited,
            MeasurementConfigured = false,
            MeasurementTriggered = false
        };
    }

    /// <summary>
    /// Cleans up all simulation resources
    /// </summary>
    public void Cleanup()
    {
        DisposeGpuEngines();
        SimulationEngine = null;
        GpuStats.Reset();
        OnConsoleLog?.Invoke("[Cleanup] ������� ��������� �����������.\n");
    }

    /// <summary>
    /// Creates decimated time series for compact export (max ~1000 points)
    /// </summary>
    public DecimatedDynamics GetDecimatedDynamics(int maxPoints = 1000)
    {
        int totalPoints = SeriesSteps.Count;
        if (totalPoints == 0)
            return new DecimatedDynamics();

        int stride = Math.Max(1, totalPoints / maxPoints);
        int resultCount = (totalPoints + stride - 1) / stride;

        int[] steps = new int[resultCount];
        int[] excited = new int[resultCount];
        double[] energy = new double[resultCount];
        double[] heavyMass = new double[resultCount];
        int[] largestCluster = new int[resultCount];
        int[] strongEdges = new int[resultCount];
        double[] spectralDim = new double[resultCount];
        double[] temp = new double[resultCount];

        for (int i = 0, j = 0; i < totalPoints && j < resultCount; i += stride, j++)
        {
            steps[j] = SeriesSteps[i];
            excited[j] = SeriesExcited[i];
            energy[j] = i < SeriesEnergy.Count ? SeriesEnergy[i] : 0;
            heavyMass[j] = i < SeriesHeavyMass.Count ? SeriesHeavyMass[i] : 0;
            largestCluster[j] = i < SeriesLargestCluster.Count ? SeriesLargestCluster[i] : 0;
            strongEdges[j] = i < SeriesStrongEdges.Count ? SeriesStrongEdges[i] : 0;
            spectralDim[j] = i < SeriesSpectralDimension.Count ? SeriesSpectralDimension[i] : 0;
            temp[j] = i < SeriesNetworkTemperature.Count ? SeriesNetworkTemperature[i] : 0;
        }

        return new DecimatedDynamics
        {
            TotalPoints = totalPoints,
            DecimationStride = stride,
            Steps = steps,
            Excited = excited,
            Energy = energy,
            HeavyMass = heavyMass,
            LargestCluster = largestCluster,
            StrongEdges = strongEdges,
            SpectralDimension = spectralDim,
            NetworkTemperature = temp
        };
    }

    /// <summary>
    /// Compact decimated time series for export
    /// </summary>
    public record DecimatedDynamics
    {
        public int TotalPoints { get; init; }
        public int DecimationStride { get; init; }
        public int[] Steps { get; init; } = [];
        public int[] Excited { get; init; } = [];
        public double[] Energy { get; init; } = [];
        public double[] HeavyMass { get; init; } = [];
        public int[] LargestCluster { get; init; } = [];
        public int[] StrongEdges { get; init; } = [];
        public double[] SpectralDimension { get; init; } = [];
        public double[] NetworkTemperature { get; init; } = [];
    }

    /// <summary>
    /// Creates a new simulation session
    /// </summary>
    public SimulationSession CreateSession(bool gpuEnabled, string? gpuDeviceName, DisplayFilters filters)
    {
        return new SimulationSession
        {
            SessionId = Guid.NewGuid(),
            StartedAt = DateTime.UtcNow,
            Config = LastConfig,
            GpuEnabled = gpuEnabled,
            GpuDeviceName = gpuDeviceName,
            Filters = filters
        };
    }

    /// <summary>
    /// Archives current session to history
    /// </summary>
    public void ArchiveSession(SessionEndReason reason, string consoleLog, string summaryText, List<ImportantEvent> events)
    {
        if (CurrentSession == null) return;

        double wallClockSeconds = (DateTime.UtcNow - SimulationWallClockStart).TotalSeconds;

        double finalSpectralDim = FinalSpectralDimension;
        double finalNetworkTemp = FinalNetworkTemperature;
        if (finalSpectralDim == 0 && SeriesSpectralDimension.Count > 0)
            finalSpectralDim = SeriesSpectralDimension[^1];
        if (finalNetworkTemp == 0 && SeriesNetworkTemperature.Count > 0)
            finalNetworkTemp = SeriesNetworkTemperature[^1];

        var session = CurrentSession with
        {
            EndedAt = DateTime.UtcNow,
            EndReason = reason,
            Config = LastConfig,
            Result = LastResult,
            ModernResult = ModernResult,
            SeriesSteps = [.. SeriesSteps],
            SeriesExcited = [.. SeriesExcited],
            SeriesHeavyMass = [.. SeriesHeavyMass],
            SeriesHeavyCount = [.. SeriesHeavyCount],
            SeriesLargestCluster = [.. SeriesLargestCluster],
            SeriesAvgDist = [.. SeriesAvgDist],
            SeriesDensity = [.. SeriesDensity],
            SeriesEnergy = [.. SeriesEnergy],
            SeriesCorr = [.. SeriesCorr],
            SeriesStrongEdges = [.. SeriesStrongEdges],
            SeriesQNorm = [.. SeriesQNorm],
            SeriesQEnergy = [.. SeriesQEnergy],
            SeriesEntanglement = [.. SeriesEntanglement],
            SeriesSpectralDimension = [.. SeriesSpectralDimension],
            SeriesNetworkTemperature = [.. SeriesNetworkTemperature],
            SeriesEffectiveG = [.. SeriesEffectiveG],
            SeriesAdaptiveThreshold = [.. SeriesAdaptiveThreshold],
            SynthesisData = SynthesisData?.ToList(),
            SynthesisCount = SynthesisCount,
            FissionCount = FissionCount,
            ConsoleLog = consoleLog,
            SummaryText = summaryText,
            LastStep = SeriesSteps.Count > 0 ? SeriesSteps[^1] : 0,
            TotalStepsPlanned = LastConfig?.TotalSteps ?? 0,
            FinalSpectralDimension = finalSpectralDim,
            FinalNetworkTemperature = finalNetworkTemp,
            WallClockDurationSeconds = wallClockSeconds,
            ImportantEvents = events
        };

        SessionHistory.Add(session);
    }

    // === Energy Ledger for Conservation Tracking ===
    public EnergyLedger EnergyLedger { get; } = new();

    /// <summary>
    /// RQ-COMPLIANT: Event-based simulation loop.
    /// 
    /// Uses priority queue where each node has its own proper time ?_i.
    /// Time dilation from gravity affects node update rate.
    /// This is the TRUE relational time evolution where:
    /// - There is no global "now"
    /// - Each node evolves according to local proper time
    /// - Heavy/curved regions run slower (GR time dilation)
    /// 
    /// RQ-Hypothesis Principle: Time emerges from quantum correlations,
    /// not from external parameter.
    /// </summary>
    /// <param name="ct">Cancellation token</param>
    /// <param name="totalEvents">Total events to process (mapped to UI "steps")</param>
    /// <param name="useGpu">Whether to use GPU acceleration</param>
    public void RunEventBasedLoop(CancellationToken ct, int totalEvents, bool useGpu)
    {
        if (SimulationEngine?.Graph == null)
            throw new InvalidOperationException("SimulationEngine not initialized");

        var graph = SimulationEngine.Graph;
        int nodeCount = graph.N;

        // Initialize asynchronous time system
        graph.InitAsynchronousTime();

        // Initialize energy ledger
        EnergyLedger.Initialize(graph.ComputeTotalEnergyUnified());
        OnConsoleLog?.Invoke($"[EventBased] Energy ledger initialized: E_0 = {EnergyLedger.TrackedEnergy:F4}\n");

        // Clear time series
        ClearTimeSeries();

        // Events processed counter (for UI "step" mapping)
        int eventsProcessed = 0;
        int metricsInterval = Math.Max(100, nodeCount); // Collect metrics every ~N events (one "sweep")
        int spectralDimInterval = nodeCount * 20; // Every ~20 sweeps
        int energyValidationInterval = nodeCount * 10; // Every ~10 sweeps

        // NOTE: Initialize to 0 to indicate "not yet computed"
        // UI will show 0.000 until first actual d_S computation
        double lastSpectralDimension = 0.0;
        List<int> excitedHistory = [];

        OnConsoleLog?.Invoke($"[EventBased] Starting event-driven simulation: {totalEvents} events, {nodeCount} nodes\n");
        OnConsoleLog?.Invoke($"[EventBased] Metrics interval: {metricsInterval}, Spectral: {spectralDimInterval}\n");

        while (!ct.IsCancellationRequested && eventsProcessed < totalEvents)
        {
            // Process batch of events (10 at a time for efficiency)
            int batchSize = Math.Min(10, totalEvents - eventsProcessed);
            graph.StepEventBasedBatch(batchSize);
            eventsProcessed += batchSize;

            // Map events to equivalent "step" for UI compatibility
            // One "step" ? N events (one full sweep of all nodes on average)
            int equivalentStep = eventsProcessed / Math.Max(1, nodeCount);

            // Lightweight metrics every step
            int excitedCount = CollectExcitedCount();
            excitedHistory.Add(excitedCount);

            // Full metrics collection at intervals
            if (eventsProcessed % metricsInterval == 0)
            {
                var metrics = CollectMetrics();
                double threshold = Math.Min(graph.GetAdaptiveHeavyThreshold(), RQGraph.HeavyClusterThreshold);

                // === FIX: Update NetworkTemperature by annealing schedule ===
                // In Event-Based mode, temperature must be computed and updated explicitly
                // to show proper cooling during simulation (Big Bang → Cold Universe)
                double startTemp = LiveConfig.HotStartTemperature;
                double currentTemp = ComputeAnnealingTemperature(equivalentStep, startTemp, totalEvents / Math.Max(1, nodeCount));
                graph.NetworkTemperature = currentTemp;

                // Effective G from live config
                double effectiveG = LiveConfig.GravitationalCoupling;

                StoreMetrics(equivalentStep, metrics.excited, metrics.heavyMass, metrics.heavyCount,
                    metrics.largestCluster, metrics.energy, metrics.strongEdges, metrics.correlation,
                    metrics.qNorm, metrics.entanglement, lastSpectralDimension, currentTemp, effectiveG, threshold,
                    metrics.totalClusters, metrics.avgClusterMass, metrics.maxClusterMass, metrics.avgDegree);

                // Auto-tuning (if enabled)
                string? tuneResult = PerformAutoTuning(
                    equivalentStep, lastSpectralDimension, metrics.excited, metrics.totalClusters,
                    metrics.largestCluster, metrics.heavyMass, nodeCount);

                if (tuneResult != null)
                {
                    OnConsoleLog?.Invoke($"[AutoTune] Event {eventsProcessed}: {tuneResult}\n");
                }
            }

            // Spectral dimension computation (expensive, less frequent)
            if (eventsProcessed % spectralDimInterval == 0 && eventsProcessed > 0)
            {
                double newDim = graph.ComputeSpectralDimension(t_max: 100, num_walkers: 50);

                if (newDim > 0 && !double.IsNaN(newDim))
                {
                    lastSpectralDimension = newDim;

                    // Check graph health
                    var clusters = graph.GetStrongCorrelationClusters(graph.GetAdaptiveHeavyThreshold());
                    int largestClusterSize = clusters.Count > 0 ? clusters.Max(c => c.Count) : 0;
                    var healthStatus = graph.CheckGraphHealth(newDim, largestClusterSize);

                    OnConsoleLog?.Invoke($"[d_S] Event {eventsProcessed}: {healthStatus.StatusDescription}\n");

                    // Fragmentation check
                    if (healthStatus.IsFragmented)
                    {
                        try
                        {
                            graph.CheckFragmentationTerminal(equivalentStep, newDim);
                            string recovery = graph.PerformGraphRecovery(healthStatus);
                            OnConsoleLog?.Invoke($"[FRAGMENTATION RECOVERY] {recovery}\n");
                        }
                        catch (GraphFragmentationException ex)
                        {
                            OnConsoleLog?.Invoke($"[FRAGMENTATION TERMINAL] {ex.Message}\n");
                            FinalizeSimulation(excitedHistory);
                            throw;
                        }
                    }
                }
            }

            // Energy conservation validation
            if (eventsProcessed % energyValidationInterval == 0 && eventsProcessed > 0)
            {
                try
                {
                    double currentEnergy = graph.ComputeTotalEnergyUnified();
                    EnergyLedger.ValidateConservation(currentEnergy);
                }
                catch (EnergyConservationException ex)
                {
                    OnConsoleLog?.Invoke($"[ENERGY] {ex.Message}\n");
                    // Don't throw - just log the violation
                }
            }

            // Progress logging every ~10000 events
            if (eventsProcessed % 10000 == 0 && eventsProcessed > 0)
            {
                double elapsed = (DateTime.UtcNow - SimulationWallClockStart).TotalSeconds;
                double eventsPerSec = eventsProcessed / elapsed;
                double globalTime = graph.GlobalTime;
                OnConsoleLog?.Invoke($"[EventBased] Events: {eventsProcessed}/{totalEvents}, " +
                    $"?_global={globalTime:F3}, speed={eventsPerSec:F0} ev/s, d_S={lastSpectralDimension:F2}\n");
            }
        }

        // Finalize
        FinalizeSimulation(excitedHistory);

        // Final energy check
        double finalEnergy = graph.ComputeTotalEnergyUnified();
        double energyDrift = Math.Abs(finalEnergy - EnergyLedger.TrackedEnergy);
        OnConsoleLog?.Invoke($"[EventBased] Simulation complete: {eventsProcessed} events, " +
            $"d_S={FinalSpectralDimension:F3}, E_drift={energyDrift:F6}\n");
    }

    // === Parallel Event Engine for Multi-Threaded Simulation ===
    private ParallelEventEngine? _parallelEngine;

    /// <summary>
    /// Configured CPU thread count from UI (numericUpDown1).
    /// Used for all Parallel.For operations to respect user preference.
    /// </summary>
    public int CpuThreadCount { get; set; } = Environment.ProcessorCount;

    /// <summary>
    /// Creates ParallelOptions with the configured CPU thread count.
    /// Use this for all Parallel.For calls to respect UI setting.
    /// </summary>
    public ParallelOptions CreateParallelOptions()
    {
        return new ParallelOptions
        {
            MaxDegreeOfParallelism = Math.Max(1, CpuThreadCount)
        };
    }

    /// <summary>
    /// RQ-COMPLIANT: Parallel event-based simulation loop.
    /// 
    /// Uses graph coloring to identify causally independent node groups,
    /// then processes each group in parallel with work-stealing thread pool.
    /// 
    /// Key RQ-Hypothesis insight:
    /// - No global "now" means nodes without causal connection can evolve independently
    /// - Causal independence = no shared edges AND no shared neighbors
    /// - Graph coloring partitions nodes into such independent sets
    /// 
    /// NOTE ON TIME: The "equivalentStep" and "eventsProcessed" are UI-ONLY counters
    /// for progress display. They do NOT represent physical time - each node has its
    /// own proper time ?_i tracked by the graph's asynchronous time system.
    /// 
    /// Performance benefits:
    /// - Thread pool is created once, reused across all sweeps
    /// - Work-stealing balances load between threads
    /// - Barrier sync between color groups ensures causality
    /// - Typical speedup: 2-4x on multi-core systems
    /// - GPU acceleration for gravity/curvature computation (10-50x for large graphs)
    /// </summary>
    /// <param name="ct">Cancellation token</param>
    /// <param name="totalEvents">Total events to process (UI progress counter only)</param>
    /// <param name="useParallel">Enable multi-threaded processing</param>
    /// <param name="useGpu">Enable GPU acceleration for gravity computation</param>
    public void RunParallelEventBasedLoop(CancellationToken ct, int totalEvents, bool useParallel = true, bool useGpu = false)
    {
        if (SimulationEngine?.Graph == null)
            throw new InvalidOperationException("SimulationEngine not initialized");

        var graph = SimulationEngine.Graph;
        int nodeCount = graph.N;

        // Initialize asynchronous time system (proper time per node)
        graph.InitAsynchronousTime();

        // Initialize energy ledger
        EnergyLedger.Initialize(graph.ComputeTotalEnergyUnified());
        OnConsoleLog?.Invoke($"[ParallelEvent] Energy ledger initialized: E_0 = {EnergyLedger.TrackedEnergy:F4}\n");

        // Initialize GPU engines if requested
        bool gpuActive = false;
        if (useGpu && GpuAvailable)
        {
            gpuActive = InitializeGpuEngines();
            if (gpuActive)
            {
                GpuStats.IsActive = true;
                OnConsoleLog?.Invoke($"[ParallelEvent] GPU acceleration ENABLED for gravity/curvature\n");
            }
            else
            {
                OnConsoleLog?.Invoke($"[ParallelEvent] GPU init failed, using CPU fallback\n");
            }
        }

        // Initialize parallel engine if requested
        if (useParallel)
        {
            _parallelEngine?.Dispose();
            _parallelEngine = new ParallelEventEngine(graph);
            _parallelEngine.ComputeGraphColoring();
            OnConsoleLog?.Invoke($"[ParallelEvent] Workers: {_parallelEngine.WorkerCount}, Colors: {_parallelEngine.ColorCount}\n");
        }

        // Clear time series
        ClearTimeSeries();

        // CRITICAL FIX: Set LiveTotalSteps to number of SWEEPS (not events)
        // UI timer reads LiveTotalSteps and compares to LiveStep (which is sweepCount)
        int totalSweeps = totalEvents / Math.Max(1, nodeCount);
        Dispatcher.LiveTotalSteps = totalSweeps;

        // UI progress counter (NOT physical time!)
        int eventsProcessed = 0;

        // FIXED: Collect metrics EVERY SWEEP for responsive UI
        // metricsInterval in events = nodeCount (one sweep)
        int metricsInterval = nodeCount; // Every sweep exactly
        int spectralDimInterval = Math.Max(nodeCount * 10, 2000); // Every 10 sweeps
        int energyValidationInterval = Math.Max(nodeCount * 5, 1000); // Every 5 sweeps
        int recolorInterval = Math.Max(nodeCount * 20, 5000); // Recompute coloring every 20 sweeps
        int progressLogInterval = Math.Max(nodeCount * 5, 1000); // Log progress every 5 sweeps
        int gravityInterval = 5; // Run GPU gravity every 5 sweeps (batched for efficiency)
        int gpuSyncInterval = 50; // Sync GPU weights back every 50 sweeps

        // NOTE: Initialize to 0 to indicate "not yet computed"
        // UI will show 0.000 until first actual d_S computation
        double lastSpectralDimension = 0.0;
        List<int> excitedHistory = [];
        double dt = 0.01;

        // Track sweep count for UI "step" display
        int sweepCount = 0;

        OnConsoleLog?.Invoke($"[ParallelEvent] Starting: {totalEvents} events, {nodeCount} nodes, parallel={useParallel}, gpu={gpuActive}\n");
        OnConsoleLog?.Invoke($"[ParallelEvent] Intervals: metrics={metricsInterval}, spectral={spectralDimInterval}, progress={progressLogInterval}\n");

        while (!ct.IsCancellationRequested && eventsProcessed < totalEvents)
        {
            int batchEvents;

            if (useParallel && _parallelEngine != null)
            {
                // Parallel sweep: process all nodes once, grouped by color
                batchEvents = _parallelEngine.ProcessParallelSweep(dt);
                eventsProcessed += batchEvents;
                sweepCount++;

                // Recompute coloring periodically (topology may change)
                if (eventsProcessed % recolorInterval == 0 && _parallelEngine.NeedsRecoloring)
                {
                    _parallelEngine.ComputeGraphColoring();
                    OnConsoleLog?.Invoke($"[ParallelEvent] Recolored: {_parallelEngine.ColorCount} colors\n");
                }
            }
            else
            {
                // Sequential fallback
                batchEvents = Math.Min(nodeCount, totalEvents - eventsProcessed);
                graph.StepEventBasedBatch(batchEvents);
                eventsProcessed += batchEvents;
                sweepCount++;
            }

            // === GPU GRAVITY STEP (batched for efficiency) ===
            // Run gravity computation on GPU every N sweeps
            if (gpuActive && sweepCount % gravityInterval == 0)
            {
                double effectiveG = LiveConfig.GravitationalCoupling;

                if (OptimizedGpuEngine != null)
                {
                    // Batch GPU computation (5 steps at once for reduced overhead)
                    OptimizedGpuEngine.StepGpuBatch(
                        batchSize: gravityInterval,
                        dt: (float)dt,
                        G: (float)effectiveG,
                        lambda: (float)PhysicsConstants.CosmologicalConstant,
                        degreePenalty: (float)PhysicsConstants.DegreePenaltyFactor,
                        diffusionRate: (float)PhysicsConstants.FieldDiffusionRate,
                        scalarMass: (float)PhysicsConstants.KleinGordonMass);

                    // Track GPU kernel launches (4 shaders × batchSize)
                    GpuStats.RecordKernelLaunch(4 * gravityInterval);
                }
                else if (graph.IsGpuGravityActive())
                {
                    // Use standalone GPU gravity engine
                    graph.EvolveGravityGpuBatch(
                        batchSize: gravityInterval,
                        dt: dt,
                        G: effectiveG,
                        lambda: PhysicsConstants.CosmologicalConstant,
                        degreePenaltyFactor: PhysicsConstants.DegreePenaltyFactor);
                    GpuStats.RecordKernelLaunch(2 * gravityInterval);
                }

                // Sync weights back to CPU graph periodically (expensive, do less often)
                if (sweepCount % gpuSyncInterval == 0)
                {
                    OptimizedGpuEngine?.SyncWeightsToGraph();
                    OptimizedGpuEngine?.SyncScalarFieldToGraph();
                    GpuStats.RecordWeightSync();
                }
            }
            else if (!gpuActive && sweepCount % gravityInterval == 0)
            {
                // CPU gravity fallback (less frequent for performance)
                double effectiveG = LiveConfig.GravitationalCoupling;
                if (effectiveG > 0)
                {
                    ImprovedNetworkGravity.EvolveNetworkGeometryOllivierDynamic(graph, dt, effectiveG);
                }
            }

            // UI "step" counter - NOT physical time, just for display
            // Maps sweep count to "equivalent step" for UI compatibility
            int equivalentStep = sweepCount;

            // Lightweight metrics every sweep (for excitedHistory tracking)
            int excitedCount = CollectExcitedCount();
            excitedHistory.Add(excitedCount);

            // Full metrics collection - THROTTLED to reduce lock contention
            // Store every 20 sweeps to reduce CPU overhead from CollectMetrics()
            // UI will display slightly delayed data but simulation runs faster
            int storeInterval = 20; // Reduced from 5 to 20 for performance
            bool shouldStoreMetrics = (sweepCount % storeInterval == 0) || sweepCount <= 5;

            if (shouldStoreMetrics)
            {
                // FIX: Update correlation weights less frequently (expensive O(E) operation)
                // Only update every 50 sweeps instead of every store to reduce CPU load
                bool shouldUpdateWeights = (sweepCount % 50 == 0) || sweepCount <= 5;
                if (shouldUpdateWeights)
                {
                    graph.UpdateCorrelationWeights();
                }

                var metrics = CollectMetrics();
                double threshold = Math.Min(graph.GetAdaptiveHeavyThreshold(), RQGraph.HeavyClusterThreshold);

                // === FIX: Update NetworkTemperature by annealing schedule ===
                // In Event-Based mode, temperature must be computed and updated explicitly
                // to show proper cooling during simulation (Big Bang → Cold Universe)
                double startTemp = LiveConfig.HotStartTemperature;
                double currentTemp = ComputeAnnealingTemperature(equivalentStep, startTemp, totalSweeps);
                graph.NetworkTemperature = currentTemp;

                double effectiveG = LiveConfig.GravitationalCoupling;

                // Compute topology metrics for Summary tab (edge/component counts)
                // Edge count: count undirected edges (i<j)
                int edgeCount = 0;
                for (int i = 0; i < nodeCount; i++)
                {
                    foreach (int j in graph.Neighbors(i))
                    {
                        if (j > i) edgeCount++;
                    }
                }
                // Component count via union-find over current topology
                int componentCount = 0;
                {
                    int[] parent = new int[nodeCount];
                    int[] rank = new int[nodeCount];
                    for (int i = 0; i < nodeCount; i++) { parent[i] = i; rank[i] = 0; }
                    int Find(int x)
                    {
                        while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
                        return x;
                    }
                    void Union(int a, int b)
                    {
                        int pa = Find(a), pb = Find(b);
                        if (pa == pb) return;
                        if (rank[pa] < rank[pb]) parent[pa] = pb;
                        else if (rank[pa] > rank[pb]) parent[pb] = pa;
                        else { parent[pb] = pa; rank[pa]++; }
                    }
                    for (int i = 0; i < nodeCount; i++)
                    {
                        foreach (int j in graph.Neighbors(i))
                        {
                            if (i < j) Union(i, j);
                        }
                    }
                    var seen = new HashSet<int>();
                    for (int i = 0; i < nodeCount; i++) seen.Add(Find(i));
                    componentCount = seen.Count;
                }

                // Store metrics - this updates Dispatcher for UI
                StoreMetrics(equivalentStep, metrics.excited, metrics.heavyMass, metrics.heavyCount,
                    metrics.largestCluster, metrics.energy, metrics.strongEdges, metrics.correlation,
                    metrics.qNorm, metrics.entanglement, lastSpectralDimension, currentTemp, effectiveG, threshold,
                    metrics.totalClusters, metrics.avgClusterMass, metrics.maxClusterMass, metrics.avgDegree,
                    edgeCount, componentCount);

                // Auto-tuning (if enabled)
                string? tuneResult = PerformAutoTuning(
                    equivalentStep, lastSpectralDimension, metrics.excited, metrics.totalClusters,
                    metrics.largestCluster, metrics.heavyMass, nodeCount);

                if (tuneResult != null)
                {
                    OnConsoleLog?.Invoke($"[AutoTune] Sweep {sweepCount}: {tuneResult}\n");
                }
            } // end shouldStoreMetrics

            // Spectral dimension (expensive, less frequent)
            // Use GPU if available for faster computation
            if (eventsProcessed % spectralDimInterval == 0 && eventsProcessed > 0)
            {
                double newDim;
                if (gpuActive && (GpuSpectralWalkEngine != null || GpuHeatKernelEngine != null))
                {
                    // GPU-accelerated spectral dimension with automatic method selection
                    // and topology synchronization
                    newDim = ComputeSpectralDimensionGpu(graph, enableCpuComparison: false);
                }
                else
                {
                    newDim = graph.ComputeSpectralDimension(t_max: 100, num_walkers: 50);
                }

                if (newDim > 0 && !double.IsNaN(newDim))
                {
                    lastSpectralDimension = newDim;
                    var clusters = graph.GetStrongCorrelationClusters(graph.GetAdaptiveHeavyThreshold());
                    int largestClusterSize = clusters.Count > 0 ? clusters.Max(c => c.Count) : 0;
                    var healthStatus = graph.CheckGraphHealth(newDim, largestClusterSize);
                    OnConsoleLog?.Invoke($"[d_S] Sweep {sweepCount}: {healthStatus.StatusDescription}\n");

                    if (healthStatus.IsFragmented)
                    {
                        try
                        {
                            graph.CheckFragmentationTerminal(equivalentStep, newDim);
                            string recovery = graph.PerformGraphRecovery(healthStatus);
                            OnConsoleLog?.Invoke($"[FRAGMENTATION RECOVERY] {recovery}\n");
                        }
                        catch (GraphFragmentationException ex)
                        {
                            OnConsoleLog?.Invoke($"[FRAGMENTATION TERMINAL] {ex.Message}\n");
                            FinalizeSimulation(excitedHistory);
                            throw;
                        }
                    }
                }
            }

            // Energy validation
            if (eventsProcessed % energyValidationInterval == 0 && eventsProcessed > 0)
            {
                try
                {
                    double currentEnergy = graph.ComputeTotalEnergyUnified();
                    EnergyLedger.ValidateConservation(currentEnergy);
                }
                catch (EnergyConservationException ex)
                {
                    OnConsoleLog?.Invoke($"[ENERGY] {ex.Message}\n");
                }
            }

            // Progress logging
            if (eventsProcessed % progressLogInterval == 0 && eventsProcessed > 0)
            {
                double elapsed = (DateTime.UtcNow - SimulationWallClockStart).TotalSeconds;
                double eventsPerSec = elapsed > 0 ? eventsProcessed / elapsed : 0;
                double sweepsPerSec = elapsed > 0 ? sweepCount / elapsed : 0;
                string parallelStats = _parallelEngine != null
                    ? $", parallel={_parallelEngine.ParallelUpdates / (double)(_parallelEngine.ParallelUpdates + _parallelEngine.SequentialUpdates + 1):P0}"
                    : "";
                string gpuStats = gpuActive ? $", GPU kernels={GpuStats.KernelLaunches}" : "";
                OnConsoleLog?.Invoke($"[ParallelEvent] Sweep {sweepCount}, Events: {eventsProcessed}/{totalEvents}, " +
                    $"speed={sweepsPerSec:F1} sweeps/s, d_S={lastSpectralDimension:F2}{parallelStats}{gpuStats}\n");
            }

            // CRITICAL FIX: Variable yielding to UI thread
            // Every 10 sweeps, brief sleep to allow UI to acquire _writeLock for buffer swap
            // Less frequent than before (was every 5) for better simulation speed
            if (sweepCount % 10 == 0)
            {
                Thread.Sleep(5); // 5ms - brief yield for UI
            }
            else if (sweepCount % 50 == 0)
            {
                Thread.Sleep(1); // 1ms yield every 50 sweeps
            }
        }

        // Final GPU sync before finalization
        if (gpuActive && OptimizedGpuEngine != null)
        {
            OptimizedGpuEngine.SyncWeightsToGraph();
            OptimizedGpuEngine.SyncScalarFieldToGraph();
        }

        // Finalize
        FinalizeSimulation(excitedHistory);

        // Report parallel stats
        if (_parallelEngine != null)
        {
            OnConsoleLog?.Invoke($"[ParallelEvent] {_parallelEngine.GetStatsSummary()}\n");
            _parallelEngine.Dispose();
            _parallelEngine = null;
        }

        // Report GPU stats
        if (gpuActive)
        {
            OnConsoleLog?.Invoke($"[GPU] Total kernels: {GpuStats.KernelLaunches}, syncs: {GpuStats.WeightSyncs}, rebuilds: {GpuStats.TopologyRebuilds}\n");
            DisposeGpuEngines();
        }

        double finalEnergy = graph.ComputeTotalEnergyUnified();
        double energyDrift = Math.Abs(finalEnergy - EnergyLedger.TrackedEnergy);
        OnConsoleLog?.Invoke($"[ParallelEvent] Complete: sweeps={sweepCount}, d_S={FinalSpectralDimension:F3}, E_drift={energyDrift:F6}\n");
    }


}