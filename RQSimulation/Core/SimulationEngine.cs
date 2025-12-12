using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;

namespace RQSimulation
{
    /// <summary>
    /// Configuration for simulation parameters.
    /// Used by both legacy and modern simulation modes.
    /// </summary>
    public sealed class SimulationConfig
    {
        // === Basic Graph Parameters ===
        public int NodeCount { get; set; }
        public double InitialEdgeProb { get; set; }
        public double InitialExcitedProb { get; set; }
        public int TargetDegree { get; set; }
        public double LambdaState { get; set; }
        public double Temperature { get; set; } = PhysicsConstants.InitialAnnealingTemperature;
        public double EdgeTrialProbability { get; set; }
        public double MeasurementThreshold { get; set; }
        public double DynamicMeasurementThreshold { get; set; } = 0.9;
        public int Seed { get; set; }
        public int TotalSteps { get; set; }
        public int LogEvery { get; set; } = 1;
        public int BaselineWindow { get; set; } = 100;

        // === Legacy impulse parameters (disabled in Modern mode) ===
        public int FirstImpulse { get; set; } = -1;
        public int ImpulsePeriod { get; set; } = -1;
        public int CalibrationStep { get; set; } = -1;

        public int VisualizationInterval { get; set; } = 10;
        public int MeasurementLogInterval { get; set; } = 100;
        public bool UseQuantumDrivenStates { get; set; } = true;
        public int FractalLevels { get; set; } = 2;
        public int FractalBranchFactor { get; set; } = 2;
        public int StrongEdgeThreshold { get; set; } = 25;

        // === Physics Modules ===
        public bool UseSpacetimePhysics { get; set; } = true;
        public bool UseSpinorField { get; set; } = true;
        public bool UseVacuumFluctuations { get; set; } = true;
        public bool UseBlackHolePhysics { get; set; } = true;
        public bool UseYangMillsGauge { get; set; } = true;
        public bool UseEnhancedKleinGordon { get; set; } = true;

        // === Physics Constants (should use PhysicsConstants values) ===
        public double GravitationalCoupling { get; set; } = PhysicsConstants.GravitationalCoupling;
        public double VacuumEnergyScale { get; set; } = PhysicsConstants.VacuumFluctuationBaseRate * 20;
        public bool UseInternalTime { get; set; } = true;
        public bool UseSpectralGeometry { get; set; } = true;
        public bool UseQuantumGraphity { get; set; } = true;
        public double InitialNetworkTemperature { get; set; } = PhysicsConstants.InitialAnnealingTemperature;
        public double AnnealingCoolingRate { get; set; } = 0.995;

        // === Relational Dynamics ===
        public bool UseRelationalTime { get; set; } = true;
        public bool UseRelationalYangMills { get; set; } = true;
        public bool UseNetworkGravity { get; set; } = true;
        public double DecoherenceRate { get; set; } = 0.001;

        // === RQ-Hypothesis Compliance ===
        public bool UseUnifiedPhysicsStep { get; set; } = true;
        public bool EnforceGaugeConstraints { get; set; } = true;
        public bool UseCausalRewiring { get; set; } = true;
        public bool UseAsynchronousTime { get; set; } = false;
        public bool UseTopologicalProtection { get; set; } = true;
        public bool ValidateEnergyConservation { get; set; } = true;
        public double EnergyConservationTolerance { get; set; } = 0.01;

        /// <summary>
        /// RQ-COMPLIANT: Use event-based simulation with per-node proper time.
        /// When enabled, each node evolves according to its local proper time ?_i
        /// with gravitational time dilation (heavy regions run slower).
        /// 
        /// This is the TRUE relational time model where:
        /// - No global "now" exists
        /// - Time emerges from quantum correlations
        /// - Causality is enforced via graph connectivity
        /// 
        /// Default: false - controlled by UI checkbox "Events Model (time free)"
        /// </summary>
        public bool UseEventBasedSimulation { get; set; } = true;

        // === Absolute Rigor Options ===
        public bool UseMexicanHatPotential { get; set; } = true;
        public bool UseHotStartAnnealing { get; set; } = true;
        public bool DisableManualParticleDetection { get; set; } = true;
        public bool UseSoftPotentialWalls { get; set; } = true;
        public bool UseGeometryMomenta { get; set; } = true;
        public bool UseTopologicalCensorship { get; set; } = true;
        public bool UseSpectralMass { get; set; } = false;
        public double HotStartTemperature { get; set; } = PhysicsConstants.InitialAnnealingTemperature;
    }

    /// <summary>
    /// Result container for simulation metrics.
    /// </summary>
    public sealed class SimulationResult
    {
        public double AverageExcited { get; set; }
        public int MaxExcited { get; set; }
        public List<string> TimeSeries { get; } = new();
        public List<string> AvalancheStats { get; } = new();
        public List<string> MeasurementEvents { get; } = new();
        public List<string> HeavyClusters { get; } = new();
        public bool MeasurementConfigured { get; set; }
        public bool MeasurementTriggered { get; set; }
        public List<string> DiagnosticsExport { get; } = new();
    }

    /// <summary>
    /// Event args for console logging.
    /// </summary>
    public sealed class ConsoleLogEventArgs : EventArgs
    {
        public ConsoleLogEventArgs(string message) => Message = message;
        public string Message { get; }
    }

    /// <summary>
    /// Event args for progress updates (kept for compatibility, but minimally used).
    /// </summary>
    public sealed class SimulationProgressEventArgs : EventArgs
    {
        public SimulationProgressEventArgs(int currentStep, int totalSteps, int currentOn, bool shouldRedraw)
        {
            CurrentStep = currentStep;
            TotalSteps = totalSteps;
            CurrentOn = currentOn;
            ShouldRedraw = shouldRedraw;
        }
        public int CurrentStep { get; }
        public int TotalSteps { get; }
        public int CurrentOn { get; }
        public bool ShouldRedraw { get; }
    }

    /// <summary>
    /// Modern simulation engine - acts as a factory for RQGraph initialization.
    /// The actual simulation loop is in Form_Main.RunModernAsync().
    /// 
    /// ARCHITECTURE:
    /// - SimulationEngine: Creates and configures RQGraph
    /// - Form_Main.RunModernAsync(): Runs the simulation loop with GPU support
    /// - RQGraph: Contains all physics methods
    /// </summary>
    public sealed class SimulationEngine
    {
        private readonly SimulationConfig _cfg;
        private readonly RQGraph _graph;
        private Action<RQGraph>? _customInitializer;

        /// <summary>
        /// Creates a new simulation engine and initializes the graph.
        /// </summary>
        public SimulationEngine(SimulationConfig cfg) : this(cfg, null)
        {
        }

        /// <summary>
        /// Creates a new simulation engine with optional custom graph initializer.
        /// </summary>
        /// <param name="cfg">Simulation configuration</param>
        /// <param name="customInitializer">Optional action to customize graph topology after basic initialization</param>
        public SimulationEngine(SimulationConfig cfg, Action<RQGraph>? customInitializer)
        {
            ArgumentNullException.ThrowIfNull(cfg);
            _cfg = cfg;
            _customInitializer = customInitializer;

            // Create graph with initial topology
            _graph = new RQGraph(
                cfg.NodeCount,
                cfg.InitialEdgeProb,
                cfg.InitialExcitedProb,
                cfg.TargetDegree,
                cfg.LambdaState,
                cfg.Temperature,
                cfg.EdgeTrialProbability,
                cfg.MeasurementThreshold,
                cfg.Seed);

            // Apply custom initializer if provided (e.g., for DNA linear chain)
            // This is done BEFORE fractal topology to allow complete override
            if (_customInitializer != null)
            {
                LogConsole("[SimulationEngine] Applying custom graph initializer...\n");
                _customInitializer(_graph);
            }
            else
            {
                // Initialize fractal topology only if no custom initializer
                _graph.InitFractalTopology(cfg.FractalLevels, cfg.FractalBranchFactor);
            }

            _graph.ComputeFractalLevels();

            // Initialize coordinates for visualization
            _graph.InitCoordinatesRandom(range: 1.0);
            _graph.RelaxCoordinatesFromCorrelation(0.05);

            // Initialize quantum wavefunction
            if (_cfg.UseQuantumDrivenStates)
            {
                _graph.ConfigureQuantumComponents(1);
                _graph.InitQuantumWavefunction();
            }

            // Initialize physics modules
            InitializePhysicsModules();

            // Log initialization
            LogConsole("[SimulationEngine] Graph initialized\n");
            LogConsole($"[SimulationEngine] N={_graph.N}, edges={_graph.FlatEdgesFrom?.Length ?? 0}\n");
        }

        /// <summary>
        /// Sets a custom initializer to be applied to the next graph created.
        /// Used by experiment system to set up special topologies (e.g., DNA linear chain).
        /// </summary>
        public void SetCustomInitializer(Action<RQGraph>? initializer)
        {
            _customInitializer = initializer;
        }

        /// <summary>
        /// Gets the current custom initializer, if any.
        /// </summary>
        public Action<RQGraph>? CustomInitializer => _customInitializer;

        /// <summary>
        /// Initializes all physics modules based on configuration.
        /// </summary>
        private void InitializePhysicsModules()
        {
            if (_cfg.UseSpacetimePhysics)
            {
                _graph.InitSpacetimeCoordinates();
            }

            if (_cfg.UseYangMillsGauge || _cfg.UseSpinorField)
            {
                _graph.InitSpinorField(0.01);
            }

            if (_cfg.UseVacuumFluctuations)
            {
                _graph.InitVacuumField();
            }

            if (_cfg.UseBlackHolePhysics)
            {
                _graph.InitBlackHolePhysics();
            }

            if (_cfg.UseYangMillsGauge)
            {
                _graph.InitYangMillsFields();
                _graph.InitEdgeGaugePhases();
            }

            if (_cfg.UseEnhancedKleinGordon)
            {
                _graph.InitEnhancedKleinGordon(0.01);
            }

            if (_cfg.UseInternalTime)
            {
                _graph.InitClockSubsystem(0.05);
            }

            if (_cfg.UseRelationalTime)
            {
                int clockSize = Math.Max(2, _graph.N / 20);
                _graph.InitInternalClock(clockSize);
            }

            if (_cfg.UseSpectralGeometry)
            {
                _graph.UpdateSpectralCoordinates();
                _graph.SyncCoordinatesFromSpectral();
            }

            if (_cfg.UseQuantumGraphity)
            {
                _graph.NetworkTemperature = _cfg.InitialNetworkTemperature;
            }

            if (_cfg.UseAsynchronousTime)
            {
                _graph.InitAsynchronousTime();
            }

            if (_cfg.UseUnifiedPhysicsStep)
            {
                // RQ-FIX: Asynchronous updates are now always enabled (RQ-compliant default)
                // UseAsynchronousUpdates property has been removed - async mode is the only option
                _graph.EnforceGaugeConstraintsEnabled = _cfg.EnforceGaugeConstraints;
                _graph.ValidateEnergyConservationEnabled = _cfg.ValidateEnergyConservation;
            }

            // Mexican Hat (Higgs) potential
            if (_cfg.UseMexicanHatPotential)
            {
                if (_cfg.UseHotStartAnnealing)
                {
                    _graph.InitScalarFieldHotStart(_cfg.HotStartTemperature);
                }
                else
                {
                    _graph.InitScalarFieldMexicanHat(0.01);
                }
            }

            // Geometry momenta for gravitational waves
            if (_cfg.UseGeometryMomenta)
            {
                _graph.InitGeometryMomenta();
            }
        }

        /// <summary>
        /// The initialized graph, ready for simulation.
        /// </summary>
        public RQGraph Graph => _graph;

        /// <summary>
        /// Event for console logging.
        /// </summary>
        public event EventHandler<ConsoleLogEventArgs>? ConsoleLog;


        private void LogConsole(string msg) => ConsoleLog?.Invoke(this, new ConsoleLogEventArgs(msg));

        /// <summary>
        /// Runs the event-driven simulation loop.
        /// Replaces the synchronous RunModernSync loop with a relativistic event-based approach.
        /// </summary>
        public void RunEventDrivenLoop()
        {
            // Priority queue: (Time, NodeIndex). Sorted by time.
            var eventQueue = new PriorityQueue<int, double>();

            // Initialize: each node schedules its first update
            // ProperTime is initialized to 0 for all nodes
            for (int i = 0; i < _graph.N; i++)
            {
                // Ensure ProperTime array exists and is initialized
                if (_graph.ProperTime == null || _graph.ProperTime.Length != _graph.N)
                    _graph.InitAsynchronousTime();
                    
                eventQueue.Enqueue(i, _graph.ProperTime[i]);
            }

            // Simulation loop
            // Note: _isRunning flag should be managed by the caller or added to this class
            // For now, we assume a cancellation token or external flag would be used in a real scenario
            // Here we just run for the configured TotalSteps (interpreted as events or time)
            
            int maxEvents = _cfg.TotalSteps * _graph.N; // Approximate equivalent work
            int eventsProcessed = 0;

            while (eventsProcessed < maxEvents && eventQueue.Count > 0)
            {
                // 1. Dequeue event with minimum local time
                if (!eventQueue.TryDequeue(out int nodeIndex, out double executionTime)) break;

                // 2. Update ONLY this node (locality)
                // Use the event-driven update method
                _graph.UpdateNodePhysics(nodeIndex, 0.01); // 0.01 is base dt, scaled internally

                // 3. Compute next event time (Lapse function)
                // dtau = N(x) * dt_coord
                // GetTimeDilation returns 1/N(x) roughly, or we use GetLocalLapse if available
                // The user's snippet used GetLocalLapse, but existing code has GetTimeDilation
                // We'll use GetTimeDilation which is already in RQGraph.EventDrivenExtensions.cs
                
                double timeDilation = _graph.GetTimeDilation(nodeIndex);
                double nextTime = executionTime + 0.01 * timeDilation; 
                
                _graph.ProperTime[nodeIndex] = nextTime;
                eventQueue.Enqueue(nodeIndex, nextTime);
                
                eventsProcessed++;
            }
        }

        // === Modern Physics Step State ===
        private PriorityQueue<int, double> _eventQueue = new();
        private double _globalTime = 0;
        private int _step = 0;
        private const double _dt = 0.1; // Base coordinate time step
        private const int TopologyUpdateInterval = 50;

        /// <summary>
        /// Modern physics step integrating all RQ-Hypothesis fixes.
        /// 1. Energy Conservation (Ledger)
        /// 2. Event-Driven Time (Local Lapse)
        /// 3. Local Action (Fields + Geometry)
        /// 4. Topological Protection (Bipartite)
        /// </summary>
        public void RunPhysicsStepModern()
        {
            // Initialize event queue if empty (first run)
            if (_eventQueue.Count == 0 && _step == 0)
            {
                for (int i = 0; i < _graph.N; i++)
                {
                    _eventQueue.Enqueue(i, _graph.ProperTime != null && i < _graph.ProperTime.Length ? _graph.ProperTime[i] : 0.0);
                }
                // Initialize Ledger if needed
                if (Math.Abs(_graph.Ledger.TotalTrackedEnergy) < 1e-10)
                {
                    _graph.Ledger.Initialize(_graph.CalculateTotalEnergy());
                }
            }

            // 1. Energy Integrity Check
            if (_cfg.ValidateEnergyConservation)
            {
                double currentEnergy = _graph.CalculateTotalEnergy();
                // Only validate if initialized
                if (Math.Abs(_graph.Ledger.TotalTrackedEnergy) > 1e-10)
                {
                    try 
                    {
                        _graph.Ledger.ValidateConservation(currentEnergy);
                    }
                    catch (EnergyConservationException ex)
                    {
                        LogConsole($"[ENERGY ERROR] {ex.Message}");
                        // Optionally halt or correct
                    }
                }
            }

            // 2. Asynchronous Update Loop (Event-Driven)
            // Process events up to _globalTime + _dt
            double targetTime = _globalTime + _dt;
            
            while (_eventQueue.TryPeek(out int node, out double time) && time <= targetTime)
            {
                _eventQueue.Dequeue();
                
                // a) Local Field Update (Yang-Mills, Dirac, Scalar)
                // Uses EnergyLedger internally for fluctuations
                _graph.UpdateLocalFields(node, _dt);
                
                // b) Local Geometry Update (Ollivier-Ricci + Volume Constraint)
                _graph.UpdateLocalGeometry(node, _dt);
                
                // c) Schedule Next Event
                // dtau = N(x) * dt_coord
                double lapse = _graph.GetLocalLapse(node);
                double nextTime = time + 0.01 * lapse; // 0.01 is base proper time step
                
                if (_graph.ProperTime != null && node < _graph.ProperTime.Length)
                {
                    _graph.ProperTime[node] = nextTime;
                }
                
                _eventQueue.Enqueue(node, nextTime);
            }

            // 3. Topological Updates (Quantum Graphity)
            // Includes Bipartite Check for Fermion consistency
            if (_step % TopologyUpdateInterval == 0)
            {
                _graph.MetropolisTopologyUpdate_Checked();
            }

            _globalTime += _dt;
            _step++;
        }
    }
}


