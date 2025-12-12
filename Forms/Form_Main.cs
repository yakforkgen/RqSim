using iSukces;
using iSukces.DrawingPanel;
using RQSimulation;
using RQSimulation.GPUOptimized;
using RQSimulation.Experiments;
using System.Text;
using System.IO;
using System.Linq;
using System.Windows.Forms;
using System.IO.Compression;
using System.Drawing; // added
using System.Threading; // added
using System.Threading.Tasks; // added
using System.Collections.Generic; // added
using System.Globalization; // added for parsing threshold
using System.Text.Json;
using RqSimForms.Forms.Interfaces; // added for JSON export

namespace RqSimForms
{
    public partial class Form_Main : Form
    {
        private DoubleBufferedPanel drawingPanel;
        private Bitmap? canvasBitmap;

        // === Simulation API (contains all non-UI logic) ===
        private readonly FormSimAPI _simApi = new();

        // Shortcut accessors for commonly used properties
        private SimulationEngine? _simulationEngine => _simApi.SimulationEngine;
        private List<int> _seriesSteps => _simApi.SeriesSteps;
        private List<int> _seriesExcited => _simApi.SeriesExcited;
        private List<double> _seriesHeavyMass => _simApi.SeriesHeavyMass;
        private List<int> _seriesLargestCluster => _simApi.SeriesLargestCluster;
        private List<double> _seriesEnergy => _simApi.SeriesEnergy;
        private List<int> _seriesHeavyCount => _simApi.SeriesHeavyCount;
        private List<int> _seriesStrongEdges => _simApi.SeriesStrongEdges;
        private List<double> _seriesCorr => _simApi.SeriesCorr;
        private List<double> _seriesQNorm => _simApi.SeriesQNorm;
        private List<double> _seriesEntanglement => _simApi.SeriesEntanglement;
        private List<double> _seriesSpectralDimension => _simApi.SeriesSpectralDimension;
        private List<double> _seriesNetworkTemperature => _simApi.SeriesNetworkTemperature;
        private List<double> _seriesEffectiveG => _simApi.SeriesEffectiveG;
        private List<double> _seriesAdaptiveThreshold => _simApi.SeriesAdaptiveThreshold;

        // Throttling mechanism for UI updates
        private DateTime _lastUiUpdate = DateTime.MinValue;
        private readonly TimeSpan _minUiUpdateInterval = TimeSpan.FromMilliseconds(250);

        // Separate throttle for chart updates (less frequent)
        private DateTime _lastChartUpdate = DateTime.MinValue;
        private readonly TimeSpan _minChartUpdateInterval = TimeSpan.FromMilliseconds(500);

        // Drawing lock to prevent multiple concurrent drawing tasks
        private int _isDrawing = 0;

        // UI Timer for periodic updates (runs on UI thread)
        private System.Windows.Forms.Timer? _uiUpdateTimer;
        private CancellationTokenSource? _modernCts;

        private PointF[]? _cachedNodePositions; // cached for hit test
        private bool _useDynamicCoords = false; // toggle dynamic layout; circle if false
        private double _displayWeightThreshold = 0.0; // текущий порог веса для отображения рёбер
        private bool _displayShowHeavyOnly = false; // режим показа только тяжёлых кластеров

        // === Graph drawing throttling (for Event-Based performance) ===
        private DateTime _lastGraphDrawTime = DateTime.MinValue;
        private readonly TimeSpan _graphDrawInterval = TimeSpan.FromSeconds(5); // Draw graph every 5 sec

        // === Graph Health UI Controls (RQ-Hypothesis compliance) ===
        private NumericUpDown numGiantClusterThreshold = null!;
        private NumericUpDown numEmergencyGiantClusterThreshold = null!;
        private NumericUpDown numGiantClusterDecoherenceRate = null!;
        private NumericUpDown numMaxDecoherenceEdgesFraction = null!;
        private NumericUpDown numCriticalSpectralDimension = null!;
        private NumericUpDown numWarningSpectralDimension = null!;

        private ConsoleBuffer? _consoleBuffer;

        public Form_Main()
        {
            InitializeComponent();

            drawingPanel = new DoubleBufferedPanel();
            drawingPanel.Dock = DockStyle.Fill;
            drawingPanel.BackColor = Color.White;
            drawingPanel.Paint += DrawingPanel_Paint;
            drawingPanel.MouseClick += DrawingPanel_MouseClick;

            // Добавляем drawingPanel в tabPage_GUI
            tabPage_GUI.Controls.Add(drawingPanel);

            // Chart panels paint handlers
            panelOnChart.Paint += PanelOnChart_Paint;
            panelHeavyChart.Paint += PanelHeavyChart_Paint;
            panelClusterChart.Paint += PanelClusterChart_Paint;
            panelEnergyChart.Paint += PanelEnergyChart_Paint;

            // TextBox для CSV вывода (already exists in Designer)
            // Настройка вкладки Synthesis Analysis
            InitializeSynthesisTab();

            // Инициализация GPU контролов
            InitializeGpuControls();

            // Инициализация контролов Graph Health (RQ-Hypothesis compliance)
            InitializeGraphHealthControls();

            // Инициализация пресетов симуляции
            InitializePresets();

            // Initialize experiment selection (sidebar combobox)
            InitializeExperiments();

            // Initialize Experiments tab (full-featured UI)
            InitializeExperimentsTab();

            // Setup console log callback
            _simApi.OnConsoleLog = msg => BeginInvoke(new Action(() => AppendConsole(msg)));

            // Initialize MaxFPS control to match default timer interval
            numericUpDown_MaxFPS.Value = 10; // 10 FPS = 100ms interval

            // Wire up live parameter change handlers for runtime updates
            WireLiveParameterHandlers();

            _consoleBuffer = new ConsoleBuffer(consoleTextBox, chkAutoScrollConsole);
        }

        /// <summary>
        /// Connects ValueChanged handlers to all NumericUpDown controls in grpSimParams and grpPhysicsConstants
        /// for live parameter updates during simulation run.
        /// </summary>
        private void WireLiveParameterHandlers()
        {
            // grpPhysicsConstants - NumericUpDown controls
            numInitialEdgeProb.ValueChanged += OnLiveParameterChanged;
            numGravitationalCoupling.ValueChanged += OnLiveParameterChanged;
            numVacuumEnergyScale.ValueChanged += OnLiveParameterChanged;
            numAnnealingCoolingRate.ValueChanged += OnLiveParameterChanged;
            numDecoherenceRate.ValueChanged += OnLiveParameterChanged;
            numHotStartTemperature.ValueChanged += OnLiveParameterChanged;
            numAdaptiveThresholdSigma.ValueChanged += OnLiveParameterChanged;
            numWarmupDuration.ValueChanged += OnLiveParameterChanged;
            numGravityTransitionDuration.ValueChanged += OnLiveParameterChanged;

            // grpSimParams - NumericUpDown controls
            numNodeCount.ValueChanged += OnLiveParameterChanged;
            numTargetDegree.ValueChanged += OnLiveParameterChanged;
            numInitialExcitedProb.ValueChanged += OnLiveParameterChanged;
            numLambdaState.ValueChanged += OnLiveParameterChanged;
            numTemperature.ValueChanged += OnLiveParameterChanged;
            numEdgeTrialProb.ValueChanged += OnLiveParameterChanged;
            numMeasurementThreshold.ValueChanged += OnLiveParameterChanged;
            numTotalSteps.ValueChanged += OnLiveParameterChanged;
            numFractalLevels.ValueChanged += OnLiveParameterChanged;
            numFractalBranchFactor.ValueChanged += OnLiveParameterChanged;
        }

        /// <summary>
        /// Handler for live parameter changes - updates LiveConfig for running simulation.
        /// Thread-safe: writes to fields that calculation thread reads.
        /// </summary>
        private void OnLiveParameterChanged(object? sender, EventArgs e)
        {
            var liveConfig = _simApi.LiveConfig;

            // grpPhysicsConstants values
            liveConfig.InitialEdgeProb = (double)numInitialEdgeProb.Value;
            liveConfig.GravitationalCoupling = (double)numGravitationalCoupling.Value;
            liveConfig.VacuumEnergyScale = (double)numVacuumEnergyScale.Value;
            liveConfig.AnnealingCoolingRate = (double)numAnnealingCoolingRate.Value;
            liveConfig.DecoherenceRate = (double)numDecoherenceRate.Value;
            liveConfig.HotStartTemperature = (double)numHotStartTemperature.Value;
            liveConfig.AdaptiveThresholdSigma = (double)numAdaptiveThresholdSigma.Value;
            liveConfig.WarmupDuration = (double)numWarmupDuration.Value;
            liveConfig.GravityTransitionDuration = (double)numGravityTransitionDuration.Value;

            // grpSimParams values
            liveConfig.TargetDegree = (int)numTargetDegree.Value;
            liveConfig.InitialExcitedProb = (double)numInitialExcitedProb.Value;
            liveConfig.LambdaState = (double)numLambdaState.Value;
            liveConfig.Temperature = (double)numTemperature.Value;
            liveConfig.EdgeTrialProb = (double)numEdgeTrialProb.Value;
            liveConfig.MeasurementThreshold = (double)numMeasurementThreshold.Value;
            liveConfig.FractalLevels = (int)numFractalLevels.Value;
            liveConfig.FractalBranchFactor = (int)numFractalBranchFactor.Value;

            liveConfig.MarkUpdated();

            // Log if simulation is running
            if (_isModernRunning && sender is NumericUpDown num)
            {
                AppendConsole($"[Live] {num.Name}: {num.Value}\n");
            }
        }

        /// <summary>
        /// Handles MaxFPS value change - updates UI timer interval
        /// </summary>
        private void numericUpDown_MaxFPS_ValueChanged(object? sender, EventArgs e)
        {
            int fps = (int)numericUpDown_MaxFPS.Value;
            if (fps <= 0) fps = 1;

            int intervalMs = 1000 / fps;

            if (_uiUpdateTimer != null)
            {
                _uiUpdateTimer.Interval = intervalMs;
            }

            AppendConsole($"[UI] FPS установлен: {fps} ({intervalMs}ms интервал)\n");
        }

        /// <summary>
        /// Инициализирует контролы GPU: заполняет comboBox_GPUIndex доступными устройствами
        /// </summary>
        private void InitializeGpuControls()
        {
            comboBox_GPUIndex.Items.Clear();

            try
            {
                // Проверяем доступность GPU через ComputeSharp
                var defaultDevice = ComputeSharp.GraphicsDevice.GetDefault();
                _simApi.GpuAvailable = defaultDevice != null;

                if (_simApi.GpuAvailable)
                {
                    // Добавляем устройство по умолчанию
                    comboBox_GPUIndex.Items.Add($"0: {defaultDevice.Name}");

                    // Пытаемся получить все устройства
                    int deviceIndex = 0;
                    foreach (var device in ComputeSharp.GraphicsDevice.EnumerateDevices())
                    {
                        if (deviceIndex > 0) // Первое уже добавлено
                        {
                            comboBox_GPUIndex.Items.Add($"{deviceIndex}: {device.Name}");
                        }
                        deviceIndex++;
                    }

                    comboBox_GPUIndex.SelectedIndex = 0;
                    AppendConsole($"[GPU] Обнаружено устройств: {comboBox_GPUIndex.Items.Count}\n");
                    AppendConsole($"[GPU] Устройство по умолчанию: {defaultDevice.Name}\n");
                }
                else
                {
                    comboBox_GPUIndex.Items.Add("Нет доступных GPU");
                    comboBox_GPUIndex.SelectedIndex = 0;
                    checkBox_EnableGPU.Checked = false;
                    checkBox_EnableGPU.Enabled = false;
                }
            }
            catch (Exception ex)
            {
                _simApi.GpuAvailable = false;
                comboBox_GPUIndex.Items.Add("GPU недоступен");
                comboBox_GPUIndex.SelectedIndex = 0;
                checkBox_EnableGPU.Checked = false;
                checkBox_EnableGPU.Enabled = false;
                AppendConsole($"[GPU] Ошибка инициализации: {ex.Message}\n");
            }
        }

        /// <summary>
        /// Проверяет доступность GPU и возвращает true если можно использовать GPU ускорение
        /// </summary>
        private bool CanUseGpu()
        {
            return _simApi.GpuAvailable && checkBox_EnableGPU.Checked;
        }

        /// <summary>
        /// Initializes Graph Health UI controls on the Settings tab.
        /// These controls allow runtime configuration of fragmentation detection
        /// and giant cluster decoherence parameters (RQ-Hypothesis compliance).
        /// </summary>
        private void InitializeGraphHealthControls()
        {
            // Expand tlpPhysicsConstants to add Graph Health section
            // Current row count is 10, we need 6 more rows (1 header + 6 params)
            int startRow = tlpPhysicsConstants.RowCount;
            tlpPhysicsConstants.RowCount = startRow + 7;

            // Add row styles for new rows
            for (int i = 0; i < 7; i++)
            {
                tlpPhysicsConstants.RowStyles.Add(new RowStyle(SizeType.Absolute, 28F));
            }

            // === Header Label ===
            var lblGraphHealthHeader = new Label
            {
                Text = "─── Graph Health (RQ) ───",
                AutoSize = true,
                ForeColor = Color.DarkBlue,
                Font = new Font(Font, FontStyle.Bold),
                Anchor = AnchorStyles.Left | AnchorStyles.Right
            };
            tlpPhysicsConstants.Controls.Add(lblGraphHealthHeader, 0, startRow);
            tlpPhysicsConstants.SetColumnSpan(lblGraphHealthHeader, 2);

            // === Giant Cluster Threshold ===
            var lblGiantClusterThreshold = new Label
            {
                Text = "Giant Cluster (% of N):",
                AutoSize = true,
                Anchor = AnchorStyles.Left | AnchorStyles.Right
            };
            numGiantClusterThreshold = new NumericUpDown
            {
                DecimalPlaces = 2,
                Increment = 0.05m,
                Minimum = 0.10m,
                Maximum = 0.90m,
                Value = (decimal)PhysicsConstants.GiantClusterThreshold,
                Dock = DockStyle.Fill
            };
            numGiantClusterThreshold.ValueChanged += OnGraphHealthParameterChanged;
            tlpPhysicsConstants.Controls.Add(lblGiantClusterThreshold, 0, startRow + 1);
            tlpPhysicsConstants.Controls.Add(numGiantClusterThreshold, 1, startRow + 1);

            // === Emergency Giant Cluster Threshold ===
            var lblEmergencyThreshold = new Label
            {
                Text = "Emergency Cluster (% of N):",
                AutoSize = true,
                Anchor = AnchorStyles.Left | AnchorStyles.Right
            };
            numEmergencyGiantClusterThreshold = new NumericUpDown
            {
                DecimalPlaces = 2,
                Increment = 0.05m,
                Minimum = 0.20m,
                Maximum = 0.95m,
                Value = (decimal)PhysicsConstants.EmergencyGiantClusterThreshold,
                Dock = DockStyle.Fill
            };
            numEmergencyGiantClusterThreshold.ValueChanged += OnGraphHealthParameterChanged;
            tlpPhysicsConstants.Controls.Add(lblEmergencyThreshold, 0, startRow + 2);
            tlpPhysicsConstants.Controls.Add(numEmergencyGiantClusterThreshold, 1, startRow + 2);

            // === Decoherence Rate ===
            var lblDecoherenceRate = new Label
            {
                Text = "Cluster Decoherence Rate:",
                AutoSize = true,
                Anchor = AnchorStyles.Left | AnchorStyles.Right
            };
            numGiantClusterDecoherenceRate = new NumericUpDown
            {
                DecimalPlaces = 3,
                Increment = 0.01m,
                Minimum = 0.01m,
                Maximum = 0.50m,
                Value = (decimal)PhysicsConstants.GiantClusterDecoherenceRate,
                Dock = DockStyle.Fill
            };
            numGiantClusterDecoherenceRate.ValueChanged += OnGraphHealthParameterChanged;
            tlpPhysicsConstants.Controls.Add(lblDecoherenceRate, 0, startRow + 3);
            tlpPhysicsConstants.Controls.Add(numGiantClusterDecoherenceRate, 1, startRow + 3);

            // === Max Decoherence Edges Fraction ===
            var lblMaxEdgesFraction = new Label
            {
                Text = "Max Edges Weakened (%):",
                AutoSize = true,
                Anchor = AnchorStyles.Left | AnchorStyles.Right
            };
            numMaxDecoherenceEdgesFraction = new NumericUpDown
            {
                DecimalPlaces = 2,
                Increment = 0.02m,
                Minimum = 0.02m,
                Maximum = 0.50m,
                Value = (decimal)PhysicsConstants.MaxDecoherenceEdgesFraction,
                Dock = DockStyle.Fill
            };
            numMaxDecoherenceEdgesFraction.ValueChanged += OnGraphHealthParameterChanged;
            tlpPhysicsConstants.Controls.Add(lblMaxEdgesFraction, 0, startRow + 4);
            tlpPhysicsConstants.Controls.Add(numMaxDecoherenceEdgesFraction, 1, startRow + 4);

            // === Critical Spectral Dimension ===
            var lblCriticalSpectralDim = new Label
            {
                Text = "Critical d_S (fragmentation):",
                AutoSize = true,
                Anchor = AnchorStyles.Left | AnchorStyles.Right
            };
            numCriticalSpectralDimension = new NumericUpDown
            {
                DecimalPlaces = 2,
                Increment = 0.1m,
                Minimum = 0.5m,
                Maximum = 2.0m,
                Value = (decimal)PhysicsConstants.CriticalSpectralDimension,
                Dock = DockStyle.Fill
            };
            numCriticalSpectralDimension.ValueChanged += OnGraphHealthParameterChanged;
            tlpPhysicsConstants.Controls.Add(lblCriticalSpectralDim, 0, startRow + 5);
            tlpPhysicsConstants.Controls.Add(numCriticalSpectralDimension, 1, startRow + 5);

            // === Warning Spectral Dimension ===
            var lblWarningSpectralDim = new Label
            {
                Text = "Warning d_S (correction):",
                AutoSize = true,
                Anchor = AnchorStyles.Left | AnchorStyles.Right
            };
            numWarningSpectralDimension = new NumericUpDown
            {
                DecimalPlaces = 2,
                Increment = 0.1m,
                Minimum = 1.0m,
                Maximum = 3.0m,
                Value = (decimal)PhysicsConstants.WarningSpectralDimension,
                Dock = DockStyle.Fill
            };
            numWarningSpectralDimension.ValueChanged += OnGraphHealthParameterChanged;
            tlpPhysicsConstants.Controls.Add(lblWarningSpectralDim, 0, startRow + 6);
            tlpPhysicsConstants.Controls.Add(numWarningSpectralDimension, 1, startRow + 6);
        }

        /// <summary>
        /// Handler for Graph Health parameter changes.
        /// Updates GraphHealthLive in _simApi for auto-tuning to use.
        /// </summary>
        private void OnGraphHealthParameterChanged(object? sender, EventArgs e)
        {
            // Update GraphHealthLive config for auto-tuning
            _simApi.GraphHealthLive.GiantClusterThreshold = (double)numGiantClusterThreshold.Value;
            _simApi.GraphHealthLive.EmergencyGiantClusterThreshold = (double)numEmergencyGiantClusterThreshold.Value;
            _simApi.GraphHealthLive.GiantClusterDecoherenceRate = (double)numGiantClusterDecoherenceRate.Value;
            _simApi.GraphHealthLive.MaxDecoherenceEdgesFraction = (double)numMaxDecoherenceEdgesFraction.Value;
            _simApi.GraphHealthLive.CriticalSpectralDimension = (double)numCriticalSpectralDimension.Value;
            _simApi.GraphHealthLive.WarningSpectralDimension = (double)numWarningSpectralDimension.Value;

            if (_isModernRunning && sender is NumericUpDown num)
            {
                AppendConsole($"[GraphHealth] {num.Name}: {num.Value} (live update applied)\n");
            }
        }

        /// <summary>
        /// Initializes simulation presets in comboBox_Presets.
        /// Call this in constructor after InitializeComponent.
        /// </summary>
        private void InitializePresets()
        {
            comboBox_Presets.Items.Clear();
            comboBox_Presets.Items.Add("Custom");
            comboBox_Presets.Items.Add("Quick Test (100 nodes, 1K steps)");
            comboBox_Presets.Items.Add("Small (300 nodes, 5K steps)");
            comboBox_Presets.Items.Add("Medium (500 nodes, 10K steps)");
            comboBox_Presets.Items.Add("Large (1000 nodes, 20K steps)");
            comboBox_Presets.Items.Add("XL Research (2000 nodes, 50K steps)");
            comboBox_Presets.Items.Add("Stable Clustering (optimized for d_S)");
            comboBox_Presets.Items.Add("High Energy (strong gravity)");
            comboBox_Presets.Items.Add("Quantum Dominated");
            comboBox_Presets.SelectedIndex = 0;
            comboBox_Presets.SelectedIndexChanged += ComboBox_Presets_SelectedIndexChanged;
        }

        /// <summary>
        /// Handles preset selection and applies corresponding configuration to UI controls.
        /// </summary>
        private void ComboBox_Presets_SelectedIndexChanged(object? sender, EventArgs e)
        {
            if (comboBox_Presets.SelectedIndex <= 0) return; // "Custom" selected

            string preset = comboBox_Presets.SelectedItem?.ToString() ?? "";

            // Default base values
            int nodeCount = 300;
            int totalSteps = 5000;
            int targetDegree = 8;
            double initialEdgeProb = 0.08;
            double initialExcitedProb = 0.10;
            double gravitationalCoupling = 0.025;
            double hotStartTemperature = 5.0;
            double decoherenceRate = 0.005;
            double temperature = 10.0;
            double edgeTrialProb = 0.02;
            double vacuumEnergyScale = 0.0001;
            double annealingCoolingRate = 0.995;
            double adaptiveThresholdSigma = 1.5;
            double warmupDuration = 200;
            double gravityTransitionDuration = 137;

            // Physics modules defaults
            bool useSpectralGeometry = true;
            bool useNetworkGravity = true;
            bool useQuantumDriven = true;
            bool useSpinorField = false;
            bool useVacuumFluctuations = true;
            bool useHotStartAnnealing = true;

            // Apply preset-specific values
            if (preset.Contains("Quick Test"))
            {
                nodeCount = 100;
                totalSteps = 1000;
                targetDegree = 6;
                warmupDuration = 50;
                gravityTransitionDuration = 50;
            }
            else if (preset.Contains("Small"))
            {
                nodeCount = 300;
                totalSteps = 5000;
                targetDegree = 8;
            }
            else if (preset.Contains("Medium"))
            {
                nodeCount = 500;
                totalSteps = 10000;
                targetDegree = 10;
                warmupDuration = 300;
                gravityTransitionDuration = 200;
            }
            else if (preset.Contains("Large"))
            {
                nodeCount = 1000;
                totalSteps = 20000;
                targetDegree = 12;
                warmupDuration = 500;
                gravityTransitionDuration = 300;
                gravitationalCoupling = 0.02; // Lower G for larger graphs
                decoherenceRate = 0.003;
            }
            else if (preset.Contains("XL Research"))
            {
                nodeCount = 2000;
                totalSteps = 50000;
                targetDegree = 15;
                warmupDuration = 1000;
                gravityTransitionDuration = 500;
                gravitationalCoupling = 0.015;
                decoherenceRate = 0.002;
                initialEdgeProb = 0.05; // Sparser initial graph
            }
            else if (preset.Contains("Stable Clustering"))
            {
                // Optimized for stable spectral dimension 2-4
                nodeCount = 500;
                totalSteps = 15000;
                targetDegree = 10;
                gravitationalCoupling = 0.02;
                hotStartTemperature = 3.0;
                decoherenceRate = 0.008;
                adaptiveThresholdSigma = 1.2; // Tighter threshold
                warmupDuration = 400;
                gravityTransitionDuration = 250;
                annealingCoolingRate = 0.998; // Slower cooling
            }
            else if (preset.Contains("High Energy"))
            {
                // Strong gravity, more structure formation
                nodeCount = 400;
                totalSteps = 10000;
                targetDegree = 12;
                gravitationalCoupling = 0.08; // Strong gravity
                hotStartTemperature = 8.0;
                decoherenceRate = 0.003;
                temperature = 15.0;
                useSpinorField = true;
            }
            else if (preset.Contains("Quantum Dominated"))
            {
                // Emphasize quantum effects
                nodeCount = 300;
                totalSteps = 8000;
                targetDegree = 8;
                gravitationalCoupling = 0.01; // Weak gravity
                hotStartTemperature = 2.0;
                decoherenceRate = 0.015; // Higher decoherence
                useQuantumDriven = true;
                useSpinorField = true;
                useVacuumFluctuations = true;
            }

            // Apply values to UI controls (suppress events temporarily)
            SuspendLayout();
            try
            {
                // Simulation Parameters
                numNodeCount.Value = nodeCount;
                numTotalSteps.Value = totalSteps;
                numTargetDegree.Value = targetDegree;
                numInitialExcitedProb.Value = (decimal)initialExcitedProb;
                numTemperature.Value = (decimal)temperature;
                numEdgeTrialProb.Value = (decimal)edgeTrialProb;

                // Physics Constants
                numInitialEdgeProb.Value = (decimal)initialEdgeProb;
                numGravitationalCoupling.Value = (decimal)gravitationalCoupling;
                numVacuumEnergyScale.Value = (decimal)vacuumEnergyScale;
                numAnnealingCoolingRate.Value = (decimal)annealingCoolingRate;
                numDecoherenceRate.Value = (decimal)decoherenceRate;
                numHotStartTemperature.Value = (decimal)hotStartTemperature;
                numAdaptiveThresholdSigma.Value = (decimal)adaptiveThresholdSigma;
                numWarmupDuration.Value = (decimal)warmupDuration;
                numGravityTransitionDuration.Value = (decimal)gravityTransitionDuration;

                // Physics Modules
                chkSpectralGeometry.Checked = useSpectralGeometry;
                chkNetworkGravity.Checked = useNetworkGravity;
                chkQuantumDriven.Checked = useQuantumDriven;
                chkSpinorField.Checked = useSpinorField;
                chkVacuumFluctuations.Checked = useVacuumFluctuations;
                chkHotStartAnnealing.Checked = useHotStartAnnealing;
            }
            finally
            {
                ResumeLayout();
            }

            AppendConsole($"[Preset] Применён: {preset}\n");
            AppendConsole($"  Nodes={nodeCount}, Steps={totalSteps}, G={gravitationalCoupling}, T_hot={hotStartTemperature}\n");
        }

        // === Experiment Manager ===
        private IExperiment? _selectedExperiment;

        /// <summary>
        /// Initializes experiment selection combobox with available experiments.
        /// </summary>
        private void InitializeExperiments()
        {
            comboBox_Experiments.Items.Clear();
            comboBox_Experiments.Items.Add("(Custom / Manual)");

            foreach (var experiment in ExperimentFactory.AvailableExperiments)
            {
                comboBox_Experiments.Items.Add(experiment.Name);
            }

            comboBox_Experiments.SelectedIndex = 0;
            comboBox_Experiments.SelectedIndexChanged += ComboBox_Experiments_SelectedIndexChanged;
        }

        /// <summary>
        /// Handles experiment selection change.
        /// Loads experiment configuration into UI controls.
        /// </summary>
        private void ComboBox_Experiments_SelectedIndexChanged(object? sender, EventArgs e)
        {
            if (comboBox_Experiments.SelectedIndex <= 0)
            {
                // Custom / Manual mode - clear experiment
                _selectedExperiment = null;
                _simApi.SetCustomInitializer(null);
                AppendConsole("[Experiment] Manual mode selected - using custom parameters\n");
                return;
            }

            string experimentName = comboBox_Experiments.SelectedItem?.ToString() ?? "";
            var experiment = ExperimentFactory.GetByName(experimentName);

            if (experiment != null)
            {
                LoadExperiment(experiment);
            }
        }

        /// <summary>
        /// Loads an experiment configuration into UI controls.
        /// 
        /// This method:
        /// 1. Gets the experiment's StartupConfig
        /// 2. Calls ApplyPhysicsOverrides() for any physics constant changes
        /// 3. Updates all UI fields to match the config
        /// 4. Sets up the custom initializer in the simulation API
        /// </summary>
        private void LoadExperiment(IExperiment experiment)
        {
            _selectedExperiment = experiment;

            var config = experiment.GetConfig();
            experiment.ApplyPhysicsOverrides();

            // Store custom initializer in simulation API
            _simApi.SetCustomInitializer(experiment.CustomInitializer);

            // Update UI fields with experiment configuration
            SuspendLayout();
            try
            {
                // === Simulation Parameters (grpSimParams) ===
                numNodeCount.Value = config.NodeCount;
                numTotalSteps.Value = config.TotalSteps;
                numTargetDegree.Value = config.TargetDegree;
                numInitialExcitedProb.Value = (decimal)config.InitialExcitedProb;
                numTemperature.Value = (decimal)config.Temperature;
                numEdgeTrialProb.Value = (decimal)config.EdgeTrialProbability;
                numLambdaState.Value = (decimal)config.LambdaState;
                numMeasurementThreshold.Value = (decimal)config.MeasurementThreshold;
                numFractalLevels.Value = config.FractalLevels;
                numFractalBranchFactor.Value = config.FractalBranchFactor;

                // === Physics Constants (grpPhysicsConstants) ===
                numInitialEdgeProb.Value = (decimal)config.InitialEdgeProb;
                numGravitationalCoupling.Value = (decimal)config.GravitationalCoupling;
                numVacuumEnergyScale.Value = (decimal)config.VacuumEnergyScale;
                numAnnealingCoolingRate.Value = (decimal)config.AnnealingCoolingRate;
                numDecoherenceRate.Value = (decimal)config.DecoherenceRate;
                numHotStartTemperature.Value = (decimal)config.HotStartTemperature;
                numAdaptiveThresholdSigma.Value = (decimal)config.AdaptiveThresholdSigma;
                numWarmupDuration.Value = (decimal)config.WarmupDuration;
                numGravityTransitionDuration.Value = (decimal)config.GravityTransitionDuration;

                // === Physics Modules (checkboxes) ===
                chkSpectralGeometry.Checked = config.UseSpectralGeometry;
                chkNetworkGravity.Checked = config.UseNetworkGravity;
                chkQuantumDriven.Checked = config.UseQuantumDrivenStates;
                chkSpinorField.Checked = config.UseSpinorField;
                chkVacuumFluctuations.Checked = config.UseVacuumFluctuations;
                chkHotStartAnnealing.Checked = config.UseHotStartAnnealing;
            }
            finally
            {
                ResumeLayout();
            }

            // Log experiment selection
            AppendConsole($"\n[Experiment] === {experiment.Name} ===\n");
            AppendConsole($"[Experiment] {experiment.Description}\n");
            AppendConsole($"[Experiment] Nodes={config.NodeCount}, Steps={config.TotalSteps}, " +
                         $"G={config.GravitationalCoupling}, T_hot={config.HotStartTemperature}\n");

            if (experiment.CustomInitializer != null)
            {
                AppendConsole($"[Experiment] Custom topology initializer will be applied\n");
            }

            // Reset preset selector to "Custom" since we're using experiment
            comboBox_Presets.SelectedIndex = 0;
        }


        private void DrawingPanel_Paint(object? sender, PaintEventArgs e)
        {
            if (canvasBitmap != null)
            {
                try
                {
                    e.Graphics.DrawImage(canvasBitmap, 0, 0);
                }
                catch (Exception ex)
                {
                    // Draw error on panel
                    e.Graphics.Clear(Color.Red);
                    e.Graphics.DrawString($"Paint error: {ex.Message}", SystemFonts.DefaultFont, Brushes.White, 10, 10);
                }
            }
            else
            {
                // Draw placeholder when no bitmap
                e.Graphics.Clear(Color.LightGray);
                e.Graphics.DrawString("No graph bitmap", SystemFonts.DefaultFont, Brushes.Black, 10, 10);
            }
        }

        private void DrawingPanel_MouseClick(object? sender, MouseEventArgs e)
        {
            if (_simulationEngine?.Graph == null || !_isModernRunning)
                return;

            // Находим ближайший узел к клику и возбуждаем его
            int nodeIndex = FindNodeAtPosition(e.Location);
            if (nodeIndex >= 0)
            {
                _simulationEngine.Graph.FlipNodeWithNeighbors(nodeIndex);
                AppendConsole($"Импульс в узел {nodeIndex}\n");
            }
        }


        private void BtnSnapshotImage_Click(object? sender, EventArgs e)
        {
            try
            {
                if (canvasBitmap is null)
                {
                    AppendConsole("Нет изображения для сохранения.\n");
                    return;
                }

                using SaveFileDialog dlg = new()
                {
                    Filter = "PNG Image (*.png)|*.png",
                    FileName = "rq-simulation-snapshot.png"
                };

                if (dlg.ShowDialog(this) == DialogResult.OK)
                {
                    canvasBitmap.Save(dlg.FileName);
                    AppendConsole($"Снимок сохранён: {dlg.FileName}\n");
                }
            }
            catch (Exception ex)
            {
                AppendConsole($"Ошибка сохранения снимка: {ex.Message}\n");
            }
        }

        private void button_SaveSimReult_Click(object sender, EventArgs e)
        {
            if (_simApi.IsModernRunning)
            {
                AppendConsole("Экспорт невозможен: симуляция ещё выполняется.\n");
                return;
            }

            // Check if there's any data to export
            bool hasSessionHistory = _simApi.SessionHistory.Count > 0;
            bool hasCurrentData = _simApi.LastResult != null || _seriesSteps.Count > 0;

            if (!hasSessionHistory && !hasCurrentData)
            {
                AppendConsole("Нет данных для сохранения. Запустите симуляцию.\n");
                return;
            }

            try
            {
                using SaveFileDialog dlg = new()
                {
                    Filter = "JSON file (*.json)|*.json",
                    FileName = $"rq-simulation-export-{DateTime.Now:yyyyMMdd-HHmmss}.json"
                };
                if (dlg.ShowDialog(this) != DialogResult.OK) return;

                // Build complete export including all session history
                var export = new
                {
                    meta = new
                    {
                        exportedAt = DateTime.UtcNow,
                        application = "RqSimForms",
                        version = Application.ProductVersion,
                        sessionCount = _simApi.SessionHistory.Count,
                        hasUnsavedCurrentSession = hasCurrentData && _simApi.CurrentSession == null,
                        currentFilters = new { weightThreshold = _displayWeightThreshold, heavyOnly = _displayShowHeavyOnly, dynamicLayout = _useDynamicCoords },
                        synthesis = new { synthesisCount = _simApi.SynthesisCount, fissionCount = _simApi.FissionCount }
                    },
                    // Export PhysicsConstants snapshot for reproducibility
                    physicsConstants = new
                    {
                        // Fundamental
                        fineStructureConstant = PhysicsConstants.FineStructureConstant,
                        strongCouplingConstant = PhysicsConstants.StrongCouplingConstant,
                        weakMixingAngle = PhysicsConstants.WeakMixingAngle,
                        // Gravity
                        gravitationalCoupling = PhysicsConstants.GravitationalCoupling,
                        warmupGravitationalCoupling = PhysicsConstants.WarmupGravitationalCoupling,
                        warmupDuration = PhysicsConstants.WarmupDuration,
                        gravityTransitionDuration = PhysicsConstants.GravityTransitionDuration,
                        cosmologicalConstant = PhysicsConstants.CosmologicalConstant,
                        // Annealing
                        initialAnnealingTemperature = PhysicsConstants.InitialAnnealingTemperature,
                        finalAnnealingTemperature = PhysicsConstants.FinalAnnealingTemperature,
                        annealingTimeConstant = PhysicsConstants.AnnealingTimeConstant,
                        // Fields
                        fieldDiffusionRate = PhysicsConstants.FieldDiffusionRate,
                        kleinGordonMass = PhysicsConstants.KleinGordonMass,
                        diracCoupling = PhysicsConstants.DiracCoupling,
                        gaugeCouplingConstant = PhysicsConstants.GaugeCouplingConstant,
                        // Quantum
                        vacuumFluctuationBaseRate = PhysicsConstants.VacuumFluctuationBaseRate,
                        pairCreationEnergyThreshold = PhysicsConstants.PairCreationEnergyThreshold,
                        spinorNormalizationThreshold = PhysicsConstants.SpinorNormalizationThreshold,
                        // Topology
                        edgeCreationBarrier = PhysicsConstants.EdgeCreationBarrier,
                        edgeAnnihilationBarrier = PhysicsConstants.EdgeAnnihilationBarrier,
                        defaultHeavyClusterThreshold = PhysicsConstants.DefaultHeavyClusterThreshold,
                        adaptiveThresholdSigma = PhysicsConstants.AdaptiveThresholdSigma,
                        minimumClusterSize = PhysicsConstants.MinimumClusterSize
                    },
                    // Export all archived sessions
                    sessions = _simApi.SessionHistory.Select(s => new
                    {
                        s.SessionId,
                        s.StartedAt,
                        s.EndedAt,
                        endReason = s.EndReason.ToString(),
                        s.GpuEnabled,
                        s.GpuDeviceName,
                        s.LastStep,
                        s.TotalStepsPlanned,
                        config = s.Config == null ? null : new
                        {
                            s.Config.NodeCount,
                            s.Config.InitialEdgeProb,
                            s.Config.InitialExcitedProb,
                            s.Config.TargetDegree,
                            s.Config.LambdaState,
                            s.Config.Temperature,
                            s.Config.EdgeTrialProbability,
                            s.Config.MeasurementThreshold,
                            s.Config.DynamicMeasurementThreshold,
                            s.Config.Seed,
                            s.Config.TotalSteps,
                            s.Config.LogEvery,
                            s.Config.BaselineWindow,
                            s.Config.FirstImpulse,
                            s.Config.ImpulsePeriod,
                            s.Config.CalibrationStep,
                            s.Config.VisualizationInterval,
                            s.Config.MeasurementLogInterval,
                            s.Config.UseQuantumDrivenStates,
                            s.Config.UseSpacetimePhysics,
                            s.Config.UseSpinorField,
                            s.Config.UseVacuumFluctuations,
                            s.Config.UseBlackHolePhysics,
                            s.Config.UseYangMillsGauge,
                            s.Config.UseEnhancedKleinGordon,
                            s.Config.UseInternalTime,
                            s.Config.UseSpectralGeometry,
                            s.Config.UseQuantumGraphity,
                            s.Config.FractalLevels,
                            s.Config.FractalBranchFactor,
                            // Extended physics parameters
                            s.Config.GravitationalCoupling,
                            s.Config.VacuumEnergyScale,
                            s.Config.InitialNetworkTemperature,
                            s.Config.AnnealingCoolingRate,
                            s.Config.UseRelationalTime,
                            s.Config.UseRelationalYangMills,
                            s.Config.UseNetworkGravity,
                            s.Config.DecoherenceRate,
                            s.Config.UseUnifiedPhysicsStep,
                            s.Config.EnforceGaugeConstraints,
                            s.Config.UseCausalRewiring,
                            s.Config.UseAsynchronousTime,
                            s.Config.UseTopologicalProtection,
                            s.Config.ValidateEnergyConservation,
                            s.Config.UseMexicanHatPotential,
                            s.Config.UseHotStartAnnealing,
                            s.Config.UseGeometryMomenta,
                            s.Config.UseTopologicalCensorship,
                            s.Config.UseSpectralMass,
                            s.Config.HotStartTemperature
                        },
                        filters = new { s.Filters.WeightThreshold, s.Filters.HeavyOnly, s.Filters.DynamicLayout },
                        // Session metadata
                        s.FinalSpectralDimension,
                        s.FinalNetworkTemperature,
                        s.WallClockDurationSeconds,
                        timeSeries = new
                        {
                            steps = s.SeriesSteps.ToArray(),
                            excited = s.SeriesExcited.ToArray(),
                            heavyMass = s.SeriesHeavyMass.ToArray(),
                            heavyCount = s.SeriesHeavyCount.ToArray(),
                            largestCluster = s.SeriesLargestCluster.ToArray(),
                            avgLinkDistance = s.SeriesAvgDist.ToArray(),
                            density = s.SeriesDensity.ToArray(),
                            energy = s.SeriesEnergy.ToArray(),
                            corr = s.SeriesCorr.ToArray(),
                            strongEdges = s.SeriesStrongEdges.ToArray(),
                            qNorm = s.SeriesQNorm.ToArray(),
                            qEnergy = s.SeriesQEnergy.ToArray(),
                            entanglement = s.SeriesEntanglement.ToArray(),
                            spectralDimension = s.SeriesSpectralDimension.ToArray(),
                            networkTemperature = s.SeriesNetworkTemperature.ToArray(),
                            effectiveG = s.SeriesEffectiveG.ToArray(),
                            adaptiveThreshold = s.SeriesAdaptiveThreshold.ToArray()
                        },
                        synthesisScatter = s.SynthesisData?.Select(x => new { x.volume, x.deltaMass }).ToArray(),
                        synthesisCount = s.SynthesisCount,
                        fissionCount = s.FissionCount,
                        result = s.Result == null ? null : new
                        {
                            s.Result.AverageExcited,
                            s.Result.MaxExcited,
                            s.Result.MeasurementConfigured,
                            s.Result.MeasurementTriggered
                        },
                        modernResult = s.ModernResult == null ? null : new
                        {
                            s.ModernResult.FinalTime,
                            s.ModernResult.ExcitedCount,
                            s.ModernResult.HeavyClusterCount,
                            s.ModernResult.HeavyClusterTotalMass,
                            s.ModernResult.HeavyClusterMaxMass,
                            s.ModernResult.HeavyClusterMeanMass,
                            s.ModernResult.ScalarFieldEnergy,
                            s.ModernResult.HiggsFieldEnergy
                        },
                        csv = s.Result == null ? null : new
                        {
                            timeSeriesCsv = s.Result.TimeSeries?.ToArray(),
                            avalancheStatsCsv = s.Result.AvalancheStats?.ToArray(),
                            measurementEventsCsv = s.Result.MeasurementEvents?.ToArray(),
                            heavyClustersCsv = s.Result.HeavyClusters?.ToArray()
                        },
                        consoleLog = s.ConsoleLog,
                        summaryText = s.SummaryText,
                        importantEvents = s.ImportantEvents.Select(ev => new { ev.Step, ev.Type, ev.Detail }).ToArray()
                    }).ToArray(),
                    // Also include current UI state (for continuity)
                    currentUiState = new
                    {
                        summary = summaryTextBox?.Text,
                        console = consoleTextBox?.Text,
                        importantEvents = lvEvents?.Items.Cast<ListViewItem>().Select(i => new
                        {
                            step = i.SubItems.Count > 0 ? i.SubItems[0].Text : null,
                            type = i.SubItems.Count > 1 ? i.SubItems[1].Text : null,
                            detail = i.SubItems.Count > 2 ? i.SubItems[2].Text : null
                        }).ToArray()
                    }
                };

                var options = new JsonSerializerOptions { WriteIndented = true };
                var json = JsonSerializer.Serialize(export, options);
                File.WriteAllText(dlg.FileName, json, Encoding.UTF8);

                AppendConsole($"Данные сохранены: {dlg.FileName}\n");
                AppendConsole($"Экспортировано сессий: {_simApi.SessionHistory.Count}\n");

                // Ask if user wants to clear session history after export
                if (_simApi.SessionHistory.Count > 0)
                {
                    var clearResult = MessageBox.Show(
                        $"Экспортировано {_simApi.SessionHistory.Count} сессий.\n\nОчистить историю сессий?",
                        "Очистка истории",
                        MessageBoxButtons.YesNo,
                        MessageBoxIcon.Question);

                    if (clearResult == DialogResult.Yes)
                    {
                        _simApi.SessionHistory.Clear();
                        AppendConsole("История сессий очищена.\n");
                    }
                }
            }
            catch (Exception ex)
            {
                AppendConsole($"Ошибка сохранения JSON: {ex.Message}\n");
            }
        }

        private void ToolViewGraph_Click(object? sender, EventArgs e) => tabControl1.SelectedTab = tabPage_GUI;
        private void ToolViewConsole_Click(object? sender, EventArgs e) => tabControl1.SelectedTab = tabPage_Console;
        private void ToolViewAnalysis_Click(object? sender, EventArgs e) => tabControl1.SelectedTab = tabPage_Sythnesis;
        private void ToolViewSummary_Click(object? sender, EventArgs e) => tabControl1.SelectedTab = tabPage_Summary;

        // Shortcut accessors for _simApi state
        private bool _isModernRunning { get => _simApi.IsModernRunning; set => _simApi.IsModernRunning = value; }
        private bool _spectrumLoggingEnabled { get => _simApi.SpectrumLoggingEnabled; set => _simApi.SpectrumLoggingEnabled = value; }

        // Результат и метрики Modern для экспорта/GUI
        private RQSimulation.ExampleModernSimulation.ScenarioResult? _modernResult => _simApi.ModernResult;

        // Synthesis data shortcuts
        private List<(int volume, double deltaMass)>? synthesisData { get => _simApi.SynthesisData; set => _simApi.SynthesisData = value; }
        private int synthesisCount { get => _simApi.SynthesisCount; set => _simApi.SynthesisCount = value; }
        private int fissionCount { get => _simApi.FissionCount; set => _simApi.FissionCount = value; }



        private void checkBox_EnableGPU_CheckedChanged(object sender, EventArgs e)
        {
            // Здесь можно добавить логику для обработки изменения состояния GPU
            AppendConsole($"GPU ускорение {(checkBox_EnableGPU.Checked ? "включено" : "выключено")}.\n");
        }

        // Обработчик кнопки запуска/остановки Modern Sim
        private void button_RunSimulation_Click(object? sender, EventArgs e)
        {
            if (_isModernRunning)
            {
                // Ask for confirmation before stopping
                var result = MessageBox.Show(
                    "Симуляция запущена. Остановить и сохранить текущую сессию?\n\n" +
                    "Все данные будут сохранены в историю сессий и могут быть экспортированы в JSON.",
                    "Остановка симуляции",
                    MessageBoxButtons.YesNo,
                    MessageBoxIcon.Question);

                if (result == DialogResult.Yes)
                {
                    AppendConsole("[Modern] Остановка симуляции по запросу пользователя...\n");
                    _modernCts?.Cancel();
                }
                return;
            }

            // Start new simulation
            _isModernRunning = true;
            _simApi.SimulationComplete = false;
            _modernCts = new CancellationTokenSource();
            button_RunModernSim.Text = "Stop Modern Sim";

            // Initialize new session
            string? gpuDeviceName = comboBox_GPUIndex.SelectedItem as string;
            var filters = new DisplayFilters
            {
                WeightThreshold = _displayWeightThreshold,
                HeavyOnly = _displayShowHeavyOnly,
                DynamicLayout = _useDynamicCoords
            };
            _simApi.CurrentSession = _simApi.CreateSession(CanUseGpu(), gpuDeviceName, filters);

            // Start UI update timer (runs on UI thread, FPS from control)
            _uiUpdateTimer = new System.Windows.Forms.Timer();
            int targetFps = Math.Max(1, (int)numericUpDown_MaxFPS.Value);
            _uiUpdateTimer.Interval = 1000 / targetFps;
            _uiUpdateTimer.Tick += UiUpdateTimer_Tick;
            _uiUpdateTimer.Start();

            // Run simulation entirely in background thread (Event-Based only)
            var ct = _modernCts.Token;
            Task.Run(() =>
            {
                try
                {
                    // Get config from UI thread
                    SimulationConfig config = null!;
                    Invoke(new Action(() => config = GetConfigFromUI()));

                    // Initialize simulation state for event-based loop
                    _simApi.InitializeSimulation(config);
                    _simApi.InitializeLiveConfig(config);

                    // RQ-HYPOTHESIS: totalEvents = sweeps * N
                    // Each sweep updates all N nodes once (in proper time order)
                    // TotalSteps in UI maps to number of sweeps for event-based mode
                    int totalSweeps = Math.Max(1, config.TotalSteps);
                    int totalEvents = totalSweeps * Math.Max(1, config.NodeCount);

                    Invoke(new Action(() => AppendConsole(
                        $"[EventBased] Mode: {totalSweeps} sweeps × {config.NodeCount} nodes = {totalEvents} events\n")));

                    // Pass GPU flag to enable GPU acceleration when checkbox is checked
                    bool useGpuAcceleration = false;
                    Invoke(new Action(() => useGpuAcceleration = CanUseGpu()));

                    _simApi.RunParallelEventBasedLoop(ct, totalEvents, useParallel: true, useGpu: useGpuAcceleration);

                    // Signal completion
                    _simApi.SimulationComplete = true;
                    BeginInvoke(new Action(() => OnSimulationCompleted(SessionEndReason.Completed)));
                }
                catch (OperationCanceledException)
                {
                    _simApi.SimulationComplete = true;
                    BeginInvoke(new Action(() => OnSimulationCompleted(SessionEndReason.CancelledByUser)));
                }
                catch (Exception ex)
                {
                    _simApi.SimulationComplete = true;
                    BeginInvoke(new Action(() =>
                    {
                        AppendConsole($"[Error] {ex.Message}\n");
                        OnSimulationCompleted(SessionEndReason.Error);
                    }));
                }
            });
        }

        /// <summary>
        /// Timer tick handler - updates UI from live metrics (runs on UI thread)
        /// Non-blocking: skips frame if dispatcher is busy with calculation thread
        /// Performance: Skips invisible tab updates, uses precomputed statistics
        /// Target: 5-10 FPS for smooth visualization
        /// 
        /// OPTIMIZATION: Graph drawing is throttled to every 5 seconds in Event-Based mode
        /// to reduce CPU overhead. Charts update more frequently (every tick).
        /// </summary>
        private void UiUpdateTimer_Tick(object? sender, EventArgs e)
        {
            if (_simApi.SimulationComplete || !_isModernRunning)
            {
                _uiUpdateTimer?.Stop();
                return;
            }

            // Force get display data with timeout - ensures we get data even under contention
            // This is critical for Event-Based mode where calculation thread holds lock frequently
            var displayData = _simApi.Dispatcher.ForceGetDisplayDataImmediate(timeoutMs: 50);

            // Read volatile fields (lock-free)
            int step = _simApi.LiveStep;
            int excited = _simApi.LiveExcited;
            double heavyMass = _simApi.LiveHeavyMass;
            int largestCluster = _simApi.LiveLargestCluster;
            int strongEdges = _simApi.LiveStrongEdges;
            double qNorm = _simApi.LiveQNorm;
            double entanglement = _simApi.LiveEntanglement;
            double correlation = _simApi.LiveCorrelation;
            int totalSteps = _simApi.LiveTotalSteps;

            // Calculate statistics without LINQ (precompute in loop for performance)
            double avgExcited = 0.0;
            int maxExcited = 0;
            if (displayData.DecimatedExcited.Length > 0)
            {
                int sum = 0;
                for (int i = 0; i < displayData.DecimatedExcited.Length; i++)
                {
                    int val = displayData.DecimatedExcited[i];
                    sum += val;
                    if (val > maxExcited) maxExcited = val;
                }
                avgExcited = (double)sum / displayData.DecimatedExcited.Length;
            }

            // OPTIMIZATION: Draw graph less frequently (every 5 seconds)
            // This reduces CPU load significantly while keeping simulation fast
            var now = DateTime.UtcNow;
            if (now - _lastGraphDrawTime >= _graphDrawInterval)
            {
                _lastGraphDrawTime = now;
                DrawGraph();
            }

            // Invalidate chart panels (always - they're lightweight)
            panelOnChart?.Invalidate();
            panelHeavyChart?.Invalidate();
            panelClusterChart?.Invalidate();
            panelEnergyChart?.Invalidate();

            int graphN = _simulationEngine?.Graph?.N ?? 100;
            string phase = excited > graphN / 3 ? "Active" : (excited > graphN / 10 ? "Moderate" : "Quiet");

            // Update all dashboard components (always update, they're lightweight)
            UpdateDashboard(step, totalSteps, excited, heavyMass, largestCluster,
                strongEdges, phase, qNorm, entanglement, correlation);
            UpdateStatusBar(step, totalSteps, excited, avgExcited, heavyMass);
            UpdateRunSummary(totalSteps, step, avgExcited, maxExcited, 0, false, false);
            UpdateLiveMetrics(0.0, 0.0, strongEdges, largestCluster, heavyMass, _spectrumLoggingEnabled);
        }

        /// <summary>
        /// Called when simulation completes (on UI thread)
        /// </summary>
        private void OnSimulationCompleted(SessionEndReason reason)
        {
            // Stop timer
            _uiUpdateTimer?.Stop();
            _uiUpdateTimer?.Dispose();
            _uiUpdateTimer = null;

            // Archive session
            var events = lvEvents?.Items.Cast<ListViewItem>()
                .Select(i => new ImportantEvent
                {
                    Step = int.TryParse(i.SubItems.Count > 0 ? i.SubItems[0].Text : "0", out int s) ? s : 0,
                    Type = i.SubItems.Count > 1 ? i.SubItems[1].Text : string.Empty,
                    Detail = i.SubItems.Count > 2 ? i.SubItems[2].Text : string.Empty
                }).ToList() ?? [];

            _simApi.ArchiveSession(reason, consoleTextBox?.Text ?? "", summaryTextBox?.Text ?? "", events);

            string reasonText = reason switch
            {
                SessionEndReason.Completed => "завершена успешно",
                SessionEndReason.CancelledByUser => "отменена",
                _ => "завершилась с ошибкой"
            };
            AppendConsole($"[Session] Сессия {_simApi.CurrentSession?.SessionId} {reasonText} и добавлена в историю.\n");

            // Cleanup
            _simApi.Cleanup();
            _cachedNodePositions = null;

            _modernCts?.Dispose();
            _modernCts = null;
            _isModernRunning = false;
            _simApi.CurrentSession = null;
            button_RunModernSim.Text = "Run Modern Sim";

            // Final UI update
            DrawGraph();
            panelOnChart?.Invalidate();
            panelHeavyChart?.Invalidate();
            panelClusterChart?.Invalidate();
            panelEnergyChart?.Invalidate();

            AppendConsole($"[Session] Всего сессий в истории: {_simApi.SessionHistory.Count}\n");
        }

        private void button_ForceRedrawGraphImage_Click(object sender, EventArgs e)
        {
            // Force synchronous redraw

            // CRITICAL: Stop the UI timer temporarily so our diagnostic bitmap is visible
            bool timerWasRunning = _uiUpdateTimer?.Enabled ?? false;
            _uiUpdateTimer?.Stop();

            // Check simulation engine state
            var simEngine = _simApi.SimulationEngine;
            if (simEngine == null)
            {
                AppendConsole("[ForceRedraw] SimulationEngine is null - cannot draw graph\n");
                return;
            }

            var graph = simEngine.Graph;
            if (graph == null)
            {
                AppendConsole("[ForceRedraw] Graph is null - cannot draw graph\n");
                return;
            }

            // Force reset drawing flag to ensure we can draw
            Interlocked.Exchange(ref _isDrawing, 0);

            // SYNCHRONOUS graph drawing (not Task.Run)
            try
            {
                int panelWidth = drawingPanel.Width;
                int panelHeight = drawingPanel.Height;

                if (panelWidth <= 0 || panelHeight <= 0)
                {
                    AppendConsole($"[ForceRedraw] Invalid panel dimensions: {panelWidth}x{panelHeight}\n");
                    return;
                }

                // Build node positions
                BuildNodePositions();
                var localPositions = _cachedNodePositions;

                if (localPositions == null || localPositions.Length == 0)
                {
                    AppendConsole("[ForceRedraw] localPositions is null or empty\n");
                    return;
                }

                // Create bitmap synchronously
                var newBitmap = new Bitmap(panelWidth, panelHeight);

                using (var g = Graphics.FromImage(newBitmap))
                {
                    g.Clear(Color.White);
                    g.SmoothingMode = System.Drawing.Drawing2D.SmoothingMode.AntiAlias;

                    int n = graph.N;

                    // Draw edges
                    for (int i = 0; i < n && i < localPositions.Length; i++)
                    {
                        var p1 = localPositions[i];
                        foreach (int j in graph.Neighbors(i))
                        {
                            if (j <= i) continue;
                            if (j >= localPositions.Length) continue;

                            double w = graph.Weights[i, j];
                            if (w < _displayWeightThreshold) continue;

                            var pen = GetEdgePen(w);
                            var p2 = localPositions[j];
                            g.DrawLine(pen, p1, p2);
                        }
                    }

                    // Draw nodes
                    float nodeSize = Math.Max(4, Math.Min(11, 500f / Math.Max(1, n)));
                    for (int i = 0; i < n && i < localPositions.Length; i++)
                    {
                        var p = localPositions[i];
                        var brush = GetNodeBrush(graph, i);
                        g.FillEllipse(brush, p.X - nodeSize / 2, p.Y - nodeSize / 2, nodeSize, nodeSize);
                        g.DrawEllipse(Pens.Black, p.X - nodeSize / 2, p.Y - nodeSize / 2, nodeSize, nodeSize);
                    }
                }

                // Swap bitmaps SYNCHRONOUSLY
                var oldBitmap = canvasBitmap;
                canvasBitmap = newBitmap;
                oldBitmap?.Dispose();

                // Force repaint
                drawingPanel.Invalidate();
                drawingPanel.Update();

                AppendConsole($"[ForceRedraw] Graph redrawn: {graph.N} nodes, step {_simApi.LiveStep}\n");
            }
            catch (Exception ex)
            {
                AppendConsole($"[ForceRedraw] Error: {ex.Message}\n");
            }

            // Also force chart repaint
            panelOnChart?.Invalidate();
            panelHeavyChart?.Invalidate();
            panelClusterChart?.Invalidate();
            panelEnergyChart?.Invalidate();

            // Force refresh the entire tab
            tabPage_GUI.Refresh();

            // Restart timer if it was running (with 2 second delay)
            if (timerWasRunning && _uiUpdateTimer != null)
            {
                var restartTimer = new System.Windows.Forms.Timer { Interval = 2000 };
                restartTimer.Tick += (s, args) =>
                {
                    restartTimer.Stop();
                    restartTimer.Dispose();
                    if (_isModernRunning && !_simApi.SimulationComplete)
                    {
                        _uiUpdateTimer?.Start();
                        AppendConsole("[UI] Timer restarted after force redraw\n");
                    }
                };
                restartTimer.Start();
            }
        }

        private void btnExpornShortJson_Click(object sender, EventArgs e)
        {
            // Quick snapshot export - current state + startup parameters + decimated dynamics
            try
            {
                using SaveFileDialog dlg = new()
                {
                    Filter = "JSON file (*.json)|*.json",
                    FileName = $"rq-snapshot-{DateTime.Now:yyyyMMdd-HHmmss}.json",
                    Title = "Сохранить снимок текущего состояния"
                };

                if (dlg.ShowDialog(this) != DialogResult.OK) return;

                // Get last available step data
                int lastStep = _seriesSteps.Count > 0 ? _seriesSteps[^1] : _simApi.LiveStep;
                int lastExcited = _seriesExcited.Count > 0 ? _seriesExcited[^1] : _simApi.LiveExcited;
                double lastHeavyMass = _seriesHeavyMass.Count > 0 ? _seriesHeavyMass[^1] : _simApi.LiveHeavyMass;
                double lastSpectralDim = _seriesSpectralDimension.Count > 0 ? _seriesSpectralDimension[^1] : _simApi.LiveSpectralDim;
                double lastTemp = _seriesNetworkTemperature.Count > 0 ? _seriesNetworkTemperature[^1] : _simApi.LiveTemp;
                double lastEffectiveG = _seriesEffectiveG.Count > 0 ? _seriesEffectiveG[^1] : _simApi.LiveEffectiveG;
                double lastThreshold = _seriesAdaptiveThreshold.Count > 0 ? _seriesAdaptiveThreshold[^1] : _simApi.LiveAdaptiveThreshold;
                int lastLargestCluster = _seriesLargestCluster.Count > 0 ? _seriesLargestCluster[^1] : _simApi.LiveLargestCluster;
                int lastStrongEdges = _seriesStrongEdges.Count > 0 ? _seriesStrongEdges[^1] : _simApi.LiveStrongEdges;

                // Get decimated dynamics (max 1000 points for ~50KB JSON)
                var dynamics = _simApi.GetDecimatedDynamics(1000);

                var snapshot = new
                {
                    meta = new
                    {
                        exportedAt = DateTime.UtcNow,
                        type = "snapshot",
                        isRunning = _isModernRunning,
                        totalDataPoints = _seriesSteps.Count,
                        decimatedPoints = dynamics.Steps.Length,
                        decimationStride = dynamics.DecimationStride
                    },
                    // GPU acceleration info
                    gpuInfo = new
                    {
                        enabled = CanUseGpu(),
                        available = _simApi.GpuAvailable,
                        deviceName = comboBox_GPUIndex.SelectedItem as string,
                        stats = new
                        {
                            kernelLaunches = _simApi.GpuStats.KernelLaunches,
                            topologyRebuilds = _simApi.GpuStats.TopologyRebuilds,
                            weightSyncs = _simApi.GpuStats.WeightSyncs
                        },
                        acceleratedOps = _simApi.GpuStats.AcceleratedOperations,
                        cpuBoundOps = _simApi.GpuStats.CpuBoundOperations
                    },
                    // Startup parameters
                    startupConfig = _simApi.LastConfig == null ? null : new
                    {
                        NodeCount = _simApi.LastConfig.NodeCount,
                        TotalSteps = _simApi.LastConfig.TotalSteps,
                        InitialEdgeProb = _simApi.LastConfig.InitialEdgeProb,
                        InitialExcitedProb = _simApi.LastConfig.InitialExcitedProb,
                        TargetDegree = _simApi.LastConfig.TargetDegree,
                        Seed = _simApi.LastConfig.Seed,
                        GravitationalCoupling = _simApi.LastConfig.GravitationalCoupling,
                        HotStartTemperature = _simApi.LastConfig.HotStartTemperature,
                        AnnealingCoolingRate = _simApi.LastConfig.AnnealingCoolingRate,
                        DecoherenceRate = _simApi.LastConfig.DecoherenceRate,
                        VacuumEnergyScale = _simApi.LastConfig.VacuumEnergyScale,
                        UseSpectralGeometry = _simApi.LastConfig.UseSpectralGeometry,
                        UseNetworkGravity = _simApi.LastConfig.UseNetworkGravity,
                        UseHotStartAnnealing = _simApi.LastConfig.UseHotStartAnnealing,
                        UseQuantumDrivenStates = _simApi.LastConfig.UseQuantumDrivenStates,
                        UseSpinorField = _simApi.LastConfig.UseSpinorField,
                        UseVacuumFluctuations = _simApi.LastConfig.UseVacuumFluctuations,
                        UseYangMillsGauge = _simApi.LastConfig.UseYangMillsGauge,
                        UseRelationalTime = _simApi.LastConfig.UseRelationalTime,
                        UseTopologicalProtection = _simApi.LastConfig.UseTopologicalProtection,
                        UseCausalRewiring = _simApi.LastConfig.UseCausalRewiring,
                        UseMexicanHatPotential = _simApi.LastConfig.UseMexicanHatPotential,
                        UseGeometryMomenta = _simApi.LastConfig.UseGeometryMomenta
                    },
                    // Physics constants
                    physicsConstants = new
                    {
                        FineStructureConstant = PhysicsConstants.FineStructureConstant,
                        StrongCouplingConstant = PhysicsConstants.StrongCouplingConstant,
                        GravitationalCoupling = PhysicsConstants.GravitationalCoupling,
                        WarmupGravitationalCoupling = PhysicsConstants.WarmupGravitationalCoupling,
                        WarmupDuration = PhysicsConstants.WarmupDuration,
                        GravityTransitionDuration = PhysicsConstants.GravityTransitionDuration,
                        InitialAnnealingTemperature = PhysicsConstants.InitialAnnealingTemperature,
                        FinalAnnealingTemperature = PhysicsConstants.FinalAnnealingTemperature,
                        AnnealingTimeConstant = PhysicsConstants.AnnealingTimeConstant,
                        AdaptiveThresholdSigma = PhysicsConstants.AdaptiveThresholdSigma,
                        DefaultHeavyClusterThreshold = PhysicsConstants.DefaultHeavyClusterThreshold,
                        MinimumClusterSize = PhysicsConstants.MinimumClusterSize,
                        CosmologicalConstant = PhysicsConstants.CosmologicalConstant
                    },
                    // Current state
                    currentState = new
                    {
                        step = lastStep,
                        totalStepsPlanned = _simApi.LastConfig?.TotalSteps ?? 0,
                        excited = lastExcited,
                        heavyMass = lastHeavyMass,
                        spectralDimension = lastSpectralDim,
                        networkTemperature = lastTemp,
                        effectiveG = lastEffectiveG,
                        adaptiveThreshold = lastThreshold,
                        largestCluster = lastLargestCluster,
                        strongEdges = lastStrongEdges
                    },
                    // Final metrics
                    finalMetrics = new
                    {
                        spectralDimension = _simApi.FinalSpectralDimension,
                        networkTemperature = _simApi.FinalNetworkTemperature,
                        wallClockSeconds = (DateTime.UtcNow - _simApi.SimulationWallClockStart).TotalSeconds
                    },
                    // Decimated dynamics (compact, ~1000 points max)
                    dynamics = new
                    {
                        steps = dynamics.Steps,
                        excited = dynamics.Excited,
                        energy = dynamics.Energy,
                        heavyMass = dynamics.HeavyMass,
                        largestCluster = dynamics.LargestCluster,
                        strongEdges = dynamics.StrongEdges,
                        spectralDimension = dynamics.SpectralDimension,
                        networkTemperature = dynamics.NetworkTemperature
                    }
                };

                var options = new JsonSerializerOptions { WriteIndented = true };
                var json = JsonSerializer.Serialize(snapshot, options);
                File.WriteAllText(dlg.FileName, json, Encoding.UTF8);
            }
            catch (Exception ex)
            {
                MessageBox.Show($"Ошибка сохранения: {ex.Message}", "Ошибка", MessageBoxButtons.OK, MessageBoxIcon.Error);
            }
        }

        private void checkBox_AutoTuning_CheckedChanged(object sender, EventArgs e)
        {
            _simApi.AutoTuningEnabled = checkBox_AutoTuning.Checked;
            AppendConsole($"[AutoTuning] {(checkBox_AutoTuning.Checked ? "включен" : "выключен")}\n");

            if (checkBox_AutoTuning.Checked)
            {
                AppendConsole("[AutoTuning] Критерии: d_S, кластеры, excited ratio, giant cluster\n");
                AppendConsole("[AutoTuning] Интервал: каждые 100 шагов (быстрая реакция)\n");
            }
        }

        private void Form_Main_Load(object sender, EventArgs e)
        {

        }
    }
}