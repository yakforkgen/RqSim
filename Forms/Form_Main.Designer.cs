namespace RqSimForms
{
    partial class Form_Main
    {
        // Add new UI controls backing fields
        private StatusStrip statusStrip;
        private ToolStripStatusLabel statusLabelSteps;
        private ToolStripStatusLabel statusLabelExcited;
        private ToolStripStatusLabel statusLabelHeavyMass;
        private ToolStrip toolStrip;
        private ToolStripButton toolRunButton;
        private ToolStripButton toolStopButton;
        private ToolStripButton toolKmcRunButton;
        private ToolStripButton toolKmcStopButton;
        private ToolStripSeparator toolSep1;
        private ToolStripDropDownButton toolViewDropDown;
        private ToolStripMenuItem toolViewGraph;
        private ToolStripMenuItem toolViewConsole;
        private ToolStripMenuItem toolViewCsv;
        private ToolStripMenuItem toolViewAnalysis;
        private ToolStripMenuItem toolViewSummary;
        private SplitContainer splitMain;
        private GroupBox grpRunStats;
        private TableLayoutPanel tlpRunStats;
        private Label lblTotalSteps;
        private Label valTotalSteps;
        private Label lblCurrentStep;
        private Label valCurrentStep;
        private Label lblExcitedAvg;
        private Label valExcitedAvg;
        private Label lblExcitedMax;
        private Label valExcitedMax;
        private Label lblAvalancheCount;
        private Label valAvalancheCount;
        private Label lblMeasurementStatus;
        private Label valMeasurementStatus;
        private GroupBox grpLiveMetrics;
        private TableLayoutPanel tlpLiveMetrics;
        private Label lblGlobalNbr;
        private Label valGlobalNbr;
        private Label lblGlobalSpont;
        private Label valGlobalSpont;
        private Label lblStrongEdges;
        private Label valStrongEdges;
        private Label lblLargestCluster;
        private Label valLargestCluster;
        private Label lblHeavyMass;
        private Label valHeavyMass;
        private Label lblSpectrumInfo;
        private Label valSpectrumInfo;
        private Label lblLightSpeed;
        private Label valLightSpeed;
        private ComboBox cmbWeightThreshold;
        private CheckBox chkShowHeavyOnly;
        private CheckBox chkAutoScrollConsole;
        private Button btnExpornShortJson;
        private Button btnSnapshotImage;
        private GroupBox grpEvents;
        private ListView lvEvents;
        private ColumnHeader colEventStep;
        private ColumnHeader colEventType;
        private ColumnHeader colEventDetail;
        // New tabs and chart panels
        private TabPage tabPage_OnExcited;
        private TabPage tabPage_Heavy;
        private TabPage tabPage_Cluster;
        private TabPage tabPage_Energy;
        private TabPage tabPage_Summary;
        private Panel panelOnChart;
        private Panel panelHeavyChart;
        private Panel panelClusterChart;
        private Panel panelEnergyChart;
        private TextBox summaryTextBox;
        private TextBox modernSimTextBox;

        // Existing fields
        private System.ComponentModel.IContainer components = null;
        private TabControl tabControl1;
        private TabPage tabPage_GUI;
        private TabPage tabPage_Console;
        private TextBox consoleTextBox;
        private TabPage tabPage_Sythnesis;
        private Button button_SaveSimReult;
        private Button button_RunModernSim;
        private TabPage tabPage_Settings;

        /// <summary>
        ///  Required designer variable.
        /// </summary>


        /// <summary>
        ///  Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }

            if (disposing)
            {
                canvasBitmap?.Dispose();
            }

            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        ///  Required method for Designer support - do not modify
        ///  the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            tabControl1 = new TabControl();
            tabPage_GUI = new TabPage();
            tabPage_OnExcited = new TabPage();
            panelOnChart = new Panel();
            tabPage_Heavy = new TabPage();
            panelHeavyChart = new Panel();
            tabPage_Cluster = new TabPage();
            panelClusterChart = new Panel();
            tabPage_Energy = new TabPage();
            panelEnergyChart = new Panel();
            tabPage_Console = new TabPage();
            consoleTextBox = new TextBox();
            tabPage_Sythnesis = new TabPage();
            tabPage_Summary = new TabPage();
            grpDashboard = new GroupBox();
            tlpDashboard = new TableLayoutPanel();
            lblDashNodes = new Label();
            valDashNodes = new Label();
            lblDashTotalSteps = new Label();
            valDashTotalSteps = new Label();
            lblDashCurrentStep = new Label();
            valDashCurrentStep = new Label();
            lblDashExcited = new Label();
            valDashExcited = new Label();
            lblDashHeavyMass = new Label();
            valDashHeavyMass = new Label();
            lblDashLargestCluster = new Label();
            valDashLargestCluster = new Label();
            lblDashStrongEdges = new Label();
            valDashStrongEdges = new Label();
            lblDashPhase = new Label();
            valDashPhase = new Label();
            lblDashQNorm = new Label();
            valDashQNorm = new Label();
            lblDashEntanglement = new Label();
            valDashEntanglement = new Label();
            lblDashCorrelation = new Label();
            valDashCorrelation = new Label();
            lblDashStatus = new Label();
            valDashStatus = new Label();
            lblDashSpectralDim = new Label();
            valDashSpectralDim = new Label();
            lblDashEffectiveG = new Label();
            valDashEffectiveG = new Label();
            lblDashGSuppression = new Label();
            valDashGSuppression = new Label();
            lblDashNetworkTemp = new Label();
            valDashNetworkTemp = new Label();
            grpEvents = new GroupBox();
            summaryTextBox = new TextBox();
            lvEvents = new ListView();
            colEventStep = new ColumnHeader();
            colEventType = new ColumnHeader();
            colEventDetail = new ColumnHeader();
            tabPage_Settings = new TabPage();
            settingsMainLayout = new TableLayoutPanel();
            grpSimParams = new GroupBox();
            tlpSimParams = new TableLayoutPanel();
            lblNodeCount = new Label();
            numNodeCount = new NumericUpDown();
            lblTargetDegree = new Label();
            numTargetDegree = new NumericUpDown();
            lblInitialExcitedProb = new Label();
            numInitialExcitedProb = new NumericUpDown();
            lblLambdaState = new Label();
            numLambdaState = new NumericUpDown();
            lblTemperature = new Label();
            numTemperature = new NumericUpDown();
            lblEdgeTrialProb = new Label();
            numEdgeTrialProb = new NumericUpDown();
            lblMeasurementThreshold = new Label();
            numMeasurementThreshold = new NumericUpDown();
            lblTotalStepsSettings = new Label();
            numTotalSteps = new NumericUpDown();
            lblFractalLevels = new Label();
            numFractalLevels = new NumericUpDown();
            numFractalBranchFactor = new NumericUpDown();
            lblFractalBranchFactor = new Label();
            grpPhysicsModules = new GroupBox();
            flpPhysics = new FlowLayoutPanel();
            chkQuantumDriven = new CheckBox();
            chkSpacetimePhysics = new CheckBox();
            chkSpinorField = new CheckBox();
            chkVacuumFluctuations = new CheckBox();
            chkBlackHolePhysics = new CheckBox();
            chkYangMillsGauge = new CheckBox();
            chkEnhancedKleinGordon = new CheckBox();
            chkInternalTime = new CheckBox();
            chkSpectralGeometry = new CheckBox();
            chkQuantumGraphity = new CheckBox();
            chkRelationalTime = new CheckBox();
            chkRelationalYangMills = new CheckBox();
            chkNetworkGravity = new CheckBox();
            chkUnifiedPhysicsStep = new CheckBox();
            chkEnforceGaugeConstraints = new CheckBox();
            chkCausalRewiring = new CheckBox();
            chkTopologicalProtection = new CheckBox();
            chkValidateEnergyConservation = new CheckBox();
            chkMexicanHatPotential = new CheckBox();
            chkHotStartAnnealing = new CheckBox();
            chkGeometryMomenta = new CheckBox();
            chkTopologicalCensorship = new CheckBox();
            grpPhysicsConstants = new GroupBox();
            tlpPhysicsConstants = new TableLayoutPanel();
            lblInitialEdgeProb = new Label();
            numInitialEdgeProb = new NumericUpDown();
            lblGravitationalCoupling = new Label();
            numGravitationalCoupling = new NumericUpDown();
            lblVacuumEnergyScale = new Label();
            numVacuumEnergyScale = new NumericUpDown();
            lblAnnealingCoolingRate = new Label();
            numAnnealingCoolingRate = new NumericUpDown();
            numDecoherenceRate = new NumericUpDown();
            numHotStartTemperature = new NumericUpDown();
            lblDecoherenceRate = new Label();
            lblHotStartTemperature = new Label();
            lblAdaptiveThresholdSigma = new Label();
            numAdaptiveThresholdSigma = new NumericUpDown();
            lblWarmupDuration = new Label();
            numWarmupDuration = new NumericUpDown();
            lblGravityTransitionDuration = new Label();
            numGravityTransitionDuration = new NumericUpDown();
            valAnnealingTimeConstant = new Label();
            tabPage_Experiments = new TabPage();
            statusStrip = new StatusStrip();
            statusLabelSteps = new ToolStripStatusLabel();
            statusLabelExcited = new ToolStripStatusLabel();
            statusLabelHeavyMass = new ToolStripStatusLabel();
            toolStrip = new ToolStrip();
            toolRunButton = new ToolStripButton();
            toolStopButton = new ToolStripButton();
            toolKmcRunButton = new ToolStripButton();
            toolKmcStopButton = new ToolStripButton();
            toolSep1 = new ToolStripSeparator();
            toolViewDropDown = new ToolStripDropDownButton();
            toolViewGraph = new ToolStripMenuItem();
            toolViewConsole = new ToolStripMenuItem();
            toolViewCsv = new ToolStripMenuItem();
            toolViewAnalysis = new ToolStripMenuItem();
            toolViewSummary = new ToolStripMenuItem();
            splitMain = new SplitContainer();
            label_ParamPresets = new Label();
            label_MaxFPS = new Label();
            label_CPUThreads = new Label();
            numericUpDown1 = new NumericUpDown();
            lblExperiments = new Label();
            comboBox_Experiments = new ComboBox();
            btnExpornShortJson = new Button();
            comboBox_Presets = new ComboBox();
            chkAutoScrollConsole = new CheckBox();
            checkBox_AutoTuning = new CheckBox();
            chkShowHeavyOnly = new CheckBox();
            numericUpDown_MaxFPS = new NumericUpDown();
            btnSnapshotImage = new Button();
            button_ForceRedrawGraphImage = new Button();
            button_SaveSimReult = new Button();
            comboBox_GPUIndex = new ComboBox();
            cmbWeightThreshold = new ComboBox();
            checkBox_EnableGPU = new CheckBox();
            grpLiveMetrics = new GroupBox();
            tlpLiveMetrics = new TableLayoutPanel();
            lblGlobalNbr = new Label();
            valGlobalNbr = new Label();
            lblGlobalSpont = new Label();
            valGlobalSpont = new Label();
            lblStrongEdges = new Label();
            valStrongEdges = new Label();
            lblLargestCluster = new Label();
            valLargestCluster = new Label();
            lblHeavyMass = new Label();
            valHeavyMass = new Label();
            lblSpectrumInfo = new Label();
            valSpectrumInfo = new Label();
            lblLightSpeed = new Label();
            valLightSpeed = new Label();
            grpRunStats = new GroupBox();
            tlpRunStats = new TableLayoutPanel();
            lblExcitedAvg = new Label();
            valExcitedAvg = new Label();
            lblExcitedMax = new Label();
            valExcitedMax = new Label();
            lblAvalancheCount = new Label();
            valAvalancheCount = new Label();
            lblMeasurementStatus = new Label();
            valMeasurementStatus = new Label();
            valCurrentStep = new Label();
            valTotalSteps = new Label();
            lblTotalSteps = new Label();
            lblCurrentStep = new Label();
            modernSimTextBox = new TextBox();
            button_RunModernSim = new Button();
            tabControl1.SuspendLayout();
            tabPage_OnExcited.SuspendLayout();
            tabPage_Heavy.SuspendLayout();
            tabPage_Cluster.SuspendLayout();
            tabPage_Energy.SuspendLayout();
            tabPage_Console.SuspendLayout();
            tabPage_Summary.SuspendLayout();
            grpDashboard.SuspendLayout();
            tlpDashboard.SuspendLayout();
            tabPage_Settings.SuspendLayout();
            settingsMainLayout.SuspendLayout();
            grpSimParams.SuspendLayout();
            tlpSimParams.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)numNodeCount).BeginInit();
            ((System.ComponentModel.ISupportInitialize)numTargetDegree).BeginInit();
            ((System.ComponentModel.ISupportInitialize)numInitialExcitedProb).BeginInit();
            ((System.ComponentModel.ISupportInitialize)numLambdaState).BeginInit();
            ((System.ComponentModel.ISupportInitialize)numTemperature).BeginInit();
            ((System.ComponentModel.ISupportInitialize)numEdgeTrialProb).BeginInit();
            ((System.ComponentModel.ISupportInitialize)numMeasurementThreshold).BeginInit();
            ((System.ComponentModel.ISupportInitialize)numTotalSteps).BeginInit();
            ((System.ComponentModel.ISupportInitialize)numFractalLevels).BeginInit();
            ((System.ComponentModel.ISupportInitialize)numFractalBranchFactor).BeginInit();
            grpPhysicsModules.SuspendLayout();
            flpPhysics.SuspendLayout();
            grpPhysicsConstants.SuspendLayout();
            tlpPhysicsConstants.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)numInitialEdgeProb).BeginInit();
            ((System.ComponentModel.ISupportInitialize)numGravitationalCoupling).BeginInit();
            ((System.ComponentModel.ISupportInitialize)numVacuumEnergyScale).BeginInit();
            ((System.ComponentModel.ISupportInitialize)numAnnealingCoolingRate).BeginInit();
            ((System.ComponentModel.ISupportInitialize)numDecoherenceRate).BeginInit();
            ((System.ComponentModel.ISupportInitialize)numHotStartTemperature).BeginInit();
            ((System.ComponentModel.ISupportInitialize)numAdaptiveThresholdSigma).BeginInit();
            ((System.ComponentModel.ISupportInitialize)numWarmupDuration).BeginInit();
            ((System.ComponentModel.ISupportInitialize)numGravityTransitionDuration).BeginInit();
            statusStrip.SuspendLayout();
            toolStrip.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)splitMain).BeginInit();
            splitMain.Panel1.SuspendLayout();
            splitMain.Panel2.SuspendLayout();
            splitMain.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)numericUpDown1).BeginInit();
            ((System.ComponentModel.ISupportInitialize)numericUpDown_MaxFPS).BeginInit();
            grpLiveMetrics.SuspendLayout();
            tlpLiveMetrics.SuspendLayout();
            grpRunStats.SuspendLayout();
            tlpRunStats.SuspendLayout();
            SuspendLayout();
            // 
            // tabControl1
            // 
            tabControl1.Controls.Add(tabPage_GUI);
            tabControl1.Controls.Add(tabPage_OnExcited);
            tabControl1.Controls.Add(tabPage_Heavy);
            tabControl1.Controls.Add(tabPage_Cluster);
            tabControl1.Controls.Add(tabPage_Energy);
            tabControl1.Controls.Add(tabPage_Console);
            tabControl1.Controls.Add(tabPage_Sythnesis);
            tabControl1.Controls.Add(tabPage_Summary);
            tabControl1.Controls.Add(tabPage_Settings);
            tabControl1.Controls.Add(tabPage_Experiments);
            tabControl1.Location = new Point(-3, 61);
            tabControl1.Name = "tabControl1";
            tabControl1.SelectedIndex = 0;
            tabControl1.Size = new Size(1350, 643);
            tabControl1.TabIndex = 0;
            // 
            // tabPage_GUI
            // 
            tabPage_GUI.Location = new Point(4, 24);
            tabPage_GUI.Name = "tabPage_GUI";
            tabPage_GUI.Padding = new Padding(3);
            tabPage_GUI.Size = new Size(1342, 615);
            tabPage_GUI.TabIndex = 0;
            tabPage_GUI.Text = "Network";
            tabPage_GUI.UseVisualStyleBackColor = true;
            // 
            // tabPage_OnExcited
            // 
            tabPage_OnExcited.Controls.Add(panelOnChart);
            tabPage_OnExcited.Location = new Point(4, 24);
            tabPage_OnExcited.Name = "tabPage_OnExcited";
            tabPage_OnExcited.Size = new Size(1342, 615);
            tabPage_OnExcited.TabIndex = 8;
            tabPage_OnExcited.Text = "Excited";
            tabPage_OnExcited.UseVisualStyleBackColor = true;
            // 
            // panelOnChart
            // 
            panelOnChart.BackColor = Color.White;
            panelOnChart.Dock = DockStyle.Fill;
            panelOnChart.Location = new Point(0, 0);
            panelOnChart.Name = "panelOnChart";
            panelOnChart.Size = new Size(1342, 615);
            panelOnChart.TabIndex = 0;
            // 
            // tabPage_Heavy
            // 
            tabPage_Heavy.Controls.Add(panelHeavyChart);
            tabPage_Heavy.Location = new Point(4, 24);
            tabPage_Heavy.Name = "tabPage_Heavy";
            tabPage_Heavy.Size = new Size(1342, 615);
            tabPage_Heavy.TabIndex = 9;
            tabPage_Heavy.Text = "HeavyMass";
            tabPage_Heavy.UseVisualStyleBackColor = true;
            // 
            // panelHeavyChart
            // 
            panelHeavyChart.BackColor = Color.White;
            panelHeavyChart.Dock = DockStyle.Fill;
            panelHeavyChart.Location = new Point(0, 0);
            panelHeavyChart.Name = "panelHeavyChart";
            panelHeavyChart.Size = new Size(1342, 615);
            panelHeavyChart.TabIndex = 0;
            // 
            // tabPage_Cluster
            // 
            tabPage_Cluster.Controls.Add(panelClusterChart);
            tabPage_Cluster.Location = new Point(4, 24);
            tabPage_Cluster.Name = "tabPage_Cluster";
            tabPage_Cluster.Size = new Size(1342, 615);
            tabPage_Cluster.TabIndex = 10;
            tabPage_Cluster.Text = "Clusters";
            tabPage_Cluster.UseVisualStyleBackColor = true;
            // 
            // panelClusterChart
            // 
            panelClusterChart.BackColor = Color.White;
            panelClusterChart.Dock = DockStyle.Fill;
            panelClusterChart.Location = new Point(0, 0);
            panelClusterChart.Name = "panelClusterChart";
            panelClusterChart.Size = new Size(1342, 615);
            panelClusterChart.TabIndex = 0;
            // 
            // tabPage_Energy
            // 
            tabPage_Energy.Controls.Add(panelEnergyChart);
            tabPage_Energy.Location = new Point(4, 24);
            tabPage_Energy.Name = "tabPage_Energy";
            tabPage_Energy.Size = new Size(1342, 615);
            tabPage_Energy.TabIndex = 11;
            tabPage_Energy.Text = "Energy";
            tabPage_Energy.UseVisualStyleBackColor = true;
            // 
            // panelEnergyChart
            // 
            panelEnergyChart.BackColor = Color.White;
            panelEnergyChart.Dock = DockStyle.Fill;
            panelEnergyChart.Location = new Point(0, 0);
            panelEnergyChart.Name = "panelEnergyChart";
            panelEnergyChart.Size = new Size(1342, 615);
            panelEnergyChart.TabIndex = 0;
            // 
            // tabPage_Console
            // 
            tabPage_Console.Controls.Add(consoleTextBox);
            tabPage_Console.Location = new Point(4, 24);
            tabPage_Console.Name = "tabPage_Console";
            tabPage_Console.Padding = new Padding(3);
            tabPage_Console.Size = new Size(1342, 615);
            tabPage_Console.TabIndex = 1;
            tabPage_Console.Text = "Console";
            tabPage_Console.UseVisualStyleBackColor = true;
            // 
            // consoleTextBox
            // 
            consoleTextBox.BackColor = Color.Black;
            consoleTextBox.Dock = DockStyle.Fill;
            consoleTextBox.Font = new Font("Consolas", 9F);
            consoleTextBox.ForeColor = Color.Lime;
            consoleTextBox.Location = new Point(3, 3);
            consoleTextBox.Multiline = true;
            consoleTextBox.Name = "consoleTextBox";
            consoleTextBox.ReadOnly = true;
            consoleTextBox.ScrollBars = ScrollBars.Both;
            consoleTextBox.Size = new Size(1336, 609);
            consoleTextBox.TabIndex = 0;
            consoleTextBox.WordWrap = false;
            // 
            // tabPage_Sythnesis
            // 
            tabPage_Sythnesis.Location = new Point(4, 24);
            tabPage_Sythnesis.Name = "tabPage_Sythnesis";
            tabPage_Sythnesis.Size = new Size(1342, 615);
            tabPage_Sythnesis.TabIndex = 3;
            tabPage_Sythnesis.Text = "Synthesis";
            tabPage_Sythnesis.UseVisualStyleBackColor = true;
            // 
            // tabPage_Summary
            // 
            tabPage_Summary.Controls.Add(grpDashboard);
            tabPage_Summary.Controls.Add(grpEvents);
            tabPage_Summary.Controls.Add(summaryTextBox);
            tabPage_Summary.Controls.Add(lvEvents);
            tabPage_Summary.Location = new Point(4, 24);
            tabPage_Summary.Name = "tabPage_Summary";
            tabPage_Summary.Size = new Size(1342, 615);
            tabPage_Summary.TabIndex = 12;
            tabPage_Summary.Text = "Summary";
            tabPage_Summary.UseVisualStyleBackColor = true;
            // 
            // grpDashboard
            // 
            grpDashboard.Anchor = AnchorStyles.None;
            grpDashboard.Controls.Add(tlpDashboard);
            grpDashboard.Location = new Point(140, 29);
            grpDashboard.Margin = new Padding(5);
            grpDashboard.Name = "grpDashboard";
            grpDashboard.Size = new Size(791, 360);
            grpDashboard.TabIndex = 9;
            grpDashboard.TabStop = false;
            grpDashboard.Text = "Real-Time Dashboard";
            // 
            // tlpDashboard
            // 
            tlpDashboard.ColumnCount = 2;
            tlpDashboard.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 60F));
            tlpDashboard.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 40F));
            tlpDashboard.Controls.Add(lblDashNodes, 0, 0);
            tlpDashboard.Controls.Add(valDashNodes, 1, 0);
            tlpDashboard.Controls.Add(lblDashTotalSteps, 0, 1);
            tlpDashboard.Controls.Add(valDashTotalSteps, 1, 1);
            tlpDashboard.Controls.Add(lblDashCurrentStep, 0, 2);
            tlpDashboard.Controls.Add(valDashCurrentStep, 1, 2);
            tlpDashboard.Controls.Add(lblDashExcited, 0, 3);
            tlpDashboard.Controls.Add(valDashExcited, 1, 3);
            tlpDashboard.Controls.Add(lblDashHeavyMass, 0, 4);
            tlpDashboard.Controls.Add(valDashHeavyMass, 1, 4);
            tlpDashboard.Controls.Add(lblDashLargestCluster, 0, 5);
            tlpDashboard.Controls.Add(valDashLargestCluster, 1, 5);
            tlpDashboard.Controls.Add(lblDashStrongEdges, 0, 6);
            tlpDashboard.Controls.Add(valDashStrongEdges, 1, 6);
            tlpDashboard.Controls.Add(lblDashPhase, 0, 7);
            tlpDashboard.Controls.Add(valDashPhase, 1, 7);
            tlpDashboard.Controls.Add(lblDashQNorm, 0, 8);
            tlpDashboard.Controls.Add(valDashQNorm, 1, 8);
            tlpDashboard.Controls.Add(lblDashEntanglement, 0, 9);
            tlpDashboard.Controls.Add(valDashEntanglement, 1, 9);
            tlpDashboard.Controls.Add(lblDashCorrelation, 0, 10);
            tlpDashboard.Controls.Add(valDashCorrelation, 1, 10);
            tlpDashboard.Controls.Add(lblDashStatus, 0, 11);
            tlpDashboard.Controls.Add(valDashStatus, 1, 11);
            tlpDashboard.Controls.Add(lblDashSpectralDim, 0, 12);
            tlpDashboard.Controls.Add(valDashSpectralDim, 1, 12);
            tlpDashboard.Controls.Add(lblDashEffectiveG, 0, 13);
            tlpDashboard.Controls.Add(valDashEffectiveG, 1, 13);
            tlpDashboard.Controls.Add(lblDashGSuppression, 0, 14);
            tlpDashboard.Controls.Add(valDashGSuppression, 1, 14);
            tlpDashboard.Controls.Add(lblDashNetworkTemp, 0, 15);
            tlpDashboard.Controls.Add(valDashNetworkTemp, 1, 15);
            tlpDashboard.Dock = DockStyle.Fill;
            tlpDashboard.Location = new Point(3, 19);
            tlpDashboard.Name = "tlpDashboard";
            tlpDashboard.RowCount = 16;
            tlpDashboard.RowStyles.Add(new RowStyle(SizeType.Absolute, 20F));
            tlpDashboard.RowStyles.Add(new RowStyle(SizeType.Absolute, 20F));
            tlpDashboard.RowStyles.Add(new RowStyle(SizeType.Absolute, 20F));
            tlpDashboard.RowStyles.Add(new RowStyle(SizeType.Absolute, 20F));
            tlpDashboard.RowStyles.Add(new RowStyle(SizeType.Absolute, 20F));
            tlpDashboard.RowStyles.Add(new RowStyle(SizeType.Absolute, 20F));
            tlpDashboard.RowStyles.Add(new RowStyle(SizeType.Absolute, 20F));
            tlpDashboard.RowStyles.Add(new RowStyle(SizeType.Absolute, 20F));
            tlpDashboard.RowStyles.Add(new RowStyle(SizeType.Absolute, 20F));
            tlpDashboard.RowStyles.Add(new RowStyle(SizeType.Absolute, 20F));
            tlpDashboard.RowStyles.Add(new RowStyle(SizeType.Absolute, 20F));
            tlpDashboard.RowStyles.Add(new RowStyle(SizeType.Absolute, 20F));
            tlpDashboard.RowStyles.Add(new RowStyle(SizeType.Absolute, 20F));
            tlpDashboard.RowStyles.Add(new RowStyle(SizeType.Absolute, 20F));
            tlpDashboard.RowStyles.Add(new RowStyle(SizeType.Absolute, 20F));
            tlpDashboard.RowStyles.Add(new RowStyle(SizeType.Absolute, 20F));
            tlpDashboard.Size = new Size(785, 338);
            tlpDashboard.TabIndex = 0;
            // 
            // lblDashNodes
            // 
            lblDashNodes.Anchor = AnchorStyles.Left | AnchorStyles.Right;
            lblDashNodes.AutoSize = true;
            lblDashNodes.Location = new Point(3, 2);
            lblDashNodes.Name = "lblDashNodes";
            lblDashNodes.Size = new Size(465, 15);
            lblDashNodes.TabIndex = 0;
            lblDashNodes.Text = "Nodes:";
            // 
            // valDashNodes
            // 
            valDashNodes.AutoSize = true;
            valDashNodes.Location = new Point(474, 0);
            valDashNodes.Name = "valDashNodes";
            valDashNodes.Size = new Size(13, 15);
            valDashNodes.TabIndex = 1;
            valDashNodes.Text = "0";
            // 
            // lblDashTotalSteps
            // 
            lblDashTotalSteps.Anchor = AnchorStyles.Left | AnchorStyles.Right;
            lblDashTotalSteps.AutoSize = true;
            lblDashTotalSteps.Location = new Point(3, 22);
            lblDashTotalSteps.Name = "lblDashTotalSteps";
            lblDashTotalSteps.Size = new Size(465, 15);
            lblDashTotalSteps.TabIndex = 2;
            lblDashTotalSteps.Text = "Total Steps:";
            // 
            // valDashTotalSteps
            // 
            valDashTotalSteps.AutoSize = true;
            valDashTotalSteps.Location = new Point(474, 20);
            valDashTotalSteps.Name = "valDashTotalSteps";
            valDashTotalSteps.Size = new Size(13, 15);
            valDashTotalSteps.TabIndex = 3;
            valDashTotalSteps.Text = "0";
            // 
            // lblDashCurrentStep
            // 
            lblDashCurrentStep.Anchor = AnchorStyles.Left | AnchorStyles.Right;
            lblDashCurrentStep.AutoSize = true;
            lblDashCurrentStep.Location = new Point(3, 42);
            lblDashCurrentStep.Name = "lblDashCurrentStep";
            lblDashCurrentStep.Size = new Size(465, 15);
            lblDashCurrentStep.TabIndex = 4;
            lblDashCurrentStep.Text = "Current Step:";
            // 
            // valDashCurrentStep
            // 
            valDashCurrentStep.AutoSize = true;
            valDashCurrentStep.Location = new Point(474, 40);
            valDashCurrentStep.Name = "valDashCurrentStep";
            valDashCurrentStep.Size = new Size(13, 15);
            valDashCurrentStep.TabIndex = 5;
            valDashCurrentStep.Text = "0";
            // 
            // lblDashExcited
            // 
            lblDashExcited.Anchor = AnchorStyles.Left | AnchorStyles.Right;
            lblDashExcited.AutoSize = true;
            lblDashExcited.Location = new Point(3, 62);
            lblDashExcited.Name = "lblDashExcited";
            lblDashExcited.Size = new Size(465, 15);
            lblDashExcited.TabIndex = 6;
            lblDashExcited.Text = "Excited:";
            // 
            // valDashExcited
            // 
            valDashExcited.AutoSize = true;
            valDashExcited.Location = new Point(474, 60);
            valDashExcited.Name = "valDashExcited";
            valDashExcited.Size = new Size(13, 15);
            valDashExcited.TabIndex = 7;
            valDashExcited.Text = "0";
            // 
            // lblDashHeavyMass
            // 
            lblDashHeavyMass.Anchor = AnchorStyles.Left | AnchorStyles.Right;
            lblDashHeavyMass.AutoSize = true;
            lblDashHeavyMass.Location = new Point(3, 82);
            lblDashHeavyMass.Name = "lblDashHeavyMass";
            lblDashHeavyMass.Size = new Size(465, 15);
            lblDashHeavyMass.TabIndex = 8;
            lblDashHeavyMass.Text = "Heavy Mass:";
            // 
            // valDashHeavyMass
            // 
            valDashHeavyMass.AutoSize = true;
            valDashHeavyMass.Location = new Point(474, 80);
            valDashHeavyMass.Name = "valDashHeavyMass";
            valDashHeavyMass.Size = new Size(28, 15);
            valDashHeavyMass.TabIndex = 9;
            valDashHeavyMass.Text = "0.00";
            // 
            // lblDashLargestCluster
            // 
            lblDashLargestCluster.Anchor = AnchorStyles.Left | AnchorStyles.Right;
            lblDashLargestCluster.AutoSize = true;
            lblDashLargestCluster.Location = new Point(3, 102);
            lblDashLargestCluster.Name = "lblDashLargestCluster";
            lblDashLargestCluster.Size = new Size(465, 15);
            lblDashLargestCluster.TabIndex = 10;
            lblDashLargestCluster.Text = "Largest Cluster:";
            // 
            // valDashLargestCluster
            // 
            valDashLargestCluster.AutoSize = true;
            valDashLargestCluster.Location = new Point(474, 100);
            valDashLargestCluster.Name = "valDashLargestCluster";
            valDashLargestCluster.Size = new Size(13, 15);
            valDashLargestCluster.TabIndex = 11;
            valDashLargestCluster.Text = "0";
            // 
            // lblDashStrongEdges
            // 
            lblDashStrongEdges.Anchor = AnchorStyles.Left | AnchorStyles.Right;
            lblDashStrongEdges.AutoSize = true;
            lblDashStrongEdges.Location = new Point(3, 122);
            lblDashStrongEdges.Name = "lblDashStrongEdges";
            lblDashStrongEdges.Size = new Size(465, 15);
            lblDashStrongEdges.TabIndex = 12;
            lblDashStrongEdges.Text = "Strong Edges:";
            // 
            // valDashStrongEdges
            // 
            valDashStrongEdges.AutoSize = true;
            valDashStrongEdges.Location = new Point(474, 120);
            valDashStrongEdges.Name = "valDashStrongEdges";
            valDashStrongEdges.Size = new Size(13, 15);
            valDashStrongEdges.TabIndex = 13;
            valDashStrongEdges.Text = "0";
            // 
            // lblDashPhase
            // 
            lblDashPhase.Anchor = AnchorStyles.Left | AnchorStyles.Right;
            lblDashPhase.AutoSize = true;
            lblDashPhase.Location = new Point(3, 142);
            lblDashPhase.Name = "lblDashPhase";
            lblDashPhase.Size = new Size(465, 15);
            lblDashPhase.TabIndex = 14;
            lblDashPhase.Text = "Phase:";
            // 
            // valDashPhase
            // 
            valDashPhase.AutoSize = true;
            valDashPhase.Location = new Point(474, 140);
            valDashPhase.Name = "valDashPhase";
            valDashPhase.Size = new Size(36, 15);
            valDashPhase.TabIndex = 15;
            valDashPhase.Text = "Quiet";
            // 
            // lblDashQNorm
            // 
            lblDashQNorm.Anchor = AnchorStyles.Left | AnchorStyles.Right;
            lblDashQNorm.AutoSize = true;
            lblDashQNorm.Location = new Point(3, 162);
            lblDashQNorm.Name = "lblDashQNorm";
            lblDashQNorm.Size = new Size(465, 15);
            lblDashQNorm.TabIndex = 16;
            lblDashQNorm.Text = "Q-Norm:";
            // 
            // valDashQNorm
            // 
            valDashQNorm.AutoSize = true;
            valDashQNorm.Location = new Point(474, 160);
            valDashQNorm.Name = "valDashQNorm";
            valDashQNorm.Size = new Size(34, 15);
            valDashQNorm.TabIndex = 17;
            valDashQNorm.Text = "0.000";
            // 
            // lblDashEntanglement
            // 
            lblDashEntanglement.Anchor = AnchorStyles.Left | AnchorStyles.Right;
            lblDashEntanglement.AutoSize = true;
            lblDashEntanglement.Location = new Point(3, 182);
            lblDashEntanglement.Name = "lblDashEntanglement";
            lblDashEntanglement.Size = new Size(465, 15);
            lblDashEntanglement.TabIndex = 18;
            lblDashEntanglement.Text = "Entanglement:";
            // 
            // valDashEntanglement
            // 
            valDashEntanglement.AutoSize = true;
            valDashEntanglement.Location = new Point(474, 180);
            valDashEntanglement.Name = "valDashEntanglement";
            valDashEntanglement.Size = new Size(34, 15);
            valDashEntanglement.TabIndex = 19;
            valDashEntanglement.Text = "0.000";
            // 
            // lblDashCorrelation
            // 
            lblDashCorrelation.Anchor = AnchorStyles.Left | AnchorStyles.Right;
            lblDashCorrelation.AutoSize = true;
            lblDashCorrelation.Location = new Point(3, 202);
            lblDashCorrelation.Name = "lblDashCorrelation";
            lblDashCorrelation.Size = new Size(465, 15);
            lblDashCorrelation.TabIndex = 20;
            lblDashCorrelation.Text = "Correlation:";
            // 
            // valDashCorrelation
            // 
            valDashCorrelation.AutoSize = true;
            valDashCorrelation.Location = new Point(474, 200);
            valDashCorrelation.Name = "valDashCorrelation";
            valDashCorrelation.Size = new Size(34, 15);
            valDashCorrelation.TabIndex = 21;
            valDashCorrelation.Text = "0.000";
            // 
            // lblDashStatus
            // 
            lblDashStatus.Anchor = AnchorStyles.Left | AnchorStyles.Right;
            lblDashStatus.AutoSize = true;
            lblDashStatus.Location = new Point(3, 222);
            lblDashStatus.Name = "lblDashStatus";
            lblDashStatus.Size = new Size(465, 15);
            lblDashStatus.TabIndex = 22;
            lblDashStatus.Text = "Status:";
            // 
            // valDashStatus
            // 
            valDashStatus.AutoSize = true;
            valDashStatus.Location = new Point(474, 220);
            valDashStatus.Name = "valDashStatus";
            valDashStatus.Size = new Size(39, 15);
            valDashStatus.TabIndex = 23;
            valDashStatus.Text = "Ready";
            // 
            // lblDashSpectralDim
            // 
            lblDashSpectralDim.Location = new Point(3, 240);
            lblDashSpectralDim.Name = "lblDashSpectralDim";
            lblDashSpectralDim.Size = new Size(100, 20);
            lblDashSpectralDim.TabIndex = 24;
            lblDashSpectralDim.Text = "SpectralDim";
            // 
            // valDashSpectralDim
            // 
            valDashSpectralDim.Location = new Point(474, 240);
            valDashSpectralDim.Name = "valDashSpectralDim";
            valDashSpectralDim.Size = new Size(100, 20);
            valDashSpectralDim.TabIndex = 25;
            // 
            // lblDashEffectiveG
            // 
            lblDashEffectiveG.Location = new Point(3, 260);
            lblDashEffectiveG.Name = "lblDashEffectiveG";
            lblDashEffectiveG.Size = new Size(100, 20);
            lblDashEffectiveG.TabIndex = 26;
            lblDashEffectiveG.Text = "EffectiveG";
            // 
            // valDashEffectiveG
            // 
            valDashEffectiveG.Location = new Point(474, 260);
            valDashEffectiveG.Name = "valDashEffectiveG";
            valDashEffectiveG.Size = new Size(100, 20);
            valDashEffectiveG.TabIndex = 27;
            // 
            // lblDashGSuppression
            // 
            lblDashGSuppression.Location = new Point(3, 280);
            lblDashGSuppression.Name = "lblDashGSuppression";
            lblDashGSuppression.Size = new Size(100, 20);
            lblDashGSuppression.TabIndex = 28;
            lblDashGSuppression.Text = "GSuppression";
            // 
            // valDashGSuppression
            // 
            valDashGSuppression.Location = new Point(474, 280);
            valDashGSuppression.Name = "valDashGSuppression";
            valDashGSuppression.Size = new Size(100, 20);
            valDashGSuppression.TabIndex = 29;
            // 
            // lblDashNetworkTemp
            // 
            lblDashNetworkTemp.Location = new Point(3, 300);
            lblDashNetworkTemp.Name = "lblDashNetworkTemp";
            lblDashNetworkTemp.Size = new Size(100, 23);
            lblDashNetworkTemp.TabIndex = 30;
            lblDashNetworkTemp.Text = "NetworkTemp";
            // 
            // valDashNetworkTemp
            // 
            valDashNetworkTemp.Location = new Point(474, 300);
            valDashNetworkTemp.Name = "valDashNetworkTemp";
            valDashNetworkTemp.Size = new Size(100, 23);
            valDashNetworkTemp.TabIndex = 31;
            // 
            // grpEvents
            // 
            grpEvents.Location = new Point(8, 396);
            grpEvents.Name = "grpEvents";
            grpEvents.Size = new Size(1044, 175);
            grpEvents.TabIndex = 8;
            grpEvents.TabStop = false;
            grpEvents.Text = "Important Events";
            // 
            // summaryTextBox
            // 
            summaryTextBox.BackColor = Color.White;
            summaryTextBox.Dock = DockStyle.Fill;
            summaryTextBox.Font = new Font("Consolas", 9F);
            summaryTextBox.ForeColor = Color.Black;
            summaryTextBox.Location = new Point(0, 0);
            summaryTextBox.Multiline = true;
            summaryTextBox.Name = "summaryTextBox";
            summaryTextBox.ReadOnly = true;
            summaryTextBox.ScrollBars = ScrollBars.Both;
            summaryTextBox.Size = new Size(1342, 615);
            summaryTextBox.TabIndex = 0;
            summaryTextBox.WordWrap = false;
            // 
            // lvEvents
            // 
            lvEvents.Anchor = AnchorStyles.None;
            lvEvents.Columns.AddRange(new ColumnHeader[] { colEventStep, colEventType, colEventDetail });
            lvEvents.FullRowSelect = true;
            lvEvents.Location = new Point(126, 17);
            lvEvents.Name = "lvEvents";
            lvEvents.Size = new Size(835, 594);
            lvEvents.TabIndex = 0;
            lvEvents.UseCompatibleStateImageBehavior = false;
            lvEvents.View = View.Details;
            // 
            // colEventStep
            // 
            colEventStep.Text = "Step";
            // 
            // colEventType
            // 
            colEventType.Text = "Type";
            colEventType.Width = 80;
            // 
            // colEventDetail
            // 
            colEventDetail.Text = "Detail";
            colEventDetail.Width = 130;
            // 
            // tabPage_Settings
            // 
            tabPage_Settings.AutoScroll = true;
            tabPage_Settings.Controls.Add(settingsMainLayout);
            tabPage_Settings.Location = new Point(4, 24);
            tabPage_Settings.Name = "tabPage_Settings";
            tabPage_Settings.Size = new Size(1342, 615);
            tabPage_Settings.TabIndex = 14;
            tabPage_Settings.Text = "Settings";
            tabPage_Settings.UseVisualStyleBackColor = true;
            // 
            // settingsMainLayout
            // 
            settingsMainLayout.Anchor = AnchorStyles.None;
            settingsMainLayout.ColumnCount = 3;
            settingsMainLayout.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 33F));
            settingsMainLayout.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 34F));
            settingsMainLayout.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 33F));
            settingsMainLayout.Controls.Add(grpSimParams, 0, 0);
            settingsMainLayout.Controls.Add(grpPhysicsModules, 1, 0);
            settingsMainLayout.Controls.Add(grpPhysicsConstants, 2, 0);
            settingsMainLayout.Location = new Point(216, 0);
            settingsMainLayout.Name = "settingsMainLayout";
            settingsMainLayout.RowCount = 1;
            settingsMainLayout.RowStyles.Add(new RowStyle(SizeType.Percent, 100F));
            settingsMainLayout.Size = new Size(994, 580);
            settingsMainLayout.TabIndex = 0;
            // 
            // grpSimParams
            // 
            grpSimParams.Controls.Add(tlpSimParams);
            grpSimParams.Dock = DockStyle.Fill;
            grpSimParams.Location = new Point(5, 5);
            grpSimParams.Margin = new Padding(5);
            grpSimParams.Name = "grpSimParams";
            grpSimParams.Size = new Size(318, 570);
            grpSimParams.TabIndex = 0;
            grpSimParams.TabStop = false;
            grpSimParams.Text = "Simulation Parameters";
            // 
            // tlpSimParams
            // 
            tlpSimParams.AutoSize = true;
            tlpSimParams.ColumnCount = 2;
            tlpSimParams.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 60F));
            tlpSimParams.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 40F));
            tlpSimParams.Controls.Add(lblNodeCount, 0, 0);
            tlpSimParams.Controls.Add(numNodeCount, 1, 0);
            tlpSimParams.Controls.Add(lblTargetDegree, 0, 1);
            tlpSimParams.Controls.Add(numTargetDegree, 1, 1);
            tlpSimParams.Controls.Add(lblInitialExcitedProb, 0, 2);
            tlpSimParams.Controls.Add(numInitialExcitedProb, 1, 2);
            tlpSimParams.Controls.Add(lblLambdaState, 0, 3);
            tlpSimParams.Controls.Add(numLambdaState, 1, 3);
            tlpSimParams.Controls.Add(lblTemperature, 0, 4);
            tlpSimParams.Controls.Add(numTemperature, 1, 4);
            tlpSimParams.Controls.Add(lblEdgeTrialProb, 0, 5);
            tlpSimParams.Controls.Add(numEdgeTrialProb, 1, 5);
            tlpSimParams.Controls.Add(lblMeasurementThreshold, 0, 6);
            tlpSimParams.Controls.Add(numMeasurementThreshold, 1, 6);
            tlpSimParams.Controls.Add(lblTotalStepsSettings, 0, 7);
            tlpSimParams.Controls.Add(numTotalSteps, 1, 7);
            tlpSimParams.Controls.Add(lblFractalLevels, 0, 8);
            tlpSimParams.Controls.Add(numFractalLevels, 1, 8);
            tlpSimParams.Controls.Add(numFractalBranchFactor, 1, 9);
            tlpSimParams.Controls.Add(lblFractalBranchFactor, 0, 9);
            tlpSimParams.Location = new Point(3, 19);
            tlpSimParams.Margin = new Padding(4, 3, 3, 3);
            tlpSimParams.Name = "tlpSimParams";
            tlpSimParams.RowCount = 10;
            tlpSimParams.RowStyles.Add(new RowStyle(SizeType.Absolute, 20F));
            tlpSimParams.RowStyles.Add(new RowStyle(SizeType.Absolute, 20F));
            tlpSimParams.RowStyles.Add(new RowStyle(SizeType.Absolute, 20F));
            tlpSimParams.RowStyles.Add(new RowStyle(SizeType.Absolute, 20F));
            tlpSimParams.RowStyles.Add(new RowStyle(SizeType.Absolute, 20F));
            tlpSimParams.RowStyles.Add(new RowStyle(SizeType.Absolute, 20F));
            tlpSimParams.RowStyles.Add(new RowStyle(SizeType.Absolute, 20F));
            tlpSimParams.RowStyles.Add(new RowStyle(SizeType.Absolute, 20F));
            tlpSimParams.RowStyles.Add(new RowStyle(SizeType.Absolute, 20F));
            tlpSimParams.RowStyles.Add(new RowStyle(SizeType.Absolute, 20F));
            tlpSimParams.Size = new Size(312, 548);
            tlpSimParams.TabIndex = 0;
            // 
            // lblNodeCount
            // 
            lblNodeCount.Anchor = AnchorStyles.Left | AnchorStyles.Right;
            lblNodeCount.AutoSize = true;
            lblNodeCount.Location = new Point(3, 2);
            lblNodeCount.Name = "lblNodeCount";
            lblNodeCount.Size = new Size(181, 15);
            lblNodeCount.TabIndex = 0;
            lblNodeCount.Text = "Node Count:";
            // 
            // numNodeCount
            // 
            numNodeCount.Dock = DockStyle.Fill;
            numNodeCount.Location = new Point(190, 3);
            numNodeCount.Maximum = new decimal(new int[] { 1000000, 0, 0, 0 });
            numNodeCount.Minimum = new decimal(new int[] { 1, 0, 0, 0 });
            numNodeCount.Name = "numNodeCount";
            numNodeCount.Size = new Size(119, 23);
            numNodeCount.TabIndex = 1;
            numNodeCount.Value = new decimal(new int[] { 250, 0, 0, 0 });
            // 
            // lblTargetDegree
            // 
            lblTargetDegree.Anchor = AnchorStyles.Left | AnchorStyles.Right;
            lblTargetDegree.AutoSize = true;
            lblTargetDegree.Location = new Point(3, 22);
            lblTargetDegree.Name = "lblTargetDegree";
            lblTargetDegree.Size = new Size(181, 15);
            lblTargetDegree.TabIndex = 2;
            lblTargetDegree.Text = "Target Degree:";
            // 
            // numTargetDegree
            // 
            numTargetDegree.Dock = DockStyle.Fill;
            numTargetDegree.Location = new Point(190, 23);
            numTargetDegree.Maximum = new decimal(new int[] { 20, 0, 0, 0 });
            numTargetDegree.Minimum = new decimal(new int[] { 2, 0, 0, 0 });
            numTargetDegree.Name = "numTargetDegree";
            numTargetDegree.Size = new Size(119, 23);
            numTargetDegree.TabIndex = 3;
            numTargetDegree.Value = new decimal(new int[] { 8, 0, 0, 0 });
            // 
            // lblInitialExcitedProb
            // 
            lblInitialExcitedProb.Anchor = AnchorStyles.Left | AnchorStyles.Right;
            lblInitialExcitedProb.AutoSize = true;
            lblInitialExcitedProb.Location = new Point(3, 42);
            lblInitialExcitedProb.Name = "lblInitialExcitedProb";
            lblInitialExcitedProb.Size = new Size(181, 15);
            lblInitialExcitedProb.TabIndex = 4;
            lblInitialExcitedProb.Text = "Initial Excited Prob:";
            // 
            // numInitialExcitedProb
            // 
            numInitialExcitedProb.DecimalPlaces = 2;
            numInitialExcitedProb.Dock = DockStyle.Fill;
            numInitialExcitedProb.Increment = new decimal(new int[] { 5, 0, 0, 131072 });
            numInitialExcitedProb.Location = new Point(190, 43);
            numInitialExcitedProb.Maximum = new decimal(new int[] { 1, 0, 0, 0 });
            numInitialExcitedProb.Name = "numInitialExcitedProb";
            numInitialExcitedProb.Size = new Size(119, 23);
            numInitialExcitedProb.TabIndex = 5;
            numInitialExcitedProb.Value = new decimal(new int[] { 10, 0, 0, 131072 });
            // 
            // lblLambdaState
            // 
            lblLambdaState.Anchor = AnchorStyles.Left | AnchorStyles.Right;
            lblLambdaState.AutoSize = true;
            lblLambdaState.Location = new Point(3, 62);
            lblLambdaState.Name = "lblLambdaState";
            lblLambdaState.Size = new Size(181, 15);
            lblLambdaState.TabIndex = 6;
            lblLambdaState.Text = "Lambda State:";
            // 
            // numLambdaState
            // 
            numLambdaState.DecimalPlaces = 2;
            numLambdaState.Dock = DockStyle.Fill;
            numLambdaState.Increment = new decimal(new int[] { 1, 0, 0, 65536 });
            numLambdaState.Location = new Point(190, 63);
            numLambdaState.Maximum = new decimal(new int[] { 2, 0, 0, 0 });
            numLambdaState.Name = "numLambdaState";
            numLambdaState.Size = new Size(119, 23);
            numLambdaState.TabIndex = 7;
            numLambdaState.Value = new decimal(new int[] { 5, 0, 0, 65536 });
            // 
            // lblTemperature
            // 
            lblTemperature.Anchor = AnchorStyles.Left | AnchorStyles.Right;
            lblTemperature.AutoSize = true;
            lblTemperature.Location = new Point(3, 82);
            lblTemperature.Name = "lblTemperature";
            lblTemperature.Size = new Size(181, 15);
            lblTemperature.TabIndex = 8;
            lblTemperature.Text = "Temperature:";
            // 
            // numTemperature
            // 
            numTemperature.DecimalPlaces = 2;
            numTemperature.Dock = DockStyle.Fill;
            numTemperature.Increment = new decimal(new int[] { 1, 0, 0, 65536 });
            numTemperature.Location = new Point(190, 83);
            numTemperature.Name = "numTemperature";
            numTemperature.Size = new Size(119, 23);
            numTemperature.TabIndex = 9;
            numTemperature.Value = new decimal(new int[] { 100, 0, 0, 65536 });
            // 
            // lblEdgeTrialProb
            // 
            lblEdgeTrialProb.Anchor = AnchorStyles.Left | AnchorStyles.Right;
            lblEdgeTrialProb.AutoSize = true;
            lblEdgeTrialProb.Location = new Point(3, 102);
            lblEdgeTrialProb.Name = "lblEdgeTrialProb";
            lblEdgeTrialProb.Size = new Size(181, 15);
            lblEdgeTrialProb.TabIndex = 10;
            lblEdgeTrialProb.Text = "Edge Trial Prob:";
            // 
            // numEdgeTrialProb
            // 
            numEdgeTrialProb.DecimalPlaces = 3;
            numEdgeTrialProb.Dock = DockStyle.Fill;
            numEdgeTrialProb.Increment = new decimal(new int[] { 1, 0, 0, 131072 });
            numEdgeTrialProb.Location = new Point(190, 103);
            numEdgeTrialProb.Maximum = new decimal(new int[] { 1, 0, 0, 0 });
            numEdgeTrialProb.Name = "numEdgeTrialProb";
            numEdgeTrialProb.Size = new Size(119, 23);
            numEdgeTrialProb.TabIndex = 11;
            numEdgeTrialProb.Value = new decimal(new int[] { 2, 0, 0, 131072 });
            // 
            // lblMeasurementThreshold
            // 
            lblMeasurementThreshold.Anchor = AnchorStyles.Left | AnchorStyles.Right;
            lblMeasurementThreshold.AutoSize = true;
            lblMeasurementThreshold.Location = new Point(3, 122);
            lblMeasurementThreshold.Name = "lblMeasurementThreshold";
            lblMeasurementThreshold.Size = new Size(181, 15);
            lblMeasurementThreshold.TabIndex = 12;
            lblMeasurementThreshold.Text = "Measurement Threshold:";
            // 
            // numMeasurementThreshold
            // 
            numMeasurementThreshold.DecimalPlaces = 3;
            numMeasurementThreshold.Dock = DockStyle.Fill;
            numMeasurementThreshold.Increment = new decimal(new int[] { 1, 0, 0, 131072 });
            numMeasurementThreshold.Location = new Point(190, 123);
            numMeasurementThreshold.Maximum = new decimal(new int[] { 1, 0, 0, 0 });
            numMeasurementThreshold.Name = "numMeasurementThreshold";
            numMeasurementThreshold.Size = new Size(119, 23);
            numMeasurementThreshold.TabIndex = 13;
            numMeasurementThreshold.Value = new decimal(new int[] { 30, 0, 0, 131072 });
            // 
            // lblTotalStepsSettings
            // 
            lblTotalStepsSettings.Anchor = AnchorStyles.Left | AnchorStyles.Right;
            lblTotalStepsSettings.AutoSize = true;
            lblTotalStepsSettings.Location = new Point(3, 142);
            lblTotalStepsSettings.Name = "lblTotalStepsSettings";
            lblTotalStepsSettings.Size = new Size(181, 15);
            lblTotalStepsSettings.TabIndex = 14;
            lblTotalStepsSettings.Text = "Total Steps:";
            // 
            // numTotalSteps
            // 
            numTotalSteps.Dock = DockStyle.Fill;
            numTotalSteps.Increment = new decimal(new int[] { 100, 0, 0, 0 });
            numTotalSteps.Location = new Point(190, 143);
            numTotalSteps.Maximum = new decimal(new int[] { 10000000, 0, 0, 0 });
            numTotalSteps.Minimum = new decimal(new int[] { 100, 0, 0, 0 });
            numTotalSteps.Name = "numTotalSteps";
            numTotalSteps.Size = new Size(119, 23);
            numTotalSteps.TabIndex = 15;
            numTotalSteps.Value = new decimal(new int[] { 500000, 0, 0, 0 });
            // 
            // lblFractalLevels
            // 
            lblFractalLevels.Anchor = AnchorStyles.Left | AnchorStyles.Right;
            lblFractalLevels.AutoSize = true;
            lblFractalLevels.Location = new Point(4, 162);
            lblFractalLevels.Margin = new Padding(4, 0, 4, 0);
            lblFractalLevels.Name = "lblFractalLevels";
            lblFractalLevels.Size = new Size(179, 15);
            lblFractalLevels.TabIndex = 16;
            lblFractalLevels.Text = "Fractal Levels:";
            // 
            // numFractalLevels
            // 
            numFractalLevels.Dock = DockStyle.Fill;
            numFractalLevels.Location = new Point(190, 163);
            numFractalLevels.Maximum = new decimal(new int[] { 10, 0, 0, 0 });
            numFractalLevels.Name = "numFractalLevels";
            numFractalLevels.Size = new Size(119, 23);
            numFractalLevels.TabIndex = 17;
            // 
            // numFractalBranchFactor
            // 
            numFractalBranchFactor.Dock = DockStyle.Fill;
            numFractalBranchFactor.Location = new Point(190, 183);
            numFractalBranchFactor.Maximum = new decimal(new int[] { 10, 0, 0, 0 });
            numFractalBranchFactor.Name = "numFractalBranchFactor";
            numFractalBranchFactor.Size = new Size(119, 23);
            numFractalBranchFactor.TabIndex = 19;
            // 
            // lblFractalBranchFactor
            // 
            lblFractalBranchFactor.AutoSize = true;
            lblFractalBranchFactor.Location = new Point(3, 180);
            lblFractalBranchFactor.Name = "lblFractalBranchFactor";
            lblFractalBranchFactor.Size = new Size(121, 15);
            lblFractalBranchFactor.TabIndex = 18;
            lblFractalBranchFactor.Text = "Fractal Branch Factor:";
            // 
            // grpPhysicsModules
            // 
            grpPhysicsModules.Controls.Add(flpPhysics);
            grpPhysicsModules.Dock = DockStyle.Fill;
            grpPhysicsModules.Location = new Point(333, 5);
            grpPhysicsModules.Margin = new Padding(5);
            grpPhysicsModules.Name = "grpPhysicsModules";
            grpPhysicsModules.Size = new Size(327, 570);
            grpPhysicsModules.TabIndex = 1;
            grpPhysicsModules.TabStop = false;
            grpPhysicsModules.Text = "Physics Modules";
            // 
            // flpPhysics
            // 
            flpPhysics.AutoSize = true;
            flpPhysics.Controls.Add(chkQuantumDriven);
            flpPhysics.Controls.Add(chkSpacetimePhysics);
            flpPhysics.Controls.Add(chkSpinorField);
            flpPhysics.Controls.Add(chkVacuumFluctuations);
            flpPhysics.Controls.Add(chkBlackHolePhysics);
            flpPhysics.Controls.Add(chkYangMillsGauge);
            flpPhysics.Controls.Add(chkEnhancedKleinGordon);
            flpPhysics.Controls.Add(chkInternalTime);
            flpPhysics.Controls.Add(chkSpectralGeometry);
            flpPhysics.Controls.Add(chkQuantumGraphity);
            flpPhysics.Controls.Add(chkRelationalTime);
            flpPhysics.Controls.Add(chkRelationalYangMills);
            flpPhysics.Controls.Add(chkNetworkGravity);
            flpPhysics.Controls.Add(chkUnifiedPhysicsStep);
            flpPhysics.Controls.Add(chkEnforceGaugeConstraints);
            flpPhysics.Controls.Add(chkCausalRewiring);
            flpPhysics.Controls.Add(chkTopologicalProtection);
            flpPhysics.Controls.Add(chkValidateEnergyConservation);
            flpPhysics.Controls.Add(chkMexicanHatPotential);
            flpPhysics.Controls.Add(chkHotStartAnnealing);
            flpPhysics.Controls.Add(chkGeometryMomenta);
            flpPhysics.Controls.Add(chkTopologicalCensorship);
            flpPhysics.Dock = DockStyle.Fill;
            flpPhysics.FlowDirection = FlowDirection.TopDown;
            flpPhysics.Location = new Point(3, 19);
            flpPhysics.Name = "flpPhysics";
            flpPhysics.Size = new Size(321, 548);
            flpPhysics.TabIndex = 0;
            // 
            // chkQuantumDriven
            // 
            chkQuantumDriven.AutoSize = true;
            chkQuantumDriven.Checked = true;
            chkQuantumDriven.CheckState = CheckState.Checked;
            chkQuantumDriven.Location = new Point(3, 3);
            chkQuantumDriven.Name = "chkQuantumDriven";
            chkQuantumDriven.Size = new Size(148, 19);
            chkQuantumDriven.TabIndex = 0;
            chkQuantumDriven.Text = "Quantum Driven States";
            // 
            // chkSpacetimePhysics
            // 
            chkSpacetimePhysics.AutoSize = true;
            chkSpacetimePhysics.Checked = true;
            chkSpacetimePhysics.CheckState = CheckState.Checked;
            chkSpacetimePhysics.Location = new Point(3, 28);
            chkSpacetimePhysics.Name = "chkSpacetimePhysics";
            chkSpacetimePhysics.Size = new Size(123, 19);
            chkSpacetimePhysics.TabIndex = 1;
            chkSpacetimePhysics.Text = "Spacetime Physics";
            // 
            // chkSpinorField
            // 
            chkSpinorField.AutoSize = true;
            chkSpinorField.Checked = true;
            chkSpinorField.CheckState = CheckState.Checked;
            chkSpinorField.Location = new Point(3, 53);
            chkSpinorField.Name = "chkSpinorField";
            chkSpinorField.Size = new Size(88, 19);
            chkSpinorField.TabIndex = 2;
            chkSpinorField.Text = "Spinor Field";
            // 
            // chkVacuumFluctuations
            // 
            chkVacuumFluctuations.AutoSize = true;
            chkVacuumFluctuations.Checked = true;
            chkVacuumFluctuations.CheckState = CheckState.Checked;
            chkVacuumFluctuations.Location = new Point(3, 78);
            chkVacuumFluctuations.Name = "chkVacuumFluctuations";
            chkVacuumFluctuations.Size = new Size(137, 19);
            chkVacuumFluctuations.TabIndex = 3;
            chkVacuumFluctuations.Text = "Vacuum Fluctuations";
            // 
            // chkBlackHolePhysics
            // 
            chkBlackHolePhysics.AutoSize = true;
            chkBlackHolePhysics.Checked = true;
            chkBlackHolePhysics.CheckState = CheckState.Checked;
            chkBlackHolePhysics.Location = new Point(3, 103);
            chkBlackHolePhysics.Name = "chkBlackHolePhysics";
            chkBlackHolePhysics.Size = new Size(124, 19);
            chkBlackHolePhysics.TabIndex = 4;
            chkBlackHolePhysics.Text = "Black Hole Physics";
            // 
            // chkYangMillsGauge
            // 
            chkYangMillsGauge.AutoSize = true;
            chkYangMillsGauge.Checked = true;
            chkYangMillsGauge.CheckState = CheckState.Checked;
            chkYangMillsGauge.Location = new Point(3, 128);
            chkYangMillsGauge.Name = "chkYangMillsGauge";
            chkYangMillsGauge.Size = new Size(119, 19);
            chkYangMillsGauge.TabIndex = 5;
            chkYangMillsGauge.Text = "Yang-Mills Gauge";
            // 
            // chkEnhancedKleinGordon
            // 
            chkEnhancedKleinGordon.AutoSize = true;
            chkEnhancedKleinGordon.Checked = true;
            chkEnhancedKleinGordon.CheckState = CheckState.Checked;
            chkEnhancedKleinGordon.Location = new Point(3, 153);
            chkEnhancedKleinGordon.Name = "chkEnhancedKleinGordon";
            chkEnhancedKleinGordon.Size = new Size(152, 19);
            chkEnhancedKleinGordon.TabIndex = 6;
            chkEnhancedKleinGordon.Text = "Enhanced Klein-Gordon";
            // 
            // chkInternalTime
            // 
            chkInternalTime.AutoSize = true;
            chkInternalTime.Checked = true;
            chkInternalTime.CheckState = CheckState.Checked;
            chkInternalTime.Location = new Point(3, 178);
            chkInternalTime.Name = "chkInternalTime";
            chkInternalTime.Size = new Size(186, 19);
            chkInternalTime.TabIndex = 7;
            chkInternalTime.Text = "Internal Time (Page-Wootters)";
            // 
            // chkSpectralGeometry
            // 
            chkSpectralGeometry.AutoSize = true;
            chkSpectralGeometry.Checked = true;
            chkSpectralGeometry.CheckState = CheckState.Checked;
            chkSpectralGeometry.Location = new Point(3, 203);
            chkSpectralGeometry.Name = "chkSpectralGeometry";
            chkSpectralGeometry.Size = new Size(123, 19);
            chkSpectralGeometry.TabIndex = 8;
            chkSpectralGeometry.Text = "Spectral Geometry";
            // 
            // chkQuantumGraphity
            // 
            chkQuantumGraphity.AutoSize = true;
            chkQuantumGraphity.Checked = true;
            chkQuantumGraphity.CheckState = CheckState.Checked;
            chkQuantumGraphity.Location = new Point(3, 228);
            chkQuantumGraphity.Name = "chkQuantumGraphity";
            chkQuantumGraphity.Size = new Size(125, 19);
            chkQuantumGraphity.TabIndex = 9;
            chkQuantumGraphity.Text = "Quantum Graphity";
            // 
            // chkRelationalTime
            // 
            chkRelationalTime.AutoSize = true;
            chkRelationalTime.Checked = true;
            chkRelationalTime.CheckState = CheckState.Checked;
            chkRelationalTime.Location = new Point(3, 253);
            chkRelationalTime.Name = "chkRelationalTime";
            chkRelationalTime.Size = new Size(108, 19);
            chkRelationalTime.TabIndex = 10;
            chkRelationalTime.Text = "Relational Time";
            // 
            // chkRelationalYangMills
            // 
            chkRelationalYangMills.AutoSize = true;
            chkRelationalYangMills.Checked = true;
            chkRelationalYangMills.CheckState = CheckState.Checked;
            chkRelationalYangMills.Location = new Point(3, 278);
            chkRelationalYangMills.Name = "chkRelationalYangMills";
            chkRelationalYangMills.Size = new Size(137, 19);
            chkRelationalYangMills.TabIndex = 11;
            chkRelationalYangMills.Text = "Relational Yang-Mills";
            // 
            // chkNetworkGravity
            // 
            chkNetworkGravity.AutoSize = true;
            chkNetworkGravity.Checked = true;
            chkNetworkGravity.CheckState = CheckState.Checked;
            chkNetworkGravity.Location = new Point(3, 303);
            chkNetworkGravity.Name = "chkNetworkGravity";
            chkNetworkGravity.Size = new Size(111, 19);
            chkNetworkGravity.TabIndex = 12;
            chkNetworkGravity.Text = "Network Gravity";
            // 
            // chkUnifiedPhysicsStep
            // 
            chkUnifiedPhysicsStep.AutoSize = true;
            chkUnifiedPhysicsStep.Checked = true;
            chkUnifiedPhysicsStep.CheckState = CheckState.Checked;
            chkUnifiedPhysicsStep.Location = new Point(3, 328);
            chkUnifiedPhysicsStep.Name = "chkUnifiedPhysicsStep";
            chkUnifiedPhysicsStep.Size = new Size(132, 19);
            chkUnifiedPhysicsStep.TabIndex = 13;
            chkUnifiedPhysicsStep.Text = "Unified Physics Step";
            // 
            // chkEnforceGaugeConstraints
            // 
            chkEnforceGaugeConstraints.AutoSize = true;
            chkEnforceGaugeConstraints.Checked = true;
            chkEnforceGaugeConstraints.CheckState = CheckState.Checked;
            chkEnforceGaugeConstraints.Location = new Point(3, 353);
            chkEnforceGaugeConstraints.Name = "chkEnforceGaugeConstraints";
            chkEnforceGaugeConstraints.Size = new Size(166, 19);
            chkEnforceGaugeConstraints.TabIndex = 14;
            chkEnforceGaugeConstraints.Text = "Enforce Gauge Constraints";
            // 
            // chkCausalRewiring
            // 
            chkCausalRewiring.AutoSize = true;
            chkCausalRewiring.Checked = true;
            chkCausalRewiring.CheckState = CheckState.Checked;
            chkCausalRewiring.Location = new Point(3, 378);
            chkCausalRewiring.Name = "chkCausalRewiring";
            chkCausalRewiring.Size = new Size(110, 19);
            chkCausalRewiring.TabIndex = 15;
            chkCausalRewiring.Text = "Causal Rewiring";
            // 
            // chkTopologicalProtection
            // 
            chkTopologicalProtection.AutoSize = true;
            chkTopologicalProtection.Checked = true;
            chkTopologicalProtection.CheckState = CheckState.Checked;
            chkTopologicalProtection.Location = new Point(3, 403);
            chkTopologicalProtection.Name = "chkTopologicalProtection";
            chkTopologicalProtection.Size = new Size(146, 19);
            chkTopologicalProtection.TabIndex = 16;
            chkTopologicalProtection.Text = "Topological Protection";
            // 
            // chkValidateEnergyConservation
            // 
            chkValidateEnergyConservation.AutoSize = true;
            chkValidateEnergyConservation.Checked = true;
            chkValidateEnergyConservation.CheckState = CheckState.Checked;
            chkValidateEnergyConservation.Location = new Point(3, 428);
            chkValidateEnergyConservation.Name = "chkValidateEnergyConservation";
            chkValidateEnergyConservation.Size = new Size(179, 19);
            chkValidateEnergyConservation.TabIndex = 17;
            chkValidateEnergyConservation.Text = "Validate Energy Conservation";
            // 
            // chkMexicanHatPotential
            // 
            chkMexicanHatPotential.AutoSize = true;
            chkMexicanHatPotential.Checked = true;
            chkMexicanHatPotential.CheckState = CheckState.Checked;
            chkMexicanHatPotential.Location = new Point(3, 453);
            chkMexicanHatPotential.Name = "chkMexicanHatPotential";
            chkMexicanHatPotential.Size = new Size(142, 19);
            chkMexicanHatPotential.TabIndex = 18;
            chkMexicanHatPotential.Text = "Mexican Hat Potential";
            // 
            // chkHotStartAnnealing
            // 
            chkHotStartAnnealing.AutoSize = true;
            chkHotStartAnnealing.Checked = true;
            chkHotStartAnnealing.CheckState = CheckState.Checked;
            chkHotStartAnnealing.Location = new Point(3, 478);
            chkHotStartAnnealing.Name = "chkHotStartAnnealing";
            chkHotStartAnnealing.Size = new Size(130, 19);
            chkHotStartAnnealing.TabIndex = 19;
            chkHotStartAnnealing.Text = "Hot Start Annealing";
            // 
            // chkGeometryMomenta
            // 
            chkGeometryMomenta.AutoSize = true;
            chkGeometryMomenta.Checked = true;
            chkGeometryMomenta.CheckState = CheckState.Checked;
            chkGeometryMomenta.Location = new Point(3, 503);
            chkGeometryMomenta.Name = "chkGeometryMomenta";
            chkGeometryMomenta.Size = new Size(133, 19);
            chkGeometryMomenta.TabIndex = 20;
            chkGeometryMomenta.Text = "Geometry Momenta";
            // 
            // chkTopologicalCensorship
            // 
            chkTopologicalCensorship.AutoSize = true;
            chkTopologicalCensorship.Checked = true;
            chkTopologicalCensorship.CheckState = CheckState.Checked;
            chkTopologicalCensorship.Location = new Point(195, 3);
            chkTopologicalCensorship.Name = "chkTopologicalCensorship";
            chkTopologicalCensorship.Size = new Size(150, 19);
            chkTopologicalCensorship.TabIndex = 21;
            chkTopologicalCensorship.Text = "Topological Censorship";
            // 
            // grpPhysicsConstants
            // 
            grpPhysicsConstants.Controls.Add(tlpPhysicsConstants);
            grpPhysicsConstants.Dock = DockStyle.Fill;
            grpPhysicsConstants.Location = new Point(670, 5);
            grpPhysicsConstants.Margin = new Padding(5);
            grpPhysicsConstants.Name = "grpPhysicsConstants";
            grpPhysicsConstants.Size = new Size(319, 570);
            grpPhysicsConstants.TabIndex = 2;
            grpPhysicsConstants.TabStop = false;
            grpPhysicsConstants.Text = "Physics Constants";
            // 
            // tlpPhysicsConstants
            // 
            tlpPhysicsConstants.AutoSize = true;
            tlpPhysicsConstants.ColumnCount = 2;
            tlpPhysicsConstants.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 73.3333359F));
            tlpPhysicsConstants.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 26.666666F));
            tlpPhysicsConstants.Controls.Add(lblInitialEdgeProb, 0, 0);
            tlpPhysicsConstants.Controls.Add(numInitialEdgeProb, 1, 0);
            tlpPhysicsConstants.Controls.Add(lblGravitationalCoupling, 0, 1);
            tlpPhysicsConstants.Controls.Add(numGravitationalCoupling, 1, 1);
            tlpPhysicsConstants.Controls.Add(lblVacuumEnergyScale, 0, 2);
            tlpPhysicsConstants.Controls.Add(numVacuumEnergyScale, 1, 2);
            tlpPhysicsConstants.Controls.Add(lblAnnealingCoolingRate, 0, 3);
            tlpPhysicsConstants.Controls.Add(numAnnealingCoolingRate, 1, 3);
            tlpPhysicsConstants.Controls.Add(numDecoherenceRate, 1, 4);
            tlpPhysicsConstants.Controls.Add(numHotStartTemperature, 1, 5);
            tlpPhysicsConstants.Controls.Add(lblDecoherenceRate, 0, 4);
            tlpPhysicsConstants.Controls.Add(lblHotStartTemperature, 0, 5);
            tlpPhysicsConstants.Controls.Add(lblAdaptiveThresholdSigma, 0, 6);
            tlpPhysicsConstants.Controls.Add(numAdaptiveThresholdSigma, 1, 6);
            tlpPhysicsConstants.Controls.Add(lblWarmupDuration, 0, 7);
            tlpPhysicsConstants.Controls.Add(numWarmupDuration, 1, 7);
            tlpPhysicsConstants.Controls.Add(lblGravityTransitionDuration, 0, 8);
            tlpPhysicsConstants.Controls.Add(numGravityTransitionDuration, 1, 8);
            tlpPhysicsConstants.Controls.Add(valAnnealingTimeConstant, 0, 9);
            tlpPhysicsConstants.Dock = DockStyle.Fill;
            tlpPhysicsConstants.Location = new Point(3, 19);
            tlpPhysicsConstants.Name = "tlpPhysicsConstants";
            tlpPhysicsConstants.RowCount = 10;
            tlpPhysicsConstants.RowStyles.Add(new RowStyle(SizeType.Absolute, 28F));
            tlpPhysicsConstants.RowStyles.Add(new RowStyle(SizeType.Absolute, 28F));
            tlpPhysicsConstants.RowStyles.Add(new RowStyle(SizeType.Absolute, 28F));
            tlpPhysicsConstants.RowStyles.Add(new RowStyle(SizeType.Absolute, 28F));
            tlpPhysicsConstants.RowStyles.Add(new RowStyle(SizeType.Absolute, 28F));
            tlpPhysicsConstants.RowStyles.Add(new RowStyle(SizeType.Absolute, 28F));
            tlpPhysicsConstants.RowStyles.Add(new RowStyle(SizeType.Absolute, 28F));
            tlpPhysicsConstants.RowStyles.Add(new RowStyle(SizeType.Absolute, 28F));
            tlpPhysicsConstants.RowStyles.Add(new RowStyle(SizeType.Absolute, 28F));
            tlpPhysicsConstants.RowStyles.Add(new RowStyle(SizeType.Absolute, 28F));
            tlpPhysicsConstants.Size = new Size(313, 548);
            tlpPhysicsConstants.TabIndex = 0;
            // 
            // lblInitialEdgeProb
            // 
            lblInitialEdgeProb.Anchor = AnchorStyles.Left | AnchorStyles.Right;
            lblInitialEdgeProb.AutoSize = true;
            lblInitialEdgeProb.Location = new Point(3, 6);
            lblInitialEdgeProb.Name = "lblInitialEdgeProb";
            lblInitialEdgeProb.Size = new Size(223, 15);
            lblInitialEdgeProb.TabIndex = 0;
            lblInitialEdgeProb.Text = "Initial Edge Prob:";
            // 
            // numInitialEdgeProb
            // 
            numInitialEdgeProb.DecimalPlaces = 4;
            numInitialEdgeProb.Dock = DockStyle.Fill;
            numInitialEdgeProb.Increment = new decimal(new int[] { 1, 0, 0, 131072 });
            numInitialEdgeProb.Location = new Point(232, 3);
            numInitialEdgeProb.Maximum = new decimal(new int[] { 1, 0, 0, 0 });
            numInitialEdgeProb.Name = "numInitialEdgeProb";
            numInitialEdgeProb.Size = new Size(78, 23);
            numInitialEdgeProb.TabIndex = 1;
            numInitialEdgeProb.Value = new decimal(new int[] { 35, 0, 0, 196608 });
            // 
            // lblGravitationalCoupling
            // 
            lblGravitationalCoupling.Anchor = AnchorStyles.Left | AnchorStyles.Right;
            lblGravitationalCoupling.AutoSize = true;
            lblGravitationalCoupling.Location = new Point(3, 34);
            lblGravitationalCoupling.Name = "lblGravitationalCoupling";
            lblGravitationalCoupling.Size = new Size(223, 15);
            lblGravitationalCoupling.TabIndex = 2;
            lblGravitationalCoupling.Text = "Gravitational Coupling (G):";
            // 
            // numGravitationalCoupling
            // 
            numGravitationalCoupling.DecimalPlaces = 4;
            numGravitationalCoupling.Dock = DockStyle.Fill;
            numGravitationalCoupling.Increment = new decimal(new int[] { 1, 0, 0, 131072 });
            numGravitationalCoupling.Location = new Point(232, 31);
            numGravitationalCoupling.Name = "numGravitationalCoupling";
            numGravitationalCoupling.Size = new Size(78, 23);
            numGravitationalCoupling.TabIndex = 3;
            numGravitationalCoupling.Value = new decimal(new int[] { 10, 0, 0, 196608 });
            // 
            // lblVacuumEnergyScale
            // 
            lblVacuumEnergyScale.Anchor = AnchorStyles.Left | AnchorStyles.Right;
            lblVacuumEnergyScale.AutoSize = true;
            lblVacuumEnergyScale.Location = new Point(3, 62);
            lblVacuumEnergyScale.Name = "lblVacuumEnergyScale";
            lblVacuumEnergyScale.Size = new Size(223, 15);
            lblVacuumEnergyScale.TabIndex = 4;
            lblVacuumEnergyScale.Text = "Vacuum Energy Scale:";
            // 
            // numVacuumEnergyScale
            // 
            numVacuumEnergyScale.DecimalPlaces = 4;
            numVacuumEnergyScale.Dock = DockStyle.Fill;
            numVacuumEnergyScale.Increment = new decimal(new int[] { 1, 0, 0, 131072 });
            numVacuumEnergyScale.Location = new Point(232, 59);
            numVacuumEnergyScale.Maximum = new decimal(new int[] { 1, 0, 0, 0 });
            numVacuumEnergyScale.Name = "numVacuumEnergyScale";
            numVacuumEnergyScale.Size = new Size(78, 23);
            numVacuumEnergyScale.TabIndex = 5;
            numVacuumEnergyScale.Value = new decimal(new int[] { 5, 0, 0, 327680 });
            // 
            // lblAnnealingCoolingRate
            // 
            lblAnnealingCoolingRate.Anchor = AnchorStyles.Left | AnchorStyles.Right;
            lblAnnealingCoolingRate.AutoSize = true;
            lblAnnealingCoolingRate.Location = new Point(3, 90);
            lblAnnealingCoolingRate.Name = "lblAnnealingCoolingRate";
            lblAnnealingCoolingRate.Size = new Size(223, 15);
            lblAnnealingCoolingRate.TabIndex = 6;
            lblAnnealingCoolingRate.Text = "Annealing Cooling Rate:";
            // 
            // numAnnealingCoolingRate
            // 
            numAnnealingCoolingRate.DecimalPlaces = 4;
            numAnnealingCoolingRate.Dock = DockStyle.Fill;
            numAnnealingCoolingRate.Increment = new decimal(new int[] { 1, 0, 0, 262144 });
            numAnnealingCoolingRate.Location = new Point(232, 87);
            numAnnealingCoolingRate.Maximum = new decimal(new int[] { 1, 0, 0, 0 });
            numAnnealingCoolingRate.Name = "numAnnealingCoolingRate";
            numAnnealingCoolingRate.Size = new Size(78, 23);
            numAnnealingCoolingRate.TabIndex = 7;
            numAnnealingCoolingRate.Value = new decimal(new int[] { 995, 0, 0, 196608 });
            // 
            // numDecoherenceRate
            // 
            numDecoherenceRate.DecimalPlaces = 4;
            numDecoherenceRate.Dock = DockStyle.Fill;
            numDecoherenceRate.Increment = new decimal(new int[] { 1, 0, 0, 262144 });
            numDecoherenceRate.Location = new Point(232, 115);
            numDecoherenceRate.Maximum = new decimal(new int[] { 1, 0, 0, 65536 });
            numDecoherenceRate.Name = "numDecoherenceRate";
            numDecoherenceRate.Size = new Size(78, 23);
            numDecoherenceRate.TabIndex = 9;
            numDecoherenceRate.Value = new decimal(new int[] { 5, 0, 0, 196608 });
            // 
            // numHotStartTemperature
            // 
            numHotStartTemperature.DecimalPlaces = 1;
            numHotStartTemperature.Dock = DockStyle.Fill;
            numHotStartTemperature.Increment = new decimal(new int[] { 5, 0, 0, 65536 });
            numHotStartTemperature.Location = new Point(232, 143);
            numHotStartTemperature.Name = "numHotStartTemperature";
            numHotStartTemperature.Size = new Size(78, 23);
            numHotStartTemperature.TabIndex = 11;
            numHotStartTemperature.Value = new decimal(new int[] { 6, 0, 0, 0 });
            // 
            // lblDecoherenceRate
            // 
            lblDecoherenceRate.Anchor = AnchorStyles.Left;
            lblDecoherenceRate.AutoSize = true;
            lblDecoherenceRate.Location = new Point(3, 118);
            lblDecoherenceRate.Name = "lblDecoherenceRate";
            lblDecoherenceRate.Size = new Size(105, 15);
            lblDecoherenceRate.TabIndex = 8;
            lblDecoherenceRate.Text = "Decoherence Rate:";
            // 
            // lblHotStartTemperature
            // 
            lblHotStartTemperature.Anchor = AnchorStyles.Left;
            lblHotStartTemperature.AutoSize = true;
            lblHotStartTemperature.Location = new Point(3, 146);
            lblHotStartTemperature.Name = "lblHotStartTemperature";
            lblHotStartTemperature.Size = new Size(127, 15);
            lblHotStartTemperature.TabIndex = 10;
            lblHotStartTemperature.Text = "Hot Start Temperature:";
            // 
            // lblAdaptiveThresholdSigma
            // 
            lblAdaptiveThresholdSigma.Anchor = AnchorStyles.Left | AnchorStyles.Right;
            lblAdaptiveThresholdSigma.AutoSize = true;
            lblAdaptiveThresholdSigma.Location = new Point(3, 174);
            lblAdaptiveThresholdSigma.Name = "lblAdaptiveThresholdSigma";
            lblAdaptiveThresholdSigma.Size = new Size(223, 15);
            lblAdaptiveThresholdSigma.TabIndex = 12;
            lblAdaptiveThresholdSigma.Text = "Adaptive Threshold σ:";
            // 
            // numAdaptiveThresholdSigma
            // 
            numAdaptiveThresholdSigma.DecimalPlaces = 2;
            numAdaptiveThresholdSigma.Dock = DockStyle.Fill;
            numAdaptiveThresholdSigma.Increment = new decimal(new int[] { 1, 0, 0, 65536 });
            numAdaptiveThresholdSigma.Location = new Point(232, 171);
            numAdaptiveThresholdSigma.Maximum = new decimal(new int[] { 5, 0, 0, 0 });
            numAdaptiveThresholdSigma.Name = "numAdaptiveThresholdSigma";
            numAdaptiveThresholdSigma.Size = new Size(78, 23);
            numAdaptiveThresholdSigma.TabIndex = 13;
            numAdaptiveThresholdSigma.Value = new decimal(new int[] { 15, 0, 0, 65536 });
            // 
            // lblWarmupDuration
            // 
            lblWarmupDuration.Anchor = AnchorStyles.Left | AnchorStyles.Right;
            lblWarmupDuration.AutoSize = true;
            lblWarmupDuration.Location = new Point(3, 202);
            lblWarmupDuration.Name = "lblWarmupDuration";
            lblWarmupDuration.Size = new Size(223, 15);
            lblWarmupDuration.TabIndex = 14;
            lblWarmupDuration.Text = "Warmup Duration:";
            // 
            // numWarmupDuration
            // 
            numWarmupDuration.Dock = DockStyle.Fill;
            numWarmupDuration.Increment = new decimal(new int[] { 10, 0, 0, 0 });
            numWarmupDuration.Location = new Point(232, 199);
            numWarmupDuration.Maximum = new decimal(new int[] { 1000, 0, 0, 0 });
            numWarmupDuration.Name = "numWarmupDuration";
            numWarmupDuration.Size = new Size(78, 23);
            numWarmupDuration.TabIndex = 15;
            numWarmupDuration.Value = new decimal(new int[] { 200, 0, 0, 0 });
            // 
            // lblGravityTransitionDuration
            // 
            lblGravityTransitionDuration.Anchor = AnchorStyles.Left | AnchorStyles.Right;
            lblGravityTransitionDuration.AutoSize = true;
            lblGravityTransitionDuration.Location = new Point(3, 230);
            lblGravityTransitionDuration.Name = "lblGravityTransitionDuration";
            lblGravityTransitionDuration.Size = new Size(223, 15);
            lblGravityTransitionDuration.TabIndex = 16;
            lblGravityTransitionDuration.Text = "Gravity Transition (1/α):";
            // 
            // numGravityTransitionDuration
            // 
            numGravityTransitionDuration.DecimalPlaces = 1;
            numGravityTransitionDuration.Dock = DockStyle.Fill;
            numGravityTransitionDuration.Increment = new decimal(new int[] { 10, 0, 0, 0 });
            numGravityTransitionDuration.Location = new Point(232, 227);
            numGravityTransitionDuration.Maximum = new decimal(new int[] { 500, 0, 0, 0 });
            numGravityTransitionDuration.Name = "numGravityTransitionDuration";
            numGravityTransitionDuration.Size = new Size(78, 23);
            numGravityTransitionDuration.TabIndex = 17;
            numGravityTransitionDuration.Value = new decimal(new int[] { 137, 0, 0, 0 });
            // 
            // valAnnealingTimeConstant
            // 
            valAnnealingTimeConstant.Anchor = AnchorStyles.Top | AnchorStyles.Left | AnchorStyles.Right;
            valAnnealingTimeConstant.AutoSize = true;
            valAnnealingTimeConstant.Location = new Point(3, 252);
            valAnnealingTimeConstant.Name = "valAnnealingTimeConstant";
            valAnnealingTimeConstant.Size = new Size(223, 15);
            valAnnealingTimeConstant.TabIndex = 18;
            valAnnealingTimeConstant.Text = "τ_anneal = (computed)";
            // 
            // tabPage_Experiments
            // 
            tabPage_Experiments.Location = new Point(4, 24);
            tabPage_Experiments.Name = "tabPage_Experiments";
            tabPage_Experiments.Size = new Size(1342, 615);
            tabPage_Experiments.TabIndex = 15;
            tabPage_Experiments.Text = "Experiments";
            tabPage_Experiments.UseVisualStyleBackColor = true;
            // 
            // statusStrip
            // 
            statusStrip.Items.AddRange(new ToolStripItem[] { statusLabelSteps, statusLabelExcited, statusLabelHeavyMass });
            statusStrip.Location = new Point(0, 715);
            statusStrip.Name = "statusStrip";
            statusStrip.Size = new Size(1663, 22);
            statusStrip.TabIndex = 4;
            // 
            // statusLabelSteps
            // 
            statusLabelSteps.Name = "statusLabelSteps";
            statusLabelSteps.Size = new Size(53, 17);
            statusLabelSteps.Text = "Step: 0/0";
            // 
            // statusLabelExcited
            // 
            statusLabelExcited.Margin = new Padding(20, 3, 0, 2);
            statusLabelExcited.Name = "statusLabelExcited";
            statusLabelExcited.Size = new Size(110, 17);
            statusLabelExcited.Text = "Excited: 0 (avg 0.00)";
            // 
            // statusLabelHeavyMass
            // 
            statusLabelHeavyMass.Name = "statusLabelHeavyMass";
            statusLabelHeavyMass.Size = new Size(94, 17);
            statusLabelHeavyMass.Text = "HeavyMass: 0.00";
            // 
            // toolStrip
            // 
            toolStrip.Items.AddRange(new ToolStripItem[] { toolRunButton, toolStopButton, toolKmcRunButton, toolKmcStopButton, toolSep1, toolViewDropDown });
            toolStrip.Location = new Point(0, 0);
            toolStrip.Name = "toolStrip";
            toolStrip.Size = new Size(1663, 25);
            toolStrip.TabIndex = 5;
            // 
            // toolRunButton
            // 
            toolRunButton.Name = "toolRunButton";
            toolRunButton.Size = new Size(23, 22);
            // 
            // toolStopButton
            // 
            toolStopButton.DisplayStyle = ToolStripItemDisplayStyle.Text;
            toolStopButton.Name = "toolStopButton";
            toolStopButton.Size = new Size(35, 22);
            toolStopButton.Text = "Stop";
            // 
            // toolKmcRunButton
            // 
            toolKmcRunButton.Name = "toolKmcRunButton";
            toolKmcRunButton.Size = new Size(23, 22);
            // 
            // toolKmcStopButton
            // 
            toolKmcStopButton.Name = "toolKmcStopButton";
            toolKmcStopButton.Size = new Size(23, 22);
            // 
            // toolSep1
            // 
            toolSep1.Name = "toolSep1";
            toolSep1.Size = new Size(6, 25);
            // 
            // toolViewDropDown
            // 
            toolViewDropDown.DisplayStyle = ToolStripItemDisplayStyle.Text;
            toolViewDropDown.DropDownItems.AddRange(new ToolStripItem[] { toolViewGraph, toolViewConsole, toolViewCsv, toolViewAnalysis, toolViewSummary });
            toolViewDropDown.Name = "toolViewDropDown";
            toolViewDropDown.Size = new Size(45, 22);
            toolViewDropDown.Text = "View";
            // 
            // toolViewGraph
            // 
            toolViewGraph.Name = "toolViewGraph";
            toolViewGraph.Size = new Size(125, 22);
            toolViewGraph.Text = "Network";
            toolViewGraph.Click += ToolViewGraph_Click;
            // 
            // toolViewConsole
            // 
            toolViewConsole.Name = "toolViewConsole";
            toolViewConsole.Size = new Size(125, 22);
            toolViewConsole.Text = "Console";
            toolViewConsole.Click += ToolViewConsole_Click;
            // 
            // toolViewCsv
            // 
            toolViewCsv.Name = "toolViewCsv";
            toolViewCsv.Size = new Size(125, 22);
            toolViewCsv.Text = "CSV";
            // 
            // toolViewAnalysis
            // 
            toolViewAnalysis.Name = "toolViewAnalysis";
            toolViewAnalysis.Size = new Size(125, 22);
            toolViewAnalysis.Text = "Analysis";
            toolViewAnalysis.Click += ToolViewAnalysis_Click;
            // 
            // toolViewSummary
            // 
            toolViewSummary.Name = "toolViewSummary";
            toolViewSummary.Size = new Size(125, 22);
            toolViewSummary.Text = "Summary";
            toolViewSummary.Click += ToolViewSummary_Click;
            // 
            // splitMain
            // 
            splitMain.Location = new Point(1353, 184);
            splitMain.Name = "splitMain";
            splitMain.Orientation = Orientation.Horizontal;
            // 
            // splitMain.Panel1
            // 
            splitMain.Panel1.AutoScroll = true;
            splitMain.Panel1.BackgroundImageLayout = ImageLayout.None;
            splitMain.Panel1.Controls.Add(label_ParamPresets);
            splitMain.Panel1.Controls.Add(label_MaxFPS);
            splitMain.Panel1.Controls.Add(label_CPUThreads);
            splitMain.Panel1.Controls.Add(numericUpDown1);
            splitMain.Panel1.Controls.Add(lblExperiments);
            splitMain.Panel1.Controls.Add(comboBox_Experiments);
            splitMain.Panel1.Controls.Add(btnExpornShortJson);
            splitMain.Panel1.Controls.Add(comboBox_Presets);
            splitMain.Panel1.Controls.Add(chkAutoScrollConsole);
            splitMain.Panel1.Controls.Add(checkBox_AutoTuning);
            splitMain.Panel1.Controls.Add(chkShowHeavyOnly);
            splitMain.Panel1.Controls.Add(numericUpDown_MaxFPS);
            splitMain.Panel1.Controls.Add(btnSnapshotImage);
            splitMain.Panel1.Controls.Add(button_ForceRedrawGraphImage);
            splitMain.Panel1.Controls.Add(button_SaveSimReult);
            splitMain.Panel1.Controls.Add(comboBox_GPUIndex);
            splitMain.Panel1.Controls.Add(cmbWeightThreshold);
            splitMain.Panel1.Controls.Add(checkBox_EnableGPU);
            // 
            // splitMain.Panel2
            // 
            splitMain.Panel2.Controls.Add(grpLiveMetrics);
            splitMain.Size = new Size(295, 518);
            splitMain.SplitterDistance = 280;
            splitMain.TabIndex = 6;
            // 
            // label_ParamPresets
            // 
            label_ParamPresets.AutoSize = true;
            label_ParamPresets.Location = new Point(5, 164);
            label_ParamPresets.Name = "label_ParamPresets";
            label_ParamPresets.Size = new Size(47, 15);
            label_ParamPresets.TabIndex = 27;
            label_ParamPresets.Text = "Presets:";
            // 
            // label_MaxFPS
            // 
            label_MaxFPS.AutoSize = true;
            label_MaxFPS.Location = new Point(173, 260);
            label_MaxFPS.Name = "label_MaxFPS";
            label_MaxFPS.Size = new Size(54, 15);
            label_MaxFPS.TabIndex = 26;
            label_MaxFPS.Text = "Max FPS:";
            // 
            // label_CPUThreads
            // 
            label_CPUThreads.AutoSize = true;
            label_CPUThreads.Location = new Point(6, 260);
            label_CPUThreads.Name = "label_CPUThreads";
            label_CPUThreads.Size = new Size(75, 15);
            label_CPUThreads.TabIndex = 25;
            label_CPUThreads.Text = "CPU threads:";
            // 
            // numericUpDown1
            // 
            numericUpDown1.Location = new Point(84, 251);
            numericUpDown1.Maximum = new decimal(new int[] { 1024, 0, 0, 0 });
            numericUpDown1.Minimum = new decimal(new int[] { 8, 0, 0, 0 });
            numericUpDown1.Name = "numericUpDown1";
            numericUpDown1.Size = new Size(56, 23);
            numericUpDown1.TabIndex = 24;
            numericUpDown1.Value = new decimal(new int[] { 8, 0, 0, 0 });
            // 
            // lblExperiments
            // 
            lblExperiments.AutoSize = true;
            lblExperiments.Location = new Point(4, 207);
            lblExperiments.Name = "lblExperiments";
            lblExperiments.Size = new Size(69, 15);
            lblExperiments.TabIndex = 23;
            lblExperiments.Text = "Experiment:";
            // 
            // comboBox_Experiments
            // 
            comboBox_Experiments.DropDownStyle = ComboBoxStyle.DropDownList;
            comboBox_Experiments.Location = new Point(7, 224);
            comboBox_Experiments.Name = "comboBox_Experiments";
            comboBox_Experiments.Size = new Size(269, 23);
            comboBox_Experiments.TabIndex = 22;
            // 
            // btnExpornShortJson
            // 
            btnExpornShortJson.AutoSize = true;
            btnExpornShortJson.Location = new Point(6, 133);
            btnExpornShortJson.Name = "btnExpornShortJson";
            btnExpornShortJson.Size = new Size(132, 25);
            btnExpornShortJson.TabIndex = 3;
            btnExpornShortJson.Text = "Export Short Now";
            btnExpornShortJson.Click += btnExpornShortJson_Click;
            // 
            // comboBox_Presets
            // 
            comboBox_Presets.DropDownStyle = ComboBoxStyle.DropDownList;
            comboBox_Presets.Items.AddRange(new object[] { "Nodes 100", "Nodes 200", "Nodes 300", "Nodes 400", "Nodes 500", "Nodes 800", "Nodes 1000" });
            comboBox_Presets.Location = new Point(7, 180);
            comboBox_Presets.Name = "comboBox_Presets";
            comboBox_Presets.Size = new Size(269, 23);
            comboBox_Presets.TabIndex = 20;
            // 
            // chkAutoScrollConsole
            // 
            chkAutoScrollConsole.AutoSize = true;
            chkAutoScrollConsole.Checked = true;
            chkAutoScrollConsole.CheckState = CheckState.Checked;
            chkAutoScrollConsole.Location = new Point(8, 45);
            chkAutoScrollConsole.Name = "chkAutoScrollConsole";
            chkAutoScrollConsole.Size = new Size(129, 19);
            chkAutoScrollConsole.TabIndex = 2;
            chkAutoScrollConsole.Text = "Auto-scroll console";
            // 
            // checkBox_AutoTuning
            // 
            checkBox_AutoTuning.AutoSize = true;
            checkBox_AutoTuning.Location = new Point(8, 17);
            checkBox_AutoTuning.Name = "checkBox_AutoTuning";
            checkBox_AutoTuning.Size = new Size(134, 19);
            checkBox_AutoTuning.TabIndex = 19;
            checkBox_AutoTuning.Text = "Auto-tuning params";
            checkBox_AutoTuning.CheckedChanged += checkBox_AutoTuning_CheckedChanged;
            // 
            // chkShowHeavyOnly
            // 
            chkShowHeavyOnly.AutoSize = true;
            chkShowHeavyOnly.Location = new Point(156, 17);
            chkShowHeavyOnly.Name = "chkShowHeavyOnly";
            chkShowHeavyOnly.Size = new Size(128, 19);
            chkShowHeavyOnly.TabIndex = 1;
            chkShowHeavyOnly.Text = "Heavy clusters only";
            // 
            // numericUpDown_MaxFPS
            // 
            numericUpDown_MaxFPS.Location = new Point(233, 251);
            numericUpDown_MaxFPS.Maximum = new decimal(new int[] { 30, 0, 0, 0 });
            numericUpDown_MaxFPS.Minimum = new decimal(new int[] { 1, 0, 0, 0 });
            numericUpDown_MaxFPS.Name = "numericUpDown_MaxFPS";
            numericUpDown_MaxFPS.Size = new Size(43, 23);
            numericUpDown_MaxFPS.TabIndex = 18;
            numericUpDown_MaxFPS.Value = new decimal(new int[] { 10, 0, 0, 0 });
            numericUpDown_MaxFPS.ValueChanged += numericUpDown_MaxFPS_ValueChanged;
            // 
            // btnSnapshotImage
            // 
            btnSnapshotImage.AutoSize = true;
            btnSnapshotImage.Location = new Point(144, 102);
            btnSnapshotImage.Name = "btnSnapshotImage";
            btnSnapshotImage.Size = new Size(132, 25);
            btnSnapshotImage.TabIndex = 4;
            btnSnapshotImage.Text = "Snapshot";
            btnSnapshotImage.Click += BtnSnapshotImage_Click;
            // 
            // button_ForceRedrawGraphImage
            // 
            button_ForceRedrawGraphImage.AutoSize = true;
            button_ForceRedrawGraphImage.Location = new Point(6, 102);
            button_ForceRedrawGraphImage.Name = "button_ForceRedrawGraphImage";
            button_ForceRedrawGraphImage.Size = new Size(132, 25);
            button_ForceRedrawGraphImage.TabIndex = 17;
            button_ForceRedrawGraphImage.Text = "Force Redraw Graph";
            button_ForceRedrawGraphImage.Click += button_ForceRedrawGraphImage_Click;
            // 
            // button_SaveSimReult
            // 
            button_SaveSimReult.AutoSize = true;
            button_SaveSimReult.Location = new Point(144, 133);
            button_SaveSimReult.Name = "button_SaveSimReult";
            button_SaveSimReult.Size = new Size(132, 25);
            button_SaveSimReult.TabIndex = 8;
            button_SaveSimReult.Text = "Save Sim Results";
            button_SaveSimReult.Click += button_SaveSimReult_Click;
            // 
            // comboBox_GPUIndex
            // 
            comboBox_GPUIndex.DropDownStyle = ComboBoxStyle.DropDownList;
            comboBox_GPUIndex.Items.AddRange(new object[] { "All (w>0.0)", "w>0.15", "w>0.3", "w>0.5", "w>0.7", "w>0.9" });
            comboBox_GPUIndex.Location = new Point(92, 70);
            comboBox_GPUIndex.Name = "comboBox_GPUIndex";
            comboBox_GPUIndex.Size = new Size(184, 23);
            comboBox_GPUIndex.TabIndex = 15;
            // 
            // cmbWeightThreshold
            // 
            cmbWeightThreshold.DropDownStyle = ComboBoxStyle.DropDownList;
            cmbWeightThreshold.Items.AddRange(new object[] { "All (w>0.0)", "w>0.15", "w>0.3", "w>0.5", "w>0.7", "w>0.9" });
            cmbWeightThreshold.Location = new Point(156, 38);
            cmbWeightThreshold.Name = "cmbWeightThreshold";
            cmbWeightThreshold.Size = new Size(120, 23);
            cmbWeightThreshold.TabIndex = 0;
            // 
            // checkBox_EnableGPU
            // 
            checkBox_EnableGPU.AutoSize = true;
            checkBox_EnableGPU.Checked = true;
            checkBox_EnableGPU.CheckState = CheckState.Checked;
            checkBox_EnableGPU.Location = new Point(9, 79);
            checkBox_EnableGPU.Name = "checkBox_EnableGPU";
            checkBox_EnableGPU.Size = new Size(87, 19);
            checkBox_EnableGPU.TabIndex = 14;
            checkBox_EnableGPU.Text = "Enable GPU";
            checkBox_EnableGPU.UseVisualStyleBackColor = true;
            checkBox_EnableGPU.CheckedChanged += checkBox_EnableGPU_CheckedChanged;
            // 
            // grpLiveMetrics
            // 
            grpLiveMetrics.Controls.Add(tlpLiveMetrics);
            grpLiveMetrics.Location = new Point(6, 12);
            grpLiveMetrics.Name = "grpLiveMetrics";
            grpLiveMetrics.Size = new Size(295, 128);
            grpLiveMetrics.TabIndex = 0;
            grpLiveMetrics.TabStop = false;
            grpLiveMetrics.Text = "Live Metrics";
            // 
            // tlpLiveMetrics
            // 
            tlpLiveMetrics.ColumnCount = 2;
            tlpLiveMetrics.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 60F));
            tlpLiveMetrics.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 40F));
            tlpLiveMetrics.Controls.Add(lblGlobalNbr, 0, 0);
            tlpLiveMetrics.Controls.Add(valGlobalNbr, 1, 0);
            tlpLiveMetrics.Controls.Add(lblGlobalSpont, 0, 1);
            tlpLiveMetrics.Controls.Add(valGlobalSpont, 1, 1);
            tlpLiveMetrics.Controls.Add(lblStrongEdges, 0, 2);
            tlpLiveMetrics.Controls.Add(valStrongEdges, 1, 2);
            tlpLiveMetrics.Controls.Add(lblLargestCluster, 0, 3);
            tlpLiveMetrics.Controls.Add(valLargestCluster, 1, 3);
            tlpLiveMetrics.Controls.Add(lblHeavyMass, 0, 4);
            tlpLiveMetrics.Controls.Add(valHeavyMass, 1, 4);
            tlpLiveMetrics.Controls.Add(lblSpectrumInfo, 0, 5);
            tlpLiveMetrics.Controls.Add(valSpectrumInfo, 1, 5);
            tlpLiveMetrics.Controls.Add(lblLightSpeed, 0, 6);
            tlpLiveMetrics.Controls.Add(valLightSpeed, 1, 6);
            tlpLiveMetrics.Location = new Point(3, 19);
            tlpLiveMetrics.Name = "tlpLiveMetrics";
            tlpLiveMetrics.RowCount = 7;
            tlpLiveMetrics.RowStyles.Add(new RowStyle());
            tlpLiveMetrics.RowStyles.Add(new RowStyle());
            tlpLiveMetrics.RowStyles.Add(new RowStyle());
            tlpLiveMetrics.RowStyles.Add(new RowStyle());
            tlpLiveMetrics.RowStyles.Add(new RowStyle());
            tlpLiveMetrics.RowStyles.Add(new RowStyle());
            tlpLiveMetrics.RowStyles.Add(new RowStyle(SizeType.Percent, 100F));
            tlpLiveMetrics.Size = new Size(289, 112);
            tlpLiveMetrics.TabIndex = 0;
            // 
            // lblGlobalNbr
            // 
            lblGlobalNbr.AutoSize = true;
            lblGlobalNbr.Location = new Point(3, 0);
            lblGlobalNbr.Name = "lblGlobalNbr";
            lblGlobalNbr.Size = new Size(67, 15);
            lblGlobalNbr.TabIndex = 0;
            lblGlobalNbr.Text = "Global Nbr:";
            // 
            // valGlobalNbr
            // 
            valGlobalNbr.AutoSize = true;
            valGlobalNbr.Location = new Point(176, 0);
            valGlobalNbr.Name = "valGlobalNbr";
            valGlobalNbr.Size = new Size(34, 15);
            valGlobalNbr.TabIndex = 1;
            valGlobalNbr.Text = "0.000";
            // 
            // lblGlobalSpont
            // 
            lblGlobalSpont.AutoSize = true;
            lblGlobalSpont.Location = new Point(3, 15);
            lblGlobalSpont.Name = "lblGlobalSpont";
            lblGlobalSpont.Size = new Size(78, 15);
            lblGlobalSpont.TabIndex = 2;
            lblGlobalSpont.Text = "Global Spont:";
            // 
            // valGlobalSpont
            // 
            valGlobalSpont.AutoSize = true;
            valGlobalSpont.Location = new Point(176, 15);
            valGlobalSpont.Name = "valGlobalSpont";
            valGlobalSpont.Size = new Size(34, 15);
            valGlobalSpont.TabIndex = 3;
            valGlobalSpont.Text = "0.000";
            // 
            // lblStrongEdges
            // 
            lblStrongEdges.AutoSize = true;
            lblStrongEdges.Location = new Point(3, 30);
            lblStrongEdges.Name = "lblStrongEdges";
            lblStrongEdges.Size = new Size(79, 15);
            lblStrongEdges.TabIndex = 4;
            lblStrongEdges.Text = "Strong edges:";
            // 
            // valStrongEdges
            // 
            valStrongEdges.AutoSize = true;
            valStrongEdges.Location = new Point(176, 30);
            valStrongEdges.Name = "valStrongEdges";
            valStrongEdges.Size = new Size(13, 15);
            valStrongEdges.TabIndex = 5;
            valStrongEdges.Text = "0";
            // 
            // lblLargestCluster
            // 
            lblLargestCluster.AutoSize = true;
            lblLargestCluster.Location = new Point(3, 45);
            lblLargestCluster.Name = "lblLargestCluster";
            lblLargestCluster.Size = new Size(86, 15);
            lblLargestCluster.TabIndex = 6;
            lblLargestCluster.Text = "Largest cluster:";
            // 
            // valLargestCluster
            // 
            valLargestCluster.AutoSize = true;
            valLargestCluster.Location = new Point(176, 45);
            valLargestCluster.Name = "valLargestCluster";
            valLargestCluster.Size = new Size(13, 15);
            valLargestCluster.TabIndex = 7;
            valLargestCluster.Text = "0";
            // 
            // lblHeavyMass
            // 
            lblHeavyMass.AutoSize = true;
            lblHeavyMass.Location = new Point(3, 60);
            lblHeavyMass.Name = "lblHeavyMass";
            lblHeavyMass.Size = new Size(73, 15);
            lblHeavyMass.TabIndex = 8;
            lblHeavyMass.Text = "Heavy mass:";
            // 
            // valHeavyMass
            // 
            valHeavyMass.AutoSize = true;
            valHeavyMass.Location = new Point(176, 60);
            valHeavyMass.Name = "valHeavyMass";
            valHeavyMass.Size = new Size(22, 15);
            valHeavyMass.TabIndex = 9;
            valHeavyMass.Text = "0.0";
            // 
            // lblSpectrumInfo
            // 
            lblSpectrumInfo.AutoSize = true;
            lblSpectrumInfo.Location = new Point(3, 75);
            lblSpectrumInfo.Name = "lblSpectrumInfo";
            lblSpectrumInfo.Size = new Size(86, 15);
            lblSpectrumInfo.TabIndex = 10;
            lblSpectrumInfo.Text = "Spectrum logs:";
            // 
            // valSpectrumInfo
            // 
            valSpectrumInfo.AutoSize = true;
            valSpectrumInfo.Location = new Point(176, 75);
            valSpectrumInfo.Name = "valSpectrumInfo";
            valSpectrumInfo.Size = new Size(22, 15);
            valSpectrumInfo.TabIndex = 11;
            valSpectrumInfo.Text = "off";
            // 
            // lblLightSpeed
            // 
            lblLightSpeed.AutoSize = true;
            lblLightSpeed.Location = new Point(3, 90);
            lblLightSpeed.Name = "lblLightSpeed";
            lblLightSpeed.Size = new Size(35, 15);
            lblLightSpeed.TabIndex = 12;
            lblLightSpeed.Text = "c_eff:";
            // 
            // valLightSpeed
            // 
            valLightSpeed.AutoSize = true;
            valLightSpeed.Location = new Point(176, 90);
            valLightSpeed.Name = "valLightSpeed";
            valLightSpeed.Size = new Size(13, 15);
            valLightSpeed.TabIndex = 13;
            valLightSpeed.Text = "0";
            // 
            // grpRunStats
            // 
            grpRunStats.Controls.Add(tlpRunStats);
            grpRunStats.Location = new Point(1352, 54);
            grpRunStats.Name = "grpRunStats";
            grpRunStats.Size = new Size(298, 121);
            grpRunStats.TabIndex = 0;
            grpRunStats.TabStop = false;
            grpRunStats.Text = "Run Summary";
            // 
            // tlpRunStats
            // 
            tlpRunStats.ColumnCount = 2;
            tlpRunStats.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 60F));
            tlpRunStats.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 40F));
            tlpRunStats.Controls.Add(lblExcitedAvg, 0, 2);
            tlpRunStats.Controls.Add(valExcitedAvg, 1, 2);
            tlpRunStats.Controls.Add(lblExcitedMax, 0, 3);
            tlpRunStats.Controls.Add(valExcitedMax, 1, 3);
            tlpRunStats.Controls.Add(lblAvalancheCount, 0, 4);
            tlpRunStats.Controls.Add(valAvalancheCount, 1, 4);
            tlpRunStats.Controls.Add(lblMeasurementStatus, 0, 5);
            tlpRunStats.Controls.Add(valMeasurementStatus, 1, 5);
            tlpRunStats.Controls.Add(valCurrentStep, 1, 0);
            tlpRunStats.Controls.Add(valTotalSteps, 1, 1);
            tlpRunStats.Controls.Add(lblTotalSteps, 0, 1);
            tlpRunStats.Controls.Add(lblCurrentStep, 0, 0);
            tlpRunStats.Dock = DockStyle.Fill;
            tlpRunStats.Location = new Point(3, 19);
            tlpRunStats.Name = "tlpRunStats";
            tlpRunStats.RowCount = 6;
            tlpRunStats.RowStyles.Add(new RowStyle());
            tlpRunStats.RowStyles.Add(new RowStyle());
            tlpRunStats.RowStyles.Add(new RowStyle());
            tlpRunStats.RowStyles.Add(new RowStyle());
            tlpRunStats.RowStyles.Add(new RowStyle());
            tlpRunStats.RowStyles.Add(new RowStyle());
            tlpRunStats.Size = new Size(292, 99);
            tlpRunStats.TabIndex = 0;
            // 
            // lblExcitedAvg
            // 
            lblExcitedAvg.AutoSize = true;
            lblExcitedAvg.Location = new Point(3, 30);
            lblExcitedAvg.Name = "lblExcitedAvg";
            lblExcitedAvg.Size = new Size(71, 15);
            lblExcitedAvg.TabIndex = 4;
            lblExcitedAvg.Text = "Avg excited:";
            // 
            // valExcitedAvg
            // 
            valExcitedAvg.AutoSize = true;
            valExcitedAvg.Location = new Point(178, 30);
            valExcitedAvg.Name = "valExcitedAvg";
            valExcitedAvg.Size = new Size(28, 15);
            valExcitedAvg.TabIndex = 5;
            valExcitedAvg.Text = "0.00";
            // 
            // lblExcitedMax
            // 
            lblExcitedMax.AutoSize = true;
            lblExcitedMax.Location = new Point(3, 45);
            lblExcitedMax.Name = "lblExcitedMax";
            lblExcitedMax.Size = new Size(72, 15);
            lblExcitedMax.TabIndex = 6;
            lblExcitedMax.Text = "Max excited:";
            // 
            // valExcitedMax
            // 
            valExcitedMax.AutoSize = true;
            valExcitedMax.Location = new Point(178, 45);
            valExcitedMax.Name = "valExcitedMax";
            valExcitedMax.Size = new Size(13, 15);
            valExcitedMax.TabIndex = 7;
            valExcitedMax.Text = "0";
            // 
            // lblAvalancheCount
            // 
            lblAvalancheCount.AutoSize = true;
            lblAvalancheCount.Location = new Point(3, 60);
            lblAvalancheCount.Name = "lblAvalancheCount";
            lblAvalancheCount.Size = new Size(70, 15);
            lblAvalancheCount.TabIndex = 8;
            lblAvalancheCount.Text = "Avalanches:";
            // 
            // valAvalancheCount
            // 
            valAvalancheCount.AutoSize = true;
            valAvalancheCount.Location = new Point(178, 60);
            valAvalancheCount.Name = "valAvalancheCount";
            valAvalancheCount.Size = new Size(13, 15);
            valAvalancheCount.TabIndex = 9;
            valAvalancheCount.Text = "0";
            // 
            // lblMeasurementStatus
            // 
            lblMeasurementStatus.AutoSize = true;
            lblMeasurementStatus.Location = new Point(3, 75);
            lblMeasurementStatus.Name = "lblMeasurementStatus";
            lblMeasurementStatus.Size = new Size(83, 15);
            lblMeasurementStatus.TabIndex = 10;
            lblMeasurementStatus.Text = "Measurement:";
            // 
            // valMeasurementStatus
            // 
            valMeasurementStatus.AutoSize = true;
            valMeasurementStatus.Location = new Point(178, 75);
            valMeasurementStatus.Name = "valMeasurementStatus";
            valMeasurementStatus.Size = new Size(29, 15);
            valMeasurementStatus.TabIndex = 11;
            valMeasurementStatus.Text = "N/A";
            // 
            // valCurrentStep
            // 
            valCurrentStep.AutoSize = true;
            valCurrentStep.Location = new Point(178, 0);
            valCurrentStep.Name = "valCurrentStep";
            valCurrentStep.Size = new Size(13, 15);
            valCurrentStep.TabIndex = 3;
            valCurrentStep.Text = "0";
            // 
            // valTotalSteps
            // 
            valTotalSteps.AutoSize = true;
            valTotalSteps.Location = new Point(178, 15);
            valTotalSteps.Name = "valTotalSteps";
            valTotalSteps.Size = new Size(13, 15);
            valTotalSteps.TabIndex = 1;
            valTotalSteps.Text = "0";
            // 
            // lblTotalSteps
            // 
            lblTotalSteps.AutoSize = true;
            lblTotalSteps.Location = new Point(3, 15);
            lblTotalSteps.Name = "lblTotalSteps";
            lblTotalSteps.Size = new Size(66, 15);
            lblTotalSteps.TabIndex = 0;
            lblTotalSteps.Text = "Total steps:";
            // 
            // lblCurrentStep
            // 
            lblCurrentStep.AutoSize = true;
            lblCurrentStep.Location = new Point(3, 0);
            lblCurrentStep.Name = "lblCurrentStep";
            lblCurrentStep.Size = new Size(75, 15);
            lblCurrentStep.TabIndex = 2;
            lblCurrentStep.Text = "Current step:";
            // 
            // modernSimTextBox
            // 
            modernSimTextBox.Location = new Point(0, 0);
            modernSimTextBox.Name = "modernSimTextBox";
            modernSimTextBox.Size = new Size(100, 23);
            modernSimTextBox.TabIndex = 0;
            // 
            // button_RunModernSim
            // 
            button_RunModernSim.AutoSize = true;
            button_RunModernSim.Location = new Point(1, 30);
            button_RunModernSim.Name = "button_RunModernSim";
            button_RunModernSim.Size = new Size(132, 25);
            button_RunModernSim.TabIndex = 9;
            button_RunModernSim.Text = "Run simulation";
            button_RunModernSim.Click += button_RunSimulation_Click;
            // 
            // Form_Main
            // 
            AutoScaleDimensions = new SizeF(7F, 15F);
            AutoScaleMode = AutoScaleMode.Font;
            ClientSize = new Size(1663, 737);
            Controls.Add(grpRunStats);
            Controls.Add(button_RunModernSim);
            Controls.Add(splitMain);
            Controls.Add(toolStrip);
            Controls.Add(statusStrip);
            Controls.Add(tabControl1);
            Name = "Form_Main";
            Text = "RQ-Sim";
            Load += Form_Main_Load;
            tabControl1.ResumeLayout(false);
            tabPage_OnExcited.ResumeLayout(false);
            tabPage_Heavy.ResumeLayout(false);
            tabPage_Cluster.ResumeLayout(false);
            tabPage_Energy.ResumeLayout(false);
            tabPage_Console.ResumeLayout(false);
            tabPage_Console.PerformLayout();
            tabPage_Summary.ResumeLayout(false);
            tabPage_Summary.PerformLayout();
            grpDashboard.ResumeLayout(false);
            tlpDashboard.ResumeLayout(false);
            tlpDashboard.PerformLayout();
            tabPage_Settings.ResumeLayout(false);
            settingsMainLayout.ResumeLayout(false);
            grpSimParams.ResumeLayout(false);
            grpSimParams.PerformLayout();
            tlpSimParams.ResumeLayout(false);
            tlpSimParams.PerformLayout();
            ((System.ComponentModel.ISupportInitialize)numNodeCount).EndInit();
            ((System.ComponentModel.ISupportInitialize)numTargetDegree).EndInit();
            ((System.ComponentModel.ISupportInitialize)numInitialExcitedProb).EndInit();
            ((System.ComponentModel.ISupportInitialize)numLambdaState).EndInit();
            ((System.ComponentModel.ISupportInitialize)numTemperature).EndInit();
            ((System.ComponentModel.ISupportInitialize)numEdgeTrialProb).EndInit();
            ((System.ComponentModel.ISupportInitialize)numMeasurementThreshold).EndInit();
            ((System.ComponentModel.ISupportInitialize)numTotalSteps).EndInit();
            ((System.ComponentModel.ISupportInitialize)numFractalLevels).EndInit();
            ((System.ComponentModel.ISupportInitialize)numFractalBranchFactor).EndInit();
            grpPhysicsModules.ResumeLayout(false);
            grpPhysicsModules.PerformLayout();
            flpPhysics.ResumeLayout(false);
            flpPhysics.PerformLayout();
            grpPhysicsConstants.ResumeLayout(false);
            grpPhysicsConstants.PerformLayout();
            tlpPhysicsConstants.ResumeLayout(false);
            tlpPhysicsConstants.PerformLayout();
            ((System.ComponentModel.ISupportInitialize)numInitialEdgeProb).EndInit();
            ((System.ComponentModel.ISupportInitialize)numGravitationalCoupling).EndInit();
            ((System.ComponentModel.ISupportInitialize)numVacuumEnergyScale).EndInit();
            ((System.ComponentModel.ISupportInitialize)numAnnealingCoolingRate).EndInit();
            ((System.ComponentModel.ISupportInitialize)numDecoherenceRate).EndInit();
            ((System.ComponentModel.ISupportInitialize)numHotStartTemperature).EndInit();
            ((System.ComponentModel.ISupportInitialize)numAdaptiveThresholdSigma).EndInit();
            ((System.ComponentModel.ISupportInitialize)numWarmupDuration).EndInit();
            ((System.ComponentModel.ISupportInitialize)numGravityTransitionDuration).EndInit();
            statusStrip.ResumeLayout(false);
            statusStrip.PerformLayout();
            toolStrip.ResumeLayout(false);
            toolStrip.PerformLayout();
            splitMain.Panel1.ResumeLayout(false);
            splitMain.Panel1.PerformLayout();
            splitMain.Panel2.ResumeLayout(false);
            ((System.ComponentModel.ISupportInitialize)splitMain).EndInit();
            splitMain.ResumeLayout(false);
            ((System.ComponentModel.ISupportInitialize)numericUpDown1).EndInit();
            ((System.ComponentModel.ISupportInitialize)numericUpDown_MaxFPS).EndInit();
            grpLiveMetrics.ResumeLayout(false);
            tlpLiveMetrics.ResumeLayout(false);
            tlpLiveMetrics.PerformLayout();
            grpRunStats.ResumeLayout(false);
            tlpRunStats.ResumeLayout(false);
            tlpRunStats.PerformLayout();
            ResumeLayout(false);
            PerformLayout();
        }

        #endregion

        private CancellationTokenSource? _kmcCts;


        private NumericUpDown numNodeCount;
        private NumericUpDown numTargetDegree;
        private NumericUpDown numInitialExcitedProb;
        private NumericUpDown numLambdaState;
        private NumericUpDown numTemperature;
        private NumericUpDown numEdgeTrialProb;
        private NumericUpDown numMeasurementThreshold;
        private NumericUpDown numTotalSteps;
        private NumericUpDown numFractalLevels;
        private NumericUpDown numFractalBranchFactor;
        private CheckBox chkQuantumDriven;
        private CheckBox chkSpacetimePhysics;
        private CheckBox chkSpinorField;
        private CheckBox chkVacuumFluctuations;
        private CheckBox chkBlackHolePhysics;
        private CheckBox chkYangMillsGauge;
        private CheckBox chkEnhancedKleinGordon;
        private CheckBox chkInternalTime;
        private CheckBox chkSpectralGeometry;
        private CheckBox chkQuantumGraphity;
        private TableLayoutPanel settingsMainLayout;
        private GroupBox grpSimParams;
        private TableLayoutPanel tlpSimParams;
        private Label lblNodeCount;
        private Label lblTargetDegree;
        private Label lblInitialExcitedProb;
        private Label lblLambdaState;
        private Label lblTemperature;
        private Label lblEdgeTrialProb;
        private Label lblMeasurementThreshold;
        private Label lblTotalStepsSettings;
        private Label lblFractalLevels;
        private Label lblFractalBranchFactor;
        private GroupBox grpPhysicsModules;
        private FlowLayoutPanel flpPhysics;
        private CheckBox checkBox_EnableGPU;
        private ComboBox comboBox_GPUIndex;
        // Physics Constants controls
        private GroupBox grpPhysicsConstants;
        private TableLayoutPanel tlpPhysicsConstants;
        private Label lblInitialEdgeProb;
        private NumericUpDown numInitialEdgeProb;
        private Label lblGravitationalCoupling;
        private NumericUpDown numGravitationalCoupling;
        private Label lblVacuumEnergyScale;
        private NumericUpDown numVacuumEnergyScale;
        private Label lblAnnealingCoolingRate;
        private NumericUpDown numAnnealingCoolingRate;
        private Label lblDecoherenceRate;
        private NumericUpDown numDecoherenceRate;
        private Label lblHotStartTemperature;
        private NumericUpDown numHotStartTemperature;
        // Additional physics constants
        private Label lblAdaptiveThresholdSigma;
        private NumericUpDown numAdaptiveThresholdSigma;
        private Label lblWarmupDuration;
        private NumericUpDown numWarmupDuration;
        private Label lblGravityTransitionDuration;
        private NumericUpDown numGravityTransitionDuration;
        private Label valAnnealingTimeConstant;
        // Extended physics checkboxes
        private CheckBox chkRelationalTime;
        private CheckBox chkRelationalYangMills;
        private CheckBox chkNetworkGravity;
        private CheckBox chkUnifiedPhysicsStep;
        private CheckBox chkEnforceGaugeConstraints;
        private CheckBox chkCausalRewiring;
        private CheckBox chkTopologicalProtection;
        private CheckBox chkValidateEnergyConservation;
        private CheckBox chkMexicanHatPotential;
        private CheckBox chkHotStartAnnealing;
        private CheckBox chkGeometryMomenta;
        private CheckBox chkTopologicalCensorship;
        private GroupBox grpDashboard;
        private TableLayoutPanel tlpDashboard;
        private Label lblDashNodes;
        private Label valDashNodes;
        private Label lblDashTotalSteps;
        private Label valDashTotalSteps;
        private Label lblDashCurrentStep;
        private Label valDashCurrentStep;
        private Label lblDashExcited;
        private Label valDashExcited;
        private Label lblDashHeavyMass;
        private Label valDashHeavyMass;
        private Label lblDashLargestCluster;
        private Label valDashLargestCluster;
        private Label lblDashStrongEdges;
        private Label valDashStrongEdges;
        private Label lblDashPhase;
        private Label valDashPhase;
        private Label lblDashQNorm;
        private Label valDashQNorm;
        private Label lblDashEntanglement;
        private Label valDashEntanglement;
        private Label lblDashCorrelation;
        private Label valDashCorrelation;
        private Label lblDashStatus;
        private Label valDashStatus;
        private Label lblDashSpectralDim;
        private Label valDashSpectralDim;
        private Label lblDashEffectiveG;
        private Label valDashEffectiveG;
        private Label lblDashGSuppression;
        private Label valDashGSuppression;
        private Label lblDashNetworkTemp;
        private Label valDashNetworkTemp;
        private Button button_ForceRedrawGraphImage;
        private NumericUpDown numericUpDown_MaxFPS;
        private CheckBox checkBox_AutoTuning;
        private ComboBox comboBox_Presets;
        private ComboBox comboBox_Experiments;
        private Label lblExperiments;
        private TabPage tabPage_Experiments;
        private NumericUpDown numericUpDown1;
        private Label label_MaxFPS;
        private Label label_CPUThreads;
        private Label label_ParamPresets;
    }
}