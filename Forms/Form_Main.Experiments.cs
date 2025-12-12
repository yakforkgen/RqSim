// Forms/Form_Main.Experiments.cs
// Partial class for Experiments tab UI and logic

using System;
using System.Collections.Generic;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.Json;
using System.Windows.Forms;
using RQSimulation.Experiments;

namespace RqSimForms;

public partial class Form_Main
{
    // === Experiment Tab Controls ===
    private TableLayoutPanel _expMainLayout = null!;
    
    // GroupBox 1: Experiment Selection & Info
    private GroupBox _grpExpSelection = null!;
    private ComboBox _cmbExpSelect = null!;
    private TextBox _txtExpDescription = null!;
    private TextBox _txtExpNotes = null!;
    private Button _btnExpNew = null!;
    private Button _btnExpLoad = null!;
    private Button _btnExpSave = null!;
    
    // GroupBox 2: Expected Results
    private GroupBox _grpExpExpected = null!;
    private TableLayoutPanel _tlpExpExpected = null!;
    private NumericUpDown _numExpSpectralDimTarget = null!;
    private NumericUpDown _numExpSpectralDimMin = null!;
    private NumericUpDown _numExpSpectralDimMax = null!;
    private NumericUpDown _numExpHeavyClusterMin = null!;
    private NumericUpDown _numExpHeavyClusterMax = null!;
    private NumericUpDown _numExpHeavyMassMin = null!;
    private NumericUpDown _numExpHeavyMassMax = null!;
    private NumericUpDown _numExpLargestClusterMax = null!;
    private NumericUpDown _numExpFinalTempMax = null!;
    private CheckBox _chkExpPhaseTransition = null!;
    private CheckBox _chkExpStabilization = null!;
    
    // GroupBox 3: Validation Summary
    private GroupBox _grpExpValidation = null!;
    private Label _lblExpValidationStatus = null!;
    private ListView _lvExpCriteria = null!;
    private TextBox _txtExpReport = null!;
    private Button _btnExpValidate = null!;
    private Button _btnExpApplyToSim = null!;
    
    // Current experiment state
    private ExperimentDefinition? _currentExperiment;
    
    /// <summary>
    /// Initializes the Experiments tab with all controls.
    /// Call this in Form_Main constructor after InitializeComponent.
    /// </summary>
    private void InitializeExperimentsTab()
    {
        // Main layout: 3 columns
        _expMainLayout = new TableLayoutPanel
        {
            Dock = DockStyle.Fill,
            ColumnCount = 3,
            RowCount = 1,
            Padding = new Padding(5)
        };
        _expMainLayout.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 30F));
        _expMainLayout.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 35F));
        _expMainLayout.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 35F));
        _expMainLayout.RowStyles.Add(new RowStyle(SizeType.Percent, 100F));
        
        // === GroupBox 1: Experiment Selection ===
        InitializeExpSelectionGroup();
        
        // === GroupBox 2: Expected Results ===
        InitializeExpExpectedGroup();
        
        // === GroupBox 3: Validation Summary ===
        InitializeExpValidationGroup();
        
        // Add to main layout
        _expMainLayout.Controls.Add(_grpExpSelection, 0, 0);
        _expMainLayout.Controls.Add(_grpExpExpected, 1, 0);
        _expMainLayout.Controls.Add(_grpExpValidation, 2, 0);
        
        tabPage_Experiments.Controls.Add(_expMainLayout);
        
        // Load built-in experiments
        LoadBuiltInExperiments();
    }
    
    private void InitializeExpSelectionGroup()
    {
        _grpExpSelection = new GroupBox
        {
            Text = "Experiment Selection",
            Dock = DockStyle.Fill,
            Padding = new Padding(5)
        };
        
        var layout = new TableLayoutPanel
        {
            Dock = DockStyle.Fill,
            ColumnCount = 2,
            RowCount = 7
        };
        layout.ColumnStyles.Add(new ColumnStyle(SizeType.AutoSize));
        layout.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 100F));
        layout.RowStyles.Add(new RowStyle(SizeType.AutoSize)); // ComboBox
        layout.RowStyles.Add(new RowStyle(SizeType.AutoSize)); // Buttons
        layout.RowStyles.Add(new RowStyle(SizeType.AutoSize)); // Label
        layout.RowStyles.Add(new RowStyle(SizeType.Percent, 50F)); // Description
        layout.RowStyles.Add(new RowStyle(SizeType.AutoSize)); // Notes label
        layout.RowStyles.Add(new RowStyle(SizeType.Percent, 50F)); // Notes
        layout.RowStyles.Add(new RowStyle(SizeType.AutoSize)); // Apply button
        
        // Row 0: Experiment selector
        var lblSelect = new Label { Text = "Experiment:", AutoSize = true, Anchor = AnchorStyles.Left | AnchorStyles.Right };
        _cmbExpSelect = new ComboBox { Dock = DockStyle.Fill, DropDownStyle = ComboBoxStyle.DropDownList };
        _cmbExpSelect.SelectedIndexChanged += CmbExpSelect_SelectedIndexChanged;
        layout.Controls.Add(lblSelect, 0, 0);
        layout.Controls.Add(_cmbExpSelect, 1, 0);
        
        // Row 1: Buttons
        var btnPanel = new FlowLayoutPanel { Dock = DockStyle.Fill, FlowDirection = FlowDirection.LeftToRight };
        _btnExpNew = new Button { Text = "New", Width = 60 };
        _btnExpNew.Click += BtnExpNew_Click;
        _btnExpLoad = new Button { Text = "Load JSON", Width = 80 };
        _btnExpLoad.Click += BtnExpLoad_Click;
        _btnExpSave = new Button { Text = "Save JSON", Width = 80 };
        _btnExpSave.Click += BtnExpSave_Click;
        btnPanel.Controls.AddRange(new Control[] { _btnExpNew, _btnExpLoad, _btnExpSave });
        layout.Controls.Add(btnPanel, 0, 1);
        layout.SetColumnSpan(btnPanel, 2);
        
        // Row 2-3: Description
        var lblDesc = new Label { Text = "Description:", AutoSize = true };
        layout.Controls.Add(lblDesc, 0, 2);
        layout.SetColumnSpan(lblDesc, 2);
        
        _txtExpDescription = new TextBox
        {
            Dock = DockStyle.Fill,
            Multiline = true,
            ReadOnly = true,
            ScrollBars = ScrollBars.Vertical,
            BackColor = SystemColors.Info
        };
        layout.Controls.Add(_txtExpDescription, 0, 3);
        layout.SetColumnSpan(_txtExpDescription, 2);
        
        // Row 4-5: Notes
        var lblNotes = new Label { Text = "Notes (editable):", AutoSize = true };
        layout.Controls.Add(lblNotes, 0, 4);
        layout.SetColumnSpan(lblNotes, 2);
        
        _txtExpNotes = new TextBox
        {
            Dock = DockStyle.Fill,
            Multiline = true,
            ScrollBars = ScrollBars.Vertical
        };
        layout.Controls.Add(_txtExpNotes, 0, 5);
        layout.SetColumnSpan(_txtExpNotes, 2);
        
        // Row 6: Apply button
        _btnExpApplyToSim = new Button 
        { 
            Text = "Apply to Simulation ?",
            Dock = DockStyle.Fill,
            Height = 30,
            Font = new Font(Font, FontStyle.Bold)
        };
        _btnExpApplyToSim.Click += BtnExpApplyToSim_Click;
        layout.Controls.Add(_btnExpApplyToSim, 0, 6);
        layout.SetColumnSpan(_btnExpApplyToSim, 2);
        
        _grpExpSelection.Controls.Add(layout);
    }
    
    private void InitializeExpExpectedGroup()
    {
        _grpExpExpected = new GroupBox
        {
            Text = "Expected Results (Predictions)",
            Dock = DockStyle.Fill,
            Padding = new Padding(5)
        };
        
        _tlpExpExpected = new TableLayoutPanel
        {
            Dock = DockStyle.Fill,
            ColumnCount = 2,
            RowCount = 12,
            AutoScroll = true
        };
        _tlpExpExpected.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 60F));
        _tlpExpExpected.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 40F));
        
        for (int i = 0; i < 12; i++)
            _tlpExpExpected.RowStyles.Add(new RowStyle(SizeType.AutoSize));
        
        int row = 0;
        
        // Spectral Dimension Target
        AddExpectedRow("d_S Target:", ref _numExpSpectralDimTarget, 4.0, 0.0, 10.0, 1, row++);
        AddExpectedRow("d_S Min:", ref _numExpSpectralDimMin, 2.0, 0.0, 10.0, 1, row++);
        AddExpectedRow("d_S Max:", ref _numExpSpectralDimMax, 6.0, 0.0, 10.0, 1, row++);
        
        // Heavy Clusters
        AddExpectedRowInt("Heavy Cluster Min:", ref _numExpHeavyClusterMin, 1, 0, 1000, row++);
        AddExpectedRowInt("Heavy Cluster Max:", ref _numExpHeavyClusterMax, 50, 0, 1000, row++);
        
        // Heavy Mass
        AddExpectedRow("Heavy Mass Min:", ref _numExpHeavyMassMin, 0.0, 0.0, 1000.0, 1, row++);
        AddExpectedRow("Heavy Mass Max:", ref _numExpHeavyMassMax, 100.0, 0.0, 1000.0, 1, row++);
        
        // Largest Cluster (as fraction)
        AddExpectedRow("Largest Cluster Max (%):", ref _numExpLargestClusterMax, 30.0, 0.0, 100.0, 0, row++);
        
        // Final Temperature
        AddExpectedRow("Final Temp Max:", ref _numExpFinalTempMax, 1.0, 0.0, 100.0, 2, row++);
        
        // Checkboxes
        _chkExpPhaseTransition = new CheckBox { Text = "Expect Phase Transition", Checked = true, AutoSize = true };
        _tlpExpExpected.Controls.Add(_chkExpPhaseTransition, 0, row);
        _tlpExpExpected.SetColumnSpan(_chkExpPhaseTransition, 2);
        row++;
        
        _chkExpStabilization = new CheckBox { Text = "Expect Stabilization", Checked = true, AutoSize = true };
        _tlpExpExpected.Controls.Add(_chkExpStabilization, 0, row);
        _tlpExpExpected.SetColumnSpan(_chkExpStabilization, 2);
        
        _grpExpExpected.Controls.Add(_tlpExpExpected);
    }
    
    private void AddExpectedRow(string label, ref NumericUpDown num, double value, double min, double max, int decimals, int row)
    {
        var lbl = new Label { Text = label, AutoSize = true, Anchor = AnchorStyles.Left | AnchorStyles.Right };
        num = new NumericUpDown
        {
            Dock = DockStyle.Fill,
            DecimalPlaces = decimals,
            Minimum = (decimal)min,
            Maximum = (decimal)max,
            Value = (decimal)value,
            Increment = (decimal)Math.Pow(10, -decimals)
        };
        _tlpExpExpected.Controls.Add(lbl, 0, row);
        _tlpExpExpected.Controls.Add(num, 1, row);
    }
    
    private void AddExpectedRowInt(string label, ref NumericUpDown num, int value, int min, int max, int row)
    {
        var lbl = new Label { Text = label, AutoSize = true, Anchor = AnchorStyles.Left | AnchorStyles.Right };
        num = new NumericUpDown
        {
            Dock = DockStyle.Fill,
            DecimalPlaces = 0,
            Minimum = min,
            Maximum = max,
            Value = value
        };
        _tlpExpExpected.Controls.Add(lbl, 0, row);
        _tlpExpExpected.Controls.Add(num, 1, row);
    }
    
    private void InitializeExpValidationGroup()
    {
        _grpExpValidation = new GroupBox
        {
            Text = "Validation Summary",
            Dock = DockStyle.Fill,
            Padding = new Padding(5)
        };
        
        var layout = new TableLayoutPanel
        {
            Dock = DockStyle.Fill,
            ColumnCount = 1,
            RowCount = 4
        };
        layout.ColumnStyles.Add(new ColumnStyle(SizeType.Percent, 100F));
        layout.RowStyles.Add(new RowStyle(SizeType.AutoSize)); // Status label
        layout.RowStyles.Add(new RowStyle(SizeType.AutoSize)); // Validate button
        layout.RowStyles.Add(new RowStyle(SizeType.Percent, 50F)); // Criteria list
        layout.RowStyles.Add(new RowStyle(SizeType.Percent, 50F)); // Report
        
        // Status label
        _lblExpValidationStatus = new Label
        {
            Text = "Not validated yet",
            Dock = DockStyle.Fill,
            Font = new Font(Font.FontFamily, 12, FontStyle.Bold),
            TextAlign = ContentAlignment.MiddleCenter,
            Height = 40,
            BackColor = SystemColors.ControlLight
        };
        layout.Controls.Add(_lblExpValidationStatus, 0, 0);
        
        // Validate button
        _btnExpValidate = new Button
        {
            Text = "Validate Results",
            Dock = DockStyle.Fill,
            Height = 30
        };
        _btnExpValidate.Click += BtnExpValidate_Click;
        layout.Controls.Add(_btnExpValidate, 0, 1);
        
        // Criteria ListView
        _lvExpCriteria = new ListView
        {
            Dock = DockStyle.Fill,
            View = View.Details,
            FullRowSelect = true,
            GridLines = true
        };
        _lvExpCriteria.Columns.Add("Criterion", 150);
        _lvExpCriteria.Columns.Add("Status", 60);
        _lvExpCriteria.Columns.Add("Expected", 100);
        _lvExpCriteria.Columns.Add("Actual", 100);
        layout.Controls.Add(_lvExpCriteria, 0, 2);
        
        // Report TextBox
        _txtExpReport = new TextBox
        {
            Dock = DockStyle.Fill,
            Multiline = true,
            ReadOnly = true,
            ScrollBars = ScrollBars.Both,
            Font = new Font("Consolas", 9F),
            WordWrap = false
        };
        layout.Controls.Add(_txtExpReport, 0, 3);
        
        _grpExpValidation.Controls.Add(layout);
    }
    
    private void LoadBuiltInExperiments()
    {
        _cmbExpSelect.Items.Clear();
        _cmbExpSelect.Items.Add("(New Custom Experiment)");
        
        foreach (var exp in ExperimentFactory.AvailableExperiments)
        {
            _cmbExpSelect.Items.Add(exp.Name);
        }
        
        _cmbExpSelect.SelectedIndex = 0;
    }
    
    #region Event Handlers
    
    private void CmbExpSelect_SelectedIndexChanged(object? sender, EventArgs e)
    {
        if (_cmbExpSelect.SelectedIndex <= 0)
        {
            // New custom experiment
            _currentExperiment = new ExperimentDefinition
            {
                Name = "Custom Experiment",
                Description = "Create a new custom experiment with your own parameters and expected results."
            };
            LoadExperimentToUI(_currentExperiment);
            return;
        }
        
        string expName = _cmbExpSelect.SelectedItem?.ToString() ?? "";
        var builtIn = ExperimentFactory.GetByName(expName);
        
        if (builtIn != null)
        {
            _currentExperiment = CreateDefinitionFromBuiltIn(builtIn);
            LoadExperimentToUI(_currentExperiment);
        }
    }
    
    private void BtnExpNew_Click(object? sender, EventArgs e)
    {
        _currentExperiment = new ExperimentDefinition
        {
            Name = $"Experiment_{DateTime.Now:yyyyMMdd_HHmmss}",
            Description = "New custom experiment"
        };
        
        // Copy current Settings tab values to experiment config
        _currentExperiment.Config = GetCurrentConfigFromUI();
        
        _cmbExpSelect.SelectedIndex = 0;
        LoadExperimentToUI(_currentExperiment);
        
        AppendConsole("[Experiment] New experiment created from current settings\n");
    }
    
    private void BtnExpLoad_Click(object? sender, EventArgs e)
    {
        using var dlg = new OpenFileDialog
        {
            Filter = "JSON files (*.json)|*.json|All files (*.*)|*.*",
            Title = "Load Experiment Definition"
        };
        
        if (dlg.ShowDialog(this) != DialogResult.OK)
            return;
        
        try
        {
            string json = File.ReadAllText(dlg.FileName);
            var options = new JsonSerializerOptions { PropertyNameCaseInsensitive = true };
            _currentExperiment = JsonSerializer.Deserialize<ExperimentDefinition>(json, options);
            
            if (_currentExperiment != null)
            {
                LoadExperimentToUI(_currentExperiment);
                _cmbExpSelect.SelectedIndex = 0; // Custom
                AppendConsole($"[Experiment] Loaded: {_currentExperiment.Name}\n");
            }
        }
        catch (Exception ex)
        {
            MessageBox.Show($"Error loading experiment: {ex.Message}", "Error", MessageBoxButtons.OK, MessageBoxIcon.Error);
        }
    }
    
    private void BtnExpSave_Click(object? sender, EventArgs e)
    {
        if (_currentExperiment == null)
        {
            MessageBox.Show("No experiment to save.", "Warning", MessageBoxButtons.OK, MessageBoxIcon.Warning);
            return;
        }
        
        // Update experiment from UI
        UpdateExperimentFromUI();
        
        using var dlg = new SaveFileDialog
        {
            Filter = "JSON files (*.json)|*.json",
            FileName = $"experiment-{_currentExperiment.Name.Replace(" ", "-").ToLowerInvariant()}-{DateTime.Now:yyyyMMdd}.json",
            Title = "Save Experiment Definition"
        };
        
        if (dlg.ShowDialog(this) != DialogResult.OK)
            return;
        
        try
        {
            var options = new JsonSerializerOptions { WriteIndented = true };
            string json = JsonSerializer.Serialize(_currentExperiment, options);
            File.WriteAllText(dlg.FileName, json, Encoding.UTF8);
            AppendConsole($"[Experiment] Saved: {dlg.FileName}\n");
        }
        catch (Exception ex)
        {
            MessageBox.Show($"Error saving experiment: {ex.Message}", "Error", MessageBoxButtons.OK, MessageBoxIcon.Error);
        }
    }
    
    private void BtnExpApplyToSim_Click(object? sender, EventArgs e)
    {
        if (_currentExperiment == null)
        {
            MessageBox.Show("No experiment selected.", "Warning", MessageBoxButtons.OK, MessageBoxIcon.Warning);
            return;
        }
        
        // Update experiment from UI first
        UpdateExperimentFromUI();
        
        // Apply config to Settings tab
        ApplyExperimentConfigToSettingsTab(_currentExperiment.Config);
        
        // Switch to Settings tab
        tabControl1.SelectedTab = tabPage_Settings;
        
        AppendConsole($"[Experiment] Applied '{_currentExperiment.Name}' config to Settings tab\n");
    }
    
    private void BtnExpValidate_Click(object? sender, EventArgs e)
    {
        if (_currentExperiment == null)
        {
            MessageBox.Show("No experiment selected.", "Warning", MessageBoxButtons.OK, MessageBoxIcon.Warning);
            return;
        }
        
        // Update expected values from UI
        UpdateExperimentFromUI();
        
        // Collect actual results from simulation
        var actual = CollectActualResultsFromSimulation();
        _currentExperiment.Actual = actual;
        _currentExperiment.LastRunAt = DateTime.UtcNow;
        
        // Validate
        var validation = ExperimentValidator.Validate(_currentExperiment.Expected, actual);
        _currentExperiment.Validation = validation;
        
        // Update UI
        UpdateValidationUI(validation);
        
        // Generate report
        _txtExpReport.Text = ExperimentValidator.GenerateReport(_currentExperiment);
        
        AppendConsole($"[Experiment] Validation: {validation.Summary}\n");
    }
    
    #endregion
    
    #region Helper Methods
    
    private ExperimentDefinition CreateDefinitionFromBuiltIn(IExperiment builtIn)
    {
        var config = builtIn.GetConfig();
        
        return new ExperimentDefinition
        {
            Name = builtIn.Name,
            Description = builtIn.Description,
            Config = ExperimentConfig.FromStartupConfig(config),
            Expected = GetDefaultExpectedForExperiment(builtIn.Name),
            Notes = $"Built-in experiment: {builtIn.Name}"
        };
    }
    
    private ExpectedResults GetDefaultExpectedForExperiment(string name)
    {
        // Set sensible defaults based on experiment type
        return name switch
        {
            "Vacuum Genesis" => new ExpectedResults
            {
                SpectralDimensionTarget = 4.0,
                SpectralDimensionMin = 3.0,
                SpectralDimensionMax = 5.0,
                HeavyClusterCountMin = 0,
                HeavyClusterCountMax = 10,
                LargestClusterFractionMax = 0.20,
                ExpectPhaseTransition = true,
                ExpectStabilization = true
            },
            "Mass Nucleation" => new ExpectedResults
            {
                SpectralDimensionTarget = 3.5,
                SpectralDimensionMin = 2.5,
                SpectralDimensionMax = 4.5,
                HeavyClusterCountMin = 5,
                HeavyClusterCountMax = 50,
                HeavyMassMin = 10.0,
                HeavyMassMax = 200.0,
                LargestClusterFractionMax = 0.25,
                ExpectPhaseTransition = true,
                ExpectStabilization = true
            },
            "Bio-Folding (DNA Hairpin)" => new ExpectedResults
            {
                SpectralDimensionTarget = 2.5,
                SpectralDimensionMin = 1.5,
                SpectralDimensionMax = 3.5,
                HeavyClusterCountMin = 1,
                HeavyClusterCountMax = 3,
                LargestClusterFractionMax = 0.50,
                ExpectPhaseTransition = true,
                ExpectStabilization = true
            },
            _ => new ExpectedResults()
        };
    }
    
    private void LoadExperimentToUI(ExperimentDefinition exp)
    {
        _txtExpDescription.Text = exp.Description;
        _txtExpNotes.Text = exp.Notes;
        
        // Load expected values
        _numExpSpectralDimTarget.Value = (decimal)exp.Expected.SpectralDimensionTarget;
        _numExpSpectralDimMin.Value = (decimal)exp.Expected.SpectralDimensionMin;
        _numExpSpectralDimMax.Value = (decimal)exp.Expected.SpectralDimensionMax;
        _numExpHeavyClusterMin.Value = exp.Expected.HeavyClusterCountMin;
        _numExpHeavyClusterMax.Value = exp.Expected.HeavyClusterCountMax;
        _numExpHeavyMassMin.Value = (decimal)exp.Expected.HeavyMassMin;
        _numExpHeavyMassMax.Value = (decimal)exp.Expected.HeavyMassMax;
        _numExpLargestClusterMax.Value = (decimal)(exp.Expected.LargestClusterFractionMax * 100);
        _numExpFinalTempMax.Value = (decimal)exp.Expected.FinalTemperatureMax;
        _chkExpPhaseTransition.Checked = exp.Expected.ExpectPhaseTransition;
        _chkExpStabilization.Checked = exp.Expected.ExpectStabilization;
        
        // Reset validation display
        _lblExpValidationStatus.Text = exp.Validation != null 
            ? exp.Validation.Summary 
            : "Not validated yet";
        _lblExpValidationStatus.BackColor = exp.Validation?.Passed == true 
            ? Color.LightGreen 
            : exp.Validation?.Passed == false 
                ? Color.LightCoral 
                : SystemColors.ControlLight;
        
        if (exp.Validation != null)
        {
            UpdateValidationUI(exp.Validation);
        }
        else
        {
            _lvExpCriteria.Items.Clear();
            _txtExpReport.Clear();
        }
    }
    
    private void UpdateExperimentFromUI()
    {
        if (_currentExperiment == null) return;
        
        _currentExperiment.Notes = _txtExpNotes.Text;
        
        _currentExperiment.Expected.SpectralDimensionTarget = (double)_numExpSpectralDimTarget.Value;
        _currentExperiment.Expected.SpectralDimensionMin = (double)_numExpSpectralDimMin.Value;
        _currentExperiment.Expected.SpectralDimensionMax = (double)_numExpSpectralDimMax.Value;
        _currentExperiment.Expected.HeavyClusterCountMin = (int)_numExpHeavyClusterMin.Value;
        _currentExperiment.Expected.HeavyClusterCountMax = (int)_numExpHeavyClusterMax.Value;
        _currentExperiment.Expected.HeavyMassMin = (double)_numExpHeavyMassMin.Value;
        _currentExperiment.Expected.HeavyMassMax = (double)_numExpHeavyMassMax.Value;
        _currentExperiment.Expected.LargestClusterFractionMax = (double)_numExpLargestClusterMax.Value / 100.0;
        _currentExperiment.Expected.FinalTemperatureMax = (double)_numExpFinalTempMax.Value;
        _currentExperiment.Expected.ExpectPhaseTransition = _chkExpPhaseTransition.Checked;
        _currentExperiment.Expected.ExpectStabilization = _chkExpStabilization.Checked;
    }
    
    private ExperimentConfig GetCurrentConfigFromUI()
    {
        return new ExperimentConfig
        {
            NodeCount = (int)numNodeCount.Value,
            TotalSteps = (int)numTotalSteps.Value,
            TargetDegree = (int)numTargetDegree.Value,
            InitialEdgeProb = (double)numInitialEdgeProb.Value,
            InitialExcitedProb = (double)numInitialExcitedProb.Value,
            GravitationalCoupling = (double)numGravitationalCoupling.Value,
            VacuumEnergyScale = (double)numVacuumEnergyScale.Value,
            HotStartTemperature = (double)numHotStartTemperature.Value,
            AnnealingCoolingRate = (double)numAnnealingCoolingRate.Value,
            DecoherenceRate = (double)numDecoherenceRate.Value,
            AdaptiveThresholdSigma = (double)numAdaptiveThresholdSigma.Value,
            WarmupDuration = (double)numWarmupDuration.Value,
            GravityTransitionDuration = (double)numGravityTransitionDuration.Value,
            Temperature = (double)numTemperature.Value,
            LambdaState = (double)numLambdaState.Value,
            EdgeTrialProbability = (double)numEdgeTrialProb.Value,
            MeasurementThreshold = (double)numMeasurementThreshold.Value,
            UseSpectralGeometry = chkSpectralGeometry.Checked,
            UseNetworkGravity = chkNetworkGravity.Checked,
            UseQuantumDrivenStates = chkQuantumDriven.Checked,
            UseSpinorField = chkSpinorField.Checked,
            UseVacuumFluctuations = chkVacuumFluctuations.Checked,
            UseHotStartAnnealing = chkHotStartAnnealing.Checked,
            FractalLevels = (int)numFractalLevels.Value,
            FractalBranchFactor = (int)numFractalBranchFactor.Value
        };
    }
    
    private void ApplyExperimentConfigToSettingsTab(ExperimentConfig config)
    {
        numNodeCount.Value = config.NodeCount;
        numTotalSteps.Value = config.TotalSteps;
        numTargetDegree.Value = config.TargetDegree;
        numInitialEdgeProb.Value = (decimal)config.InitialEdgeProb;
        numInitialExcitedProb.Value = (decimal)config.InitialExcitedProb;
        numGravitationalCoupling.Value = (decimal)config.GravitationalCoupling;
        numVacuumEnergyScale.Value = (decimal)config.VacuumEnergyScale;
        numHotStartTemperature.Value = (decimal)config.HotStartTemperature;
        numAnnealingCoolingRate.Value = (decimal)config.AnnealingCoolingRate;
        numDecoherenceRate.Value = (decimal)config.DecoherenceRate;
        numAdaptiveThresholdSigma.Value = (decimal)config.AdaptiveThresholdSigma;
        numWarmupDuration.Value = (decimal)config.WarmupDuration;
        numGravityTransitionDuration.Value = (decimal)config.GravityTransitionDuration;
        numTemperature.Value = (decimal)config.Temperature;
        numLambdaState.Value = (decimal)config.LambdaState;
        numEdgeTrialProb.Value = (decimal)config.EdgeTrialProbability;
        numMeasurementThreshold.Value = (decimal)config.MeasurementThreshold;
        chkSpectralGeometry.Checked = config.UseSpectralGeometry;
        chkNetworkGravity.Checked = config.UseNetworkGravity;
        chkQuantumDriven.Checked = config.UseQuantumDrivenStates;
        chkSpinorField.Checked = config.UseSpinorField;
        chkVacuumFluctuations.Checked = config.UseVacuumFluctuations;
        chkHotStartAnnealing.Checked = config.UseHotStartAnnealing;
        numFractalLevels.Value = config.FractalLevels;
        numFractalBranchFactor.Value = config.FractalBranchFactor;
    }
    
    private ActualResults CollectActualResultsFromSimulation()
    {
        // Get final values from time series
        int finalStep = _seriesSteps.Count > 0 ? _seriesSteps[^1] : _simApi.LiveStep;
        int totalSteps = _simApi.LastConfig?.TotalSteps ?? (int)numTotalSteps.Value;
        int nodeCount = _simApi.LastConfig?.NodeCount ?? (int)numNodeCount.Value;
        
        double finalSpectralDim = _seriesSpectralDimension.Count > 0 ? _seriesSpectralDimension[^1] : _simApi.LiveSpectralDim;
        double initialSpectralDim = _seriesSpectralDimension.Count > 0 ? _seriesSpectralDimension[0] : finalSpectralDim;
        
        int heavyCount = _seriesHeavyCount.Count > 0 ? _seriesHeavyCount[^1] : 0;
        double heavyMass = _seriesHeavyMass.Count > 0 ? _seriesHeavyMass[^1] : _simApi.LiveHeavyMass;
        int largestCluster = _seriesLargestCluster.Count > 0 ? _seriesLargestCluster[^1] : _simApi.LiveLargestCluster;
        
        double finalTemp = _seriesNetworkTemperature.Count > 0 ? _seriesNetworkTemperature[^1] : _simApi.LiveTemp;
        double effectiveG = _seriesEffectiveG.Count > 0 ? _seriesEffectiveG[^1] : _simApi.LiveEffectiveG;
        
        int excitedCount = _seriesExcited.Count > 0 ? _seriesExcited[^1] : _simApi.LiveExcited;
        int strongEdges = _seriesStrongEdges.Count > 0 ? _seriesStrongEdges[^1] : _simApi.LiveStrongEdges;
        double correlation = _seriesCorr.Count > 0 ? _seriesCorr[^1] : _simApi.LiveCorrelation;
        double qNorm = _seriesQNorm.Count > 0 ? _seriesQNorm[^1] : _simApi.LiveQNorm;
        double entanglement = _seriesEntanglement.Count > 0 ? _seriesEntanglement[^1] : _simApi.LiveEntanglement;
        
        // Calculate variance in last 20%
        double spectralDimVariance = CalculateVarianceLast20Percent(_seriesSpectralDimension);
        double heavyMassVariance = CalculateVarianceLast20Percent(_seriesHeavyMass);
        
        double wallClock = (DateTime.UtcNow - _simApi.SimulationWallClockStart).TotalSeconds;
        
        return ExperimentValidator.CollectResults(
            finalStep, totalSteps, !_isModernRunning, null,
            finalSpectralDim, initialSpectralDim,
            heavyCount, heavyMass, largestCluster, nodeCount,
            finalTemp, effectiveG, excitedCount, strongEdges,
            correlation, qNorm, entanglement,
            spectralDimVariance, heavyMassVariance, wallClock);
    }
    
    private double CalculateVarianceLast20Percent(List<double> series)
    {
        if (series.Count < 5) return 0;
        
        int start = (int)(series.Count * 0.8);
        var last20 = series.Skip(start).ToList();
        
        if (last20.Count == 0) return 0;
        
        double mean = last20.Average();
        double variance = last20.Sum(x => (x - mean) * (x - mean)) / last20.Count;
        return variance;
    }
    
    private double CalculateVarianceLast20Percent(List<int> series)
    {
        return CalculateVarianceLast20Percent(series.Select(x => (double)x).ToList());
    }
    
    private void UpdateValidationUI(ValidationSummary validation)
    {
        // Update status label
        _lblExpValidationStatus.Text = validation.Summary;
        _lblExpValidationStatus.BackColor = validation.Passed ? Color.LightGreen : Color.LightCoral;
        
        // Update criteria list
        _lvExpCriteria.Items.Clear();
        
        foreach (var criterion in validation.Criteria)
        {
            var item = new ListViewItem(criterion.Name);
            item.SubItems.Add(criterion.Passed ? "? PASS" : "? FAIL");
            item.SubItems.Add(criterion.Expected);
            item.SubItems.Add(criterion.Actual);
            
            item.BackColor = criterion.Passed 
                ? Color.LightGreen 
                : criterion.Level == ValidationLevel.Critical 
                    ? Color.LightCoral 
                    : Color.LightYellow;
            
            _lvExpCriteria.Items.Add(item);
        }
        
        // Auto-size columns
        foreach (ColumnHeader col in _lvExpCriteria.Columns)
        {
            col.Width = -2; // Auto-size to content
        }
    }
    
    #endregion
}
