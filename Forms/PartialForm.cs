using RQSimulation;
using RQSimulation.GPUOptimized;
using System.Text;
namespace RqSimForms;

public partial class Form_Main
{


    private void InitializeSynthesisTab()
    {
        // ????????? ??????? ?? ??? ?????: ????? ??????, ?????? ?????
        TableLayoutPanel layoutPanel = new TableLayoutPanel
        {
            Dock = DockStyle.Fill,
            RowCount = 2,
            ColumnCount = 1
        };
        layoutPanel.RowStyles.Add(new RowStyle(SizeType.Percent, 40F)); // ????????? ?????
        layoutPanel.RowStyles.Add(new RowStyle(SizeType.Percent, 60F)); // ??????

        // TextBox ??? ?????????? ???????
        TextBox synthesisTextBox = new TextBox
        {
            Dock = DockStyle.Fill,
            Multiline = true,
            ReadOnly = true,
            ScrollBars = ScrollBars.Both,
            WordWrap = false,
            Font = new Font("Consolas", 9F),
            BackColor = Color.White,
            ForeColor = Color.Black,
            Name = "synthesisTextBox"
        };

        // Panel ??? ???????
        Panel chartPanel = new Panel
        {
            Dock = DockStyle.Fill,
            BackColor = Color.White,
            Name = "synthesisChartPanel"
        };
        chartPanel.Paint += SynthesisChartPanel_Paint;

        layoutPanel.Controls.Add(synthesisTextBox, 0, 0);
        layoutPanel.Controls.Add(chartPanel, 0, 1);

        tabPage_Sythnesis.Controls.Add(layoutPanel);
    }
    private int FindNodeAtPosition(Point clickPos)
    {
        if (_simulationEngine?.Graph == null) return -1;
        var graph = _simulationEngine.Graph;
        int n = graph.N;
        if (n <= 0) return -1;
        // ensure positions exist
        if (_cachedNodePositions == null || _cachedNodePositions.Length != n)
            BuildNodePositions();
        float hitRadius = Math.Max(12f, 500f / Math.Max(50, n));
        float hitRadiusSq = hitRadius * hitRadius;
        int best = -1; float bestDist = float.MaxValue;
        for (int i = 0; i < n; i++)
        {
            var p = _cachedNodePositions[i];
            float dx = clickPos.X - p.X; float dy = clickPos.Y - p.Y; float d2 = dx * dx + dy * dy;
            if (d2 < hitRadiusSq && d2 < bestDist) { bestDist = d2; best = i; }
        }
        return best;
    }

    private void BuildNodePositions()
    {
        if (_simulationEngine?.Graph == null) return;
        var graph = _simulationEngine.Graph;
        int n = graph.N; if (n <= 0) { _cachedNodePositions = null; return; }
        _cachedNodePositions = new PointF[n];
        int w = Math.Max(drawingPanel.Width, 1); int h = Math.Max(drawingPanel.Height, 1);
        float margin = 20f;
        if (_useDynamicCoords && graph.Coordinates != null && graph.Coordinates.Length == n)
        {
            double minX = graph.Coordinates.Min(c => c.X);
            double maxX = graph.Coordinates.Max(c => c.X);
            double minY = graph.Coordinates.Min(c => c.Y);
            double maxY = graph.Coordinates.Max(c => c.Y);
            double centerX = (minX + maxX) * 0.5; double centerY = (minY + maxY) * 0.5;
            // uniform scale using max radial distance to avoid flattening into a line
            double maxR = 0.0;
            for (int i = 0; i < n; i++) { var (cx, cy) = graph.Coordinates[i]; double dx = cx - centerX; double dy = cy - centerY; double r = Math.Sqrt(dx * dx + dy * dy); if (r > maxR) maxR = r; }
            if (maxR < 1e-9) maxR = 1.0;
            float radiusPixels = (Math.Min(w, h) - 2 * margin) * 0.5f;
            double scale = radiusPixels / maxR;
            for (int i = 0; i < n; i++)
            {
                var (cx, cy) = graph.Coordinates[i];
                double dx = cx - centerX; double dy = cy - centerY;
                float x = (float)(w / 2.0 + dx * scale);
                float y = (float)(h / 2.0 + dy * scale);
                _cachedNodePositions[i] = new PointF(x, y);
            }
        }
        else
        {
            // stable circle layout
            float radius = (Math.Min(w, h) - 2 * margin) * 0.5f;
            for (int i = 0; i < n; i++)
            {
                double ang = 2 * Math.PI * i / n;
                float x = w / 2f + radius * (float)Math.Cos(ang);
                float y = h / 2f + radius * (float)Math.Sin(ang);
                _cachedNodePositions[i] = new PointF(x, y);
            }
        }
    }


    private void AnalyzeSynthesis(List<string> avalancheStats)
    {
        if (InvokeRequired)
        {
            Invoke(() => AnalyzeSynthesis(avalancheStats));
            return;
        }

        var synthesisTextBox = tabPage_Sythnesis.Controls
            .OfType<TableLayoutPanel>().FirstOrDefault()?
            .Controls.OfType<TextBox>().FirstOrDefault();

        if (synthesisTextBox == null || avalancheStats.Count < 2)
            return;

        // ?????? CSV ??????
        string[] headers = avalancheStats[0].Split(',');
        int volumeIdx = Array.IndexOf(headers, "volume");
        int deltaMassIdx = Array.IndexOf(headers, "deltaHeavyMass");

        if (volumeIdx < 0 || deltaMassIdx < 0)
        {
            synthesisTextBox.Text = "??????: ?? ??????? ??????? volume ??? deltaHeavyMass ? CSV";
            return;
        }

        synthesisData = new List<(int, double)>();
        synthesisCount = 0;
        fissionCount = 0;
        int neutralCount = 0;

        double totalSynthesisMass = 0.0;
        double totalFissionMass = 0.0;
        double maxSynthesis = double.MinValue;
        double maxFission = double.MaxValue;

        for (int i = 1; i < avalancheStats.Count; i++)
        {
            string[] fields = avalancheStats[i].Split(',');

            if (fields.Length <= Math.Max(volumeIdx, deltaMassIdx))
                continue;

            if (int.TryParse(fields[volumeIdx], out int volume) &&
                double.TryParse(fields[deltaMassIdx], out double deltaMass))
            {
                synthesisData.Add((volume, deltaMass));

                if (deltaMass > 0.001)
                {
                    synthesisCount++;
                    totalSynthesisMass += deltaMass;
                    if (deltaMass > maxSynthesis) maxSynthesis = deltaMass;
                }
                else if (deltaMass < -0.001)
                {
                    fissionCount++;
                    totalFissionMass += Math.Abs(deltaMass);
                    if (deltaMass < maxFission) maxFission = deltaMass;
                }
                else
                {
                    neutralCount++;
                }
            }
        }

        // ????????? ????????? ?????
        var sb = new StringBuilder();
        sb.AppendLine("=== ?????? ???????/??????? ?????? ????????? ===\n");
        sb.AppendLine($"????? ????? ????????????????: {avalancheStats.Count - 1}");
        sb.AppendLine();

        sb.AppendLine("--- ????????????? ????? ---");
        sb.AppendLine($"? ?????? (deltaHeavyMass > 0):      {synthesisCount,4} ????? ({100.0 * synthesisCount / (avalancheStats.Count - 1):F1}%)");
        sb.AppendLine($"? ?????? (deltaHeavyMass < 0):      {fissionCount,4} ????? ({100.0 * fissionCount / (avalancheStats.Count - 1):F1}%)");
        sb.AppendLine($"? ??????????? (|deltaMass| < 0.001): {neutralCount,4} ????? ({100.0 * neutralCount / (avalancheStats.Count - 1):F1}%)");
        sb.AppendLine();

        sb.AppendLine("--- ?????????? ?? ????? ---");
        sb.AppendLine($"????????? ????? ?????????????:   {totalSynthesisMass,8:F2}");
        sb.AppendLine($"????????? ????? ?????????:       {totalFissionMass,8:F2}");
        sb.AppendLine($"?????? ??????? ?????:            {totalSynthesisMass - totalFissionMass,8:F2}");
        sb.AppendLine();

        if (synthesisCount > 0)
        {
            sb.AppendLine($"???????????? ?????? ?? ??????:   {maxSynthesis,8:F2}");
            sb.AppendLine($"??????? ?????? (??? ?????????):  {totalSynthesisMass / synthesisCount,8:F2}");
        }

        if (fissionCount > 0)
        {
            sb.AppendLine($"???????????? ?????? ?? ??????:   {maxFission,8:F2}");
            sb.AppendLine($"??????? ?????? (??? ?????????):  {-totalFissionMass / fissionCount,8:F2}");
        }

        sb.AppendLine();
        sb.AppendLine("--- ???????? ? ??????? ---");
        sb.AppendLine("?????? > ?????? ? ???????????? (???????????? ?????? ?????????)");
        sb.AppendLine("?????? > ?????? ? ???????????? (?????????? ????????)");
        sb.AppendLine();
        sb.AppendLine("?????? ???? ?????????? ??????????? deltaHeavyMass ?? ?????? ??????.");
        sb.AppendLine("??????? ????? (y=0) ????????? ?????? ? ??????.");

        synthesisTextBox.Text = sb.ToString();

        // ?????????????? ??????
        var chartPanel = tabPage_Sythnesis.Controls
            .OfType<TableLayoutPanel>().FirstOrDefault()?
            .Controls.OfType<Panel>().FirstOrDefault();
        chartPanel?.Invalidate();
    }


    /// <summary>
    /// Optimized graph drawing with thread-safe GDI resources.
    /// FIX: Creates fresh bitmap AND fresh GDI objects each time to avoid race conditions.
    /// GDI objects (Pen, Brush) are NOT thread-safe - must create new ones in background thread!
    /// 
    /// FIX 2: All early returns inside Task.Run now go through finally block to reset _isDrawing flag.
    /// </summary>
    private void DrawGraph()
    {
        var graph = _simulationEngine?.Graph;
        if (graph == null) return;

        // Non-blocking lock - if already drawing, skip this frame
        if (Interlocked.CompareExchange(ref _isDrawing, 1, 0) != 0) return;

        int panelWidth = drawingPanel.Width;
        int panelHeight = drawingPanel.Height;

        // Validate panel dimensions before starting async work
        if (panelWidth <= 0 || panelHeight <= 0)
        {
            Interlocked.Exchange(ref _isDrawing, 0);
            return;
        }

        // Capture display settings (read on UI thread)
        bool useDynamic = _useDynamicCoords;
        bool showHeavyOnly = _displayShowHeavyOnly;
        double weightThreshold = _displayWeightThreshold;

        Task.Run(() =>
        {
            Bitmap? produced = null;
            Pen[]? localPens = null;
            SolidBrush[]? localBrushes = null;
            Pen? blackPen = null;

            try
            {
                // Re-check graph availability (could be disposed during async transition)
                var localGraph = _simulationEngine?.Graph;
                if (localGraph == null)
                {
                    // Early exit - finally will reset _isDrawing
                    return;
                }

                if (useDynamic) NormalizeToCircle();
                BuildNodePositions();
                var localPositions = _cachedNodePositions;
                if (localPositions == null || localPositions.Length == 0)
                {
                    // Early exit - finally will reset _isDrawing
                    return;
                }

                HashSet<int>? heavyNodes = null;
                if (showHeavyOnly)
                {
                    try
                    {
                        var clusters = localGraph.GetStrongCorrelationClusters(localGraph.GetAdaptiveHeavyThreshold());
                        heavyNodes = new HashSet<int>(clusters.SelectMany(c => c));
                        if (heavyNodes.Count == 0) heavyNodes = null;
                    }
                    catch { heavyNodes = null; }
                }

                double effectiveWeightThreshold = weightThreshold;
                if (effectiveWeightThreshold >= 0.95) effectiveWeightThreshold = 0.0;

                // Create thread-local GDI objects (CRITICAL for thread safety)
                localPens = CreateEdgePenSet();
                localBrushes = CreateNodeBrushSet();
                blackPen = new Pen(Color.Black, 1f);

                // Create fresh bitmap for this frame
                produced = new Bitmap(panelWidth, panelHeight);
                using (var g = Graphics.FromImage(produced))
                {
                    g.Clear(Color.White);
                    g.SmoothingMode = System.Drawing.Drawing2D.SmoothingMode.AntiAlias;
                    int n = localGraph.N;

                    // Draw edges using thread-local Pens
                    for (int i = 0; i < n; i++)
                    {
                        if (heavyNodes != null && !heavyNodes.Contains(i)) continue;
                        if (i >= localPositions.Length) continue;
                        var p1 = localPositions[i];
                        foreach (int j in localGraph.Neighbors(i))
                        {
                            if (j <= i) continue;
                            if (heavyNodes != null && !heavyNodes.Contains(j)) continue;
                            if (j >= localPositions.Length) continue;
                            double w = localGraph.Weights[i, j];
                            if (w < effectiveWeightThreshold) continue;

                            var pen = GetEdgePenFromSet(localPens, w);
                            var p2 = localPositions[j];
                            g.DrawLine(pen, p1, p2);
                        }
                    }

                    // Draw nodes using thread-local Brushes
                    int nodeCountToDraw = heavyNodes == null ? localGraph.N : heavyNodes.Count;
                    float nodeSize = Math.Max(4, Math.Min(11, 500f / Math.Max(1, nodeCountToDraw)));
                    for (int i = 0; i < localGraph.N; i++)
                    {
                        if (heavyNodes != null && !heavyNodes.Contains(i)) continue;
                        if (i >= localPositions.Length) continue;
                        var p = localPositions[i];

                        var brush = GetNodeBrushFromSet(localBrushes, localGraph, i);
                        g.FillEllipse(brush, p.X - nodeSize / 2, p.Y - nodeSize / 2, nodeSize, nodeSize);
                        g.DrawEllipse(blackPen, p.X - nodeSize / 2, p.Y - nodeSize / 2, nodeSize, nodeSize);
                    }
                }

                // Pass bitmap to UI thread for display
                if (!IsDisposed && IsHandleCreated)
                {
                    var bitmapToDisplay = produced;
                    produced = null; // Transfer ownership to UI thread

                    BeginInvoke(new Action(() =>
                    {
                        if (IsDisposed)
                        {
                            bitmapToDisplay?.Dispose();
                            return;
                        }

                        // Swap bitmaps
                        var oldBitmap = canvasBitmap;
                        canvasBitmap = bitmapToDisplay;
                        oldBitmap?.Dispose();

                        // Force immediate repaint
                        drawingPanel.Invalidate();
                        drawingPanel.Update();
                    }));
                }
            }
            catch (Exception ex)
            {
                if (!IsDisposed && IsHandleCreated)
                {
                    try
                    {
                        BeginInvoke(new Action(() => AppendConsole($"[DrawGraph] error: {ex.Message}\n")));
                    }
                    catch { /* Ignore if form is closing */ }
                }
            }
            finally
            {
                // Dispose thread-local GDI resources
                if (localPens != null)
                {
                    foreach (var pen in localPens) pen?.Dispose();
                }
                if (localBrushes != null)
                {
                    foreach (var brush in localBrushes) brush?.Dispose();
                }
                blackPen?.Dispose();

                // Dispose bitmap if not transferred to UI thread
                produced?.Dispose();

                // CRITICAL: Always reset drawing flag - this was the bug!
                Interlocked.Exchange(ref _isDrawing, 0);
            }
        });
    }

    private Pen[] CreateEdgePenSet()
    {
        // Thread-safe pen set for background drawing (matches DrawingOptimizations colors/weights)
        return
        [
            new Pen(Color.FromArgb(220, 255, 0, 0), 3f),      // VeryStrong > 0.7
            new Pen(Color.FromArgb(200, 255, 69, 0), 2.5f),   // Strong > 0.5
            new Pen(Color.FromArgb(180, 255, 140, 0), 2f),    // Medium > 0.3
            new Pen(Color.FromArgb(120, 218, 165, 32), 1.2f), // Weak > 0.15
            new Pen(Color.FromArgb(80, 160, 160, 160), 0.6f)  // VeryWeak
        ];
    }

    private SolidBrush[] CreateNodeBrushSet()
    {
        // Thread-safe brush set for background drawing (matches DrawingOptimizations colors)
        // Index: [0]=Excited, [1]=Refractory, [2]=Rest, [3]=Composite, [4]=Clock
        return
        [
            new SolidBrush(Color.Red),       // Excited
            new SolidBrush(Color.Orange),    // Refractory
            new SolidBrush(Color.LightBlue), // Rest
            new SolidBrush(Color.DarkRed),   // Composite
            new SolidBrush(Color.Blue)       // Clock
        ];
    }

    private static Pen GetEdgePenFromSet(Pen[] localPens, double w)
    {
        // Same thresholds as GetEdgePen in DrawingOptimizations.cs
        return w switch
        {
            > 0.7 => localPens[0],  // VeryStrong
            > 0.5 => localPens[1],  // Strong
            > 0.3 => localPens[2],  // Medium
            > 0.15 => localPens[3], // Weak
            _ => localPens[4]       // VeryWeak
        };
    }

    private static SolidBrush GetNodeBrushFromSet(SolidBrush[] localBrushes, RQGraph localGraph, int i)
    {
        // Same logic as GetNodeBrush in DrawingOptimizations.cs
        if (localGraph.PhysicsProperties[i].IsClock)
            return localBrushes[4]; // Clock = Blue
        if (localGraph.PhysicsProperties[i].Type == ParticleType.Composite)
            return localBrushes[3]; // Composite = DarkRed

        return localGraph.State[i] switch
        {
            NodeState.Excited => localBrushes[0],    // Red
            NodeState.Refractory => localBrushes[1], // Orange
            _ => localBrushes[2]                     // Rest = LightBlue
        };
    }

    private void PanelOnChart_Paint(object? sender, PaintEventArgs e)
    {
        var data = _simApi.Dispatcher.ForceGetDisplayDataImmediate(timeoutMs: 20);
        if (data.DecimatedSteps.Length == 0)
        {
            e.Graphics.Clear(Color.White);
            e.Graphics.DrawString("Нет данных", new Font("Arial", 10), Brushes.Gray, 10, 10);
            return;
        }
        DrawSimpleLineChartFast(e.Graphics, panelOnChart, data.DecimatedSteps, data.DecimatedExcited, "Excited nodes", Color.Red);
    }

    private void PanelHeavyChart_Paint(object? sender, PaintEventArgs e)
    {
        var data = _simApi.Dispatcher.ForceGetDisplayDataImmediate(timeoutMs: 20);
        if (data.DecimatedSteps.Length == 0)
        {
            e.Graphics.Clear(Color.White);
            e.Graphics.DrawString("Нет данных", new Font("Arial", 10), Brushes.Gray, 10, 10);
            return;
        }
        DrawSimpleLineChartFast(e.Graphics, panelHeavyChart, data.DecimatedSteps, data.DecimatedHeavyMass, "Heavy mass", Color.DarkOrange);
    }

    private void PanelClusterChart_Paint(object? sender, PaintEventArgs e)
    {
        var data = _simApi.Dispatcher.ForceGetDisplayDataImmediate(timeoutMs: 20);
        if (data.DecimatedSteps.Length == 0)
        {
            e.Graphics.Clear(Color.White);
            e.Graphics.DrawString("Нет данных", new Font("Arial", 10), Brushes.Gray, 10, 10);
            return;
        }
        DrawSimpleLineChartFast(e.Graphics, panelClusterChart, data.DecimatedSteps, data.DecimatedLargestCluster, "Largest cluster size", Color.Blue);
    }

    private void PanelEnergyChart_Paint(object? sender, PaintEventArgs e)
    {
        var data = _simApi.Dispatcher.ForceGetDisplayDataImmediate(timeoutMs: 20);
        if (data.DecimatedSteps.Length == 0)
        {
            e.Graphics.Clear(Color.White);
            e.Graphics.DrawString("Нет данных", new Font("Arial", 10), Brushes.Gray, 10, 10);
            return;
        }

        if (data.DecimatedNetworkTemp.Length == 0)
        {
            DrawSimpleLineChartFast(e.Graphics, panelEnergyChart, data.DecimatedSteps, data.DecimatedEnergy, "Total energy", Color.Green);
            return;
        }

        DrawDualLineChartFast(
            e.Graphics,
            panelEnergyChart,
            data.DecimatedSteps,
            data.DecimatedEnergy,
            data.DecimatedNetworkTemp,
            "Energy vs Network Temp",
            Color.ForestGreen,
            Color.MediumSlateBlue);
    }

    /// <summary>
    /// Fast chart drawing using pre-decimated array data (no LINQ, no List conversions)
    /// </summary>
    private void DrawSimpleLineChartFast(Graphics g, Panel panel, int[] xSteps, int[] values, string title, Color color)
    {
        g.Clear(Color.White);
        if (xSteps.Length == 0 || values.Length == 0)
        {
            g.DrawString("Нет данных", new Font("Arial", 10), Brushes.Gray, 10, 10);
            return;
        }

        int width = panel.Width;
        int height = panel.Height;
        int marginLeft = 50;
        int marginRight = 20;
        int marginTop = 30;
        int marginBottom = 40;
        int plotWidth = width - marginLeft - marginRight;
        int plotHeight = height - marginTop - marginBottom;

        if (plotWidth <= 10 || plotHeight <= 10)
            return;

        int minX = xSteps[0];
        int maxX = xSteps[^1];
        int minY = int.MaxValue, maxY = int.MinValue;
        for (int i = 0; i < values.Length; i++)
        {
            if (values[i] < minY) minY = values[i];
            if (values[i] > maxY) maxY = values[i];
        }

        if (maxY == minY)
        {
            maxY += 1;
            minY -= 1;
        }

        // Draw axes
        using Pen axisPen = new(Color.Black, 1);
        g.DrawLine(axisPen, marginLeft, marginTop, marginLeft, marginTop + plotHeight);
        g.DrawLine(axisPen, marginLeft, marginTop + plotHeight, marginLeft + plotWidth, marginTop + plotHeight);

        using Font f = new("Consolas", 9f);
        g.DrawString(title, f, Brushes.Black, marginLeft, 5);
        g.DrawString(maxY.ToString(), f, Brushes.Black, 5, marginTop - 5);
        g.DrawString(minY.ToString(), f, Brushes.Black, 5, marginTop + plotHeight - 12);
        g.DrawString(minX.ToString(), f, Brushes.Black, marginLeft, marginTop + plotHeight + 5);
        g.DrawString(maxX.ToString(), f, Brushes.Black, marginLeft + plotWidth - 40, marginTop + plotHeight + 5);

        // Draw line
        using Pen linePen = new(color, 1.6f);
        double xRange = maxX - minX;
        double yRange = maxY - minY;
        if (xRange <= 0) xRange = 1;

        float prevX = 0, prevY = 0;
        for (int i = 0; i < xSteps.Length && i < values.Length; i++)
        {
            float x = marginLeft + (float)((xSteps[i] - minX) / xRange * plotWidth);
            float y = marginTop + (float)((maxY - values[i]) / yRange * plotHeight);
            if (i > 0)
                g.DrawLine(linePen, prevX, prevY, x, y);
            prevX = x;
            prevY = y;
        }
    }

    /// <summary>
    /// Fast chart drawing for double arrays
    /// </summary>
    private void DrawSimpleLineChartFast(Graphics g, Panel panel, int[] xSteps, double[] values, string title, Color color)
    {
        g.Clear(Color.White);
        if (xSteps.Length == 0 || values.Length == 0)
        {
            g.DrawString("Нет данных", new Font("Arial", 10), Brushes.Gray, 10, 10);
            return;
        }

        int width = panel.Width;
        int height = panel.Height;
        int marginLeft = 50;
        int marginRight = 20;
        int marginTop = 30;
        int marginBottom = 40;
        int plotWidth = width - marginLeft - marginRight;
        int plotHeight = height - marginTop - marginBottom;

        if (plotWidth <= 10 || plotHeight <= 10)
            return;

        int minX = xSteps[0];
        int maxX = xSteps[^1];
        double minY = double.MaxValue, maxY = double.MinValue;
        for (int i = 0; i < values.Length; i++)
        {
            if (values[i] < minY) minY = values[i];
            if (values[i] > maxY) maxY = values[i];
        }

        if (Math.Abs(maxY - minY) < 1e-9)
        {
            maxY += 1.0;
            minY -= 1.0;
        }

        // Draw axes
        using Pen axisPen = new(Color.Black, 1);
        g.DrawLine(axisPen, marginLeft, marginTop, marginLeft, marginTop + plotHeight);
        g.DrawLine(axisPen, marginLeft, marginTop + plotHeight, marginLeft + plotWidth, marginTop + plotHeight);

        using Font f = new("Consolas", 9f);
        g.DrawString(title, f, Brushes.Black, marginLeft, 5);
        g.DrawString(maxY.ToString("F2"), f, Brushes.Black, 5, marginTop - 5);
        g.DrawString(minY.ToString("F2"), f, Brushes.Black, 5, marginTop + plotHeight - 12);
        g.DrawString(minX.ToString(), f, Brushes.Black, marginLeft, marginTop + plotHeight + 5);
        g.DrawString(maxX.ToString(), f, Brushes.Black, marginLeft + plotWidth - 40, marginTop + plotHeight + 5);

        // Draw line
        using Pen linePen = new(color, 1.6f);
        double xRange = maxX - minX;
        double yRange = maxY - minY;
        if (xRange <= 0) xRange = 1;

        float prevX = 0, prevY = 0;
        for (int i = 0; i < xSteps.Length && i < values.Length; i++)
        {
            float x = marginLeft + (float)((xSteps[i] - minX) / xRange * plotWidth);
            float y = marginTop + (float)((maxY - values[i]) / yRange * plotHeight);
            if (i > 0)
                g.DrawLine(linePen, prevX, prevY, x, y);
            prevX = x;
            prevY = y;
        }
    }
    private void DrawDualLineChartFast(Graphics g, Panel panel, int[] xSteps, double[] primaryValues, double[] secondaryValues,
        string title, Color primaryColor, Color secondaryColor)
    {
        g.Clear(Color.White);
        if (xSteps.Length == 0 || primaryValues.Length == 0 || secondaryValues.Length == 0)
        {
            g.DrawString("Нет данных", new Font("Arial", 10), Brushes.Gray, 10, 10);
            return;
        }

        int width = panel.Width;
        int height = panel.Height;
        int marginLeft = 50;
        int marginRight = 50; // Increased for right axis
        int marginTop = 30;
        int marginBottom = 40;
        int plotWidth = width - marginLeft - marginRight;
        int plotHeight = height - marginTop - marginBottom;

        if (plotWidth <= 10 || plotHeight <= 10)
            return;

        int sampleCount = Math.Min(xSteps.Length, Math.Min(primaryValues.Length, secondaryValues.Length));
        if (sampleCount == 0)
        {
            g.DrawString("Нет данных", new Font("Arial", 10), Brushes.Gray, 10, 10);
            return;
        }

        int minX = xSteps[0];
        int maxX = xSteps[sampleCount - 1];

        // Calculate ranges separately
        double minY1 = double.MaxValue, maxY1 = double.MinValue;
        double minY2 = double.MaxValue, maxY2 = double.MinValue;

        for (int i = 0; i < sampleCount; i++)
        {
            double v1 = primaryValues[i];
            if (v1 < minY1) minY1 = v1;
            if (v1 > maxY1) maxY1 = v1;

            double v2 = secondaryValues[i];
            if (v2 < minY2) minY2 = v2;
            if (v2 > maxY2) maxY2 = v2;
        }

        if (Math.Abs(maxY1 - minY1) < 1e-9) { maxY1 += 1.0; minY1 -= 1.0; }
        if (Math.Abs(maxY2 - minY2) < 1e-9) { maxY2 += 1.0; minY2 -= 1.0; }

        using Pen axisPen = new(Color.Black, 1);
        // Left Y axis
        g.DrawLine(axisPen, marginLeft, marginTop, marginLeft, marginTop + plotHeight);
        // Right Y axis
        g.DrawLine(axisPen, marginLeft + plotWidth, marginTop, marginLeft + plotWidth, marginTop + plotHeight);
        // X axis
        g.DrawLine(axisPen, marginLeft, marginTop + plotHeight, marginLeft + plotWidth, marginTop + plotHeight);

        using Font labelFont = new("Consolas", 9f);
        g.DrawString(title, labelFont, Brushes.Black, marginLeft, 5);

        // Left axis labels (Primary)
        using Brush primaryBrush = new SolidBrush(primaryColor);
        g.DrawString(maxY1.ToString("F2"), labelFont, primaryBrush, 5, marginTop - 5);
        g.DrawString(minY1.ToString("F2"), labelFont, primaryBrush, 5, marginTop + plotHeight - 12);

        // Right axis labels (Secondary)
        using Brush secondaryBrush = new SolidBrush(secondaryColor);
        string maxLabel2 = maxY2.ToString("F2");
        string minLabel2 = minY2.ToString("F2");
        g.DrawString(maxLabel2, labelFont, secondaryBrush, width - marginRight + 5, marginTop - 5);
        g.DrawString(minLabel2, labelFont, secondaryBrush, width - marginRight + 5, marginTop + plotHeight - 12);

        // X axis labels
        g.DrawString(minX.ToString(), labelFont, Brushes.Black, marginLeft, marginTop + plotHeight + 5);
        g.DrawString(maxX.ToString(), labelFont, Brushes.Black, marginLeft + plotWidth - 40, marginTop + plotHeight + 5);

        double xRange = maxX - minX;
        double yRange1 = maxY1 - minY1;
        double yRange2 = maxY2 - minY2;
        if (xRange <= 0) xRange = 1;

        using Pen primaryPen = new(primaryColor, 1.6f);
        using Pen secondaryPen = new(secondaryColor, 1.6f);

        // Draw Primary
        float prevX = 0, prevY = 0;
        for (int i = 0; i < sampleCount; i++)
        {
            float x = marginLeft + (float)((xSteps[i] - minX) / xRange * plotWidth);
            float y = marginTop + (float)((maxY1 - primaryValues[i]) / yRange1 * plotHeight);
            if (i > 0) g.DrawLine(primaryPen, prevX, prevY, x, y);
            prevX = x; prevY = y;
        }

        // Draw Secondary
        prevX = 0; prevY = 0;
        for (int i = 0; i < sampleCount; i++)
        {
            float x = marginLeft + (float)((xSteps[i] - minX) / xRange * plotWidth);
            float y = marginTop + (float)((maxY2 - secondaryValues[i]) / yRange2 * plotHeight);
            if (i > 0) g.DrawLine(secondaryPen, prevX, prevY, x, y);
            prevX = x; prevY = y;
        }

        // Legend
        int legendX = marginLeft + 10;
        int legendY = marginTop + 10;
        using Font legendFont = new("Consolas", 8.5f);
        
        g.FillRectangle(primaryBrush, legendX, legendY, 12, 12);
        g.DrawString("Energy", legendFont, Brushes.Black, legendX + 18, legendY - 1);
        
        legendY += 18;
        g.FillRectangle(secondaryBrush, legendX, legendY, 12, 12);
        g.DrawString("Network Temp", legendFont, Brushes.Black, legendX + 18, legendY - 1);
    }



    // Legacy method kept for compatibility (used by Synthesis chart)
    private void DrawSimpleLineChart(Graphics g, Panel panel, List<int> xSteps, List<double> values, string title, Color color)
    {
        g.Clear(Color.White);
        if (xSteps.Count == 0 || values.Count == 0)
        {
            g.DrawString("??? ??????", new Font("Arial", 10), Brushes.Gray, 10, 10);
            return;
        }
        int width = panel.Width;
        int height = panel.Height;
        int marginLeft = 50;
        int marginRight = 20;
        int marginTop = 30;
        int marginBottom = 40;
        int plotWidth = width - marginLeft - marginRight;
        int plotHeight = height - marginTop - marginBottom;
        if (plotWidth <= 10 || plotHeight <= 10)
            return;
        int minX = xSteps.First();
        int maxX = xSteps.Last();
        double minY = values.Min();
        double maxY = values.Max();
        if (Math.Abs(maxY - minY) < 1e-9)
        {
            maxY += 1.0;
            minY -= 1.0;
        }
        using Pen axisPen = new Pen(Color.Black, 1);
        g.DrawLine(axisPen, marginLeft, marginTop, marginLeft, marginTop + plotHeight);
        g.DrawLine(axisPen, marginLeft, marginTop + plotHeight, marginLeft + plotWidth, marginTop + plotHeight);
        using Font f = new("Consolas", 9f);
        g.DrawString(title, f, Brushes.Black, marginLeft, 5);
        g.DrawString(maxY.ToString("F2"), f, Brushes.Black, 5, marginTop - 5);
        g.DrawString(minY.ToString("F2"), f, Brushes.Black, 5, marginTop + plotHeight - 12);
        g.DrawString(minX.ToString(), f, Brushes.Black, marginLeft, marginTop + plotHeight + 5);
        g.DrawString(maxX.ToString(), f, Brushes.Black, marginLeft + plotWidth - 40, marginTop + plotHeight + 5);
        PointF? prev = null;
        using Pen linePen = new(color, 1.6f);
        int count = Math.Min(xSteps.Count, values.Count);
        for (int i = 0; i < count; i++)
        {
            float x = marginLeft + (float)((xSteps[i] - minX) / (double)(maxX - minX)) * plotWidth;
            float y = marginTop + (float)((maxY - values[i]) / (maxY - minY)) * plotHeight;
            if (prev != null)
                g.DrawLine(linePen, prev.Value, new PointF(x, y));
            prev = new PointF(x, y);
        }
    }

    private void SynthesisChartPanel_Paint(object? sender, PaintEventArgs e)
    {
        if (synthesisData == null || synthesisData.Count == 0)
        {
            // ?????????? ?????????, ???? ?????? ???
            e.Graphics.DrawString("?????? ????? ???????? ????? ?????????? ?????????",
                new Font("Arial", 12), Brushes.Gray, 10, 10);
            return;
        }

        Graphics g = e.Graphics;
        g.SmoothingMode = System.Drawing.Drawing2D.SmoothingMode.AntiAlias;
        g.Clear(Color.White);

        Panel panel = (Panel)sender!;
        int width = panel.Width;
        int height = panel.Height;

        // ???????
        int marginLeft = 60;
        int marginRight = 20;
        int marginTop = 20;
        int marginBottom = 50;

        int plotWidth = width - marginLeft - marginRight;
        int plotHeight = height - marginTop - marginBottom;

        if (plotWidth <= 0 || plotHeight <= 0)
            return;

        // ??????? ????????? ??????
        int minVolume = synthesisData.Min(d => d.volume);
        int maxVolume = synthesisData.Max(d => d.volume);
        double minDelta = synthesisData.Min(d => d.deltaMass);
        double maxDelta = synthesisData.Max(d => d.deltaMass);

        // ????????? ????????? ?????
        double deltaRange = Math.Max(Math.Abs(maxDelta), Math.Abs(minDelta)) * 1.1;
        minDelta = -deltaRange;
        maxDelta = deltaRange;

        int volumeRange = maxVolume - minVolume;
        if (volumeRange == 0) volumeRange = 1;

        // ?????? ???
        using (Pen axisPen = new Pen(Color.Black, 2))
        {
            // Y-???
            g.DrawLine(axisPen, marginLeft, marginTop, marginLeft, height - marginBottom);
            // X-???
            g.DrawLine(axisPen, marginLeft, height - marginBottom, width - marginRight, height - marginBottom);
        }

        // ??????? ????? y=0 (????????? ?????? ? ??????)
        int zeroY = marginTop + (int)(plotHeight * (maxDelta / (maxDelta - minDelta)));
        using (Pen zeroPen = new Pen(Color.Red, 2))
        {
            zeroPen.DashStyle = System.Drawing.Drawing2D.DashStyle.Dash;
            g.DrawLine(zeroPen, marginLeft, zeroY, width - marginRight, zeroY);
        }

        // ??????? ????
        using (Font axisFont = new Font("Arial", 10))
        {
            // Y-???
            g.DrawString("Delta Heavy Mass", axisFont, Brushes.Black, 5, marginTop + plotHeight / 2 - 30);
            g.DrawString($"{maxDelta:F1}", axisFont, Brushes.Black, 5, marginTop);
            g.DrawString("0", axisFont, Brushes.Red, 5, zeroY - 7);
            g.DrawString($"{minDelta:F1}", axisFont, Brushes.Black, 5, height - marginBottom - 15);

            // X-???
            StringFormat sf = new StringFormat { Alignment = StringAlignment.Center };
            g.DrawString("Avalanche Volume", axisFont, Brushes.Black,
                marginLeft + plotWidth / 2, height - marginBottom + 30, sf);
            g.DrawString($"{minVolume}", axisFont, Brushes.Black, marginLeft, height - marginBottom + 5);
            g.DrawString($"{maxVolume}", axisFont, Brushes.Black, width - marginRight - 30, height - marginBottom + 5);
        }

        // ?????? ?????
        foreach (var (volume, deltaMass) in synthesisData)
        {
            int x = marginLeft + (int)((volume - minVolume) * plotWidth / (double)volumeRange);
            int y = marginTop + (int)((maxDelta - deltaMass) * plotHeight / (maxDelta - minDelta));

            Color pointColor = deltaMass > 0 ? Color.Blue : (deltaMass < 0 ? Color.OrangeRed : Color.Gray);
            int pointSize = 4;

            using (Brush pointBrush = new SolidBrush(pointColor))
            {
                g.FillEllipse(pointBrush, x - pointSize / 2, y - pointSize / 2, pointSize, pointSize);
            }
        }

        // ???????
        int legendX = width - marginRight - 150;
        int legendY = marginTop + 10;

        using (Font legendFont = new Font("Arial", 9))
        {
            g.FillRectangle(Brushes.Blue, legendX, legendY, 10, 10);
            g.DrawString("?????? (? > 0)", legendFont, Brushes.Black, legendX + 15, legendY - 2);

            g.FillRectangle(Brushes.OrangeRed, legendX, legendY + 20, 10, 10);
            g.DrawString("?????? (? < 0)", legendFont, Brushes.Black, legendX + 15, legendY + 18);
        }
    }


    private void AppendSummary(string text)
    {
        if (summaryTextBox == null) return;
        if (summaryTextBox.InvokeRequired)
        {
            summaryTextBox.Invoke(() => AppendSummary(text));
            return;
        }
        summaryTextBox.AppendText(text + (text.EndsWith("\n") ? string.Empty : "\n"));
    }

    private void AppendConsole(string text)
    {
        // Guard against disposed controls/form
        if (IsDisposed || consoleTextBox == null || consoleTextBox.IsDisposed)
            return;

        if (consoleTextBox.InvokeRequired)
        {
            if (!consoleTextBox.IsHandleCreated) return;
            consoleTextBox.BeginInvoke(new Action(() => AppendConsole(text)));
            return;
        }

        consoleTextBox.AppendText(text);
        if (chkAutoScrollConsole != null && !chkAutoScrollConsole.IsDisposed && chkAutoScrollConsole.Checked)
        {
            consoleTextBox.SelectionStart = consoleTextBox.Text.Length;
            consoleTextBox.ScrollToCaret();
        }
    }

    protected override void OnFormClosing(FormClosingEventArgs e)
    {
        // Cancel modern simulation if running
        _modernCts?.Cancel();
        _simApi.Cleanup();

        // Dispose cached GDI resources to prevent memory leaks
        DisposeGdiResources();

        base.OnFormClosing(e);
    }

    private void NormalizeToCircle()
    {
        if (_simulationEngine?.Graph == null || _simulationEngine.Graph.Coordinates == null) return;
        int n = _simulationEngine.Graph.N;
        var coords = _simulationEngine.Graph.Coordinates;
        // ?????????? ?? ???????? ????, ????? ????????? ????????????? ???????
        var indexed = coords.Select((c, i) => (i, angle: Math.Atan2(c.Y, c.X))).OrderBy(x => x.angle).ToArray();
        for (int k = 0; k < indexed.Length; k++)
        {
            double ang = 2 * Math.PI * k / n;
            coords[indexed[k].i] = (Math.Cos(ang), Math.Sin(ang));
        }
    }



    /// <summary>
    /// Reads simulation configuration from UI controls with clamped ranges
    /// </summary>
    private SimulationConfig GetConfigFromUI()
    {
        var config = new SimulationConfig
        {
            NodeCount = (int)numNodeCount.Value,
            InitialEdgeProb = (double)numInitialEdgeProb.Value, // From UI: connectivity at Big Bang
            InitialExcitedProb = (double)numInitialExcitedProb.Value,
            TargetDegree = (int)numTargetDegree.Value,
            LambdaState = (double)numLambdaState.Value,
            Temperature = (double)numTemperature.Value,
            EdgeTrialProbability = (double)numEdgeTrialProb.Value,
            MeasurementThreshold = (double)numMeasurementThreshold.Value,
            Seed = 42,
            TotalSteps = Math.Max(1, (int)numTotalSteps.Value),
            LogEvery = 1,
            BaselineWindow = 50,
            FirstImpulse = -1, // Disabled - use hot start annealing
            ImpulsePeriod = -1, // Disabled
            CalibrationStep = -1, // Disabled
            VisualizationInterval = 1,
            MeasurementLogInterval = 50,
            FractalLevels = (int)numFractalLevels.Value,
            FractalBranchFactor = (int)numFractalBranchFactor.Value,

            UseQuantumDrivenStates = chkQuantumDriven.Checked,
            UseSpacetimePhysics = chkSpacetimePhysics.Checked,
            UseSpinorField = chkSpinorField.Checked,
            UseVacuumFluctuations = chkVacuumFluctuations.Checked,
            UseBlackHolePhysics = chkBlackHolePhysics.Checked,
            UseYangMillsGauge = chkYangMillsGauge.Checked,
            UseEnhancedKleinGordon = chkEnhancedKleinGordon.Checked,
            UseInternalTime = chkInternalTime.Checked,
            UseSpectralGeometry = chkSpectralGeometry.Checked,
            UseQuantumGraphity = chkQuantumGraphity.Checked,

            // === Physics Constants (from new panel) ===
            GravitationalCoupling = (double)numGravitationalCoupling.Value,
            VacuumEnergyScale = (double)numVacuumEnergyScale.Value,
            AnnealingCoolingRate = (double)numAnnealingCoolingRate.Value,
            DecoherenceRate = (double)numDecoherenceRate.Value,
            HotStartTemperature = (double)numHotStartTemperature.Value,
            InitialNetworkTemperature = (double)numHotStartTemperature.Value,

            // === Extended Physics Flags ===
            UseRelationalTime = chkRelationalTime.Checked,
            UseRelationalYangMills = chkRelationalYangMills.Checked,
            UseNetworkGravity = chkNetworkGravity.Checked,
            UseUnifiedPhysicsStep = chkUnifiedPhysicsStep.Checked,
            EnforceGaugeConstraints = chkEnforceGaugeConstraints.Checked,
            UseCausalRewiring = chkCausalRewiring.Checked,
            UseTopologicalProtection = chkTopologicalProtection.Checked,
            ValidateEnergyConservation = chkValidateEnergyConservation.Checked,
            UseMexicanHatPotential = chkMexicanHatPotential.Checked,
            UseHotStartAnnealing = chkHotStartAnnealing.Checked,
            UseGeometryMomenta = chkGeometryMomenta.Checked,
            UseTopologicalCensorship = chkTopologicalCensorship.Checked,

            UseEventBasedSimulation = true
        };

        // Clamp ranges for robustness
        if (config.InitialExcitedProb < 0) config.InitialExcitedProb = 0;
        if (config.InitialExcitedProb > 1) config.InitialExcitedProb = 1;
        if (config.EdgeTrialProbability < 0) config.EdgeTrialProbability = 0;
        if (config.EdgeTrialProbability > 1) config.EdgeTrialProbability = 1;
        if (config.MeasurementThreshold < 0) config.MeasurementThreshold = 0;
        if (config.VisualizationInterval < 1) config.VisualizationInterval = 1;
        if (config.MeasurementLogInterval < 1) config.MeasurementLogInterval = 1;
        if (config.InitialEdgeProb < 0) config.InitialEdgeProb = 0;
        if (config.InitialEdgeProb > 1) config.InitialEdgeProb = 1;
        if (config.GravitationalCoupling < 0) config.GravitationalCoupling = 0;
        if (config.DecoherenceRate < 0) config.DecoherenceRate = 0;

        return config;
    }

    /// <summary>
    /// Updates dashboard metrics with current simulation state
    /// </summary>
    private void UpdateDashboard(int step, int totalSteps, int excited, double heavyMass,
        int largestCluster, int strongEdges, string phase, double qNorm,
        double entanglement, double correlation)
    {
        if (InvokeRequired)
        {
            Invoke(() => UpdateDashboard(step, totalSteps, excited, heavyMass, largestCluster,
                strongEdges, phase, qNorm, entanglement, correlation));
            return;
        }

        valDashNodes.Text = _simulationEngine?.Graph?.N.ToString() ?? "0";
        valDashTotalSteps.Text = totalSteps.ToString();
        valDashCurrentStep.Text = step.ToString();
        valDashExcited.Text = excited.ToString();
        valDashHeavyMass.Text = heavyMass.ToString("F2");
        valDashLargestCluster.Text = largestCluster.ToString();
        valDashStrongEdges.Text = strongEdges.ToString();
        valDashPhase.Text = phase;
        valDashQNorm.Text = qNorm.ToString("F6");
        valDashEntanglement.Text = entanglement.ToString("F6");
        valDashCorrelation.Text = correlation.ToString("F6");
        valDashStatus.Text = _isModernRunning ? "Running..." : "Ready";

        // New spectral metrics
        double spectralDim = _simApi.LiveSpectralDim;
        double effectiveG = _simApi.LiveEffectiveG;
        double networkTemp = _simApi.LiveTemp;

        // Calculate gSuppression = effectiveG / targetG (after warmup+transition)
        double targetG = _simApi.LastConfig?.GravitationalCoupling ?? 0.2;
        double gSuppression = targetG > 0 ? effectiveG / targetG : 1.0;
        gSuppression = Math.Clamp(gSuppression, 0.0, 2.0);

        valDashSpectralDim.Text = spectralDim.ToString("F3");
        valDashEffectiveG.Text = effectiveG.ToString("F4");
        valDashGSuppression.Text = gSuppression.ToString("F3");
        valDashNetworkTemp.Text = networkTemp.ToString("F3");

        // Color-code spectral dimension for quick visual feedback
        if (spectralDim < 1.5)
            valDashSpectralDim.ForeColor = Color.Red;
        else if (spectralDim > 4.0)
            valDashSpectralDim.ForeColor = Color.DarkOrange;
        else
            valDashSpectralDim.ForeColor = Color.Green;

        // Color-code gSuppression
        if (gSuppression < 0.3)
            valDashGSuppression.ForeColor = Color.Red;
        else if (gSuppression < 0.7)
            valDashGSuppression.ForeColor = Color.DarkOrange;
        else
            valDashGSuppression.ForeColor = Color.Black;

        // === Extended Live Metrics from Dispatcher ===
        double clusterRatio = _simApi.Dispatcher.LiveClusterRatio;
        double avgDegree = _simApi.Dispatcher.LiveAvgDegree;
        int edgeCount = _simApi.Dispatcher.LiveEdgeCount;
        int componentCount = _simApi.Dispatcher.LiveComponentCount;

        // Color-code LargestCluster by ratio (giant cluster warning)
        if (clusterRatio >= 0.7)
            valDashLargestCluster.ForeColor = Color.Red;
        else if (clusterRatio >= 0.5)
            valDashLargestCluster.ForeColor = Color.DarkOrange;
        else
            valDashLargestCluster.ForeColor = Color.Green;

        // Update LargestCluster text to show ratio
        valDashLargestCluster.Text = $"{largestCluster} ({clusterRatio:P0})";
    }



    // Helpers to update newly added UI controls
    private void UpdateStatusBar(int currentStep, int totalSteps, int currentOn, double avgExcited, double? heavyMass)
    {
        if (statusLabelSteps is null) return;

        statusLabelSteps.Text = $"Step: {currentStep}/{totalSteps}";
        statusLabelExcited.Text = $"Excited: {currentOn} (avg {avgExcited:F2})";

        // Extended status bar with key topology metrics
        double clusterRatio = _simApi.Dispatcher.LiveClusterRatio;
        double avgDegree = _simApi.Dispatcher.LiveAvgDegree;
        int edgeCount = _simApi.Dispatcher.LiveEdgeCount;
        int componentCount = _simApi.Dispatcher.LiveComponentCount;

        statusLabelHeavyMass.Text = $"Giant:{clusterRatio:P0} | E:{edgeCount} | <k>:{avgDegree:F1} | Comp:{componentCount}";
    }

    private void UpdateRunSummary(int totalSteps, int currentStep, double avgExcited, int maxExcited, int avalanches, bool measurementConfigured, bool measurementTriggered)
    {
        if (valTotalSteps is null) return;
        valTotalSteps.Text = totalSteps.ToString();
        valCurrentStep.Text = currentStep.ToString();
        valExcitedAvg.Text = avgExcited.ToString("F2");
        valExcitedMax.Text = maxExcited.ToString();
        valAvalancheCount.Text = avalanches.ToString();
        valMeasurementStatus.Text = measurementConfigured
            ? (measurementTriggered ? "TRIGGERED" : "READY")
            : "NOT CONFIGURED";
    }

    private void UpdateLiveMetrics(double globalNbr, double globalSpont, int strongEdges, int largestCluster, double heavyMass, bool spectrumOn)
    {
        if (valGlobalNbr is null) return;
        valGlobalNbr.Text = globalNbr.ToString("F3");
        valGlobalSpont.Text = globalSpont.ToString("F3");
        valStrongEdges.Text = strongEdges.ToString();
        valLargestCluster.Text = largestCluster.ToString();
        valHeavyMass.Text = heavyMass.ToString("F2");
        valSpectrumInfo.Text = spectrumOn ? "on" : "off";
    }

    private void AddImportantEvent(int step, string type, string detail)
    {
        if (lvEvents is null) return;
        var item = new ListViewItem(new[] { step.ToString(), type, detail });
        lvEvents.Items.Add(item);
        if (lvEvents.Items.Count > 1000)
            lvEvents.Items.RemoveAt(0);
    }

}


