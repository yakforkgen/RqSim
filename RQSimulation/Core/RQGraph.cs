using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Runtime.CompilerServices;

namespace RQSimulation
{
    // Core type declarations (single source)
    public enum ParticleType { Vacuum, Fermion, Boson, Composite }
    public enum NodeState { Rest = 0, Excited = 1, Refractory = 2 }
    public struct NodePhysics
    {
        public ParticleType Type;
        public double Mass;
        public double Charge;
        public double Spin;
        public int Occupation;
        public int MaxOccupation;
        public bool IsClock;
    }

    public partial class RQGraph
    {




        public const double HeavyClusterThreshold = PhysicsConstants.DefaultHeavyClusterThreshold;
        public const int HeavyClusterMinSize = 3;

        public int N { get; private set; }
        public bool[,] Edges { get; private set; }
        public double[,] Weights { get; private set; }
        public NodeState[] State { get; private set; }

        /// <summary>
        /// Incremented each time graph topology changes (edges added/removed).
        /// Used by ParallelEventEngine to detect when coloring needs recomputation.
        /// </summary>
        public int TopologyVersion { get; private set; }
        private NodeState[] _nextState;
        private bool[,] _measurementBond;
        private int[] _degree;
        private int[] _charges;
        private int[] _refractoryCounter;
        private readonly Random _rng;
        private readonly int _targetDegree;
        private readonly double _lambdaState;
        private double _temperature;
        private readonly double _edgeTrialProbability;
        private double _measurementThreshold;
        private int[] _targetDegreePerNode;
        public NodePhysics[] PhysicsProperties;

        /// <summary>
        /// DEPRECATED: Use graph-based distances (ShortestPathDistance, GetGraphDistance) for physics calculations.
        /// Coordinates are only for visualization purposes and should not affect physical dynamics.
        /// </summary>
        [Obsolete("Use ONLY for rendering, not for physics! Use ShortestPathDistance() or GetGraphDistance() instead.")]
        public (double X, double Y)[] Coordinates;
        private double[] _nodeEnergy;
        public double[] NodeEnergy => _nodeEnergy;
        private double _adaptiveHeavyThreshold = HeavyClusterThreshold;
        public double AdaptiveHeavyThreshold => _adaptiveHeavyThreshold;
        public double[] ProperTime;
        public double[,] EdgeDelay;
        public sbyte[,] EdgeDirection;
        private double[,] _targetDistance;
        private const double EpsDist = 1e-6;
        public double EnergyBreakThreshold { get; set; } = 5.0;
        public double EnergySewThreshold { get; set; } = 1.0;
        private double _hebbRate = 0.25;
        private double _decayRate = 0.002;
        private bool _hasHeavyClusters;

        // Add missing public properties/fields required by other partials and examples
        public int QuantumComponents { get; private set; } = 1;
        public double CorrelationTime { get; private set; } = 0.0;
        public double EdgeDecayLength { get; set; } = 0.2;
        public double PropagationLength { get; set; } = 0.1;
        public double EnergyDiffusionRate { get; set; } = 0.05;

        public int[] FractalLevel; // fractal level per node
        public double[] LocalPotential; // local potential per node
        public double[] StoredEnergy; // checklist energy accumulator
        public bool[] ParticleTag;    // checklist particle birth marker
        public double[,] PathWeight;  // checklist dynamic corridor weights
        // Removed StructuralMass as per RQ-Hypothesis unification
        public int[] Domain; // checklist item 3 large-scale domains

        /// <summary>
        /// GPU Gravity Engine for accelerated network geometry evolution.
        /// When non-null, EvolveNetworkGeometry will use GPU computation.
        /// 
        /// NOTE: Currently prepared but not active due to pre-existing GPU shader
        /// compilation issues (missing ComputeSharp descriptors). Once shaders are
        /// fixed, GPU computation will activate automatically for 10-100x speedup.
        /// </summary>
        public RQSimulation.GPUOptimized.GpuGravityEngine? GpuGravity { get; set; }

        // Diagnostic fields for spectral dimension analysis
        private double _lastSpectralSlope;
        private string _lastSpectralMethod = "none";

        // RQ-FIX: Spectral dimension smoothing to prevent method-switching jumps

        private const double SpectralDimSmoothingAlpha = 0.3; // EMA smoothing factor

        /// <summary>
        /// Last computed slope from eigenvalue or return probability analysis.
        /// Useful for debugging d_S computation.
        /// </summary>
        public double LastSpectralSlope => _lastSpectralSlope;

        /// <summary>
        /// Method used for last spectral dimension computation.
        /// Values: "HeatKernel", "RandomWalk", "Laplacian", "fallback"
        /// </summary>
        public string LastSpectralMethod => _lastSpectralMethod;

        /// <summary>
        /// Smoothed spectral dimension (EMA filtered to prevent jumps).
        /// </summary>
        public double SmoothedSpectralDimension => _smoothedSpectralDimension;

        /// <summary>
        /// Cache the spectral dimension with EMA smoothing.
        /// Called internally by ComputeSpectralDimension.
        /// NOTE: This method updates _smoothedSpectralDimension - the main implementation
        /// is also in RQGraph.SpectralGeometry.cs but works with _cachedSpectralDim.
        /// </summary>
        private void ApplySpectralDimensionSmoothing(double rawValue)
        {
            if (rawValue > 0 && !double.IsNaN(rawValue))
            {
                // Exponential moving average to smooth jumps between methods
                _smoothedSpectralDimension = SpectralDimSmoothingAlpha * rawValue
                                           + (1.0 - SpectralDimSmoothingAlpha) * _smoothedSpectralDimension;
            }
        }

        /// <summary>
        /// Вычисляет Спектральную Размерность графа.
        /// 
        /// RQ-гипотеза предсказывает:
        /// - d_S ≈ 2 на Планковских масштабах (UV)
        /// - d_S ≈ 4 на макроскопических масштабах (IR)
        /// - Переход между режимами при t ~ 1/m_Planck
        /// 
        /// Используем три метода с автоматическим выбором:
        /// 1. Heat Kernel trace: Tr(e^{-tL}) ~ t^{-d_S/2} — наиболее надёжный
        /// 2. Random Walk: P(t) ~ t^{-d_S/2} — физически мотивированный
        /// 3. Laplacian eigenvalues: λ_k ~ k^{2/d_S} — для очень плотных графов
        /// </summary>
        public double ComputeSpectralDimension(int t_max = 100, int num_walkers = 500)
        {
            if (N < 10)
            {
                _lastSpectralMethod = "fallback";
                _lastSpectralSlope = 0;
                return 2.0;
            }

            // Вычисляем характеристики графа
            int edgeCount = 0;
            for (int i = 0; i < N; i++)
                edgeCount += _degree[i];
            edgeCount /= 2;

            int maxEdges = N * (N - 1) / 2;
            double density = maxEdges > 0 ? (double)edgeCount / maxEdges : 0;
            double avgDegree = N > 0 ? 2.0 * edgeCount / N : 0;

            // === RQ-FIX: Выбор метода на основе физических соображений ===

            // Для очень плотных графов (density > 20%) — Laplacian работает лучше
            if (density > 0.20)
            {
                return ComputeSpectralDimensionLaplacian();
            }

            // Для графов с низкой связностью — Random Walk
            // Для умеренных графов — Heat Kernel (более стабильный)
            double result;
            if (avgDegree < 4.0)
            {
                // Разреженный граф — Random Walk
                result = ComputeSpectralDimensionRandomWalk(t_max, num_walkers);
            }
            else
            {
                // Умеренно связанный — Heat Kernel как основной
                double d_hk = ComputeSpectralDimensionHeatKernel(t_max);

                // Если Heat Kernel даёт разумный результат, используем его
                if (d_hk > 1.0 && d_hk < 8.0)
                {
                    result = d_hk;
                }
                else
                {
                    // Иначе пробуем Random Walk
                    double d_rw = ComputeSpectralDimensionRandomWalk(t_max, num_walkers);
                    if (d_rw > 1.0 && d_rw < 8.0)
                    {
                        result = d_rw;
                    }
                    else
                    {
                        // Fallback к Laplacian
                        result = ComputeSpectralDimensionLaplacian();
                    }
                }
            }

            // Кэшируем результат и применяем сглаживание
            ApplySpectralDimensionSmoothing(result);
            return result;
        }

        /// <summary>
        /// Heat Kernel метод для спектральной размерности.
        /// 
        /// Tr(e^{-tL}) ~ t^{-d_S/2} для малых t
        /// 
        /// Это наиболее физически мотивированный метод:
        /// Heat kernel — это пропагатор диффузии на графе,
        /// его асимптотика определяет эффективную размерность.
        /// 
        /// Алгоритм: вычисляем trace через стохастическую оценку
        /// с использованием random vectors (Hutchinson's trick).
        /// </summary>
        private double ComputeSpectralDimensionHeatKernel(int t_max)
        {
            _lastSpectralMethod = "HeatKernel";

            // Строим нормализованный Лапласиан
            var degree = new double[N];
            for (int i = 0; i < N; i++)
            {
                double deg = 0;
                foreach (int j in Neighbors(i))
                    deg += Weights[i, j];
                degree[i] = deg;
            }

            // Используем стохастическую оценку trace через random vectors
            int numVectors = Math.Min(50, N); // Hutchinson estimator
            var traces = new double[t_max + 1];

            var rng = new Random(42); // Детерминистический seed для воспроизводимости

            for (int v = 0; v < numVectors; v++)
            {
                // Random ±1 vector (Rademacher)
                var z = new double[N];
                for (int i = 0; i < N; i++)
                    z[i] = rng.NextDouble() < 0.5 ? -1.0 : 1.0;

                // Эволюция: y(t+1) = (I - dt*L_norm) y(t)
                // где L_norm = D^{-1/2} L D^{-1/2}
                var y = (double[])z.Clone();
                double dt = 0.1; // Шаг времени для diffusion

                for (int t = 1; t <= t_max; t++)
                {
                    var yNew = new double[N];

                    for (int i = 0; i < N; i++)
                    {
                        double sum = y[i]; // Diagonal term (1 - dt*1) = 1 - dt

                        if (degree[i] > 1e-10)
                        {
                            double dInvSqrt_i = 1.0 / Math.Sqrt(degree[i]);

                            foreach (int j in Neighbors(i))
                            {
                                if (degree[j] > 1e-10)
                                {
                                    double dInvSqrt_j = 1.0 / Math.Sqrt(degree[j]);
                                    double L_ij = -Weights[i, j] * dInvSqrt_i * dInvSqrt_j;
                                    sum -= dt * L_ij * y[j];
                                }
                            }

                            // Diagonal: L_ii = 1 (for normalized Laplacian)
                            sum -= dt * y[i];
                        }

                        yNew[i] = sum;
                    }

                    y = yNew;

                    // Trace estimate: z^T e^{-tL} z ≈ Tr(e^{-tL}) для random z
                    double traceEst = 0;
                    for (int i = 0; i < N; i++)
                        traceEst += z[i] * y[i];

                    traces[t] += traceEst / numVectors;
                }
            }

            // Fit: log(Tr) = -d_S/2 * log(t) + const
            // Используем диапазон t где trace ещё не затухло
            var logT = new List<double>();
            var logTrace = new List<double>();

            for (int t = 5; t <= Math.Min(t_max, 50); t++)
            {
                if (traces[t] > 1e-10)
                {
                    logT.Add(Math.Log(t * 0.1)); // Effective time = t * dt
                    logTrace.Add(Math.Log(traces[t]));
                }
            }

            if (logT.Count < 5)
            {
                _lastSpectralSlope = 0;
                return 2.0; // Недостаточно данных
            }

            // Linear regression
            double meanLogT = logT.Average();
            double meanLogTrace = logTrace.Average();

            double num = 0, den = 0;
            for (int i = 0; i < logT.Count; i++)
            {
                double dT = logT[i] - meanLogT;
                double dTr = logTrace[i] - meanLogTrace;
                num += dT * dTr;
                den += dT * dT;
            }

            if (Math.Abs(den) < 1e-10)
            {
                _lastSpectralSlope = 0;
                return 2.0;
            }

            double slope = num / den;
            _lastSpectralSlope = slope;

            // Tr(e^{-tL}) ~ t^{-d_S/2} => slope = -d_S/2 => d_S = -2*slope
            double d_S = -2.0 * slope;

            return Math.Clamp(d_S, 1.0, 10.0);
        }

        /// <summary>
        /// Вычисляет долю узлов в наибольшей связной компоненте.
        /// </summary>
        private double ComputeLargestComponentFraction()
        {
            if (N == 0) return 0.0;

            // Union-Find для быстрого определения компонент
            int[] parent = new int[N];
            int[] rank = new int[N];
            for (int i = 0; i < N; i++) { parent[i] = i; rank[i] = 0; }

            int Find(int x)
            {
                if (parent[x] != x) parent[x] = Find(parent[x]);
                return parent[x];
            }

            void Union(int x, int y)
            {
                int px = Find(x), py = Find(y);
                if (px == py) return;
                if (rank[px] < rank[py]) parent[px] = py;
                else if (rank[px] > rank[py]) parent[py] = px;
                else { parent[py] = px; rank[px]++; }
            }

            // Объединяем по рёбрам
            for (int i = 0; i < N; i++)
            {
                foreach (int j in Neighbors(i))
                {
                    if (i < j) Union(i, j);
                }
            }

            // Считаем размер каждой компоненты
            var componentSizes = new Dictionary<int, int>();
            for (int i = 0; i < N; i++)
            {
                int root = Find(i);
                componentSizes.TryGetValue(root, out int size);
                componentSizes[root] = size + 1;
            }

            int maxSize = componentSizes.Values.Max();
            return (double)maxSize / N;
        }

        /// <summary>
        /// Laplacian-based spectral dimension calculation.
        /// Uses eigenvalue scaling: λ_k ~ k^(2/d_S) for small k.
        /// 
        /// Для RQ-гипотезы: этот метод лучше работает для плотных графов
        /// где random walk быстро термализуется.
        /// 
        /// Улучшения:
        /// - Больше eigenvalues (до 30)
        /// - Взвешенная регрессия (меньший вес для больших k)
        /// - Outlier rejection
        /// </summary>
        public double ComputeSpectralDimensionLaplacian()
        {
            _lastSpectralMethod = "Laplacian";

            if (N < 20)
            {
                _lastSpectralSlope = 0;
                return 2.0;
            }

            try
            {
                // Compute Laplacian spectrum (k smallest eigenvalues)
                int k = Math.Min(30, N / 3);
                var (eigenvalues, _) = ComputeLaplacianSpectrum(k);

                if (eigenvalues == null || eigenvalues.Length < 5)
                {
                    _lastSpectralSlope = 0;
                    return 2.0;
                }

                // Для нормализованного Лапласиана eigenvalues ∈ [0, 2]
                // λ_0 = 0 (тривиальный), λ_1 > 0 (spectral gap)

                var logK = new List<double>();
                var logLambda = new List<double>();
                var weights = new List<double>();

                // Используем eigenvalues начиная с λ_1
                // Вес обратно пропорционален k (первые eigenvalues важнее)
                for (int i = 1; i < eigenvalues.Length && i <= 25; i++)
                {
                    double lambda = eigenvalues[i];

                    // Фильтруем некорректные eigenvalues
                    if (lambda > 1e-6 && lambda < 1.99)
                    {
                        logK.Add(Math.Log(i));
                        logLambda.Add(Math.Log(lambda));
                        weights.Add(1.0 / i); // Weighted regression
                    }
                }

                if (logK.Count < 4)
                {
                    _lastSpectralSlope = 0;
                    return 2.0;
                }

                // Weighted linear regression
                double sumW = weights.Sum();
                double meanLogK = 0, meanLogLambda = 0;
                for (int i = 0; i < logK.Count; i++)
                {
                    meanLogK += weights[i] * logK[i];
                    meanLogLambda += weights[i] * logLambda[i];
                }
                meanLogK /= sumW;
                meanLogLambda /= sumW;

                double num = 0, den = 0;
                for (int i = 0; i < logK.Count; i++)
                {
                    double dK = logK[i] - meanLogK;
                    double dL = logLambda[i] - meanLogLambda;
                    num += weights[i] * dK * dL;
                    den += weights[i] * dK * dK;
                }

                if (Math.Abs(den) < 1e-10)
                {
                    _lastSpectralSlope = 0;
                    return 2.0;
                }

                double slope = num / den;
                _lastSpectralSlope = slope;

                // λ_k ~ k^(2/d_S) => log(λ) = (2/d_S) * log(k)
                // slope = 2/d_S => d_S = 2/slope

                double d_S;
                if (slope > 0.1)
                {
                    d_S = 2.0 / slope;
                }
                else if (slope > 0.01)
                {
                    // Очень малый slope — высокоразмерная геометрия
                    d_S = 2.0 / slope;
                }
                else
                {
                    // Degenerate case — все eigenvalues примерно равны
                    // Это означает expander graph / mean-field topology
                    d_S = 4.0; // Default to 4D для RQ
                }

                return Math.Clamp(d_S, 1.0, 10.0);
            }
            catch
            {
                _lastSpectralSlope = 0;
                _lastSpectralMethod = "fallback";
                return 2.0;
            }
        }

        /// <summary>
        /// Random Walk Return Probability метод.
        /// P(t) ~ t^(-d_S/2) => d_S = -2 * d(ln P)/d(ln t)
        /// 
        /// Физика: random walker на d-мерном пространстве возвращается
        /// в начальную точку с вероятностью P(t) ~ t^{-d/2}.
        /// 
        /// Улучшения для RQ:
        /// - Weighted random walk (учитывает Weights[i,j])
        /// - Больше walkers для статистики
        /// - Правильный диапазон t (не слишком малый, не слишком большой)
        /// - Robust regression
        /// 
        /// BUG FIX: Изолированные узлы (degree=0) не должны считаться как returns,
        /// так как walker фактически не двигался. Это согласовано с GPU реализацией.
        /// </summary>
        private double ComputeSpectralDimensionRandomWalk(int t_max, int num_walkers)
        {
            _lastSpectralMethod = "RandomWalk";

            // Масштабируем количество walkers с размером графа
            int effectiveWalkers = Math.Max(num_walkers, N * 20);

            // Ограничиваем t_max — слишком большое t даёт шум
            int effectiveTmax = Math.Min(t_max, 100);

            var returnCounts = new int[effectiveTmax + 1];
            object lockObj = new();

            // Параллельный запуск walkers
            Parallel.For(0, effectiveWalkers, () => new Random(Guid.NewGuid().GetHashCode()),
                (w, state, localRng) =>
                {
                    int startNode = localRng.Next(N);
                    int currentNode = startNode;
                    var localCounts = new int[effectiveTmax + 1];

                    for (int t = 1; t <= effectiveTmax; t++)
                    {
                        // BUG FIX: Track whether walker actually moved
                        int prevNode = currentNode;
                        currentNode = RandomWalkStep(currentNode, localRng);
                        
                        // Only count return if walker actually moved this step
                        // Isolated nodes (staying in place) should NOT count as returns
                        bool didMove = (currentNode != prevNode);
                        if (didMove && currentNode == startNode)
                            localCounts[t]++;
                    }

                    lock (lockObj)
                    {
                        for (int t = 1; t <= effectiveTmax; t++)
                            returnCounts[t] += localCounts[t];
                    }

                    return localRng;
                },
                _ => { });

            // Вычисляем вероятности
            var returnProb = new double[effectiveTmax + 1];
            for (int t = 1; t <= effectiveTmax; t++)
                returnProb[t] = (double)returnCounts[t] / effectiveWalkers;

            // Robust linear regression с отбрасыванием нулей и выбросов
            // Используем диапазон t ∈ [10, t_max/2] — средние времена
            var logT = new List<double>();
            var logP = new List<double>();

            int t_start = Math.Max(5, effectiveTmax / 10);
            int t_end = effectiveTmax * 2 / 3;

            for (int t = t_start; t <= t_end; t++)
            {
                if (returnProb[t] > 1e-8)
                {
                    logT.Add(Math.Log(t));
                    logP.Add(Math.Log(returnProb[t]));
                }
            }

            if (logT.Count < 5)
            {
                _lastSpectralSlope = 0;
                // Мало возвратов — граф слишком большой или фрагментированный
                // Пробуем Laplacian метод
                return ComputeSpectralDimensionLaplacian();
            }

            // Linear regression
            double meanLogT = logT.Average();
            double meanLogP = logP.Average();

            double num = 0, den = 0;
            for (int i = 0; i < logT.Count; i++)
            {
                double dT = logT[i] - meanLogT;
                double dP = logP[i] - meanLogP;
                num += dT * dP;
                den += dT * dT;
            }

            if (Math.Abs(den) < 1e-10)
            {
                _lastSpectralSlope = 0;
                return 2.0;
            }

            double slope = num / den;
            _lastSpectralSlope = slope;

            // P(t) ~ t^{-d_S/2} => log(P) = (-d_S/2) * log(t)
            // slope = -d_S/2 => d_S = -2 * slope

            // Для правильной физики slope должен быть отрицательным
            // (вероятность возврата убывает со временем)
            if (slope >= 0)
            {
                // Walker застревает (компактный граф или ловушки)
                // Fallback к Laplacian
                return ComputeSpectralDimensionLaplacian();
            }

            double d_S = -2.0 * slope;

            return Math.Clamp(d_S, 1.0, 10.0);
        }

        private int RandomWalkStep(int node, Random rng)
        {
            double totalWeight = 0;
            foreach (int nb in Neighbors(node)) totalWeight += Weights[node, nb];
            if (totalWeight <= 0) return node;
            double r = rng.NextDouble() * totalWeight;
            double sum = 0;
            foreach (int nb in Neighbors(node))
            {
                sum += Weights[node, nb];
                if (r <= sum) return nb;
            }
            return node;
        }
        /// <summary>
        /// Вычисляет энтропию запутанности между кластером (Region A) и остальным миром (Region B).
        /// S = -Tr(rho_A * ln rho_A), где rho_A - редуцированная матрица плотности.
        /// </summary>
        public double ComputeEntanglementEntropy(List<int> regionA)
        {
            if (_waveMulti == null) return 0;

            // Упрощенная оценка для чистых состояний на графе:
            // Энтропия ~ Сумма квадратов амплитуд на границе разреза (Cut)
            // Это аппроксимация для графовых состояний (Graph States).

            double boundaryEntropy = 0.0;

            foreach (int i in regionA)
            {
                foreach (int j in Neighbors(i))
                {
                    // Если ребро пересекает границу (j не в Region A)
                    if (!regionA.Contains(j))
                    {
                        // Вклад в запутанность пропорционален весу связи и корреляции состояний
                        double w = Weights[i, j];
                        // S_link ~ -w * ln(w)
                        if (w > 0.001)
                            boundaryEntropy += -w * Math.Log(w);
                    }
                }
            }

            return boundaryEntropy;
        }

        /// <summary>
        /// Проверка Area Law: Возвращает отношение Энтропии к "Площади" границы.
        /// Если ratio ~ const при росте кластера, значит пространство голографично.
        /// </summary>
        public double CheckHolographicPrinciple(List<int> cluster)
        {
            double entropy = ComputeEntanglementEntropy(cluster);

            // "Площадь" в графе = Сумма весов разрезанных ребер
            double area = 0;
            foreach (int i in cluster)
            {
                foreach (int j in Neighbors(i))
                {
                    if (!cluster.Contains(j)) area += Weights[i, j];
                }
            }

            return area > 0 ? entropy / area : 0;
        }


        private void InitProperTime()
        {
            ProperTime = new double[N];
            for (int i = 0; i < N; i++) ProperTime[i] = 0.0;
        }
        private void InitEdgeData()
        {
            EdgeDelay = new double[N, N];
            EdgeDirection = new sbyte[N, N];
        }

        /// <summary>
        /// DEPRECATED: This method uses external coordinates which violates RQ-hypothesis.
        /// Use GetGraphDistanceWeighted() for physics calculations instead.
        /// </summary>
        [Obsolete("Use GetGraphDistanceWeighted() for physics, this method is for rendering only.")]
        public double GetPhysicalDistance(int a, int b)
        {
#pragma warning disable CS0618 // Suppress obsolete warning for Coordinates access
            if (Coordinates == null || a < 0 || b < 0 || a >= N || b >= N) return 0.0;
            var (x1, y1) = Coordinates[a]; var (x2, y2) = Coordinates[b];
#pragma warning restore CS0618
            double dx = x1 - x2; double dy = y1 - y2; return Math.Sqrt(dx * dx + dy * dy);
        }

        /// <summary>
        /// Compute graph distance using shortest path with metric d = -ln(w).
        /// This is the RQ-compliant distance measure based purely on graph topology.
        /// Implements checklist item 1.1: Background-independent distance.
        /// </summary>
        /// <param name="startNode">Starting node index</param>
        /// <param name="endNode">Ending node index</param>
        /// <returns>Topological distance based on weighted shortest path</returns>
        public double GetGraphDistanceWeighted(int startNode, int endNode)
        {
            // Delegate to existing Dijkstra implementation
            return ShortestPathDistance(startNode, endNode);
        }

        /// <summary>
        /// Updates target distances from edge weights.
        /// Called after weight updates to keep distance matrix synchronized.
        /// Uses formula: d_ij = -log(w_ij + ε)
        /// </summary>
        public void UpdateTargetDistancesFromWeights()
        {
            if (_targetDistance == null || _targetDistance.GetLength(0) != N || _targetDistance.GetLength(1) != N)
                _targetDistance = new double[N, N];
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    if (!Edges[i, j] || i == j) { _targetDistance[i, j] = 0.0; continue; }
                    double w = Weights[i, j]; _targetDistance[i, j] = w <= 0.0 ? 0.0 : -Math.Log(w + EpsDist);
                }
            }
        }
        private void SetupEdgeMetric(int i, int j)
        {
            double d = GetPhysicalDistance(i, j);
            double avgDeg = 0.0; for (int k = 0; k < N; k++) avgDeg += _degree[k];
            avgDeg = N > 0 ? avgDeg / N : 1.0;
            double maxSignalSpeed = 1.0 / Math.Max(1.0, avgDeg);
            double delay = d / maxSignalSpeed;
            EdgeDelay[i, j] = delay; EdgeDelay[j, i] = delay;
            EdgeDirection[i, j] = 1; EdgeDirection[j, i] = -1;
        }
        private void UpdateEdgeDelaysFromDistances()
        {
            if (EdgeDelay == null || EdgeDelay.GetLength(0) != N || EdgeDelay.GetLength(1) != N)
                EdgeDelay = new double[N, N];
            if (_targetDistance == null || _targetDistance.GetLength(0) != N || _targetDistance.GetLength(1) != N)
                UpdateTargetDistancesFromWeights();
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    if (!Edges[i, j] || i == j) { EdgeDelay[i, j] = 0.0; continue; }
                    EdgeDelay[i, j] = _targetDistance[i, j];
                }
            }
        }

        // Remove duplicate AddEdge here (defined in CoreHelpers partial) – keep only one implementation
        // public void AddEdge(int i, int j) { }


        public void RelaxCoordinatesFromCorrelation(double eta)
        {
            if (Coordinates == null || Coordinates.Length != N) return; if (eta <= 0) return;
            UpdateTargetDistancesFromWeights();
            var dX = new double[N]; var dY = new double[N];
            for (int i = 0; i < N; i++)
            {
                double fx = 0.0, fy = 0.0; double xi = Coordinates[i].X; double yi = Coordinates[i].Y;
                for (int j = 0; j < N; j++)
                {
                    if (!Edges[i, j] || i == j) continue;
                    double dx = Coordinates[j].X - xi; double dy = Coordinates[j].Y - yi; double r = Math.Sqrt(dx * dx + dy * dy) + 1e-9;
                    double target = _targetDistance[i, j]; if (target <= 0.0) continue;
                    double diff = r - target; double k = diff / r; fx += k * dx; fy += k * dy;
                }
                dX[i] = -eta * fx; dY[i] = -eta * fy;
            }
            for (int i = 0; i < N; i++) Coordinates[i] = (Coordinates[i].X + dX[i], Coordinates[i].Y + dY[i]);
        }

        public double ComputeEffectiveLightSpeed()
        {
            if (EdgeDelay == null || Coordinates == null) return 0.0;
            double sumDist = 0.0, sumDelay = 0.0; int edgeCount = 0;
            for (int i = 0; i < N; i++)
            {
                for (int j = i + 1; j < N; j++)
                {
                    if (!Edges[i, j]) continue;
                    double d = GetPhysicalDistance(i, j); double delay = EdgeDelay[i, j]; if (delay <= 0) continue;
                    sumDist += d; sumDelay += delay; edgeCount++;
                }
            }
            if (edgeCount == 0 || sumDelay <= 0) return 0.0;
            double cEff = (sumDist / edgeCount) / (sumDelay / edgeCount);
            return double.IsFinite(cEff) ? cEff : 0.0;
        }

        // ComputeTotalEnergy with string field energy contribution per checklist
        public double ComputeTotalEnergy()
        {
            double energy = 0.0; bool useLocal = _targetDegreePerNode != null && _targetDegreePerNode.Length == N;
            for (int i = 0; i < N; i++)
            {
                int target = useLocal ? _targetDegreePerNode[i] : _targetDegree;
                double diff = _degree[i] - target;
                energy += diff * diff;
            }
            for (int i = 0; i < N; i++)
            {
                foreach (int nb in Neighbors(i)) if (i < nb) energy -= Weights[i, nb] * Weights[i, nb];
            }
            if (_stringEnergy != null && _stringEnergy.GetLength(0) == N)
            {
                double stringE = 0.0;
                for (int i = 0; i < N; i++) for (int j = i + 1; j < N; j++) stringE += _stringEnergy[i, j];
                energy += stringE;
            }
            return energy;
        }

        // Per-node correlation mass distribution (used by gravity wrapper)
        public double[] ComputePerNodeCorrelationMass(double wThreshold = PhysicsConstants.DefaultHeavyClusterThreshold, int minSize = HeavyClusterMinSize)
        {
            var nodeCorrMass = new double[N]; var clusters = GetStrongCorrelationClusters(wThreshold);
            foreach (var cluster in clusters)
            {
                int size = cluster.Count; if (size < minSize) continue;
                double mass = 0.0;
                for (int a = 0; a < cluster.Count; a++)
                {
                    int v = cluster[a];
                    for (int b = a + 1; b < cluster.Count; b++)
                    {
                        int u = cluster[b]; if (!Edges[v, u]) continue;
                        double w = Weights[v, u]; if (w <= wThreshold) continue; mass += (w - wThreshold);
                    }
                }
                if (mass <= 0.0) continue;
                double perNode = mass / size; foreach (int v in cluster) nodeCorrMass[v] += perNode;
            }
            return nodeCorrMass;
        }

        private void InitClocks(double fraction = 0.05)
        {
            _clockNodes.Clear();
            var heavy = GetStrongCorrelationClusters(AdaptiveHeavyThreshold);
            var largest = heavy.OrderByDescending(c => c.Count).FirstOrDefault();
            if (largest != null && largest.Count > 0)
            {
                foreach (int idx in largest)
                {
                    if (PhysicsProperties != null && PhysicsProperties.Length == N)
                        PhysicsProperties[idx].IsClock = true;
                    _clockNodes.Add(idx);
                }
                return;
            }
            int count = Math.Max(1, (int)(N * fraction));
            for (int i = 0; i < count; i++)
            {
                int idx = _rng.Next(N);
                if (PhysicsProperties != null && PhysicsProperties.Length == N)
                    PhysicsProperties[idx].IsClock = true;
                _clockNodes.Add(idx);
            }
        }

        // Measurement back-action helper referenced in Measurement partial
        private void ApplyMeasurementBackAction(List<int> system, List<int> apparatus)
        {
            if (system == null || apparatus == null) return;
            foreach (int i in system)
            {
                foreach (int j in system)
                {
                    if (i == j) continue;
                    if (Edges[i, j])
                    {
                        double w = Weights[i, j] * 0.8; Weights[i, j] = w; Weights[j, i] = w;
                    }
                }
                if (_nodeEnergy != null && i < _nodeEnergy.Length) _nodeEnergy[i] += 1.0;
            }
            foreach (int i in apparatus)
            {
                foreach (int j in apparatus)
                {
                    if (i == j) continue;
                    if (Edges[i, j])
                    {
                        double w = Weights[i, j] * 0.8; Weights[i, j] = w; Weights[j, i] = w;
                    }
                }
                if (_nodeEnergy != null && i < _nodeEnergy.Length) _nodeEnergy[i] += 1.0;
            }
            foreach (int i in system)
            {
                for (int j = 0; j < N; j++)
                {
                    if (system.Contains(j) || apparatus.Contains(j)) continue;
                    if (Edges[i, j]) { double w = Weights[i, j] * 0.9; Weights[i, j] = w; Weights[j, i] = w; }
                }
            }
            foreach (int i in apparatus)
            {
                for (int j = 0; j < N; j++)
                {
                    if (system.Contains(j) || apparatus.Contains(j)) continue;
                    if (Edges[i, j]) { double w = Weights[i, j] * 0.9; Weights[i, j] = w; Weights[j, i] = w; }
                }
            }
            UpdateEdgeDelaysFromDistances();
        }

        // Local curvature helper referenced by Updates partial
        public double GetLocalCurvature(int node)
        {
            double totalDegree = 0; for (int i = 0; i < N; i++) totalDegree += _degree[i];
            double avgDeg = N > 0 ? totalDegree / N : 0.0; double localDeg = _degree[node];
            return avgDeg == 0.0 ? 0.0 : (localDeg - avgDeg) / avgDeg;
        }
        public double GetLocalCorrelationDensity(int node)
        {
            double sumWeights = 0.0; foreach (var nb in Neighbors(node)) sumWeights += Math.Abs(Weights[node, nb]);
            double effEdges = sumWeights; double cellArea = 1.0; return cellArea > 0 ? effEdges / cellArea : 0.0;
        }

        // Missing private fields referenced in InitClocks
        private List<int> _clockNodes = new();

        // Heavy mass / correlation helpers required by CoreHelpers
        private void EnforcePlanckCorrelationLimit()
        {
            for (int i = 0; i < N; i++)
            {
                double local = GetLocalCorrelationDensity(i);
                if (local <= 1.0) continue;
                foreach (var nb in Neighbors(i))
                {
                    double w = Weights[i, nb];
                    double newW = w * 0.9;
                    Weights[i, nb] = newW; Weights[nb, i] = newW;
                }
            }
        }

        // Heavy mass delta overloads required by CoreHelpers Step
        public void ApplyHeavyMassDelta(double deltaMass)
        {
            if (N <= 0 || _nodeEnergy == null) return;
            double perNode = deltaMass / N;
            for (int i = 0; i < _nodeEnergy.Length; i++) _nodeEnergy[i] -= perNode;
        }
        public void ApplyHeavyMassDelta(double deltaMass, IEnumerable<int> cluster)
        {
            if (cluster == null) { ApplyHeavyMassDelta(deltaMass); return; }
            var list = cluster.Where(i => i >= 0 && i < N).ToList();
            if (list.Count == 0) { ApplyHeavyMassDelta(deltaMass); return; }
            double share = Math.Abs(deltaMass) / list.Count;
            foreach (int v in list)
            {
                foreach (int nb in Neighbors(v))
                {
                    _nodeEnergy[nb] -= share;
                    if (!Edges[v, nb] && _nodeEnergy[nb] < EnergySewThreshold)
                    {
                        AddEdge(v, nb); Weights[v, nb] = 0.1; Weights[nb, v] = 0.1;
                    }
                }
            }
        }

        /// <summary>
        /// RQ-HYPOTHESIS: Initializes random coordinates FOR VISUALIZATION PURPOSES ONLY.
        /// 
        /// WARNING: These coordinates are external/background coordinates that violate
        /// the relational principle of RQ-Hypothesis. They MUST NOT be used in:
        /// - Distance calculations (use ShortestPathDistance or GetGraphDistanceWeighted)
        /// - Mass calculations (use spectral coordinates from Laplacian eigenvectors)
        /// - Physics evolution (all physics must be topology-based)
        /// 
        /// These coordinates are used ONLY by the UI drawing system for visualization.
        /// The actual "positions" in RQ-Hypothesis emerge from correlation structure,
        /// computed via spectral geometry (SpectralCoordinates property).
        /// </summary>
        [System.ComponentModel.EditorBrowsable(System.ComponentModel.EditorBrowsableState.Advanced)]
        public void InitCoordinatesRandom(double range = 1.0)
        {
            if (N <= 0) return;
#pragma warning disable CS0618 // Suppress obsolete warning - this method IS the UI initializer
            if (Coordinates == null || Coordinates.Length != N) Coordinates = new (double X, double Y)[N];
            for (int i = 0; i < N; i++)
            {
                double x = (_rng.NextDouble() * 2.0 - 1.0) * range;
                double y = (_rng.NextDouble() * 2.0 - 1.0) * range;
                Coordinates[i] = (x, y);
            }
#pragma warning restore CS0618
        }

        // Einstein dynamics stub referenced by examples (simplified)
        public void ApplyEinsteinDynamics()
        {
            if (_nodeEnergy == null || _nodeEnergy.Length != N)
            {
                _nodeEnergy = new double[N];
            }
            for (int i = 0; i < N; i++)
            {
                double baseE = 0.0;
                foreach (var nb in Neighbors(i)) baseE += Weights[i, nb];
                _nodeEnergy[i] = baseE;
            }
        }

        public void ComputeFractalLevels()
        {
            if (N <= 0) { FractalLevel = Array.Empty<int>(); return; }
            FractalLevel = new int[N];
            for (int i = 0; i < N; i++)
            {
                int deg = _degree != null && i < _degree.Length ? _degree[i] : 0;
                FractalLevel[i] = (int)Math.Log(Math.Max(1, deg), 2);
            }
        }

        public void InitFractalTopology(int levels, int branchFactor)
        {
            if (N <= 0) return;
            int[] nodeIndices = Enumerable.Range(0, N).ToArray();
            // initialise domain & structural mass if needed
            if (Domain == null || Domain.Length != N)
            {
                Domain = new int[N];
                for (int i = 0; i < N; i++) Domain[i] = _rng.Next(levels * 4);
            }
            // RQ-FIX: Removed StructuralMass initialization
            for (int lvl = 1; lvl <= levels; lvl++)
            {
                int groupCount = (int)Math.Pow(branchFactor, lvl);
                int size = Math.Max(1, N / groupCount);
                for (int g = 0; g < groupCount; g++)
                {
                    var group = nodeIndices.Skip(g * size).Take(size).ToList();
                    foreach (int i in group)
                        foreach (int j in group)
                        {
                            if (i < j && _rng.NextDouble() < 0.5)
                            {
                                AddEdge(i, j);
                                double w = 0.15 + 0.45 * _rng.NextDouble();
                                Weights[i, j] = w; Weights[j, i] = w;
                            }
                        }
                }
                for (int g1 = 0; g1 < groupCount; g1++)
                {
                    for (int g2 = g1 + 1; g2 < groupCount; g2++)
                    {
                        if (_rng.NextDouble() < 0.1)
                        {
                            int a = g1 * size;
                            int b = g2 * size;
                            if (a < N && b < N)
                            {
                                AddEdge(nodeIndices[a], nodeIndices[b]);
                                double w = 0.15 + 0.25 * _rng.NextDouble();
                                Weights[nodeIndices[a], nodeIndices[b]] = w;
                                Weights[nodeIndices[b], nodeIndices[a]] = w;
                            }
                        }
                    }
                }
            }
            UpdateEdgeDelaysFromDistances();
            // initialize local potential after topology defined
            if (LocalPotential == null || LocalPotential.Length != N)
            {
                LocalPotential = new double[N];
                for (int i = 0; i < N; i++) LocalPotential[i] = _rng.NextDouble() * 0.1;
            }
        }
        public enum LocalPhase { Cold, Warm, Hot }
        public LocalPhase[] NodePhase { get; set; }
        public double GlobalExcitationRegulator { get; set; } = 1.0;
        public void UpdateStoredEnergyAndPaths()
        {
            if (LocalPotential == null || StoredEnergy == null) return;
            for (int i = 0; i < N; i++)
            {
                double lp = LocalPotential[i];
                if (lp > 0.7) StoredEnergy[i] += 0.02 * lp; else StoredEnergy[i] *= 0.98;
                if (StoredEnergy[i] < 0) StoredEnergy[i] = 0;
            }
            if (PathWeight != null)
            {
                for (int i = 0; i < N; i++)
                {
                    foreach (int nb in Neighbors(i))
                    {
                        double diff = Math.Abs(LocalPotential[i] - LocalPotential[nb]);
                        double w = Math.Exp(-diff);
                        PathWeight[i, nb] = w;
                        PathWeight[nb, i] = w;
                    }
                }
            }
        }

        public enum GraphPhase { Quiet, MetaStable, Active }
        public GraphPhase CurrentPhase { get; set; } = GraphPhase.Quiet;

        /// <summary>
        /// Global energy ledger for strict conservation.
        /// </summary>
        private readonly EnergyLedger _ledger = new EnergyLedger();
        public EnergyLedger Ledger => _ledger;
    }
}
