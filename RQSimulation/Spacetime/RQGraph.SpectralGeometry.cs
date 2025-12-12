using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;

namespace RQSimulation
{
    /// <summary>
    /// Spectral embedding for emergent geometry.
    /// Coordinates emerge from graph Laplacian eigenvectors rather than being fixed.
    /// </summary>
    public partial class RQGraph
    {
        // Spectral coordinates (emergent from Laplacian)
        private double[] _spectralX;
        private double[] _spectralY;
        private double[] _spectralZ;

        // Graph Laplacian eigenvalues (for spectral dimension analysis)
        private double[] _laplacianEigenvalues;

        // Cached spectral dimension from last computation
        private double _cachedSpectralDimension = 2.0;
        private int _spectralDimCacheStep = -1;

        // Сглаженное спектральное измерение
        private double _smoothedSpectralDimension = double.NaN;
        private const double SpectralDimensionEMA = 0.3;  // Smoothing factor

        /// <summary>
        /// Numerical tolerance for eigenvalue computations and matrix operations.
        /// Used to avoid division by near-zero values in spectral decomposition.
        /// </summary>
        private const double SpectralTolerance = 1e-10;

        /// <summary>
        /// Number of power iterations for eigenvalue computation.
        /// Increased from 100 to 300 for better convergence on dense graphs.
        /// </summary>
        private const int PowerIterations = 300;

        /// <summary>
        /// Convergence threshold for power iteration.
        /// Iteration stops when change is below this value.
        /// </summary>
        private const double ConvergenceThreshold = 1e-8;

        /// <summary>
        /// Spectral X coordinates derived from first non-trivial eigenvector
        /// </summary>
        public double[] SpectralX => _spectralX ?? Array.Empty<double>();

        /// <summary>
        /// Spectral Y coordinates derived from second non-trivial eigenvector  
        /// </summary>
        public double[] SpectralY => _spectralY ?? Array.Empty<double>();

        /// <summary>
        /// Spectral Z coordinates derived from third non-trivial eigenvector
        /// </summary>
        public double[] SpectralZ => _spectralZ ?? Array.Empty<double>();

        /// <summary>
        /// Compute the graph Laplacian matrix L = D - W
        /// where D is the degree matrix and W is the correlation weight matrix.
        /// </summary>
        public double[,] ComputeGraphLaplacian()
        {
            var L = new double[N, N];

            for (int i = 0; i < N; i++)
            {
                double degree = 0.0;
                for (int j = 0; j < N; j++)
                {
                    if (i != j && Edges[i, j])
                    {
                        double w = Math.Max(0.0, Weights[i, j]);
                        L[i, j] = -w;
                        degree += w;
                    }
                }
                L[i, i] = degree;
            }

            return L;
        }

        /// <summary>
        /// Compute the normalized graph Laplacian L_norm = D^(-1/2) L D^(-1/2)
        /// Better for spectral embedding as it handles varying node degrees.
        /// </summary>
        public double[,] ComputeNormalizedLaplacian()
        {
            var L = ComputeGraphLaplacian();
            var Lnorm = new double[N, N];

            // Compute D^(-1/2) diagonal entries
            var dInvSqrt = new double[N];
            for (int i = 0; i < N; i++)
            {
                double deg = L[i, i];
                dInvSqrt[i] = deg > SpectralTolerance ? 1.0 / Math.Sqrt(deg) : 0.0;
            }

            // L_norm = D^(-1/2) L D^(-1/2)
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    Lnorm[i, j] = L[i, j] * dInvSqrt[i] * dInvSqrt[j];
                }
            }

            return Lnorm;
        }

        /// <summary>
        /// Compute eigenvalues and eigenvectors of the graph Laplacian using improved power iteration.
        /// Returns the k smallest non-zero eigenvalues and their eigenvectors.
        /// 
        /// RQ-FIX: Increased iterations (300), added convergence check, better orthogonalization.
        /// </summary>
        public (double[] eigenvalues, double[,] eigenvectors) ComputeLaplacianSpectrum(int k = 4)
        {
            k = Math.Min(k, N);
            var L = ComputeNormalizedLaplacian();

            var eigenvalues = new double[k];
            var eigenvectors = new double[N, k];

            // Store found eigenvectors for orthogonalization
            var foundVectors = new List<double[]>();

            for (int ev = 0; ev < k; ev++)
            {
                // Power iteration with deflation
                var v = new double[N];
                // Use deterministic initialization for reproducibility
                for (int i = 0; i < N; i++)
                    v[i] = Math.Sin((i + 1) * (ev + 1) * 0.1) + 0.5;

                // Normalize
                double norm = Math.Sqrt(v.Sum(x => x * x));
                if (norm > SpectralTolerance)
                    for (int i = 0; i < N; i++) v[i] /= norm;

                // Inverse power iteration to find smallest eigenvalues
                // We use (σI - L) to convert to largest eigenvalue problem
                var shiftedL = new double[N, N];
                double shift = 2.0; // Larger than max eigenvalue of normalized Laplacian
                for (int i = 0; i < N; i++)
                {
                    for (int j = 0; j < N; j++)
                    {
                        shiftedL[i, j] = (i == j ? shift : 0.0) - L[i, j];
                    }
                }

                double prevNorm = 0.0;

                // Power iteration on shifted matrix with convergence check
                for (int iter = 0; iter < PowerIterations; iter++)
                {
                    var vNew = new double[N];

                    // Matrix-vector multiplication
                    for (int i = 0; i < N; i++)
                    {
                        double sum = 0.0;
                        for (int j = 0; j < N; j++)
                            sum += shiftedL[i, j] * v[j];
                        vNew[i] = sum;
                    }

                    // Orthogonalize against ALL previously found eigenvectors (Gram-Schmidt)
                    for (int p = 0; p < foundVectors.Count; p++)
                    {
                        var prev = foundVectors[p];
                        double dot = 0.0;
                        for (int i = 0; i < N; i++) dot += vNew[i] * prev[i];
                        for (int i = 0; i < N; i++) vNew[i] -= dot * prev[i];
                    }

                    // Normalize
                    norm = Math.Sqrt(vNew.Sum(x => x * x));
                    if (norm > SpectralTolerance)
                    {
                        for (int i = 0; i < N; i++) v[i] = vNew[i] / norm;

                        // Convergence check: if norm change is small, we've converged
                        if (iter > 10 && Math.Abs(norm - prevNorm) < ConvergenceThreshold * norm)
                            break;

                        prevNorm = norm;
                    }
                    else
                    {
                        break;
                    }
                }

                // Compute Rayleigh quotient for eigenvalue
                double numerator = 0.0, denominator = 0.0;
                var Lv = new double[N];
                for (int i = 0; i < N; i++)
                {
                    double sum = 0.0;
                    for (int j = 0; j < N; j++)
                        sum += L[i, j] * v[j];
                    Lv[i] = sum;
                }
                for (int i = 0; i < N; i++)
                {
                    numerator += v[i] * Lv[i];
                    denominator += v[i] * v[i];
                }

                eigenvalues[ev] = denominator > SpectralTolerance ? numerator / denominator : 0.0;

                // Store eigenvector
                var evCopy = new double[N];
                for (int i = 0; i < N; i++)
                {
                    eigenvectors[i, ev] = v[i];
                    evCopy[i] = v[i];
                }
                foundVectors.Add(evCopy);
            }

            // Sort by eigenvalue
            var indices = Enumerable.Range(0, k).OrderBy(i => eigenvalues[i]).ToArray();
            var sortedEigenvalues = indices.Select(i => eigenvalues[i]).ToArray();
            var sortedEigenvectors = new double[N, k];
            for (int col = 0; col < k; col++)
            {
                for (int row = 0; row < N; row++)
                {
                    sortedEigenvectors[row, col] = eigenvectors[row, indices[col]];
                }
            }

            _laplacianEigenvalues = sortedEigenvalues;
            return (sortedEigenvalues, sortedEigenvectors);
        }

        /// <summary>
        /// Update spectral coordinates from graph Laplacian eigenvectors.
        /// This implements emergent geometry from correlation structure.
        /// </summary>
        public void UpdateSpectralCoordinates()
        {
            if (N <= 3) return;

            var (eigenvalues, eigenvectors) = ComputeLaplacianSpectrum(4);

            // Initialize arrays
            if (_spectralX == null || _spectralX.Length != N) _spectralX = new double[N];
            if (_spectralY == null || _spectralY.Length != N) _spectralY = new double[N];
            if (_spectralZ == null || _spectralZ.Length != N) _spectralZ = new double[N];

            // Use eigenvectors 1, 2, 3 (skip first trivial one for connected graph)
            // Scale by inverse eigenvalue for meaningful distances
            int xIdx = eigenvalues[1] > SpectralTolerance ? 1 : 0;
            int yIdx = Math.Min(xIdx + 1, eigenvalues.Length - 1);
            int zIdx = Math.Min(yIdx + 1, eigenvalues.Length - 1);

            double scaleX = eigenvalues[xIdx] > SpectralTolerance ? 1.0 / Math.Sqrt(eigenvalues[xIdx]) : 1.0;
            double scaleY = eigenvalues[yIdx] > SpectralTolerance ? 1.0 / Math.Sqrt(eigenvalues[yIdx]) : 1.0;
            double scaleZ = eigenvalues[zIdx] > SpectralTolerance ? 1.0 / Math.Sqrt(eigenvalues[zIdx]) : 1.0;

            for (int i = 0; i < N; i++)
            {
                _spectralX[i] = eigenvectors[i, xIdx] * scaleX;
                _spectralY[i] = eigenvectors[i, yIdx] * scaleY;
                _spectralZ[i] = eigenvectors[i, zIdx] * scaleZ;
            }
        }

        /// <summary>
        /// Get emergent distance between two nodes using spectral coordinates.
        /// This replaces fixed coordinate-based distance.
        /// </summary>
        public double GetSpectralDistance(int i, int j)
        {
            if (_spectralX == null || _spectralY == null || i < 0 || j < 0 || i >= N || j >= N)
                return GetGraphDistance(i, j); // Fallback

            double dx = _spectralX[i] - _spectralX[j];
            double dy = _spectralY[i] - _spectralY[j];
            double dz = (_spectralZ != null && i < _spectralZ.Length && j < _spectralZ.Length)
                ? _spectralZ[i] - _spectralZ[j] : 0.0;

            return Math.Sqrt(dx * dx + dy * dy + dz * dz);
        }

        /// <summary>
        /// Get graph-theoretic distance (shortest path in terms of hops).
        /// Used as fallback and for pure relational geometry.
        /// </summary>
        public int GetGraphDistance(int i, int j)
        {
            if (i == j) return 0;
            if (i < 0 || j < 0 || i >= N || j >= N) return int.MaxValue;

            // BFS for shortest path
            var dist = new int[N];
            for (int k = 0; k < N; k++) dist[k] = -1;

            var queue = new Queue<int>();
            queue.Enqueue(i);
            dist[i] = 0;

            while (queue.Count > 0)
            {
                int v = queue.Dequeue();
                if (v == j) return dist[j];

                foreach (int u in Neighbors(v))
                {
                    if (dist[u] == -1)
                    {
                        dist[u] = dist[v] + 1;
                        queue.Enqueue(u);
                    }
                }
            }

            return dist[j] >= 0 ? dist[j] : int.MaxValue;
        }

        /// <summary>
        /// Get effective distance using correlation weights: d = 1/w
        /// Strong correlations mean short effective distance.
        /// </summary>
        public double GetCorrelationDistance(int i, int j)
        {
            if (i == j) return 0.0;
            if (!Edges[i, j]) return double.PositiveInfinity;

            double w = Weights[i, j];
            return w > SpectralTolerance ? 1.0 / w : double.PositiveInfinity;
        }

        /// <summary>
        /// Sync visualization coordinates from spectral embedding.
        /// Keeps Coordinates array for Form_Main drawing while deriving from emergent geometry.
        /// </summary>
        public void SyncCoordinatesFromSpectral()
        {
            if (_spectralX == null || _spectralY == null) return;
            if (Coordinates == null || Coordinates.Length != N)
                Coordinates = new (double X, double Y)[N];

            // Find bounds for normalization
            double minX = _spectralX.Min();
            double maxX = _spectralX.Max();
            double minY = _spectralY.Min();
            double maxY = _spectralY.Max();

            double rangeX = maxX - minX;
            double rangeY = maxY - minY;
            double maxRange = Math.Max(rangeX, rangeY);

            if (maxRange < SpectralTolerance) maxRange = 1.0;

            for (int i = 0; i < N; i++)
            {
                // Normalize to [-1, 1] range
                double x = rangeX > SpectralTolerance ? 2.0 * (_spectralX[i] - minX) / rangeX - 1.0 : 0.0;
                double y = rangeY > SpectralTolerance ? 2.0 * (_spectralY[i] - minY) / rangeY - 1.0 : 0.0;

                // Apply uniform scaling to preserve aspect ratio
                Coordinates[i] = (x, y);
            }
        }

        /// <summary>
        /// Get cached spectral dimension without recomputation.
        /// Useful for UI updates where fresh computation is too expensive.
        /// </summary>
        public double GetCachedSpectralDimension() => _cachedSpectralDimension;

        /// <summary>
        /// Get the Laplacian eigenvalues array for external analysis.
        /// </summary>
        public double[] ComputeLaplacianEigenvalues()
        {
            if (_laplacianEigenvalues == null || _laplacianEigenvalues.Length == 0)
            {
                ComputeLaplacianSpectrum(Math.Min(10, N));
            }
            return _laplacianEigenvalues ?? Array.Empty<double>();
        }
        
        // ================================================================
        // CHECKLIST ITEM 6: SPECTRAL MASS FROM LAPLACIAN EIGENVALUES
        // ================================================================
        
        /// <summary>
        /// Compute mass of a cluster from spectral structure (λ₂ eigenvalue).
        /// 
        /// CHECKLIST ITEM 6: Topological mass from graph Laplacian.
        /// 
        /// In the RQ-hypothesis, mass emerges from the graph's spectral properties.
        /// The Fiedler eigenvalue λ₂ (smallest non-zero eigenvalue of Laplacian)
        /// characterizes the "algebraic connectivity" of a subgraph.
        /// 
        /// Physical interpretation:
        /// - λ₂ ≈ 0 → loosely connected (easy to separate) → low inertia → low mass
        /// - λ₂ large → tightly connected (hard to separate) → high inertia → high mass
        /// 
        /// Mass formula: m ∝ MassScale / λ₂
        /// 
        /// This relates to:
        /// - Higgs mechanism: mass from broken symmetry (here: broken translational symmetry)
        /// - Gravity: massive objects curve spacetime (here: mass affects graph geometry)
        /// </summary>
        /// <param name="clusterNodes">List of node indices in the cluster</param>
        /// <param name="massScale">Scale factor for mass (default from PhysicsConstants)</param>
        /// <returns>Computed spectral mass</returns>
        public double ComputeSpectralMass(List<int> clusterNodes, double massScale = 1.0)
        {
            if (clusterNodes == null || clusterNodes.Count < 2)
                return 0.0;
                
            // Build subgraph Laplacian for the cluster
            int n = clusterNodes.Count;
            var L = BuildSubgraphLaplacian(clusterNodes);
            
            // Find λ₂ using power iteration on (λ_max * I - L)
            double lambda2 = ComputeSmallestNonzeroEigenvalue(L, n);
            
            // Avoid division by zero
            if (lambda2 < SpectralTolerance)
            {
                // Disconnected cluster: effectively infinite mass (or interpret as separate masses)
                return massScale * n; // Use node count as proxy for disconnected case
            }
            
            // Mass inversely proportional to algebraic connectivity
            // More connected = smaller λ₂ = higher mass? No, opposite:
            // Higher λ₂ = harder to separate = more "rigid" = higher effective mass
            // 
            // Actually, for physical mass interpretation:
            // m ∝ 1/λ₂ makes sense: tighter binding (larger λ₂) → smaller wavelength → higher mass (de Broglie)
            // But cluster stability suggests: m ∝ n/λ₂ (size matters too)
            
            double mass = massScale * n / lambda2;
            
            return mass;
        }
        
        /// <summary>
        /// Build the Laplacian matrix for a subgraph defined by a node list.
        /// </summary>
        /// <param name="nodes">List of node indices in the subgraph</param>
        /// <returns>n×n Laplacian matrix where n = nodes.Count</returns>
        public double[,] BuildSubgraphLaplacian(List<int> nodes)
        {
            int n = nodes.Count;
            var L = new double[n, n];
            
            // Map global indices to local indices
            var localIndex = new Dictionary<int, int>();
            for (int i = 0; i < n; i++)
            {
                localIndex[nodes[i]] = i;
            }
            
            for (int i = 0; i < n; i++)
            {
                int globalI = nodes[i];
                double degree = 0.0;
                
                for (int j = 0; j < n; j++)
                {
                    if (i == j) continue;
                    
                    int globalJ = nodes[j];
                    
                    if (Edges[globalI, globalJ])
                    {
                        double w = Math.Max(0.0, Weights[globalI, globalJ]);
                        L[i, j] = -w;
                        degree += w;
                    }
                }
                
                L[i, i] = degree;
            }
            
            return L;
        }
        
        /// <summary>
        /// Compute the smallest non-zero eigenvalue (λ₂, Fiedler value) of a Laplacian matrix.
        /// 
        /// Uses inverse power iteration to find smallest eigenvalue, then deflation.
        /// </summary>
        /// <param name="L">Laplacian matrix</param>
        /// <param name="n">Matrix dimension</param>
        /// <returns>λ₂ (Fiedler eigenvalue)</returns>
        public double ComputeSmallestNonzeroEigenvalue(double[,] L, int n)
        {
            if (n < 2) return 0.0;
            
            // First eigenvalue of Laplacian is always 0 (constant eigenvector)
            // We need λ₂, the first non-trivial eigenvalue
            
            // Shift matrix: L_shifted = λ_max * I - L
            // This converts smallest eigenvalue problem to largest eigenvalue problem
            double lambdaMax = EstimateLargestEigenvalue(L, n);
            
            var Lshifted = new double[n, n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    Lshifted[i, j] = (i == j ? lambdaMax : 0.0) - L[i, j];
                }
            }
            
            // Trivial eigenvector (constant): v₀ = [1,1,...,1]/√n
            var v0 = new double[n];
            double norm0 = 1.0 / Math.Sqrt(n);
            for (int i = 0; i < n; i++) v0[i] = norm0;
            
            // Power iteration to find largest eigenvector of L_shifted
            var v = new double[n];
            var random = new Random(42);
            for (int i = 0; i < n; i++) v[i] = random.NextDouble() - 0.5;
            
            // Orthogonalize against trivial eigenvector
            double dot0 = 0.0;
            for (int i = 0; i < n; i++) dot0 += v[i] * v0[i];
            for (int i = 0; i < n; i++) v[i] -= dot0 * v0[i];
            
            // Normalize
            double norm = Math.Sqrt(v.Sum(x => x * x));
            if (norm > SpectralTolerance)
                for (int i = 0; i < n; i++) v[i] /= norm;
            
            // Power iteration
            for (int iter = 0; iter < PowerIterations; iter++)
            {
                var vNew = new double[n];
                
                // Matrix-vector multiply
                for (int i = 0; i < n; i++)
                {
                    double sum = 0.0;
                    for (int j = 0; j < n; j++)
                        sum += Lshifted[i, j] * v[j];
                    vNew[i] = sum;
                }
                
                // Orthogonalize against trivial eigenvector
                dot0 = 0.0;
                for (int i = 0; i < n; i++) dot0 += vNew[i] * v0[i];
                for (int i = 0; i < n; i++) vNew[i] -= dot0 * v0[i];
                
                // Normalize
                norm = Math.Sqrt(vNew.Sum(x => x * x));
                if (norm > SpectralTolerance)
                {
                    for (int i = 0; i < n; i++) v[i] = vNew[i] / norm;
                }
                else
                {
                    break;
                }
            }
            
            // Compute Rayleigh quotient for λ₂
            var Lv = new double[n];
            for (int i = 0; i < n; i++)
            {
                double sum = 0.0;
                for (int j = 0; j < n; j++)
                    sum += L[i, j] * v[j];
                Lv[i] = sum;
            }
            
            double numerator = 0.0, denominator = 0.0;
            for (int i = 0; i < n; i++)
            {
                numerator += v[i] * Lv[i];
                denominator += v[i] * v[i];
            }
            
            return denominator > SpectralTolerance ? numerator / denominator : 0.0;
        }
        
        /// <summary>
        /// Estimate the largest eigenvalue of a matrix using power iteration.
        /// Used for shifting in inverse iteration.
        /// </summary>
        private double EstimateLargestEigenvalue(double[,] M, int n)
        {
            var v = new double[n];
            for (int i = 0; i < n; i++) v[i] = 1.0 / Math.Sqrt(n);
            
            double lambda = 0.0;
            
            for (int iter = 0; iter < 50; iter++)
            {
                var vNew = new double[n];
                
                for (int i = 0; i < n; i++)
                {
                    double sum = 0.0;
                    for (int j = 0; j < n; j++)
                        sum += M[i, j] * v[j];
                    vNew[i] = sum;
                }
                
                double norm = Math.Sqrt(vNew.Sum(x => x * x));
                if (norm > SpectralTolerance)
                {
                    lambda = norm;
                    for (int i = 0; i < n; i++) v[i] = vNew[i] / norm;
                }
                else
                {
                    break;
                }
            }
            
            // Add small buffer to ensure shift is above all eigenvalues
            return lambda * 1.1 + 1.0;
        }
        
        /// <summary>
        /// Compute spectral mass for all heavy clusters in the graph.
        /// 
        /// Returns a dictionary mapping cluster index to computed mass.
        /// Updates the StructuralMass field for nodes in clusters.
        /// </summary>
        public Dictionary<int, double> ComputeAllClusterSpectralMasses(double threshold = -1)
        {
            if (threshold < 0)
                threshold = GetAdaptiveHeavyThreshold();
                
            var clusters = GetStrongCorrelationClusters(threshold);
            var masses = new Dictionary<int, double>();
            
            for (int c = 0; c < clusters.Count; c++)
            {
                var cluster = clusters[c];
                if (cluster.Count < PhysicsConstants.MinimumClusterSize)
                    continue;
                    
                double mass = ComputeSpectralMass(cluster);
                masses[c] = mass;
                
                // RQ-FIX: Removed StructuralMass update. 
                // Mass is now derived purely from correlations or spectral properties on demand.
            }
            
            return masses;
        }
        
        /// <summary>
        /// Get the spectral gap (difference between λ₂ and λ₁=0).
        /// 
        /// The spectral gap indicates:
        /// - How quickly random walks mix (larger gap = faster mixing)
        /// - Connectivity robustness (larger gap = more connected)
        /// - In RQ-hypothesis: relates to information propagation speed
        /// </summary>
        public double GetSpectralGap()
        {
            var eigenvalues = ComputeLaplacianEigenvalues();
            
            if (eigenvalues.Length < 2)
                return 0.0;
                
            // λ₁ should be ~0, λ₂ is the spectral gap
            // Find first eigenvalue > threshold (skip numerical zero)
            for (int i = 0; i < eigenvalues.Length; i++)
            {
                if (eigenvalues[i] > SpectralTolerance)
                    return eigenvalues[i];
            }
            
            return 0.0;
        }
    }
}
