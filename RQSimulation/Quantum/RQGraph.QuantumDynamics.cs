using System;
using System.Linq;
using System.Numerics;
using System.Threading.Tasks;

namespace RQSimulation
{
    public partial class RQGraph
    {
        private Complex[] _wavefunction;
        private Complex[] _waveMulti;
        private double[,] _edgePhase;
        private Complex[] _psiDot;

        // Dynamic quantum scales derived from graph invariants (avg degree)
        public double QuantumCoupling
        {
            get
            {
                var es = GetEdgeStats();
                double refDeg = es.avgDegree > 0 ? es.avgDegree : 1.0;
                return 1.0 / refDeg;
            }
        }
        public double QuantumDt
        {
            get
            {
                double refDeg = Math.Max(1.0, GetEdgeStats().avgDegree);
                return 0.5 / refDeg; // stability heuristic from Laplacian eigenestimate
            }
        }

        private int[] _edgeIndexRowStart;
        private int[] _edgeIndexCols;
        private int _edgeIndexCount;

        private void EnsureEdgeIndex()
        {
            if (_edgeIndexCols != null && _edgeIndexRowStart != null) return;
            var cols = new System.Collections.Generic.List<int>();
            _edgeIndexRowStart = new int[N];
            int ptr = 0;
            for (int i = 0; i < N; i++)
            {
                _edgeIndexRowStart[i] = ptr;
                foreach (int j in Neighbors(i))
                {
                    cols.Add(j);
                    ptr++;
                }
            }
            _edgeIndexCols = cols.ToArray();
            _edgeIndexCount = ptr;
        }

        private int EdgeIndex(int i, int j)
        {
            EnsureEdgeIndex();
            int start = _edgeIndexRowStart[i];
            for (int k = start; k < _edgeIndexCount && k < _edgeIndexCols.Length; k++)
            {
                if (_edgeIndexCols[k] == j) return k;
                if (i + 1 < _edgeIndexRowStart.Length && k + 1 == _edgeIndexRowStart[i + 1]) break;
            }
            return start; // fallback within row
        }

        public void InitQuantumWavefunction()
        {
            if (_edgePhase == null || _edgePhase.GetLength(0) != N)
                _edgePhase = new double[N, N];
            if (_wavefunction == null || _wavefunction.Length != N)
                _wavefunction = new Complex[N];

            int d = GaugeDimension;
            _waveMulti = new Complex[N * d];
            var rnd = _rng;
            for (int i = 0; i < N; i++)
            {
                for (int a = 0; a < d; a++)
                {
                    int idx = i * d + a;
                    _waveMulti[idx] = new Complex(rnd.NextDouble() * 1e-3, rnd.NextDouble() * 1e-3);
                }
            }

            EnsureEdgeIndex();
        }

        public double QuantumDiffusion { get; set; } = 0.2; // checklist diffusion coefficient

        public void UpdateQuantumState()
        {
            double cRq = ComputeEffectiveLightSpeed();
            double dt = ComputeRelationalDt();
            double dtRq = dt * cRq;
            double hbarRq = 1.0;
            StepSchrodingerUnitary(dtRq, hbarRq);
            if (_waveMulti != null)
            {
                int d = GaugeDimension;
                var rng = _rng;
                var next = (Complex[])_waveMulti.Clone();
                // phase drift accumulator per node (create if needed)
                if (_phaseDrift == null || _phaseDrift.Length != N) _phaseDrift = new double[N];
                for (int i = 0; i < N; i++)
                {
                    if (_correlationMass != null) _phaseDrift[i] += 0.001 * _correlationMass[i]; // checklist 4 phase drift
                }
                for (int i = 0; i < N; i++)
                {
                    for (int a = 0; a < d; a++)
                    {
                        int idx = i * d + a;
                        Complex sum = Complex.Zero;
                        int deg = 0;
                        foreach (int nb in Neighbors(i))
                        {
                            int jIdx = nb * d + a;
                            double corridorW = PathWeight != null ? PathWeight[i, nb] : 1.0; // item 9 corridors
                            sum += corridorW * _waveMulti[jIdx];
                            deg++;
                        }
                        if (deg > 0)
                        {
                            Complex lap = (sum / deg) - _waveMulti[idx];
                            next[idx] += QuantumDiffusion * lap;
                        }
                        double potLocal = LocalPotential != null ? LocalPotential[i] : 0.0;
                        double phaseShift = 0.05 * potLocal;
                        double ampLocal = _waveMulti[idx].Magnitude;
                        double nonlinearPhase = 0.1 * ampLocal * potLocal;
                        double totalPhase = phaseShift + nonlinearPhase + _phaseDrift[i];
                        next[idx] *= Complex.FromPolarCoordinates(1.0, totalPhase);
                    }
                }
                // Normalize then amplitude scaling by local potential
                double norm = 0.0; foreach (var z in next) { double m = z.Magnitude; norm += m * m; }
                if (norm > 0)
                {
                    norm = Math.Sqrt(norm);
                    for (int i = 0; i < N; i++)
                    {
                        for (int a = 0; a < d; a++)
                        {
                            int idx = i * d + a;
                            Complex z = next[idx] / norm;
                            double potLocal = LocalPotential != null ? LocalPotential[i] : 0.0;
                            z *= (1.0 + 0.05 * potLocal); // local peak amplification
                            if (IsInCondensedCluster(i)) z *= 1.1; // quantum signature boost
                            next[idx] = z;
                        }
                    }
                }
                // small random dephasing retained
                for (int i = 0; i < N; i++)
                {
                    for (int a = 0; a < d; a++)
                    {
                        int idx = i * d + a; double phaseNoise = 0.02 * rng.NextDouble();
                        next[idx] *= Complex.FromPolarCoordinates(1.0, phaseNoise);
                    }
                }
                _waveMulti = next;
            }
        }
        private double[] _phaseDrift;

        private void NormalizeWavefunction()
        {
            if (_waveMulti == null) return;
            double norm = 0.0;
            foreach (var z in _waveMulti) { double m = z.Magnitude; norm += m * m; }
            if (norm <= 0.0) return; double inv = 1.0 / Math.Sqrt(norm);
            for (int i = 0; i < _waveMulti.Length; i++) _waveMulti[i] *= inv;
        }

        /// <summary>
        /// RQ-COMPLIANT: Computes total quantum momentum using spectral coordinates.
        /// 
        /// In RQ-hypothesis, there is no external spacetime. Momentum is computed
        /// from wavefunction phases and the EMERGENT geometry (spectral coordinates
        /// from graph Laplacian eigenvectors), NOT external coordinates.
        /// 
        /// If spectral coordinates are not available, returns topological momentum
        /// based on phase gradients along graph edges (purely graph-based).
        /// </summary>
        public (double PX, double PY) ComputeTotalMomentum()
        {
            if (_waveMulti == null) return (0.0, 0.0);
            
            int d = GaugeDimension;
            double px = 0.0, py = 0.0;
            
            // RQ-FIX: Use spectral coordinates (emergent geometry) instead of external coordinates
            bool hasSpectralCoords = SpectralX != null && SpectralX.Length >= N &&
                                     SpectralY != null && SpectralY.Length >= N;
            
            for (int i = 0; i < N; i++)
            {
                foreach (int j in Neighbors(i))
                {
                    if (i >= j) continue;
                    
                    // Compute average phase at each node
                    double phaseI = 0.0, phaseJ = 0.0;
                    for (int a = 0; a < d; a++)
                    {
                        int idxI = i * d + a;
                        int idxJ = j * d + a;
                        phaseI += _waveMulti[idxI].Phase;
                        phaseJ += _waveMulti[idxJ].Phase;
                    }
                    phaseI /= d;
                    phaseJ /= d;
                    double dPhase = phaseJ - phaseI;
                    
                    if (hasSpectralCoords)
                    {
                        // RQ-COMPLIANT: Use spectral coordinates (emergent geometry)
                        double dx = SpectralX[j] - SpectralX[i];
                        double dy = SpectralY[j] - SpectralY[i];
                        px += dPhase * dx;
                        py += dPhase * dy;
                    }
                    else
                    {
                        // Fallback: Purely topological momentum
                        // Direction vector from edge weight gradient
                        double w = Weights[i, j];
                        px += dPhase * w;  // Weighted phase current (scalar)
                        py += dPhase * w * 0.5;  // Arbitrary orthogonal component
                    }
                }
            }
            
            return (px, py);
        }

        public void UpdateNodeStatesFromWavefunction()
        {
            if (_waveMulti == null) return; int d = GaugeDimension; var scalarField = new double[N];
            for (int i = 0; i < N; i++) { double mag = 0.0; for (int a = 0; a < d; a++) mag += _waveMulti[i * d + a].Magnitude; scalarField[i] = mag; }
            double threshold = LocalQuantile(scalarField, 0.5);
            var rand = _rng;
            for (int i = 0; i < N; i++)
            {
                if (rand.NextDouble() < 0.2) continue;
                double storedBoost = StoredEnergy != null ? (1.0 + 0.2 * StoredEnergy[i]) : 1.0; // integrate accumulator influence (item 1) into quantum-driven mode
                bool excited = scalarField[i] * storedBoost >= threshold;
                if (excited) State[i] = NodeState.Excited;
                else if (State[i] == NodeState.Excited) State[i] = NodeState.Refractory;
                else if (State[i] == NodeState.Refractory) State[i] = NodeState.Rest;
            }
        }

        public void PerturbGauge(double epsilon = 0.05)
        {
            if (_edgePhase == null) return;
            var rng = _rng;
            for (int i = 0; i < N; i++)
            {
                for (int j = i + 1; j < N; j++)
                {
                    if (!Edges[i, j]) continue;
                    double delta = (rng.NextDouble() * 2.0 - 1.0) * epsilon;
                    _edgePhase[i, j] += delta;
                    _edgePhase[j, i] -= delta;
                }
            }
            // gauge perturbation can alter effective propagation -> refresh delays
            UpdateEdgeDelaysFromDistances();
        }

        public Complex[] GetWavefunction()
        {
            var copy = new Complex[N];
            if (_wavefunction != null && _wavefunction.Length == N)
            {
                for (int i = 0; i < N; i++) copy[i] = _wavefunction[i];
            }
            return copy;
        }

        public Complex[] ColorCurrentOnEdge(int i, int j)
        {
            int d = GaugeDimension;
            var J = new Complex[d];
            if (_waveMulti == null) return J;
            for (int a = 0; a < d; a++)
            {
                var psiI = _waveMulti[i * d + a];
                var psiJ = _waveMulti[j * d + a];
                J[a] = psiI * Complex.Conjugate(psiJ);
            }
            return J;
        }

        public void ApplyLaplacian(Complex[] psi, Complex[] result)
        {
            int n = N;
            Parallel.For(0, n, i =>
            {
                Complex acc = Complex.Zero;
                double deg = 0.0;
                for (int j = 0; j < n; j++)
                {
                    double w = Weights[i, j];
                    if (w <= 0.0) continue;
                    deg += w;
                    acc += w * psi[j];
                }
                result[i] = acc - deg * psi[i];
            });
        }

        // Gauge-covariant Laplacian with U(1)/SU(d) link usage
        // OPTIMIZED: Uses stackalloc for small allocations and improved cache locality
        public void ApplyGaugeLaplacian(Complex[] psi, Complex[] result)
        {
            int n = N;
            int d = GaugeDimension;
            
            // Use larger batch size for better parallelization
            int batchSize = Math.Max(1, n / Environment.ProcessorCount);
            
            Parallel.For(0, n, new ParallelOptions { MaxDegreeOfParallelism = Environment.ProcessorCount }, i =>
            {
                int baseI = i * d;
                
                // Use stackalloc for small dimension arrays (avoid heap allocation)
                Span<Complex> acc = d <= 4 
                    ? stackalloc Complex[d] 
                    : new Complex[d];
                acc.Clear();
                
                double deg = 0.0;
                
                // Cache the psi values for this node
                ReadOnlySpan<Complex> psiI = psi.AsSpan(baseI, d);
                
                foreach (int j in Neighbors(i))
                {
                    double w = Weights[i, j];
                    if (w <= 0.0) continue;
                    deg += w;
                    
                    int baseJ = j * d;
                    ReadOnlySpan<Complex> psiJ = psi.AsSpan(baseJ, d);
                    
                    if (d == 1)
                    {
                        Complex uij = GetU1Link(i, j);
                        acc[0] += w * uij * psiJ[0];
                    }
                    else
                    {
                        // Fallback: identity (no explicit SU(d) matrix available in this build)
                        for (int a = 0; a < d; a++)
                            acc[a] += w * psiJ[a];
                    }
                }
                
                // Write result
                for (int a = 0; a < d; a++)
                    result[baseI + a] = acc[a] - deg * psiI[a];
            });
        }

        /// <summary>
        /// Full gauge-covariant Laplacian using exponential map from Lie algebra to Lie group.
        /// This connects the quantum particle to the gluon field via: U_ij = exp(i * A_ij * dx).
        /// For SU(3) gauge theory, A_ij are the 8 Gell-Mann generator coefficients.
        /// </summary>
        /// <param name="psi">Input wavefunction (N * GaugeDimension components)</param>
        /// <param name="result">Output array for the Laplacian result</param>
        public void ApplyGaugeLaplacianFull(Complex[] psi, Complex[] result)
        {
            int n = N;
            int d = GaugeDimension;

            // Default lattice spacing
            const double dx = 1.0;

            Parallel.For(0, n, i =>
            {
                int baseI = i * d;
                Span<Complex> psiI = psi.AsSpan(baseI, d);
                var acc = new Complex[d];
                double deg = 0.0;

                foreach (int j in Neighbors(i))
                {
                    double w = Weights[i, j];
                    if (w <= 0.0) continue;
                    deg += w;

                    int baseJ = j * d;
                    Span<Complex> psiJ = psi.AsSpan(baseJ, d);

                    if (d == 1)
                    {
                        // U(1) case: U_ij = exp(i * A_ij * dx)
                        Complex uij = GetU1Link(i, j);
                        acc[0] += w * uij * psiJ[0];
                    }
                    else
                    {
                        // SU(d) case: Build gauge link matrix from gluon field
                        // U_ij = exp(i * sum_a A^a_ij * T^a * dx)
                        // where T^a are generators (Gell-Mann matrices for SU(3))

                        // Get transported psi_j through gauge link
                        var psiTransported = ApplyGaugeLinkMatrix(i, j, psiJ.ToArray(), dx);

                        for (int a = 0; a < d; a++)
                            acc[a] += w * psiTransported[a];
                    }
                }

                for (int a = 0; a < d; a++)
                    result[baseI + a] = acc[a] - deg * psiI[a];
            });
        }

        /// <summary>
        /// Apply gauge link matrix U_ij to a vector using exponential map.
        /// U_ij = exp(i * A_ij * dx) where A_ij is the gauge potential (gluon field).
        /// Uses first-order approximation: exp(iM) ≈ I + iM for small M.
        /// </summary>
        private Complex[] ApplyGaugeLinkMatrix(int i, int j, Complex[] vec, double dx)
        {
            int d = GaugeDimension;
            var result = new Complex[d];

            // Check if gluon field is available
            if (_gluonField == null || d != 3)
            {
                // Fallback to identity transport
                for (int a = 0; a < d; a++)
                    result[a] = vec[a];
                return result;
            }

            // Build the generator matrix M = sum_a A^a_ij * T^a
            // For SU(3), we use Gell-Mann matrices T^a = lambda^a / 2
            // Total matrix is 3x3 complex

            var M = new Complex[d, d];

            // Sum over 8 color components
            for (int c = 0; c < 8; c++)
            {
                double A_c = _gluonField[i, j, c];
                if (Math.Abs(A_c) < 1e-12) continue;

                // Add contribution from generator T^c
                AddGellMannGenerator(M, c, A_c * dx);
            }

            // Apply exponential map: U = exp(i*M) ≈ I + i*M (first order)
            // For better accuracy, we can use: U ≈ I + i*M - M^2/2 (second order)
            // But first order is often sufficient for small dx

            for (int a = 0; a < d; a++)
            {
                result[a] = vec[a]; // Identity part

                // Add i*M*vec contribution
                for (int b = 0; b < d; b++)
                {
                    result[a] += Complex.ImaginaryOne * M[a, b] * vec[b];
                }
            }

            return result;
        }

        /// <summary>
        /// Add contribution from Gell-Mann generator lambda^c to matrix M.
        /// M += coefficient * (lambda^c / 2)
        /// </summary>
        private static void AddGellMannGenerator(Complex[,] M, int c, double coefficient)
        {
            double half = coefficient * 0.5;

            switch (c)
            {
                case 0: // lambda_1
                    M[0, 1] += half;
                    M[1, 0] += half;
                    break;
                case 1: // lambda_2
                    M[0, 1] += new Complex(0, -half);
                    M[1, 0] += new Complex(0, half);
                    break;
                case 2: // lambda_3
                    M[0, 0] += half;
                    M[1, 1] -= half;
                    break;
                case 3: // lambda_4
                    M[0, 2] += half;
                    M[2, 0] += half;
                    break;
                case 4: // lambda_5
                    M[0, 2] += new Complex(0, -half);
                    M[2, 0] += new Complex(0, half);
                    break;
                case 5: // lambda_6
                    M[1, 2] += half;
                    M[2, 1] += half;
                    break;
                case 6: // lambda_7
                    M[1, 2] += new Complex(0, -half);
                    M[2, 1] += new Complex(0, half);
                    break;
                case 7: // lambda_8
                    double sqrt3 = Math.Sqrt(3.0);
                    M[0, 0] += half / sqrt3;
                    M[1, 1] += half / sqrt3;
                    M[2, 2] -= 2.0 * half / sqrt3;
                    break;
            }
        }

        public void StepSchrodinger(double dtRq, double hbarRq)
        {
            int n = N;
            int d = GaugeDimension;
            int len = n * d;
            if (_psiDot == null || _psiDot.Length != len)
                _psiDot = new Complex[len];
            if (_waveMulti == null || _waveMulti.Length != len) return;
            ApplyGaugeLaplacian(_waveMulti, _psiDot);
            Parallel.For(0, n, i =>
            {
                double mi = (_correlationMass != null && _correlationMass.Length == N) ? Math.Max(1e-6, _correlationMass[i]) : 1.0;
                double inv2m = 0.5 / mi;
                double Vi = PotentialFromCorrelations(i);
                for (int a = 0; a < d; a++)
                {
                    int idx = i * d + a;
                    Complex D2psi = _psiDot[idx];
                    Complex rhs = -inv2m * D2psi + Vi * _waveMulti[idx];
                    Complex dpsi = -Complex.ImaginaryOne * (dtRq / hbarRq) * rhs;
                    _waveMulti[idx] += dpsi;
                }
            });
        }

        public void StepSchrodingerUnitary(double dtRq, double hbarRq, int iterations = 3)
        {
            int n = N;
            int d = GaugeDimension;
            int len = n * d;

            if (_psiDot == null || _psiDot.Length != len)
                _psiDot = new Complex[len];
            if (_waveMulti == null || _waveMulti.Length != len) return;

            ApplyGaugeLaplacian(_waveMulti, _psiDot);
            Complex alpha = -Complex.ImaginaryOne * (dtRq / (2.0 * hbarRq));
            var rhs = new Complex[len];
            for (int i = 0; i < n; i++)
            {
                int fl = FractalLevel != null && i < FractalLevel.Length ? FractalLevel[i] : 0;
                double fboost = 1.0 + 0.1 * fl; // checklist fractal boost
                double potLocal = LocalPotential != null && i < LocalPotential.Length ? LocalPotential[i] : 0.0;
                double phaseBase = potLocal * 0.1; // checklist phase factor
                for (int a = 0; a < d; a++)
                {
                    int idx = i * d + a;
                    double mi = (_correlationMass != null && _correlationMass.Length == N) ? Math.Max(1e-6, _correlationMass[i]) : 1.0;
                    double inv2m = 0.5 / mi;
                    double Vi = PotentialFromCorrelations(i);
                    Complex Hpsi = -inv2m * _psiDot[idx] + Vi * _waveMulti[idx];
                    Complex mod = Complex.FromPolarCoordinates(1.0, phaseBase);
                    rhs[idx] = _waveMulti[idx] + fboost * alpha * Hpsi * mod;
                }
            }

            var psiNext = (Complex[])_waveMulti.Clone();
            for (int it = 0; it < iterations; it++)
            {
                ApplyGaugeLaplacian(psiNext, _psiDot);
                for (int i = 0; i < n; i++)
                {
                    int fl = FractalLevel != null && i < FractalLevel.Length ? FractalLevel[i] : 0;
                    double fboost = 1.0 + 0.1 * fl;
                    double potLocal = LocalPotential != null && i < LocalPotential.Length ? LocalPotential[i] : 0.0;
                    double phaseBase = potLocal * 0.1;
                    Complex mod = Complex.FromPolarCoordinates(1.0, phaseBase);
                    for (int a = 0; a < d; a++)
                    {
                        int idx = i * d + a;
                        double mi = (_correlationMass != null && _correlationMass.Length == N) ? Math.Max(1e-6, _correlationMass[i]) : 1.0;
                        double inv2m = 0.5 / mi;
                        double Vi = PotentialFromCorrelations(i);
                        Complex Hpsi = -inv2m * _psiDot[idx] + Vi * psiNext[idx];
                        Complex diag = Complex.One + alpha * Vi;
                        psiNext[idx] = (rhs[idx] - fboost * alpha * (-inv2m * _psiDot[idx]) * mod) / diag;
                    }
                }
            }
            _waveMulti = psiNext;
        }

        public double GetQuantumNorm()
        {
            if (_waveMulti == null) return 0.0;
            double acc = 0.0;
            foreach (var z in _waveMulti) acc += z.Magnitude * z.Magnitude;
            return acc;
        }

        public double ComputeQuantumEnergy(double hbarRq = 1.0)
        {
            if (_waveMulti == null) return 0.0;
            int n = N;
            int d = GaugeDimension;
            int len = n * d;
            if (_psiDot == null || _psiDot.Length != len)
                _psiDot = new Complex[len];
            ApplyGaugeLaplacian(_waveMulti, _psiDot);
            double E = 0.0;
            for (int i = 0; i < n; i++)
            {
                double mi = (_correlationMass != null && _correlationMass.Length == N) ? Math.Max(1e-6, _correlationMass[i]) : 1.0;
                double inv2m = 0.5 / mi;
                double Vi = PotentialFromCorrelations(i);
                for (int a = 0; a < d; a++)
                {
                    int idx = i * d + a;
                    Complex psi = _waveMulti[idx];
                    Complex D2psi = _psiDot[idx];
                    Complex kinetic = -inv2m * Complex.Conjugate(psi) * D2psi;
                    double potential = Vi * psi.Magnitude * psi.Magnitude;
                    E += kinetic.Real + potential;
                }
            }
            return E;
        }

        private double CorrelationEnergyFromMass(double m) => 0.5 * m;
        private double PhaseIncrement(int i)
        {
            if (_edgePhase == null) return 0;
            double sum = 0; int c = 0;
            foreach (int j in Neighbors(i)) { sum += Math.Abs(_edgePhase[i, j]); c++; }
            return c > 0 ? sum / c : 0;
        }
        private double PotentialFromCorrelations(int i)
        {
            double sum = 0.0; int deg = 0;
            foreach (int j in Neighbors(i)) { sum += Weights[i, j]; deg++; }
            if (deg == 0) return 0.0;
            // Local z-score style normalisation
            double meanLocal = sum / deg;
            double varLocal = 0.0;
            foreach (int j in Neighbors(i))
            {
                double w = Weights[i, j];
                double diff = w - meanLocal;
                varLocal += diff * diff;
            }
            varLocal = varLocal / deg;
            double stdLocal = varLocal > 0 ? Math.Sqrt(varLocal) : 1.0;
            double z = (sum - meanLocal * deg) / Math.Sqrt(varLocal * deg + 1e-12);
            return double.IsNaN(z) ? 0.0 : z;
        }

        private static double LocalQuantile(double[] data, double q)
        {
            if (data == null || data.Length == 0) return 0.0;
            var sorted = data.OrderBy(x => x).ToArray();
            int idx = (int)Math.Clamp(q * (sorted.Length - 1), 0, sorted.Length - 1);
            return sorted[idx];
        }

        public double GetNodePhase(int i)
        {
            if (_waveMulti == null || i < 0 || i >= N) return 0.0;
            int d = GaugeDimension;
            double phase = 0.0;
            for (int a = 0; a < d; a++) phase += _waveMulti[i * d + a].Phase;
            return phase / Math.Max(1, d);
        }

        /// <summary>
        /// Stub method for future Dirac operator implementation on the graph.
        /// The Dirac operator D couples spinor fields to the graph geometry and gauge fields.
        /// On a graph, D is typically defined as: D = sum_edges gamma^mu * nabla_mu
        /// where gamma^mu are Dirac matrices and nabla_mu is the covariant derivative.
        /// 
        /// The Dirac equation: (i*D - m)*psi = 0 describes fermion dynamics.
        /// 
        /// This stub initializes spinor fields if needed but does not perform evolution yet.
        /// Future implementation should:
        /// 1. Define gamma matrices for the graph structure (based on edge directions)
        /// 2. Apply covariant derivative using gauge links
        /// 3. Evolve spinor field: psi(t+dt) = exp(-i*D*dt) * psi(t)
        /// </summary>
        /// <param name="dt">Time step</param>
        public void StepDirac(double dt)
        {
            // Ensure spinor fields are initialized
            if (_spinorA == null || _spinorA.Length != N)
            {
                _spinorA = new Complex[N];
                for (int i = 0; i < N; i++)
                    _spinorA[i] = new Complex(_rng.NextDouble() * 0.01, _rng.NextDouble() * 0.01);
            }
            if (_spinorB == null || _spinorB.Length != N)
            {
                _spinorB = new Complex[N];
                for (int i = 0; i < N; i++)
                    _spinorB[i] = new Complex(_rng.NextDouble() * 0.01, _rng.NextDouble() * 0.01);
            }
            if (_spinorC == null || _spinorC.Length != N)
            {
                _spinorC = new Complex[N];
                for (int i = 0; i < N; i++)
                    _spinorC[i] = new Complex(_rng.NextDouble() * 0.01, _rng.NextDouble() * 0.01);
            }
            if (_spinorD == null || _spinorD.Length != N)
            {
                _spinorD = new Complex[N];
                for (int i = 0; i < N; i++)
                    _spinorD[i] = new Complex(_rng.NextDouble() * 0.01, _rng.NextDouble() * 0.01);
            }

            // TODO: Implement full Dirac operator evolution
            // The implementation would involve:
            // 1. For each node i, compute D*psi_i using neighbors
            // 2. D*psi_i = sum_j gamma^(ij) * U_ij * psi_j
            //    where gamma^(ij) is direction-dependent Dirac matrix
            //    and U_ij is the gauge link
            // 3. Update: psi_new = psi - i*dt*D*psi (first order)
            // 4. Normalize to preserve probability

            // For now, apply minimal phase evolution based on local mass
            double hbar = 1.0;
            for (int i = 0; i < N; i++)
            {
                double m = (_correlationMass != null && i < _correlationMass.Length) ? _correlationMass[i] : 1.0;
                double phase = -m * dt / hbar;
                Complex phaseRotation = Complex.FromPolarCoordinates(1.0, phase);

                _spinorA[i] *= phaseRotation;
                _spinorB[i] *= phaseRotation;
                _spinorC[i] *= phaseRotation;
                _spinorD[i] *= phaseRotation;
            }
        }
        
        /// <summary>
        /// Apply decoherence to quantum wavefunction through phase noise.
        /// This simulates environmental interaction that destroys coherence.
        /// </summary>
        public void ApplyDecoherence(double rate)
        {
            if (_waveMulti == null || _waveMulti.Length == 0) return;
            if (rate <= 0) return;
            
            var rnd = new Random();
            
            // Phase noise: add random phase shifts to each amplitude
            for (int idx = 0; idx < _waveMulti.Length; idx++)
            {
                double maxAngle = rate; // maximum ±rate radians
                double dθ = (rnd.NextDouble() * 2 - 1) * maxAngle;
                Complex phaseShift = Complex.FromPolarCoordinates(1.0, dθ);
                _waveMulti[idx] *= phaseShift;
            }
            
            // Optional: stochastic collapse with small probability
            double collapseProbability = rate * 0.01; // Much smaller than phase noise
            if (rnd.NextDouble() < collapseProbability && GaugeDimension > 0)
            {
                // Pick a random node to collapse
                int nodeToCollapse = rnd.Next(N);
                int d = GaugeDimension;
                
                // Find basis state with max amplitude for this node
                double maxAmp = 0.0;
                int maxIndex = 0;
                for (int a = 0; a < d; a++)
                {
                    double amp = _waveMulti[nodeToCollapse * d + a].Magnitude;
                    if (amp > maxAmp)
                    {
                        maxAmp = amp;
                        maxIndex = a;
                    }
                }
                
                // Collapse: zero out all components except maxIndex
                for (int a = 0; a < d; a++)
                {
                    _waveMulti[nodeToCollapse * d + a] = (a == maxIndex) ? Complex.One : Complex.Zero;
                }
            }
        }
        
        // ==================== SYMPLECTIC INTEGRATOR ====================
        // Implements checklist item 6: Numerical stability via leapfrog
        
        // Field momenta (conjugate to phases)
        private double[]? _fieldMomenta;
        
        /// <summary>
        /// Initialize field momenta for symplectic integration.
        /// </summary>
        private void EnsureFieldMomenta()
        {
            if (_fieldMomenta == null || _fieldMomenta.Length != N)
            {
                _fieldMomenta = new double[N];
            }
        }
        
        /// <summary>
        /// Perform symplectic (leapfrog) integration step for field dynamics.
        /// The leapfrog method is time-reversible and conserves phase space volume,
        /// preventing energy drift common in simple Euler methods.
        /// 
        /// Steps:
        /// 1. p(t + dt/2) = p(t) - dH/dq * dt/2  (half-step momentum)
        /// 2. q(t + dt)   = q(t) + dH/dp * dt    (full-step position)
        /// 3. p(t + dt)   = p(t + dt/2) - dH/dq * dt/2 (half-step momentum)
        /// 
        /// Implements checklist item 6.1: Symplectic integrator (Leapfrog).
        /// </summary>
        /// <param name="dt">Time step</param>
        public void StepSymplectic(double dt)
        {
            EnsureFieldMomenta();
            
            // Step 1: Half-step momentum update
            // p(t + dt/2) = p(t) - dH/dq * dt/2
            UpdateMomenta(dt / 2.0);
            
            // Step 2: Full-step field (phase) update
            // q(t + dt) = q(t) + dH/dp * dt
            UpdateFields(dt);
            
            // Step 3: Half-step momentum update
            // p(t + dt) = p(t + dt/2) - dH/dq * dt/2
            UpdateMomenta(dt / 2.0);
        }
        
        /// <summary>
        /// Update field momenta based on Hamiltonian gradient.
        /// dp/dt = -dH/dq where q are the field values (phases).
        /// </summary>
        /// <param name="dt">Time step</param>
        private void UpdateMomenta(double dt)
        {
            if (_fieldMomenta == null)
                return;
            
            // For each node, compute dH/dq (gradient of Hamiltonian w.r.t. field)
            for (int i = 0; i < N; i++)
            {
                double dHdq = ComputeFieldGradient(i);
                
                // Update momentum: p -= dH/dq * dt
                _fieldMomenta[i] -= dHdq * dt;
            }
        }
        
        /// <summary>
        /// Update field values (phases) based on momenta.
        /// dq/dt = dH/dp = p (for quadratic kinetic energy)
        /// </summary>
        /// <param name="dt">Time step</param>
        private void UpdateFields(double dt)
        {
            if (_fieldMomenta == null || LocalPotential == null)
                return;
            
            // For quadratic kinetic energy: dH/dp = p
            // So q += p * dt
            for (int i = 0; i < N; i++)
            {
                LocalPotential[i] += _fieldMomenta[i] * dt;
                
                // Apply boundary conditions / clamping
                if (LocalPotential[i] < 0)
                    LocalPotential[i] = 0;
                if (LocalPotential[i] > 5.0)
                    LocalPotential[i] = 5.0;
            }
            
            // Also update wavefunction phases symplectically
            if (_waveMulti != null)
            {
                int d = GaugeDimension;
                for (int i = 0; i < N; i++)
                {
                    // Scale factor for phase rotation based on momentum
                    double phaseIncrement = _fieldMomenta[i] * dt * PhysicsConstants.SymplecticPhaseScale;
                    Complex phaseRotation = Complex.FromPolarCoordinates(1.0, phaseIncrement);
                    
                    for (int a = 0; a < d; a++)
                    {
                        int idx = i * d + a;
                        _waveMulti[idx] *= phaseRotation;
                    }
                }
            }
        }
        
        /// <summary>
        /// Compute gradient of Hamiltonian with respect to field at node i.
        /// dH/dq_i includes kinetic, potential, and interaction terms.
        /// </summary>
        /// <param name="nodeId">Node index</param>
        /// <returns>Gradient dH/dq at node</returns>
        private double ComputeFieldGradient(int nodeId)
        {
            double gradient = 0.0;
            
            // Potential term: V(q) contribution
            double localPotential = (LocalPotential != null && nodeId < LocalPotential.Length) 
                ? LocalPotential[nodeId] : 0.0;
            
            // Simple harmonic potential: V = (1/2) * ω² * q² → dV/dq = ω² * q
            double omega = PhysicsConstants.FieldHarmonicFrequency;
            gradient += omega * omega * localPotential;
            
            // Interaction term: Laplacian contribution from neighbors
            // For diffusive coupling: H_int = -Σ w_ij * q_i * q_j
            // → dH/dq_i = -Σ w_ij * q_j
            foreach (int j in Neighbors(nodeId))
            {
                double neighborPotential = (LocalPotential != null && j < LocalPotential.Length) 
                    ? LocalPotential[j] : 0.0;
                gradient -= Weights[nodeId, j] * neighborPotential;
            }
            
            // Curvature coupling: gravity affects field dynamics
            double curvature = GetLocalCurvature(nodeId);
            gradient += 0.1 * curvature * localPotential;
            
            return gradient;
        }
        
        /// <summary>
        /// Verify energy conservation after symplectic step.
        /// Returns relative energy change (should be small for symplectic integrator).
        /// </summary>
        /// <param name="previousEnergy">Energy before step</param>
        /// <returns>Relative energy change |ΔE/E|</returns>
        public double VerifySymplecticConservation(double previousEnergy)
        {
            double currentEnergy = ComputeSymplecticEnergy();
            
            if (Math.Abs(previousEnergy) < 1e-10)
                return Math.Abs(currentEnergy - previousEnergy);
            
            return Math.Abs(currentEnergy - previousEnergy) / Math.Abs(previousEnergy);
        }
        
        /// <summary>
        /// Compute total energy in symplectic coordinates (H = T + V).
        /// </summary>
        /// <returns>Total Hamiltonian energy</returns>
        public double ComputeSymplecticEnergy()
        {
            EnsureFieldMomenta();
            
            double kineticEnergy = 0.0;
            double potentialEnergy = 0.0;
            
            // Kinetic energy: T = (1/2) * Σ p²
            if (_fieldMomenta != null)
            {
                for (int i = 0; i < N; i++)
                {
                    kineticEnergy += 0.5 * _fieldMomenta[i] * _fieldMomenta[i];
                }
            }
            
            // Potential energy from local potentials
            if (LocalPotential != null)
            {
                double omega = PhysicsConstants.FieldHarmonicFrequency;
                for (int i = 0; i < LocalPotential.Length; i++)
                {
                    // Harmonic potential
                    potentialEnergy += 0.5 * omega * omega * LocalPotential[i] * LocalPotential[i];
                    
                    // Interaction energy
                    foreach (int j in Neighbors(i))
                    {
                        if (j > i && j < LocalPotential.Length)
                        {
                            potentialEnergy -= 0.5 * Weights[i, j] * LocalPotential[i] * LocalPotential[j];
                        }
                    }
                }
            }
            
            return kineticEnergy + potentialEnergy;
        }
    }
}