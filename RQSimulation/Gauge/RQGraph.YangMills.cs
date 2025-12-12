using System;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Runtime.Intrinsics;
using System.Runtime.Intrinsics.X86;

namespace RQSimulation
{
    /// <summary>
    /// Implements Yang-Mills gauge field dynamics for electroweak (SU(2)×U(1))
    /// and strong (SU(3)) interactions in the RQ framework.
    /// Gauge fields live on edges and couple to matter fields on nodes.
    /// </summary>
    public partial class RQGraph
    {
        private double[,,]? _gluonField;      // [i,j,8]
        private double[,,]? _weakField;       // [i,j,3]
        private double[,]? _hyperchargeField; // [i,j]
        private double[,,]? _gluonFieldStrength;
        private double[,,]? _weakFieldStrength;
        private double[,]? _hyperchargeFieldStrength;

        // Reusable delta buffers to avoid per-step large allocations (freeze mitigation)
        private double[,,]? _gluonDelta;
        private double[,,]? _weakDelta;
        private double[,]? _hyperDelta;

        public double StrongCoupling { get; set; } = 1.0;
        public double WeakCoupling { get; set; } = 0.65;
        public double HypergaugeCoupling { get; set; } = 0.35;
        private readonly double[] _colorWeights = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };

        [ThreadStatic] private static bool _isEvolvingYangMills;

        // Precompute node color densities (first 3 color components) reused inside evolution.
        private double[]? _cachedColorDensity;
        private double[]? _cachedWeakDensity;
        private double[]? _cachedHyperDensity;

        // SU(3) structure constants - ONLY access via GetFabcFull() to ensure initialization
        private static readonly Lazy<double[,,]> s_fabcFullCache = new Lazy<double[,,]>(() =>
        {
            var f = new double[8, 8, 8];

            // (0,1,2) family
            f[0, 1, 2] = 1.0; f[1, 2, 0] = 1.0; f[2, 0, 1] = 1.0;
            f[0, 2, 1] = -1.0; f[2, 1, 0] = -1.0; f[1, 0, 2] = -1.0;

            // (0,3,6) family
            f[0, 3, 6] = 0.5; f[3, 6, 0] = 0.5; f[6, 0, 3] = 0.5;
            f[0, 6, 3] = -0.5; f[6, 3, 0] = -0.5; f[3, 0, 6] = -0.5;

            // (0,4,5) family
            f[0, 4, 5] = -0.5; f[4, 5, 0] = -0.5; f[5, 0, 4] = -0.5;
            f[0, 5, 4] = 0.5; f[5, 4, 0] = 0.5; f[4, 0, 5] = 0.5;

            // (1,3,5) family
            f[1, 3, 5] = 0.5; f[3, 5, 1] = 0.5; f[5, 1, 3] = 0.5;
            f[1, 5, 3] = -0.5; f[5, 3, 1] = -0.5; f[3, 1, 5] = -0.5;

            // (1,4,6) family
            f[1, 4, 6] = 0.5; f[4, 6, 1] = 0.5; f[6, 1, 4] = 0.5;
            f[1, 6, 4] = -0.5; f[6, 4, 1] = -0.5; f[4, 1, 6] = -0.5;

            // (2,3,4) family
            f[2, 3, 4] = 0.5; f[3, 4, 2] = 0.5; f[4, 2, 3] = 0.5;
            f[2, 4, 3] = -0.5; f[4, 3, 2] = -0.5; f[3, 2, 4] = -0.5;

            // (2,5,6) family
            f[2, 5, 6] = -0.5; f[5, 6, 2] = -0.5; f[6, 2, 5] = -0.5;
            f[2, 6, 5] = 0.5; f[6, 5, 2] = 0.5; f[5, 2, 6] = 0.5;

            // (3,4,7) family - sqrt(3)/2
            double sqrt3_2 = 0.8660254037844386;
            f[3, 4, 7] = sqrt3_2; f[4, 7, 3] = sqrt3_2; f[7, 3, 4] = sqrt3_2;
            f[3, 7, 4] = -sqrt3_2; f[7, 4, 3] = -sqrt3_2; f[4, 3, 7] = -sqrt3_2;

            // (5,6,7) family - sqrt(3)/2
            f[5, 6, 7] = sqrt3_2; f[6, 7, 5] = sqrt3_2; f[7, 5, 6] = sqrt3_2;
            f[5, 7, 6] = -sqrt3_2; f[7, 6, 5] = -sqrt3_2; f[6, 5, 7] = -sqrt3_2;

            return f;
        }, LazyThreadSafetyMode.ExecutionAndPublication);

        [MethodImpl(MethodImplOptions.NoInlining)] // Changed from AggressiveInlining to NoInlining
        private static double[,,] GetFabcFull()
        {
            return s_fabcFullCache.Value;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double AbsSquared(Complex z)
        {
            double re = z.Real;
            double im = z.Imaginary;
            return re * re + im * im;
        }

        public void InitYangMillsFields()
        {
            // CHECKLIST ITEM 5: Check for conflicting gauge representations
            // Both _gaugeSU (matrix) and _gluonField (component) represent SU(3).
            // Using both simultaneously causes inconsistencies.
            if (_gaugeSU != null && GaugeDimension == 3)
            {
                Console.WriteLine("[WARNING] Both _gaugeSU (matrix) and _gluonField (component) active for SU(3).");
                Console.WriteLine("[WARNING] Disabling _gaugeSU to avoid conflict. Use component fields only.");
                _gaugeSU = null!;
            }
            
            // RQ-FIX: Initialize SU(3) matrices for proper lattice gauge theory
            if (GaugeDimension == 3)
            {
                ConfigureGaugeDimension(3);
            }

            _gluonField = new double[N, N, 8];
            _weakField = new double[N, N, 3];
            _hyperchargeField = new double[N, N];
            _gluonFieldStrength = new double[N, N, 8];
            _weakFieldStrength = new double[N, N, 3];
            _hyperchargeFieldStrength = new double[N, N];

            for (int i = 0; i < N; i++)
            {
                foreach (int j in Neighbors(i))
                {
                    for (int a = 0; a < 8; a++)
                        _gluonField[i, j, a] = (_rng.NextDouble() - 0.5) * 0.01;
                    for (int a = 0; a < 3; a++)
                        _weakField[i, j, a] = (_rng.NextDouble() - 0.5) * 0.01;
                    _hyperchargeField[i, j] = (_rng.NextDouble() - 0.5) * 0.01;
                }
            }
        }

        public void ComputeGluonFieldStrength()
        {
            // Redirect to optimized version in RQGraph.YangMills.Optimized.cs
            ComputeGluonFieldStrengthOptimized();
        }

        [MethodImpl(MethodImplOptions.NoInlining)] // Changed from AggressiveInlining to NoInlining
        private static bool TryGetStructureConstant(int a, int b, int c, out double fabc)
        {
            // Bounds check
            if ((uint)a >= 8u || (uint)b >= 8u || (uint)c >= 8u)
            {
                fabc = 0.0;
                return false;
            }

            // Diagonal check
            if (a == b || b == c || a == c)
            {
                fabc = 0.0;
                return false;
            }

            // Get cached structure constants
            var fabcFull = GetFabcFull();
            fabc = fabcFull[a, b, c];
            return fabc != 0.0;
        }

        public void ComputeWeakFieldStrength()
        {
            // Redirect to optimized version in RQGraph.YangMills.Optimized.cs
            ComputeWeakFieldStrengthOptimized();
        }

        public void ComputeHyperchargeFieldStrength()
        {
            if (_hyperchargeField == null || _hyperchargeFieldStrength == null) return;

            for (int i = 0; i < N; i++)
            {
                foreach (int j in Neighbors(i))
                {
                    if (j <= i) continue;

                    double curl = 0;
                    foreach (int k in Neighbors(i))
                    {
                        if (k == j) continue;
                        if (Edges[k, j])
                            curl += _hyperchargeField[i, k] - _hyperchargeField[k, j];
                    }

                    _hyperchargeFieldStrength[i, j] = curl;
                    _hyperchargeFieldStrength[j, i] = -curl;
                }
            }
        }

        public double ComputeYangMillsAction()
        {
            ComputeGluonFieldStrength();
            ComputeWeakFieldStrength();
            ComputeHyperchargeFieldStrength();

            double S = 0;

            if (_gluonFieldStrength != null)
            {
                for (int i = 0; i < N; i++)
                {
                    foreach (int j in Neighbors(i))
                    {
                        if (j <= i) continue;

                        for (int a = 0; a < 8; a++)
                        {
                            double F = _gluonFieldStrength[i, j, a];
                            S += 0.25 * F * F / (StrongCoupling * StrongCoupling);
                        }
                    }
                }
            }

            if (_weakFieldStrength != null)
            {
                for (int i = 0; i < N; i++)
                {
                    foreach (int j in Neighbors(i))
                    {
                        if (j <= i) continue;

                        for (int a = 0; a < 3; a++)
                        {
                            double W = _weakFieldStrength[i, j, a];
                            S += 0.25 * W * W / (WeakCoupling * WeakCoupling);
                        }
                    }
                }
            }

            if (_hyperchargeFieldStrength != null)
            {
                for (int i = 0; i < N; i++)
                {
                    foreach (int j in Neighbors(i))
                    {
                        if (j <= i) continue;

                        double B = _hyperchargeFieldStrength[i, j];
                        S += 0.25 * B * B / (HypergaugeCoupling * HypergaugeCoupling);
                    }
                }
            }

            return S;
        }

        public void EvolveYangMillsFields(double dt)
        {
            if (_isEvolvingYangMills) return;

            try
            {
                _isEvolvingYangMills = true;

                // RQ-FIX: Use proper SU(3) lattice evolution if available
                // This ensures unitarity via exponential map updates
                if (GaugeDimension == 3 && _gaugeSU3 != null)
                {
                    // 1. Pure gauge evolution (Wilson action)
                    UpdateNonAbelianGauge(dt);
                    
                    // 2. Matter coupling (currents)
                    UpdateGaugeFromMatter(dt);
                    
                    // 3. Enforce constraints
                    EnforceGaugeConstraintsInternal();
                    
                    // Sync back to component fields for visualization/compatibility
                    // (Optional, but good for consistency)
                    // SyncSU3ToComponents(); 
                    
                    // Continue to evolve Weak and Hypercharge fields below...
                }

                if (_gluonField == null) InitYangMillsFields();
                if (_waveMulti == null) return;

                int d = GaugeDimension;
                if (d < 3) return;

                int lenWave = _waveMulti.Length;

                // Cache color densities
                if (_cachedColorDensity == null || _cachedColorDensity.Length != N)
                    _cachedColorDensity = new double[N];

                for (int i = 0; i < N; i++)
                {
                    double rho = 0.0;
                    int baseIdx = i * d;

                    if (baseIdx < lenWave)
                    {
                        if (baseIdx + 0 < lenWave) rho += AbsSquared(_waveMulti[baseIdx + 0]);
                        if (baseIdx + 1 < lenWave) rho += AbsSquared(_waveMulti[baseIdx + 1]);
                        if (baseIdx + 2 < lenWave) rho += AbsSquared(_waveMulti[baseIdx + 2]);
                    }

                    _cachedColorDensity[i] = rho;
                }

                // Cache weak densities
                if (_spinorA != null)
                {
                    if (_cachedWeakDensity == null || _cachedWeakDensity.Length != N)
                        _cachedWeakDensity = new double[N];

                    for (int i = 0; i < N; i++)
                    {
                        double rho = 0.0;
                        double magA = _spinorA[i].Magnitude;
                        rho += magA * magA;

                        if (_spinorB != null && i < _spinorB.Length)
                        {
                            double magB = _spinorB[i].Magnitude;
                            rho += magB * magB;
                        }

                        _cachedWeakDensity[i] = rho;
                    }
                }

                // Cache hypercharge densities
                if (_spinorA != null || _spinorC != null)
                {
                    if (_cachedHyperDensity == null || _cachedHyperDensity.Length != N)
                        _cachedHyperDensity = new double[N];

                    for (int i = 0; i < N; i++)
                    {
                        double h = 0.0;

                        if (_spinorA != null && i < _spinorA.Length)
                        {
                            double magA = _spinorA[i].Magnitude;
                            h += -0.5 * magA * magA;

                            if (_spinorB != null && i < _spinorB.Length)
                            {
                                double magB = _spinorB[i].Magnitude;
                                h += -0.5 * magB * magB;
                            }
                        }

                        if (_spinorC != null && i < _spinorC.Length)
                        {
                            double magC = _spinorC[i].Magnitude;
                            h += -1.0 * magC * magC;

                            if (_spinorD != null && i < _spinorD.Length)
                            {
                                double magD = _spinorD[i].Magnitude;
                                h += -1.0 * magD * magD;
                            }
                        }

                        _cachedHyperDensity[i] = h;
                    }
                }

                // Initialize delta buffers
                if (_gluonDelta == null || _gluonDelta.GetLength(0) != N)
                    _gluonDelta = new double[N, N, 8];
                if (_weakDelta == null || _weakDelta.GetLength(0) != N)
                    _weakDelta = new double[N, N, 3];
                if (_hyperDelta == null || _hyperDelta.GetLength(0) != N)
                    _hyperDelta = new double[N, N];

                // Ensure spinor fields exist
                bool needsSpinors = (_spinorA == null || _spinorA.Length != N) ||
                                   (_spinorC == null || _spinorC.Length != N);
                if (needsSpinors)
                {
                    if (_spinorA == null || _spinorA.Length != N)
                        _spinorA = new Complex[N];
                    if (_spinorC == null || _spinorC.Length != N)
                        _spinorC = new Complex[N];
                }

                // Compute field strengths
                ComputeGluonFieldStrength();
                ComputeWeakFieldStrength();
                ComputeHyperchargeFieldStrength();

                int[] scratchI = new int[N];

                // Compute field evolution
                for (int i = 0; i < N; i++)
                {
                    var neighI = GetNeighborSpan(i, ref scratchI);

                    foreach (int j in neighI)
                    {
                        double[] Aij = new double[8]; // Changed from stackalloc to heap allocation
                        double[] Fij = new double[8]; // Changed from stackalloc to heap allocation

                        for (int t = 0; t < 8; t++)
                        {
                            Aij[t] = _gluonField![i, j, t];
                            Fij[t] = _gluonFieldStrength![i, j, t];
                        }

                        // Gluon field evolution (RQ-FIX: Multiplicative update per Checklist Item 3)
                        // Replaces additive update with unitary SU(3) matrix exponential
                        if (GaugeDimension == 3)
                        {
                            // 1. Get current link matrix U_ij
                            var U_old = GetSU3Link(i, j);
                            
                            // 2. Compute Force F (gradient of Action)
                            // F^a = div E^a + g f^{abc} A^b E^c - J^a
                            double[] force = new double[8];
                            
                            for (int a = 0; a < 8; a++)
                            {
                                double divF = 0.0;
                                for (int idx = 0; idx < neighI.Length; idx++)
                                {
                                    int k = neighI[idx];
                                    if (Edges[k, j])
                                        divF += _gluonFieldStrength![k, j, a] - Fij[a];
                                }

                                double selfInt = 0.0;
                                for (int b = 0; b < 8; b++)
                                {
                                    double Ab = Aij[b];
                                    if (Ab == 0.0) continue;

                                    for (int c = 0; c < 8; c++)
                                    {
                                        if (b == c) continue;
                                        if (!TryGetStructureConstant(a, b, c, out double fabc))
                                            continue;

                                        selfInt += StrongCoupling * fabc * Ab * Fij[c];
                                    }
                                }

                                double J = ComputeColorCurrentCached(i, j, a);
                                
                                // Force component
                                force[a] = divF + selfInt - J;
                            }
                            
                            // 3. Exponential update: U_new = exp(-i * dt * F) * U_old
                            // We use -dt because we minimize action (gradient descent)
                            for(int a=0; a<8; a++) force[a] *= -dt;
                            
                            var Update = Gauge.SU3Matrix.Exp(force);
                            var U_new = Update.Multiply(U_old);
                            
                            // 4. Update link
                            SetSU3Link(i, j, U_new);
                            
                            // Note: We do NOT update _gluonDelta here as we applied the update directly
                        }
                        else
                        {
                            // Fallback for non-SU(3) or legacy mode
                            for (int a = 0; a < 8; a++)
                            {
                                double divF = 0.0;
                                for (int idx = 0; idx < neighI.Length; idx++)
                                {
                                    int k = neighI[idx];
                                    if (Edges[k, j])
                                        divF += _gluonFieldStrength![k, j, a] - Fij[a];
                                }

                                double selfInt = 0.0;
                                for (int b = 0; b < 8; b++)
                                {
                                    double Ab = Aij[b];
                                    if (Ab == 0.0) continue;

                                    for (int c = 0; c < 8; c++)
                                    {
                                        if (b == c) continue;

                                        if (!TryGetStructureConstant(a, b, c, out double fabc))
                                            continue;

                                        selfInt += StrongCoupling * fabc * Ab * Fij[c];
                                    }
                                }

                                double J = ComputeColorCurrentCached(i, j, a);
                                _gluonDelta![i, j, a] = dt * (divF + selfInt - J);
                            }
                        }

                        // Weak field evolution
                        // CHECKLIST ITEM 2: Use explicit SU(2) structure constants
                        // instead of manual index cycling (a+1)%3, (a+2)%3
                        for (int a = 0; a < 3; a++)
                        {
                            double divW = 0.0;
                            for (int idx = 0; idx < neighI.Length; idx++)
                            {
                                int k = neighI[idx];
                                if (Edges[k, j])
                                    divW += _weakFieldStrength![k, j, a] - _weakFieldStrength[i, j, a];
                            }

                            // Use explicit Levi-Civita structure constants for SU(2)
                            // f^{abc} = ε_{abc} replaces the manual (a+1)%3, (a+2)%3 trick
                            double selfInt = 0.0;
                            for (int b = 0; b < 3; b++)
                            {
                                for (int c = 0; c < 3; c++)
                                {
                                    double fabc = Gauge.PauliMatrices.GetStructureConstant(a, b, c);
                                    if (fabc != 0.0)
                                    {
                                        selfInt += WeakCoupling * fabc * _weakField![i, j, b] * _weakFieldStrength![i, j, c];
                                    }
                                }
                            }

                            double Jw = ComputeWeakCurrent(i, j, a);
                            _weakDelta![i, j, a] = dt * (divW + selfInt - Jw);
                        }

                        // Hypercharge field evolution
                        double divB = 0.0;
                        for (int idx = 0; idx < neighI.Length; idx++)
                        {
                            int k = neighI[idx];
                            if (Edges[k, j])
                                divB += _hyperchargeFieldStrength![k, j] - _hyperchargeFieldStrength[i, j];
                        }

                        double Jh = ComputeHyperchargeCurrent(i, j);
                        _hyperDelta![i, j] = dt * (divB - Jh);
                    }
                }

                // Apply field updates
                for (int i = 0; i < N; i++)
                {
                    var neighI = GetNeighborSpan(i, ref scratchI);

                    foreach (int j in neighI)
                    {
                        // Skip gluon update if handled by SU(3) lattice engine (multiplicative update)
                        if (GaugeDimension != 3)
                        {
                            for (int a = 0; a < 8; a++)
                                _gluonField![i, j, a] += _gluonDelta![i, j, a];
                        }

                        for (int a = 0; a < 3; a++)
                            _weakField![i, j, a] += _weakDelta![i, j, a];

                        _hyperchargeField![i, j] += _hyperDelta![i, j];
                    }
                }
                
                // Checklist item 3: Enforce gauge constraints after evolution step
                // This projects phases onto the constraint surface (Gauss law)
                EnforceGaugeConstraintsInternal();
            }
            finally
            {
                _isEvolvingYangMills = false;
            }
        }
        
        /// <summary>
        /// Enforce gauge constraints after Yang-Mills evolution.
        /// 
        /// CHECKLIST ITEM 7: Extended Gauss law enforcement to non-Abelian fields.
        /// 
        /// For U(1): ∇·E = ρ (implemented in EnforceGaussLaw)
        /// For SU(2)/SU(3): ∇·E^a + g f^{abc} A^b E^c = ρ^a (covariant divergence)
        /// 
        /// This ensures gauge invariance is maintained after numerical evolution.
        /// </summary>
        private void EnforceGaugeConstraintsInternal()
        {
            // Enforce U(1) Gauss law if gauge phases are initialized
            if (_edgePhaseU1 != null)
            {
                EnforceGaussLaw();
            }
            
            // CHECKLIST ITEM 7: Enforce SU(2) Gauss law for weak field
            if (_weakField != null && _weakFieldStrength != null)
            {
                EnforceWeakGaussLaw();
            }
            
            // CHECKLIST ITEM 7: Enforce SU(3) Gauss law for gluon field
            if (_gluonField != null && _gluonFieldStrength != null)
            {
                EnforceStrongGaussLaw();
            }
        }
        
        /// <summary>
        /// Enforce Gauss law constraint for SU(2) weak field.
        /// Iteratively corrects field to satisfy ∇·E^a = ρ^a.
        /// </summary>
        private void EnforceWeakGaussLaw()
        {
            if (_weakField == null || _weakFieldStrength == null || _cachedWeakDensity == null)
                return;
                
            const int maxIterations = 5;
            const double relaxation = 0.1;
            
            for (int iter = 0; iter < maxIterations; iter++)
            {
                double maxViolation = 0.0;
                
                for (int i = 0; i < N; i++)
                {
                    // Compute covariant divergence for each SU(2) component
                    for (int a = 0; a < 3; a++)
                    {
                        double divE = 0.0;
                        double totalDegree = 0.0;
                        
                        foreach (int j in Neighbors(i))
                        {
                            divE += _weakFieldStrength[i, j, a];
                            totalDegree += 1.0;
                        }
                        
                        // Charge density (weak isospin)
                        double rho_a = _cachedWeakDensity[i] * (a == 0 ? 1.0 : 0.0); // Only T^3 component carries charge
                        
                        // Gauss law violation: ∇·E^a - ρ^a
                        double violation = divE - rho_a;
                        maxViolation = Math.Max(maxViolation, Math.Abs(violation));
                        
                        // Correction: subtract gradient of a gauge function
                        if (totalDegree > 0 && Math.Abs(violation) > 1e-10)
                        {
                            double correction = relaxation * violation / totalDegree;
                            foreach (int j in Neighbors(i))
                            {
                                _weakField![i, j, a] -= correction;
                            }
                        }
                    }
                }
                
                // Early exit if converged
                if (maxViolation < 1e-6)
                    break;
            }
        }
        
        /// <summary>
        /// Enforce Gauss law constraint for SU(3) gluon field.
        /// Iteratively corrects field to satisfy ∇·E^a + g f^{abc} A^b E^c = ρ^a.
        /// </summary>
        private void EnforceStrongGaussLaw()
        {
            if (_gluonField == null || _gluonFieldStrength == null || _cachedColorDensity == null)
                return;
                
            const int maxIterations = 5;
            const double relaxation = 0.1;
            
            for (int iter = 0; iter < maxIterations; iter++)
            {
                double maxViolation = 0.0;
                
                for (int i = 0; i < N; i++)
                {
                    // Compute covariant divergence for each SU(3) component
                    for (int a = 0; a < 8; a++)
                    {
                        double divE = 0.0;
                        double totalDegree = 0.0;
                        
                        foreach (int j in Neighbors(i))
                        {
                            double E_ija = _gluonFieldStrength[i, j, a];
                            divE += E_ija;
                            
                            // Non-Abelian correction: f^{abc} A^b E^c
                            for (int b = 0; b < 8; b++)
                            {
                                for (int c = 0; c < 8; c++)
                                {
                                    if (TryGetStructureConstant(a, b, c, out double fabc))
                                    {
                                        divE += StrongCoupling * fabc * _gluonField[i, j, b] * _gluonFieldStrength[i, j, c];
                                    }
                                }
                            }
                            
                            totalDegree += 1.0;
                        }
                        
                        // Color charge density (from wavefunction)
                        double rho_a = _cachedColorDensity[i] * (a < 3 ? 1.0 : 0.0); // First 3 generators carry primary charge
                        
                        // Gauss law violation
                        double violation = divE - rho_a;
                        maxViolation = Math.Max(maxViolation, Math.Abs(violation));
                        
                        // Correction
                        if (totalDegree > 0 && Math.Abs(violation) > 1e-10)
                        {
                            double correction = relaxation * violation / totalDegree;
                            foreach (int j in Neighbors(i))
                            {
                                _gluonField![i, j, a] -= correction;
                            }
                        }
                    }
                }
                
                // Early exit if converged
                if (maxViolation < 1e-6)
                    break;
            }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private double ComputeColorCurrentCached(int i, int j, int a)
        {
            if (_cachedColorDensity == null) return 0.0;
            if ((uint)i >= (uint)_cachedColorDensity.Length ||
                (uint)j >= (uint)_cachedColorDensity.Length) return 0.0;

            double rhoI = _cachedColorDensity[i];
            double rhoJ = _cachedColorDensity[j];

            // Checklist G.1: Background independence - use graph topology, not Coordinates
            // Current flows from high to low density along edge (gradient)
            double densityGradient = rhoJ - rhoI;  // Direction: i -> j
            double edgeWeight = Weights[i, j];

            double colorWeight = (uint)a < (uint)_colorWeights.Length ? _colorWeights[a] : 1.0;
            return StrongCoupling * colorWeight * densityGradient * edgeWeight * PhysicsConstants.GaugeCurrentScaleFactor;
        }

        [MethodImpl(MethodImplOptions.NoInlining)]
        private double ComputeWeakCurrent(int i, int j, int a)
        {
            if (_cachedWeakDensity == null) return 0.0;
            if ((uint)i >= (uint)_cachedWeakDensity.Length ||
                (uint)j >= (uint)_cachedWeakDensity.Length) return 0.0;

            double rhoI = _cachedWeakDensity[i];
            double rhoJ = _cachedWeakDensity[j];

            // Checklist G.1: Background independence - use graph topology, not Coordinates
            // Current flows from high to low density along edge
            double densityGradient = rhoJ - rhoI;  // Direction: i -> j
            double edgeWeight = Weights[i, j];

            double compWeight = 1.0 + 0.05 * a;
            return WeakCoupling * compWeight * densityGradient * edgeWeight * PhysicsConstants.GaugeCurrentScaleFactor;
        }

        [MethodImpl(MethodImplOptions.NoInlining)]
        private double ComputeHyperchargeCurrent(int i, int j)
        {
            double current = 0.0;
            
            // Spinor/fermion contribution
            // Checklist G.1: Use graph topology, NOT external Coordinates for physics
            // Background-independent current: J ∝ ρ_gradient * weight
            if (_cachedHyperDensity != null &&
                (uint)i < (uint)_cachedHyperDensity.Length &&
                (uint)j < (uint)_cachedHyperDensity.Length)
            {
                double hI = _cachedHyperDensity[i];
                double hJ = _cachedHyperDensity[j];

                // Current flows from high to low density along edge
                // Using edge weight as "conductance" - purely relational, no coordinates
                double edgeWeight = Weights[i, j];
                double densityGradient = hJ - hI;  // Direction: i -> j

                current += HypergaugeCoupling * densityGradient * edgeWeight * PhysicsConstants.GaugeCurrentScaleFactor;
            }
            
            // Scalar field back-reaction contribution (Checklist E.2)
            // The scalar field current J_μ = Im[φ* D_μ φ] contributes to Maxwell equations
            // For real scalar: J_ij = g * φ_i * φ_j * sin(θ_ij)
            double scalarCurrent = ComputeScalarFieldCurrent(i, j);
            current += scalarCurrent;
            
            return current;
        }

        public double ElectromagneticField(int i, int j)
        {
            if (_weakField == null || _hyperchargeField == null) return 0;

            double sinThW = Math.Sqrt(0.23);
            double cosThW = Math.Sqrt(1 - 0.23);
            return _weakField[i, j, 2] * sinThW + _hyperchargeField[i, j] * cosThW;
        }

        public double ZBosonField(int i, int j)
        {
            if (_weakField == null || _hyperchargeField == null) return 0;

            double sinThW = Math.Sqrt(0.23);
            double cosThW = Math.Sqrt(1 - 0.23);
            return _weakField[i, j, 2] * cosThW - _hyperchargeField[i, j] * sinThW;
        }

        public (Complex WPlus, Complex WMinus) WBosonFields(int i, int j)
        {
            if (_weakField == null) return (Complex.Zero, Complex.Zero);

            double inv = 1.0 / Math.Sqrt(2);
            Complex wPlus = inv * new Complex(_weakField[i, j, 0], -_weakField[i, j, 1]);
            Complex wMinus = inv * new Complex(_weakField[i, j, 0], _weakField[i, j, 1]);
            return (wPlus, wMinus);
        }

        /// <summary>
        /// RQ-HYPOTHESIS COMPLIANT: Confinement potential using graph-based distance.
        /// 
        /// PHYSICS FIX: Uses -log(weight) as distance metric instead of GetPhysicalDistance.
        /// In RQ-Hypothesis, the metric emerges from correlation strength:
        ///   d(i,j) = -ln(w_ij)
        /// Strong correlations (w → 1) give small distance, weak correlations give large distance.
        /// 
        /// Linear confinement: V(r) = σ·r + const (QCD string tension)
        /// where r is now the graph-derived distance, not external Euclidean distance.
        /// </summary>
        public double ConfinementPotential(int i, int j)
        {
            if (_gluonFieldStrength == null || !Edges[i, j]) return 0;

            double E2 = 0;
            for (int a = 0; a < 8; a++)
                E2 += _gluonFieldStrength[i, j, a] * _gluonFieldStrength[i, j, a];

            double stringTension = 0.2;
            
            // RQ-FIX: Use graph-based distance instead of external coordinates
            // Distance = -ln(weight), where stronger correlation = shorter distance
            // This is the emergent metric from correlation structure
            double w = Weights[i, j];
            double r = w > 1e-10 ? -Math.Log(w + 1e-10) : 10.0; // Clamp for numerical stability
            
            return stringTension * r + E2 * 0.01;
        }

        public double RunningStrongCoupling(double Q)
        {
            int Nf = 6;
            double Lambda = 0.2;

            if (Q <= Lambda) return 1.0;

            double b0 = (33 - 2 * Nf) / (12.0 * Math.PI);
            return 1.0 / (b0 * Math.Log(Q * Q / (Lambda * Lambda)));
        }
        
        /// <summary>
        /// Get U(1) link variable U_ab = exp(i * phase) for edge (nodeA, nodeB).
        /// In lattice gauge theory, links carry the gauge field.
        /// </summary>
        /// <param name="nodeA">Start node</param>
        /// <param name="nodeB">End node</param>
        /// <returns>Complex link variable exp(i*phase)</returns>
        public Complex GetLinkVariable(int nodeA, int nodeB)
        {
            if (_edgePhaseU1 == null || !Edges[nodeA, nodeB])
                return Complex.One;
            
            double phase = _edgePhaseU1[nodeA, nodeB];
            return Complex.FromPolarCoordinates(1.0, phase);
        }
        
        /// <summary>
        /// Calculate plaquette flux (Wilson loop) for triangle (nodeA, nodeB, nodeC).
        /// Returns the product of link variables U_ab * U_bc * U_ca as complex number.
        /// For vacuum, this should be ~1 (trivial holonomy).
        /// Implements RQ-hypothesis checklist item 4.1.
        /// </summary>
        /// <param name="nodeA">First vertex</param>
        /// <param name="nodeB">Second vertex</param>
        /// <param name="nodeC">Third vertex</param>
        /// <returns>Complex Wilson loop flux: U_ab * U_bc * U_ca</returns>
        public Complex CalculatePlaquetteFlux(int nodeA, int nodeB, int nodeC)
        {
            // Product of phases U_ab * U_bc * U_ca
            var U_ab = GetLinkVariable(nodeA, nodeB);
            var U_bc = GetLinkVariable(nodeB, nodeC);
            var U_ca = GetLinkVariable(nodeC, nodeA);

            return U_ab * U_bc * U_ca; // Should be ~1 for vacuum
        }
    }
}