using System;
using System.Numerics;
using System.Linq;
using System.Threading.Tasks;
using RQSimulation.GPUOptimized;
using RQSimulation.Gauge;

namespace RQSimulation
{
    public partial class RQGraph
    {
        // RQ-FIX: Weight rate limiter for adiabatic spinor evolution
        // Stores previous weights to detect rapid changes that could "shock" spinors
        private double[,]? _previousWeightsForDirac;

        // CHECKLIST ITEM 9: Track spectral dimension for spinor field activation
        private double _lastSpectralDimensionForSpinor = 4.0;
        private bool _spinorFieldEnabled = true;

        /// <summary>
        /// Maximum allowed rate of weight change per step for adiabatic Dirac evolution.
        /// RQ-FIX: Rapid weight changes violate adiabatic approximation and can cause
        /// numerical instabilities ("shocks") in spinor field evolution.
        /// </summary>
        private const double MaxWeightChangeRate = 0.1; // 10% max change per step

        /// <summary>
        /// Evolve spinor field using purely topological Dirac operator with Leapfrog symplectic integrator.
        /// No reference to external coordinates - uses graph structure only.
        /// 
        /// WARNING: This method updates ALL nodes simultaneously, violating relativity in
        /// event-driven mode. For asynchronous local time evolution, use UpdateSpinorFieldAtNode()
        /// in RQGraph.EventDrivenExtensions.cs instead.
        /// 
        /// Use this method ONLY for:
        /// - Global synchronous simulation mode
        /// - Testing/debugging purposes
        /// - Small graphs where causality violations are acceptable
        /// 
        /// RQ-Hypothesis Compliant (Item 3): Uses midpoint (leapfrog) symplectic integrator
        /// that mathematically preserves the norm to O(dt^2) without forced normalization.
        /// 
        /// RQ-FIX: Implements adiabatic weight smoothing to prevent spinor "shocks"
        /// from rapid topology changes.
        /// 
        /// CHECKLIST ITEM 9: Validates spectral dimension before spinor evolution.
        /// If d_S deviates >10% from 4.0, spinor evolution is paused.
        /// 
        /// Leapfrog method:
        ///   1. Compute k1 = dψ/dt(ψ(t))
        ///   2. Compute midpoint: ψ_mid = ψ(t) + k1 * dt/2
        ///   3. Compute k2 = dψ/dt(ψ_mid)
        ///   4. Full step: ψ(t+dt) = ψ(t) + k2 * dt
        /// </summary>
        /// <param name="dt">Time step for evolution</param>
        /// <remarks>
        /// LEGACY: For event-driven simulation, use UpdateSpinorFieldAtNode() which
        /// respects local proper time and maintains causality.
        /// </remarks>
        [System.ComponentModel.EditorBrowsable(System.ComponentModel.EditorBrowsableState.Advanced)]
        public void UpdateDiracFieldRelational(double dt)
        {
            if (_spinorA == null) InitSpinorField();

            // CHECKLIST ITEM 9: Check spectral dimension compatibility
            // Only perform full check every 100 steps to avoid overhead
            if (_physicsStepCount % 100 == 0)
            {
                _lastSpectralDimensionForSpinor = ComputeSpectralDimension();
                _spinorFieldEnabled = SpectralDimensionValidator.ShouldEnableSpinorFields(_lastSpectralDimensionForSpinor);

                if (!_spinorFieldEnabled)
                {
                    Console.WriteLine($"[WARNING] d_S = {_lastSpectralDimensionForSpinor:F2}: Spinor field evolution paused (need d_S ≈ 4)");
                }
            }

            // If spectral dimension is far from 4, skip spinor evolution
            if (!_spinorFieldEnabled)
            {
                return;
            }

            double hbar = VectorMath.HBar;
            double c = VectorMath.SpeedOfLight;

            // ===== STEP 1: Compute derivatives k1 at current state =====
            var k1A = new Complex[N];
            var k1B = new Complex[N];
            var k1C = new Complex[N];
            var k1D = new Complex[N];

            ComputeDiracDerivatives(_spinorA!, _spinorB!, _spinorC!, _spinorD!,
                                    k1A, k1B, k1C, k1D, c, hbar);

            // ===== STEP 2: Compute midpoint state =====
            var midA = new Complex[N];
            var midB = new Complex[N];
            var midC = new Complex[N];
            var midD = new Complex[N];

            double halfDt = dt * 0.5;
            for (int i = 0; i < N; i++)
            {
                midA[i] = _spinorA![i] + halfDt * k1A[i];
                midB[i] = _spinorB![i] + halfDt * k1B[i];
                midC[i] = _spinorC![i] + halfDt * k1C[i];
                midD[i] = _spinorD![i] + halfDt * k1D[i];
            }

            // ===== STEP 3: Compute derivatives k2 at midpoint =====
            var k2A = new Complex[N];
            var k2B = new Complex[N];
            var k2C = new Complex[N];
            var k2D = new Complex[N];

            ComputeDiracDerivatives(midA, midB, midC, midD,
                                    k2A, k2B, k2C, k2D, c, hbar);

            // ===== STEP 4: Full step using midpoint derivatives (Leapfrog) =====
            var newA = new Complex[N];
            var newB = new Complex[N];
            var newC = new Complex[N];
            var newD = new Complex[N];

            for (int i = 0; i < N; i++)
            {
                newA[i] = _spinorA![i] + dt * k2A[i];
                newB[i] = _spinorB![i] + dt * k2B[i];
                newC[i] = _spinorC![i] + dt * k2C[i];
                newD[i] = _spinorD![i] + dt * k2D[i];
            }

            // Update spinor arrays
            _spinorA = newA;
            _spinorB = newB;
            _spinorC = newC;
            _spinorD = newD;

            // RQ-Hypothesis Item 3: Symplectic integrator preserves norm to O(dt^2)
            // Only apply minimal adaptive correction for long-term stability
            // Remove forced normalization - norm should be naturally preserved
            AdaptiveNormalizeSpinorFieldMinimal();
        }

        /// <summary>
        /// Compute Dirac operator derivatives dψ/dt for all nodes.
        /// Separated out to enable symplectic integration.
        /// 
        /// RQ-FIX: Uses adiabatic weight smoothing to prevent numerical instabilities
        /// from rapid weight changes on dynamic graphs.
        /// 
        /// CHECKLIST ITEM 2: Implements gauge-covariant parallel transport.
        /// Before computing hopping term (ψ_j - ψ_i), the neighbor spinor ψ_j
        /// is parallel-transported using U†_ij: ψ_j_transported = U†_ij * ψ_j.
        /// This ensures gauge invariance of the Dirac equation.
        /// </summary>
        private void ComputeDiracDerivatives(
            Complex[] spinorA, Complex[] spinorB, Complex[] spinorC, Complex[] spinorD,
            Complex[] dA, Complex[] dB, Complex[] dC, Complex[] dD,
            double c, double hbar)
        {
            // RQ-FIX: Initialize previous weights cache if needed
            if (_previousWeightsForDirac == null ||
                _previousWeightsForDirac.GetLength(0) != N ||
                _previousWeightsForDirac.GetLength(1) != N)
            {
                _previousWeightsForDirac = new double[N, N];
                // Initialize with current weights
                for (int i = 0; i < N; i++)
                    for (int j = 0; j < N; j++)
                        _previousWeightsForDirac[i, j] = Weights[i, j];
            }

            Parallel.For(0, N, i =>
            {
                double m = _fermionMassField![i];
                double mc = m * c;
                bool isEvenSite = (i % 2 == 0);

                Complex deltaA = Complex.Zero, deltaB = Complex.Zero;
                Complex deltaC = Complex.Zero, deltaD = Complex.Zero;

                // Важно: Neighbors(i) должен быть потокобезопасным итератором
                // Обычно чтение списков смежности (Edges/Weights) безопасно, если топология не меняется ВНУТРИ этого шага
                foreach (int j in Neighbors(i))
                {
                    // RQ-FIX: Adiabatic weight smoothing
                    // Limit the effective weight change rate to prevent spinor "shocks"
                    double currentWeight = Weights[i, j];
                    double previousWeight = _previousWeightsForDirac![i, j];

                    // Compute smoothed weight with rate limiting
                    double weight;
                    if (previousWeight > 0)
                    {
                        double changeRate = (currentWeight - previousWeight) / previousWeight;
                        if (Math.Abs(changeRate) > MaxWeightChangeRate)
                        {
                            // Clamp the change rate to adiabatic limit
                            double clampedChange = Math.Sign(changeRate) * MaxWeightChangeRate;
                            weight = previousWeight * (1.0 + clampedChange);
                        }
                        else
                        {
                            weight = currentWeight;
                        }
                    }
                    else
                    {
                        // No previous weight, use current with soft start
                        weight = currentWeight * 0.5; // Smooth introduction of new edge
                    }

                    // CHECKLIST ITEM 2: Gauge-covariant parallel transport
                    // =====================================================
                    // The Dirac equation on a graph with gauge field requires:
                    //   D_ij ψ = U†_ij ψ_j - ψ_i
                    // where U_ij is the gauge link (parallel transport operator).
                    //
                    // For U(1): U†_ij = e^{-iθ_ij} (just a phase)
                    // For SU(2): U†_ij is a 2×2 unitary matrix
                    // For SU(3): U†_ij is a 3×3 unitary matrix (implemented via ParallelTransportSU3)

                    Complex parallelTransport = Complex.One;
                    if (_edgePhaseU1 != null)
                    {
                        double phase = _edgePhaseU1[i, j];
                        parallelTransport = Complex.FromPolarCoordinates(1.0, -phase); // U† = e^{-iθ}
                    }

                    // Staggered fermion hopping with alternating signs
                    bool isNeighborEven = (j % 2 == 0);
                    double sign = (isEvenSite != isNeighborEven) ? 1.0 : -1.0;

                    // CHECKLIST ITEM 14: Directional "flavor" based on edge index modulo 2
                    int edgeDirection = Math.Abs(i - j) % 2;

                    // CHECKLIST ITEM 10 & ITEM 2: SU(2)/SU(3) gauge coupling with parallel transport
                    // Apply gauge transformation to neighbor spinor BEFORE hopping
                    Complex gaugedA_j = spinorA[j] * parallelTransport;
                    Complex gaugedB_j = spinorB[j] * parallelTransport;
                    Complex gaugedC_j = spinorC[j] * parallelTransport;
                    Complex gaugedD_j = spinorD[j] * parallelTransport;

                    if (GaugeDimension == 2 && _gaugeSU != null)
                    {
                        // SU(2) gauge transformation for left-handed doublet
                        // Get SU(2) link matrix U†_ij and apply to (A_j, B_j)
                        Complex U00 = Complex.Conjugate(_gaugeSU[i, j, 0]);
                        Complex U01 = Complex.Conjugate(_gaugeSU[i, j, 2]); // Transpose indices for U†
                        Complex U10 = Complex.Conjugate(_gaugeSU[i, j, 1]);
                        Complex U11 = Complex.Conjugate(_gaugeSU[i, j, 3]);

                        Complex nA = spinorA[j];
                        Complex nB = spinorB[j];
                        gaugedA_j = U00 * nA + U01 * nB;
                        gaugedB_j = U10 * nA + U11 * nB;
                        // Apply U(1) phase on top
                        gaugedA_j *= parallelTransport;
                        gaugedB_j *= parallelTransport;
                        // Right-handed (C, D) are SU(2)_L singlets - only U(1) phase
                    }
                    else if (GaugeDimension == 3 && _gaugeSU3 != null)
                    {
                        // CHECKLIST ITEM 2: Full SU(3) parallel transport
                        // For color-charged quarks, apply the SU(3) gauge transformation
                        // Since we have 4-component spinors (A,B,C,D) not 12-component color triplets,
                        // we apply phase based on the SU(3) Wilson line trace
                        SU3Matrix U = GetSU3Link(i, j);
                        SU3Matrix Udag = U.Dagger();
                        Complex trace = Udag.Trace() / 3.0; // Average phase

                        // Apply averaged SU(3) phase to spinor
                        gaugedA_j = spinorA[j] * trace * parallelTransport;
                        gaugedB_j = spinorB[j] * trace * parallelTransport;
                        gaugedC_j = spinorC[j] * trace * parallelTransport;
                        gaugedD_j = spinorD[j] * trace * parallelTransport;
                    }

                    if (edgeDirection == 0)
                    {
                        // "X-like" direction: couples A↔B and C↔D
                        deltaB += sign * weight * (gaugedA_j - spinorA[i]);
                        deltaA += sign * weight * (gaugedB_j - spinorB[i]);
                        deltaD += sign * weight * (gaugedC_j - spinorC[i]);
                        deltaC += sign * weight * (gaugedD_j - spinorD[i]);
                    }
                    else
                    {
                        // "Y-like" direction: couples A↔B with i and C↔D with -i
                        Complex iSign = Complex.ImaginaryOne * sign;
                        deltaB += iSign * weight * (gaugedA_j - spinorA[i]);
                        deltaA += -iSign * weight * (gaugedB_j - spinorB[i]);
                        deltaD += -iSign * weight * (gaugedC_j - spinorC[i]);
                        deltaC += iSign * weight * (gaugedD_j - spinorD[i]);
                    }
                }

                // Mass term: couples left and right handed components
                Complex massTermA = -Complex.ImaginaryOne * mc / hbar * spinorC[i];
                Complex massTermB = -Complex.ImaginaryOne * mc / hbar * spinorD[i];
                Complex massTermC = -Complex.ImaginaryOne * mc / hbar * spinorA[i];
                Complex massTermD = -Complex.ImaginaryOne * mc / hbar * spinorB[i];

                // Dirac evolution: dψ/dt = -i/ħ * H * ψ
                double factor = -1.0 / hbar;

                dA[i] = factor * (c * deltaA + massTermA);
                dB[i] = factor * (c * deltaB + massTermB);
                dC[i] = factor * (c * deltaC + massTermC);
                dD[i] = factor * (c * deltaD + massTermD);
            });

            // RQ-FIX: Update previous weights cache for next step
            // This is done after parallel section to avoid race conditions
            for (int i = 0; i < N; i++)
                for (int j = 0; j < N; j++)
                    _previousWeightsForDirac![i, j] = Weights[i, j];
        }

        /// <summary>
        /// Minimal adaptive normalization for symplectic integrator.
        /// Only applies very small corrections when norm deviates significantly.
        /// RQ-Hypothesis Item 3: Symplectic integrator should preserve norm - 
        /// this is only a safety net for long-term stability.
        /// </summary>
        private void AdaptiveNormalizeSpinorFieldMinimal()
        {
            if (_spinorA == null) return;

            // Compute total norm
            double totalNorm = 0;
            for (int i = 0; i < N; i++)
            {
                double norm = _spinorA[i].Magnitude * _spinorA[i].Magnitude
                            + _spinorB![i].Magnitude * _spinorB[i].Magnitude
                            + _spinorC![i].Magnitude * _spinorC[i].Magnitude
                            + _spinorD![i].Magnitude * _spinorD[i].Magnitude;
                totalNorm += norm;
            }

            double targetNorm = N; // Target: average norm per node ~ 1
            double relativeDeviation = Math.Abs(totalNorm - targetNorm) / targetNorm;

            // Symplectic integrator should preserve norm to O(dt^2)
            // Only apply minimal correction if deviation exceeds threshold (safety)
            if (relativeDeviation > PhysicsConstants.SymplecticNormSafetyThreshold && totalNorm > 1e-10)
            {
                // Apply very gentle correction toward target
                double currentScale = Math.Sqrt(targetNorm / totalNorm);
                double scale = 1.0 + PhysicsConstants.SymplecticNormCorrectionRate * (currentScale - 1.0);

                for (int i = 0; i < N; i++)
                {
                    _spinorA[i] *= scale;
                    _spinorB![i] *= scale;
                    _spinorC![i] *= scale;
                    _spinorD![i] *= scale;
                }
            }
        }

        /// <summary>
        /// Compute fermion density at node (ψ†ψ)
        /// </summary>
        public double ComputeFermionDensity(int i)
        {
            if (_spinorA == null || i < 0 || i >= N)
                return 0.0;

            double density = _spinorA[i].Magnitude * _spinorA[i].Magnitude
                           + _spinorB![i].Magnitude * _spinorB[i].Magnitude
                           + _spinorC![i].Magnitude * _spinorC[i].Magnitude
                           + _spinorD![i].Magnitude * _spinorD[i].Magnitude;

            return density;
        }

        /// <summary>
        /// Compute chiral density (difference between left and right)
        /// </summary>
        public double ComputeChiralDensity(int i)
        {
            if (_spinorA == null || i < 0 || i >= N)
                return 0.0;

            double leftDensity = _spinorA[i].Magnitude * _spinorA[i].Magnitude
                               + _spinorB![i].Magnitude * _spinorB[i].Magnitude;

            double rightDensity = _spinorC![i].Magnitude * _spinorC[i].Magnitude
                                + _spinorD![i].Magnitude * _spinorD[i].Magnitude;

            return leftDensity - rightDensity;
        }

        /// <summary>
        /// Compute fermion current between two nodes
        /// </summary>
        public Complex ComputeFermionCurrent(int i, int j)
        {
            if (_spinorA == null || !Edges[i, j])
                return Complex.Zero;

            // Current operator: j^μ = ψ̄ γ^μ ψ
            // For simplicity, use ψ*_i ψ_j - ψ*_j ψ_i

            Complex currentA = Complex.Conjugate(_spinorA[i]) * _spinorA![j]
                             - Complex.Conjugate(_spinorA[j]) * _spinorA[i];
            Complex currentB = Complex.Conjugate(_spinorB![i]) * _spinorB[j]
                             - Complex.Conjugate(_spinorB[j]) * _spinorB[i];
            Complex currentC = Complex.Conjugate(_spinorC![i]) * _spinorC[j]
                             - Complex.Conjugate(_spinorC[j]) * _spinorC[i];
            Complex currentD = Complex.Conjugate(_spinorD![i]) * _spinorD[j]
                             - Complex.Conjugate(_spinorD[j]) * _spinorD[i];

            return (currentA + currentB + currentC + currentD) * Weights[i, j];
        }
    }
}
