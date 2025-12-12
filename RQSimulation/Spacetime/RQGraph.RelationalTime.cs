using System;
using System.Linq;
using System.Numerics;

namespace RQSimulation
{
    /// <summary>
    /// Relational Time (Page-Wootters Mechanism) - Extended Implementation
    /// Implements fully relational dynamics without external time parameter dt
    /// 
    /// RQ-HYPOTHESIS PHYSICS:
    /// =======================
    /// The "lapse function" N(x) in ADM formalism controls how fast proper time
    /// flows at each point relative to coordinate time:
    ///   dτ = N(x) × dt
    /// 
    /// In a relational graph:
    ///   N_i = 1 / sqrt(1 + |R_i| + m_i/⟨m⟩)
    /// 
    /// where R_i is local Ricci scalar and m_i is local mass.
    /// Higher curvature/mass → slower time (gravitational time dilation).
    /// </summary>
    public partial class RQGraph
    {
        private const int DefaultClockSizeDivisor = 20; // Clock uses 1/20 (5%) of nodes
        private const double MinRelationalDt = 0.001;
        private const double MaxRelationalDt = 0.1;
        
        private int[] _clockNodesArray = Array.Empty<int>();
        private Complex[] _lastClockStateVector = Array.Empty<Complex>();
        
        /// <summary>
        /// Cached lapse function values for performance.
        /// Updated by UpdateLapseFunctions().
        /// </summary>
        private double[]? _lapseFunction;

        /// <summary>
        /// Initialize internal clock subsystem based on connectivity
        /// </summary>
        public void InitInternalClock(int clockSize)
        {
            if (clockSize <= 0 || clockSize > N)
                clockSize = Math.Max(2, N / DefaultClockSizeDivisor);

            // Select clock nodes based on high connectivity (hub nodes)
            var nodesByDegree = Enumerable.Range(0, N)
                .OrderByDescending(i => Neighbors(i).Count())
                .Take(clockSize)
                .ToArray();

            _clockNodesArray = nodesByDegree;
            
            // Mark as clock nodes in physics properties
            foreach (int idx in _clockNodesArray)
            {
                if (PhysicsProperties != null && PhysicsProperties.Length == N)
                    PhysicsProperties[idx].IsClock = true;
            }

            // Initialize clock state vector
            _lastClockStateVector = new Complex[clockSize];
            SaveClockState();
            
            // Initialize lapse function
            _lapseFunction = new double[N];
            UpdateLapseFunctions();
        }

        /// <summary>
        /// Save current quantum state of clock nodes
        /// </summary>
        public void SaveClockState()
        {
            if (_clockNodesArray.Length == 0 || _waveMulti == null)
                return;

            int d = GaugeDimension;
            for (int k = 0; k < _clockNodesArray.Length; k++)
            {
                int i = _clockNodesArray[k];
                if (i * d < _waveMulti.Length)
                {
                    // Store the first component of the gauge field
                    _lastClockStateVector[k] = _waveMulti[i * d];
                }
            }
        }

        /// <summary>
        /// Compute relational time increment based on clock state change (Fubini-Study metric)
        /// Extended version with explicit clock state tracking
        /// </summary>
        public double ComputeRelationalDtExtended()
        {
            if (_clockNodesArray.Length == 0 || _waveMulti == null)
                return 0.01; // Fallback to small fixed dt

            int d = GaugeDimension;
            double distanceSquared = 0.0;
            double normalization = 0.0;

            // Compute Fubini-Study distance between current and last clock state
            for (int k = 0; k < _clockNodesArray.Length; k++)
            {
                int i = _clockNodesArray[k];
                if (i * d >= _waveMulti.Length)
                    continue;

                Complex current = _waveMulti[i * d];
                Complex last = _lastClockStateVector[k];

                // |ψ_current - ψ_last|^2
                Complex diff = current - last;
                double diffMagnitudeSquared = diff.Real * diff.Real + diff.Imaginary * diff.Imaginary;
                
                distanceSquared += diffMagnitudeSquared;
                
                // Normalization factor
                double currentNorm = current.Magnitude;
                double lastNorm = last.Magnitude;
                normalization += currentNorm + lastNorm;
            }

            if (normalization < 1e-10)
                return 0.01; // Fallback

            // Fubini-Study metric: arccos(|<ψ|φ>|)
            // Approximation: dt ≈ sqrt(distance) / normalization
            double dt = Math.Sqrt(distanceSquared) / (normalization + 1e-10);
            
            // Clamp to reasonable range
            dt = Math.Clamp(dt, MinRelationalDt, MaxRelationalDt);
            
            return dt;
        }

        /// <summary>
        /// Advance internal clock by updating saved state
        /// </summary>
        public void AdvanceInternalClockState()
        {
            SaveClockState();
        }
        
        // ================================================================
        // RQ-HYPOTHESIS: LOCAL LAPSE FUNCTION (ADM FORMALISM)
        // ================================================================
        
        /// <summary>
        /// Get the local lapse function N_i for node i.
        /// 
        /// PHYSICS: The lapse function controls how proper time relates to
        /// coordinate time: dτ_i = N_i × dt.
        /// 
        /// In General Relativity, N = sqrt(-g_00) which encodes gravitational
        /// time dilation. On a graph:
        /// 
        ///   N_i = 1 / sqrt(1 + |R_i|/R_scale + m_i/m_scale)
        /// 
        /// where:
        /// - R_i = local Ricci scalar (average of edge curvatures)
        /// - m_i = local mass (from correlation or fields)
        /// - R_scale, m_scale = characteristic scales for normalization
        /// 
        /// Implements RQ-Hypothesis Checklist: Local Lapse Function (Step 3).
        /// </summary>
        /// <param name="node">Node index</param>
        /// <returns>Lapse function N_i ∈ (0, 1]</returns>
        public double GetLocalLapse(int node)
        {
            if (node < 0 || node >= N)
                return 1.0;
            
            // Use cached value if available
            if (_lapseFunction != null && node < _lapseFunction.Length)
            {
                return _lapseFunction[node];
            }
            
            // Compute on-demand
            return ComputeLocalLapseUncached(node);
        }
        
        /// <summary>
        /// Compute lapse function without using cache.
        /// </summary>
        private double ComputeLocalLapseUncached(int node)
        {
            // === Local Ricci scalar ===
            double R_local = Math.Abs(GetLocalCurvature(node));
            
            // === Local mass ===
            double m_local = 0.0;
            if (_correlationMass != null && node < _correlationMass.Length)
            {
                m_local = _correlationMass[node];
            }
            
            // === Characteristic scales ===
            // R_scale = average curvature magnitude in graph
            // m_scale = average mass
            double R_scale = Math.Max(0.1, _avgCurvature);
            double m_scale = Math.Max(0.1, _avgCorrelationMass);
            
            // === Lapse function formula ===
            // N = 1 / sqrt(1 + |R|/R_0 + m/m_0)
            // This ensures N → 0 as curvature/mass → ∞ (black hole limit)
            // and N → 1 in flat vacuum
            double denominator = 1.0 + R_local / R_scale + m_local / m_scale;
            double N = 1.0 / Math.Sqrt(denominator);
            
        // Clamp to prevent numerical issues
            return Math.Clamp(N, PhysicsConstants.MinTimeDilation, PhysicsConstants.MaxTimeDilation);
        }
        
        /// <summary>
        /// Update all lapse function values.
        /// Should be called after topology or mass changes.
        /// </summary>
        public void UpdateLapseFunctions()
        {
            if (_lapseFunction == null || _lapseFunction.Length != N)
            {
                _lapseFunction = new double[N];
            }
            
            // First pass: compute average curvature and mass for scaling
            UpdateAverageCurvature();
            
            // Second pass: compute lapse for each node
            Parallel.For(0, N, i =>
            {
                _lapseFunction[i] = ComputeLocalLapseUncached(i);
            });
        }
        
        /// <summary>
        /// Cached average curvature for lapse computation.
        /// </summary>
        private double _avgCurvature = 0.1;
        
        /// <summary>
        /// Update cached average curvature.
        /// </summary>
        private void UpdateAverageCurvature()
        {
            double sum = 0.0;
            int count = 0;
            
            for (int i = 0; i < N; i++)
            {
                double R = Math.Abs(GetLocalCurvature(i));
                sum += R;
                count++;
            }
            
            _avgCurvature = count > 0 ? sum / count : 0.1;
            if (_avgCurvature < 0.01)
                _avgCurvature = 0.1; // Prevent division issues
        }
        
        // NOTE: ComputeLocalProperTime, GetTimeDilation, and UpdateNodePhysics
        // are defined in RQGraph.AsynchronousTime.cs to avoid duplication.
        // The lapse function provides the ADM-style N_i = 1/sqrt(1 + |R| + m).
    }
}
