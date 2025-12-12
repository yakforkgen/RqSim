using System;
using System.Numerics;
using System.Runtime.CompilerServices;

namespace RQSimulation
{
    /// <summary>
    /// Implements relativistic spacetime dynamics based on the RQ hypothesis.
    /// Spacetime emerges from correlation structure; curvature arises from
    /// correlation density gradients.
    /// </summary>
    public partial class RQGraph
    {
        // Spacetime coordinates for each node (t, x, y, z)
        private double[]? _nodeTimeCoord;
        private double[]? _nodeZCoord;  // 3rd spatial dimension (z)

        // 4-velocity components per node
        private double[]? _velocityT;
        private double[]? _velocityX;
        private double[]? _velocityY;
        private double[]? _velocityZ;

        // Spacetime interval accumulator
        private double[]? _properTimeAccum;

        // Causal structure tracking
        private bool[,]? _causallyConnected;
        private double[,]? _lightConeDistance;

        /// <summary>
        /// Initializes relativistic spacetime coordinates from 2D layout.
        /// </summary>
        public void InitSpacetimeCoordinates()
        {
            _nodeTimeCoord = new double[N];
            _nodeZCoord = new double[N];
            _velocityT = new double[N];
            _velocityX = new double[N];
            _velocityY = new double[N];
            _velocityZ = new double[N];
            _properTimeAccum = new double[N];

            // Initialize from existing 2D coordinates if available
            for (int i = 0; i < N; i++)
            {
                _nodeTimeCoord[i] = 0.0;  // Start at t=0
                _nodeZCoord[i] = 0.0;     // z=0 initially

                // Initial 4-velocity is purely timelike (at rest)
                _velocityT[i] = 1.0;
                _velocityX[i] = 0.0;
                _velocityY[i] = 0.0;
                _velocityZ[i] = 0.0;

                _properTimeAccum[i] = 0.0;
            }

            InitCausalStructure();
        }

        /// <summary>
        /// Initialize causal structure matrices.
        /// </summary>
        private void InitCausalStructure()
        {
            _causallyConnected = new bool[N, N];
            _lightConeDistance = new double[N, N];

            // Initially, adjacent nodes are potentially causally connected
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    if (i == j)
                    {
                        _causallyConnected[i, j] = true;
                        _lightConeDistance[i, j] = 0.0;
                    }
                    else if (Edges[i, j])
                    {
                        _causallyConnected[i, j] = true;
                        double d = GetPhysicalDistance(i, j);
                        _lightConeDistance[i, j] = d / VectorMath.SpeedOfLight;
                    }
                }
            }
        }

        /// <summary>
        /// Computes the spacetime interval squared between two nodes.
        /// Returns negative for timelike, positive for spacelike, zero for null.
        /// 
        /// RQ-Hypothesis Compliant: Uses topological graph distance instead of 
        /// coordinate-based distance. The spatial interval is computed from the
        /// weighted shortest path through the graph, respecting background independence.
        /// 
        /// ds² = -c²dt² + dx² where dx = GetGraphDistanceWeighted(i, j)
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public double SpacetimeIntervalSquared(int i, int j)
        {
            if (_nodeTimeCoord == null)
                return 0.0;

            // Validate indices
            if (i < 0 || i >= N || j < 0 || j >= N)
                return 0.0;

            // RQ-Hypothesis: Use topological graph distance instead of coordinate difference
            // This removes dependence on external Euclidean background space
            double graphDistance = GetGraphDistanceWeighted(i, j);
            double dx_squared = graphDistance * graphDistance;

            double dt = _nodeTimeCoord[j] - _nodeTimeCoord[i];

            // Minkowski metric: s² = -c²dt² + dx²
            double c = VectorMath.SpeedOfLight;
            return -c * c * dt * dt + dx_squared;
        }

        /// <summary>
        /// Checks if two nodes are causally connected (within each other's light cones).
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public bool AreCausallyConnected(int i, int j)
        {
            if (_causallyConnected == null) return true;
            return _causallyConnected[i, j];
        }

        /// <summary>
        /// Updates causal structure based on current spacetime coordinates.
        /// RQ-Hypothesis Compliant: Uses topological graph distance for causality checks.
        /// </summary>
        public void UpdateCausalStructure()
        {
            if (_nodeTimeCoord == null) return;

            // Внешний цикл параллелим
            Parallel.For(0, N, i =>
            {
                // Внутренний цикл последовательный (или тоже можно распараллелить, но overhead может быть выше)
                for (int j = 0; j < N; j++)
                {
                    if (i == j) continue;

                    // Чтение координат/времени (Read-Only)
                    double ds2 = SpacetimeIntervalSquared(i, j);
                    bool causal = ds2 <= 0;

                    // Запись в матрицу bool (атомарная операция для байта)
                    // Разные i,j гарантируют отсутствие гонок
                    _causallyConnected![i, j] = causal;

                    if (causal)
                    {
                        double dist = GetGraphDistanceWeighted(i, j); // Тяжелая операция поиска пути?
                                                                      // Если GetGraphDistanceWeighted использует кэшированные пути, то ок. 
                                                                      // Если запускает Dijkstra - это будет ОЧЕНЬ тяжело, но параллелизм тут спасет.
                        _lightConeDistance![i, j] = Math.Abs(_nodeTimeCoord[j] - _nodeTimeCoord[i]);
                    }
                }
            });
        }
        /// <summary>
        /// Computes enhanced gravitational time dilation for a node using Schwarzschild metric.
        /// This provides a more accurate version than the simple one in Statistics.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public double GetSchwarzschildTimeDilation(int node)
        {
            if (_correlationMass == null || node >= _correlationMass.Length)
                return 1.0;

            double localMass = _correlationMass[node];
            double avgMass = _avgCorrelationMass > 0 ? _avgCorrelationMass : 1.0;

            // Effective "radius" from average link distance
            double r = 1.0;
            int count = 0;
            foreach (int nb in Neighbors(node))
            {
                double d = GetPhysicalDistance(node, nb);
                if (d > 0) { r += d; count++; }
            }
            if (count > 0) r = r / count;
            if (r < 0.01) r = 0.01;

            // Schwarzschild-like time dilation
            double rs = 2.0 * localMass / (VectorMath.SpeedOfLight * VectorMath.SpeedOfLight);
            double factor = 1.0 - rs / r;
            return factor > 0 ? Math.Sqrt(factor) : 0.01;
        }

        /// <summary>
        /// Advances proper time for all nodes accounting for gravitational time dilation.
        /// </summary>
        public void AdvanceProperTime(double globalDt)
        {
            if (_properTimeAccum == null || _nodeTimeCoord == null)
                return;

            for (int i = 0; i < N; i++)
            {
                double dilation = GetGravitationalTimeDilation(i);
                double localDt = globalDt * dilation;

                _properTimeAccum[i] += localDt;
                _nodeTimeCoord[i] += globalDt;  // Coordinate time advances uniformly
            }
        }

        /// <summary>
        /// Applies a Lorentz boost to a node's 4-velocity.
        /// </summary>
        public void ApplyLorentzBoost(int node, double vx, double vy, double vz)
        {
            if (_velocityT == null || _velocityX == null || _velocityY == null || _velocityZ == null)
                return;

            var (tNew, xNew, yNew, zNew) = VectorMath.LorentzBoost(
                _velocityT[node], _velocityX[node], _velocityY[node], _velocityZ[node],
                vx, vy, vz);

            _velocityT[node] = tNew;
            _velocityX[node] = xNew;
            _velocityY[node] = yNew;
            _velocityZ[node] = zNew;

            // Renormalize to unit 4-velocity
            Normalize4Velocity(node);
        }

        /// <summary>
        /// Normalizes the 4-velocity of a node to have proper length -c^2.
        /// </summary>
        private void Normalize4Velocity(int node)
        {
            if (_velocityT == null) return;

            double c2 = VectorMath.SpeedOfLight * VectorMath.SpeedOfLight;
            double spatial = _velocityX![node] * _velocityX[node] +
                           _velocityY![node] * _velocityY[node] +
                           _velocityZ![node] * _velocityZ[node];

            // u · u = -c^2 => u_t = sqrt(c^2 + spatial)
            double newT = Math.Sqrt(c2 + spatial);
            double scale = newT / Math.Max(1e-10, Math.Abs(_velocityT[node]));

            _velocityT[node] = newT;
            // Spatial components remain unchanged in magnitude
        }

        /// <summary>
        /// Computes geodesic deviation (relative acceleration magnitude) between two nearby nodes
        /// due to spacetime curvature.
        /// RQ-Hypothesis Compliant: Uses topological graph distance instead of coordinates.
        /// Returns scalar magnitude of tidal acceleration.
        /// </summary>
        public double ComputeGeodesicDeviationScalar(int i, int j)
        {
            if (_correlationMass == null)
                return 0.0;

            // Topological distance
            double r = GetGraphDistanceWeighted(i, j);
            if (r < 1e-10) return 0.0;

            // Curvature tensor component approximation from mass gradient
            // In RQ, mass is correlation density, so gradient implies curvature gradient
            double mi = _correlationMass[i];
            double mj = _correlationMass[j];
            double massGradient = (mj - mi) / r;

            // Tidal acceleration magnitude ~ R * separation
            // a = -R * d
            // Here we approximate R ~ massGradient (Newtonian-like tidal force gradient)
            double curvature = massGradient * 0.1; 
            double a = curvature * r;

            return a;
        }

        // Tolerance for light cone constraint check (accounts for numerical precision and lattice discretization)
        private const double LightConeTolerance = 0.1;

        /// <summary>
        /// Propagates signals respecting the light cone constraint.
        /// Updates node states only if causally connected.
        /// </summary>
        public void PropagateSignalsWithCausality()
        {
            if (_causallyConnected == null) return;

            var newStates = new NodeState[N];
            Array.Copy(State, newStates, N);

            for (int i = 0; i < N; i++)
            {
                if (State[i] != NodeState.Rest) continue;

                int excitedNeighborCount = 0;
                double weightSum = 0.0;

                foreach (int j in Neighbors(i))
                {
                    // Only propagate if causally connected
                    if (!_causallyConnected[i, j]) continue;
                    if (State[j] != NodeState.Excited) continue;

                    // Check light cone constraint
                    double dt = _nodeTimeCoord != null ? Math.Abs(_nodeTimeCoord[i] - _nodeTimeCoord[j]) : 0;
                    double dr = GetPhysicalDistance(i, j);
                    if (dr > VectorMath.SpeedOfLight * dt + LightConeTolerance) continue;

                    excitedNeighborCount++;
                    weightSum += Weights[i, j];
                }

                if (excitedNeighborCount > 0)
                {
                    double prob = weightSum / (1.0 + weightSum);
                    if (_rng.NextDouble() < prob)
                    {
                        newStates[i] = NodeState.Excited;
                    }
                }
            }

            Array.Copy(newStates, State, N);
        }

        /// <summary>
        /// Computes the Ricci tensor component R_00 for a node based on local correlation structure.
        /// </summary>
        public double ComputeLocalRicci00(int node)
        {
            if (_correlationMass == null) return 0.0;

            // R_00 ≈ 4πG * (ρ + 3p/c^2) in GR
            // Here ρ is the local correlation mass density
            double localMass = _correlationMass[node];
            int deg = Degree(node);
            if (deg == 0) return 0.0;

            // Estimate local volume from average link distance
            double avgDist = 0.0;
            foreach (int nb in Neighbors(node))
            {
                avgDist += GetPhysicalDistance(node, nb);
            }
            avgDist /= deg;

            double localVolume = Math.PI * avgDist * avgDist;  // Approximate as disk
            double density = localMass / Math.Max(1e-10, localVolume);

            // Include pressure estimate from edge tensions
            double pressure = 0.0;
            foreach (int nb in Neighbors(node))
            {
                double w = Weights[node, nb];
                pressure += w * w;
            }
            pressure /= deg;

            double c2 = VectorMath.SpeedOfLight * VectorMath.SpeedOfLight;
            return 4.0 * Math.PI * (density + 3.0 * pressure / c2);
        }

        /// <summary>
        /// Updates coordinates based on Einstein's field equations (simplified).
        /// The metric responds to energy-momentum distribution.
        /// </summary>
        public void EvolveMetricFromEinstein(double dt)
        {
            if (Coordinates == null || _correlationMass == null) return;

            var newCoords = new (double X, double Y)[N];
            Array.Copy(Coordinates, newCoords, N);

            // Compute Ricci tensor and evolve metric
            for (int i = 0; i < N; i++)
            {
                double R00 = ComputeLocalRicci00(i);

                // Metric perturbation from curvature
                double h = R00 * dt * 0.01;

                // Move node towards average position weighted by curvature
                double sumX = 0, sumY = 0;
                int count = 0;
                foreach (int nb in Neighbors(i))
                {
                    double d = GetPhysicalDistance(i, nb);
                    double w = Weights[i, nb];
                    sumX += Coordinates[nb].X * w;
                    sumY += Coordinates[nb].Y * w;
                    count++;
                }

                if (count > 0)
                {
                    double avgX = sumX / count;
                    double avgY = sumY / count;

                    // Geodesic motion: move slightly towards weighted center
                    newCoords[i].X += h * (avgX - Coordinates[i].X);
                    newCoords[i].Y += h * (avgY - Coordinates[i].Y);
                }
            }

            for (int i = 0; i < N; i++)
            {
                Coordinates[i] = newCoords[i];
            }
        }
    }
}
