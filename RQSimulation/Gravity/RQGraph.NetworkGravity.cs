using System;
using System.Collections.Generic;
using System.Linq;

namespace RQSimulation
{
    /// <summary>
    /// Network Gravity: Geometry evolution from graph curvature
    /// Implements background-independent gravitational dynamics on the network
    /// </summary>
    public partial class RQGraph
    {
        // Note: Constants moved to PhysicsConstants for centralized configuration
        // - GravitationalCoupling: 2.0 (increased from 0.5)
        // - CurvatureTermScale: 0.05
        // - CosmologicalConstant: 0.0001 (reduced from 0.002)
        
        /// <summary>
        /// Calculate graph curvature for edge (i,j) using Ollivier-Ricci or Forman-Ricci approximation
        /// Positive curvature indicates clustering, negative indicates tree-like structure
        /// </summary>
        public double CalculateGraphCurvature(int i, int j)
        {
            // RQ-HYPOTHESIS FIX (Item 5): Use Ollivier-Ricci curvature by default
            // Forman-Ricci is too coarse for deformed lattices.
            // We delegate to the GPU-optimized implementation (which runs on CPU if GPU not avail).
            return GPUOptimized.OllivierRicciCurvature.ComputeOllivierRicciJaccard(this, i, j);
        }

        /// <summary>
        /// Calculate local volume (weighted degree) of a node.
        /// Used for volume constraint in gravity.
        /// </summary>
        public double GetLocalVolume(int i)
        {
            double vol = 0.0;
            foreach (int j in Neighbors(i))
            {
                vol += Weights[i, j];
            }
            return vol;
        }

        /// <summary>
        /// Evolve network geometry based on curvature and stress-energy (mass).
        /// Replaces coordinate-based EvolveMetricFromEinstein with relational dynamics.
        /// 
        /// RQ-Hypothesis Compliant (Item 4): Uses gradient descent on action S.
        /// dW/dt = -? * dS/dW where S = S_geometry + S_matter
        /// S_geometry = ? R_ij * w_ij (discrete Einstein-Hilbert)
        /// S_matter = stress-energy contribution from correlation mass
        /// 
        /// GPU Modes:
        /// - Full GPU: Curvature AND gravity computed on GPU (fastest, requires topology init)
        /// - Hybrid GPU: Curvature on CPU, gravity on GPU
        /// - CPU: All computation on CPU (fallback)
        /// </summary>
        public void EvolveNetworkGeometry(double dt)
        {
            if (Weights == null || Edges == null)
                return;

            // GPU acceleration if available
            if (GpuGravity != null)
            {
                // Check if full GPU mode is available (topology buffers initialized)
                if (GpuGravity.IsTopologyInitialized)
                {
                    // FULL GPU MODE: Curvature AND gravity computed on GPU
                    // This eliminates the CPU bottleneck entirely!
                    float[] weights = GetAllWeightsFlat();
                    float[] masses = GetNodeMasses();
                    int[] edgesFrom = FlatEdgesFrom;
                    int[] edgesTo = FlatEdgesTo;
                    
                    GpuGravity.EvolveFullGpuStep(
                        weights,
                        masses,
                        edgesFrom,
                        edgesTo,
                        (float)dt,
                        (float)PhysicsConstants.GravitationalCoupling,
                        (float)PhysicsConstants.CosmologicalConstant,
                        (float)PhysicsConstants.DegreePenaltyFactor
                    );
                    
                    UpdateWeightsFromFlat(weights);
                    UpdateTargetDistancesFromWeights();
                    return;
                }
                else
                {
                    // HYBRID GPU MODE: Curvature on CPU, gravity on GPU
                    // Falls back to this when topology buffers not initialized
                    float[] weights = GetAllWeightsFlat();
                    float[] curvatures = GetAllCurvaturesFlat(); // CPU bottleneck
                    float[] masses = GetNodeMasses();
                    int[] edgesFrom = FlatEdgesFrom;
                    int[] edgesTo = FlatEdgesTo;
                    
                    GpuGravity.EvolveGravityGpu(
                        weights,
                        curvatures,
                        masses,
                        edgesFrom,
                        edgesTo,
                        (float)dt,
                        (float)PhysicsConstants.GravitationalCoupling,
                        (float)PhysicsConstants.CosmologicalConstant
                    );
                    
                    UpdateWeightsFromFlat(weights);
                    UpdateTargetDistancesFromWeights();
                    return;
                }
            }

            // CPU implementation with updated constants from PhysicsConstants
            // RQ-HYPOTHESIS FIX: Use unified NodeMasses instead of just correlation mass
            // This ensures gravity responds to ALL field contributions (Einstein equations)
            UpdateNodeMasses();

            double[,] deltaWeights = new double[N, N];
            // Use updated GravitationalCoupling from PhysicsConstants
            double learningRate = PhysicsConstants.GravitationalCoupling * dt;

            // Parallelize outer loop, process j>i only to avoid races
            System.Threading.Tasks.Parallel.For(0, N, i =>
            {
                // RQ-FIX: Use GetNodeTotalMass from UnifiedMass.cs which includes:
                // - Fermion mass (Dirac spinors)
                // - Scalar field energy (Higgs)
                // - Gauge field energy (photons, gluons)
                // - Correlation mass (topological)
                // - Vacuum energy (cosmological)
                double massI = GetNodeTotalMass(i);
                double volI = GetLocalVolume(i);

                foreach (int j in Neighbors(i))
                {
                    if (j <= i) continue;

                    double massJ = GetNodeTotalMass(j);
                    double volJ = GetLocalVolume(j);
                    
                    double dS_geometry;
                    if (PhysicsConstants.PreferOllivierRicciCurvature)
                    {
                        dS_geometry = CalculateOllivierRicciCurvature(i, j);
                    }
                    else
                    {
                        dS_geometry = CalculateGraphCurvature(i, j);
                    }

                    // Stress-energy tensor now includes ALL field contributions
                    double stressEnergyTensor = (massI + massJ) * 0.5;
                    double dS_matter = -stressEnergyTensor * PhysicsConstants.CurvatureTermScale;
                    double dS_cosmological = PhysicsConstants.CosmologicalConstant;
                    
                    // Volume constraint (CDT-like stabilization)
                    // Penalize deviation from target volume (valence)
                    double volConstraintI = volI - PhysicsConstants.TargetVolume;
                    double volConstraintJ = volJ - PhysicsConstants.TargetVolume;
                    double dS_volume = PhysicsConstants.VolumeConstraintLambda * (volConstraintI + volConstraintJ) * 0.5;

                    double dS_total = dS_geometry + dS_matter + dS_cosmological + dS_volume;

                    double delta = -learningRate * dS_total;

                    // Safe write for unique [i,j] pair and symmetric [j,i]
                    deltaWeights[i, j] = delta;
                    deltaWeights[j, i] = delta;
                }
            });

            // Apply weight updates with bounds checking (parallelize outer loop, j>i)
            System.Threading.Tasks.Parallel.For(0, N, i =>
            {
                foreach (int j in Neighbors(i))
                {
                    if (j <= i) continue;

                    double newWeight = Weights[i, j] + deltaWeights[i, j];
                    
                    // Soft clamp for bridge edges (heuristic: if weight is low, don't let it drop to zero easily)
                    // This prevents fragmentation
                    if (newWeight < 0.05 && Weights[i, j] > 0.05)
                    {
                        // Resistance to breaking weak links
                        newWeight = Math.Max(newWeight, 0.01);
                    }

                    newWeight = Math.Clamp(newWeight, 0.0, 1.0);
                    Weights[i, j] = newWeight;
                    Weights[j, i] = newWeight;
                }
            });

            UpdateTargetDistancesFromWeights();
        }

        /// <summary>
        /// Compute average curvature of the network
        /// </summary>
        public double ComputeAverageCurvature()
        {
            if (Edges == null || Weights == null)
                return 0.0;

            double totalCurvature = 0.0;
            int edgeCount = 0;

            for (int i = 0; i < N; i++)
            {
                foreach (int j in Neighbors(i))
                {
                    if (j <= i) continue; // Count each edge once

                    totalCurvature += CalculateGraphCurvature(i, j);
                    edgeCount++;
                }
            }

            return edgeCount > 0 ? totalCurvature / edgeCount : 0.0;
        }

        /// <summary>
        /// Compute curvature scalar (sum of all edge curvatures)
        /// Analogous to Ricci scalar R in GR
        /// </summary>
        public double ComputeCurvatureScalar()
        {
            if (Edges == null || Weights == null)
                return 0.0;

            double scalar = 0.0;

            for (int i = 0; i < N; i++)
            {
                foreach (int j in Neighbors(i))
                {
                    if (j <= i) continue;

                    scalar += CalculateGraphCurvature(i, j);
                }
            }

            return scalar;
        }
        
        /// <summary>
        /// Evolve geometry for a single edge using gradient descent on Hamiltonian.
        /// The weight change is driven by local curvature: dw/dt ? R_ij
        /// Positive curvature increases weight (contracts space), negative decreases.
        /// Implements checklist item 3.2: Geodesic geodinamics via gradient descent.
        /// </summary>
        /// <param name="i">First node of edge</param>
        /// <param name="j">Second node of edge</param>
        /// <param name="learningRate">Step size for gradient descent</param>
        public void EvolveGeometry(int i, int j, double learningRate)
        {
            if (!Edges[i, j] || i == j)
                return;
            
            // dH/dw ? -R_ij (curvature dictates metric change)
            // Positive curvature ? increase weight (contract space, strengthen link)
            // Negative curvature ? decrease weight (expand space, weaken link)
            // RQ-FIX: Use CalculateGraphCurvature which delegates to Ollivier-Ricci if configured
            double curvature = CalculateGraphCurvature(i, j);
            
            // Weight update: w += ? * R (gradient ascent on curvature)
            double newWeight = Weights[i, j] + learningRate * curvature;
            
            // Clamp to valid range [0, 1]
            newWeight = Math.Clamp(newWeight, 0.0, 1.0);
            
            Weights[i, j] = newWeight;
            Weights[j, i] = newWeight; // Symmetric
        }
        
        /// <summary>
        /// Evolve entire graph geometry using curvature-driven gradient descent.
        /// All edges are updated according to their local Forman-Ricci curvature.
        /// Implements checklist item 3.2: Full geometrodynamics step.
        /// </summary>
        /// <param name="learningRate">Step size for gradient descent</param>
        public void EvolveGeometryFull(double learningRate)
        {
            if (Weights == null || Edges == null)
                return;
            
            // Buffer for weight updates (to avoid order-dependent artifacts)
            double[,] deltaWeights = new double[N, N];
            
            for (int i = 0; i < N; i++)
            {
                foreach (int j in Neighbors(i))
                {
                    if (j <= i) continue; // Each edge once
                    
                    double curvature = ComputeFormanRicciCurvature(i, j);
                    deltaWeights[i, j] = learningRate * curvature;
                    deltaWeights[j, i] = deltaWeights[i, j];
                }
            }
            
            // Apply updates
            for (int i = 0; i < N; i++)
            {
                foreach (int j in Neighbors(i))
                {
                    if (j <= i) continue;
                    
                    double newWeight = Weights[i, j] + deltaWeights[i, j];
                    newWeight = Math.Clamp(newWeight, 0.0, 1.0);
                    Weights[i, j] = newWeight;
                    Weights[j, i] = newWeight;
                }
            }
            
            // Update derived quantities
            UpdateTargetDistancesFromWeights();
        }
        
        /// <summary>
        /// Calculate Ollivier-Ricci curvature for edge (i,j)
        /// This is the NEW implementation (CHECKLIST ITEM 4)
        /// More sensitive to geometry than Forman-Ricci
        /// </summary>
        public double CalculateOllivierRicciCurvature(int i, int j)
        {
            // Delegate to GPU-optimized implementation
            return GPUOptimized.OllivierRicciCurvature.ComputeOllivierRicciJaccard(this, i, j);
        }
        
        // === Checklist B.1 & F.1: Canonical Momenta for Geometry ===
        // Edge weight momenta ?_ij for Hamiltonian dynamics
        private double[,]? _geometryMomenta;
        
        /// <summary>
        /// Public accessor for edge momenta ?_ij (Checklist F.1).
        /// These are the canonical conjugate momenta to edge weights.
        /// K = ? ??/(2M) is the geometric kinetic energy.
        /// </summary>
        public double[,]? EdgeMomenta => _geometryMomenta;
        
        /// <summary>
        /// Initialize canonical momenta for edge weights.
        /// This enables inertial geometry dynamics and gravitational waves.
        /// Implements RQ-Hypothesis Checklist B.1 and F.1.
        /// </summary>
        public void InitGeometryMomenta()
        {
            if (N <= 0) return;
            
            _geometryMomenta = new double[N, N];
            // Initialize to zero (static geometry initially)
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    _geometryMomenta[i, j] = 0.0;
                }
            }
        }
        
        /// <summary>
        /// Check if geometry momenta are initialized.
        /// </summary>
        public bool HasGeometryMomenta => _geometryMomenta != null && _geometryMomenta.GetLength(0) == N;
        
        /// <summary>
        /// Evolve geometry with canonical momenta using Velocity Verlet integrator.
        /// This is a symplectic integrator that preserves the Hamiltonian structure
        /// and provides better energy conservation than Euler methods.
        /// 
        /// Velocity Verlet (Checklist F.2):
        ///   Step 1: w(t+?dt) = w(t) + ?(t) * dt/(2M)
        ///   Step 2: F = -?H_total/?w  (forces from curvature and matter)
        ///   Step 3: ?(t+dt) = ?(t) + F * dt
        ///   Step 4: w(t+dt) = w(t+?dt) + ?(t+dt) * dt/(2M)
        /// 
        /// This gives geometry inertia and enables gravitational waves.
        /// Implements RQ-Hypothesis Checklist B.1 and F.2.
        /// </summary>
        /// <param name="dt">Time step</param>
        public void EvolveGeometryHamiltonian(double dt)
        {
            if (!HasGeometryMomenta || Weights == null || Edges == null)
            {
                // Fall back to non-Hamiltonian evolution
                EvolveNetworkGeometry(dt);
                return;
            }
            
            // Compute correlation mass (stress-energy) if not available
            if (_correlationMass == null || _correlationMass.Length != N)
            {
                _correlationMass = ComputePerNodeCorrelationMass();
            }
            
            double M = PhysicsConstants.GeometryMomentumMass;  // Inertial mass for geometry
            double damping = PhysicsConstants.GeometryDamping;
            double halfDtOverM = 0.5 * dt / M;
            
            // ===== VELOCITY VERLET STEP 1: Half-step position update =====
            // w(t+?dt) = w(t) + ?(t) * dt/(2M)
            for (int i = 0; i < N; i++)
            {
                foreach (int j in Neighbors(i))
                {
                    if (j <= i) continue; // Process each edge once
                    
                    double pi_ij = _geometryMomenta![i, j];
                    double halfStepW = Weights[i, j] + pi_ij * halfDtOverM;
                    
                    // Apply soft walls (Checklist A.4)
                    halfStepW = ApplySoftWalls(halfStepW);
                    
                    Weights[i, j] = halfStepW;
                    Weights[j, i] = halfStepW;
                }
            }
            
            // ===== VELOCITY VERLET STEP 2: Calculate forces =====
            // F = -?H_total/?w = curvature force + matter force
            // RQ-FIX: Use unified mass from GetNodeTotalMass for matter force calculation
            
            double[,] forces = new double[N, N];
            for (int i = 0; i < N; i++)
            {
                foreach (int j in Neighbors(i))
                {
                    if (j <= i) continue;
                    
                    // Curvature contribution: positive curvature increases weight
                    double curvature = ComputeFormanRicciCurvature(i, j);
                    
                    // Matter contribution: mass attracts, increases weight
                    // RQ-FIX: Use GetNodeTotalMass from UnifiedMass.cs
                    double massI = GetNodeTotalMass(i);
                    double massJ = GetNodeTotalMass(j);
                    double matterForce = PhysicsConstants.GravitationalCoupling * (massI + massJ) * 0.5;
                    
                    // Total force = curvature force + matter force - damping
                    double force = curvature + matterForce;
                    force -= damping * _geometryMomenta![i, j];
                    
                    forces[i, j] = force;
                    forces[j, i] = force;
                }
            }
            
            // ===== VELOCITY VERLET STEP 3: Full-step momentum update =====
            // ?(t+dt) = ?(t) + F * dt
            for (int i = 0; i < N; i++)
            {
                foreach (int j in Neighbors(i))
                {
                    if (j <= i) continue;
                    
                    _geometryMomenta![i, j] += forces[i, j] * dt;
                    _geometryMomenta[j, i] = _geometryMomenta[i, j]; // Symmetric
                }
            }
            
            // ===== VELOCITY VERLET STEP 4: Second half-step position update =====
            // w(t+dt) = w(t+?dt) + ?(t+dt) * dt/(2M)
            for (int i = 0; i < N; i++)
            {
                foreach (int j in Neighbors(i))
                {
                    if (j <= i) continue;
                    
                    double pi_ij = _geometryMomenta![i, j];
                    double newWeight = Weights[i, j] + pi_ij * halfDtOverM;
                    
                    // Apply soft walls (Checklist A.4)
                    newWeight = ApplySoftWalls(newWeight);
                    
                    Weights[i, j] = newWeight;
                    Weights[j, i] = newWeight;
                }
            }
            
            UpdateTargetDistancesFromWeights();
        }
        
        /// <summary>
        /// Apply soft potential walls to edge weight to prevent hard clipping.
        /// Uses smooth tanh/exp transitions at boundaries (Checklist A.4).
        /// </summary>
        private double ApplySoftWalls(double weight)
        {
            double upperWall = PhysicsConstants.WeightUpperSoftWall;
            double lowerWall = PhysicsConstants.WeightLowerSoftWall;
            
            if (weight > upperWall)
            {
                weight = upperWall + (1.0 - upperWall) * Math.Tanh(weight - upperWall);
            }
            if (weight < lowerWall)
            {
                weight = lowerWall * Math.Exp(weight / lowerWall - 1);
            }
            
            // Safety bounds
            return Math.Max(PhysicsConstants.WeightAbsoluteMinimum, 
                           Math.Min(PhysicsConstants.WeightAbsoluteMaximum, weight));
        }
        
        /// <summary>
        /// Compute total kinetic energy in geometry momenta (gravitational wave energy).
        /// K = (1/2M) ? ?_ij? (Checklist F.3)
        /// </summary>
        public double ComputeGeometryKineticEnergy()
        {
            if (!HasGeometryMomenta) return 0.0;
            
            double energy = 0.0;
            double mass = PhysicsConstants.GeometryMomentumMass;
            
            for (int i = 0; i < N; i++)
            {
                foreach (int j in Neighbors(i))
                {
                    if (j <= i) continue;
                    double p = _geometryMomenta[i, j];
                    energy += 0.5 * p * p / mass;
                }
            }
            
            return energy;
        }
        
        // === GPU Acceleration Helper Methods ===
        
        /// <summary>
        /// Cached flat arrays for GPU edge indices.
        /// These are built once and reused as long as topology doesn't change.
        /// </summary>
        private int[]? _flatEdgesFrom;
        private int[]? _flatEdgesTo;
        private int _flatEdgeCacheVersion = -1; // Track topology changes
        
        /// <summary>
        /// Public accessors for flat edge arrays (required by GPU engine).
        /// Arrays are lazily computed and cached based on _topologyVersion.
        /// When topology changes (via InvalidateTopologyCache), arrays are automatically rebuilt on next access.
        /// </summary>
        public int[] FlatEdgesFrom
        {
            get
            {
                if (_flatEdgesFrom == null || _flatEdgeCacheVersion != _topologyVersion)
                    RebuildFlatEdgeArrays();
                return _flatEdgesFrom!;
            }
        }
        
        public int[] FlatEdgesTo
        {
            get
            {
                if (_flatEdgesTo == null || _flatEdgeCacheVersion != _topologyVersion)
                    RebuildFlatEdgeArrays();
                return _flatEdgesTo!;
            }
        }
        
        // Topology version counter (incremented when edges change)
        private int _topologyVersion = 0;
        
        /// <summary>
        /// Call this method when edges are added or removed from the graph.
        /// Invalidates cached flat edge arrays by incrementing the topology version.
        /// The arrays will be automatically rebuilt on next access to FlatEdgesFrom/FlatEdgesTo.
        /// 
        /// Should be called from:
        /// - AddEdge() methods
        /// - RemoveEdge() methods
        /// - Any topology modification operations
        /// </summary>
        public void InvalidateTopologyCache()
        {
            _topologyVersion++;
        }
        
        /// <summary>
        /// Build flat arrays of edge indices for GPU processing.
        /// Each edge (i,j) where i less than j is stored as one entry.
        /// </summary>
        private void RebuildFlatEdgeArrays()
        {
            // Count edges
            int edgeCount = 0;
            for (int i = 0; i < N; i++)
            {
                foreach (int j in Neighbors(i))
                {
                    if (j > i) edgeCount++;
                }
            }
            
            _flatEdgesFrom = new int[edgeCount];
            _flatEdgesTo = new int[edgeCount];
            
            int idx = 0;
            for (int i = 0; i < N; i++)
            {
                foreach (int j in Neighbors(i))
                {
                    if (j > i)
                    {
                        _flatEdgesFrom[idx] = i;
                        _flatEdgesTo[idx] = j;
                        idx++;
                    }
                }
            }
            
            _flatEdgeCacheVersion = _topologyVersion;
        }
        
        /// <summary>
        /// Get all edge weights as a flat array for GPU processing.
        /// Order matches FlatEdgesFrom/FlatEdgesTo arrays.
        /// </summary>
        public float[] GetAllWeightsFlat()
        {
            var edgesFrom = FlatEdgesFrom; // Ensure arrays are built
            var edgesTo = FlatEdgesTo;
            int edgeCount = edgesFrom.Length;
            
            float[] weights = new float[edgeCount];
            for (int e = 0; e < edgeCount; e++)
            {
                int i = edgesFrom[e];
                int j = edgesTo[e];
                weights[e] = (float)Weights[i, j];
            }
            
            return weights;
        }
        
        /// <summary>
        /// Get all edge curvatures as a flat array for GPU processing.
        /// Uses Forman-Ricci curvature (computed on CPU for now).
        /// Order matches FlatEdgesFrom/FlatEdgesTo arrays.
        /// </summary>
        public float[] GetAllCurvaturesFlat()
        {
            var edgesFrom = FlatEdgesFrom;
            var edgesTo = FlatEdgesTo;
            int edgeCount = edgesFrom.Length;
            
            float[] curvatures = new float[edgeCount];
            for (int e = 0; e < edgeCount; e++)
            {
                int i = edgesFrom[e];
                int j = edgesTo[e];
                // Use ComputeFormanRicciCurvature from RQGraph.UnifiedEnergy.cs
                curvatures[e] = (float)ComputeFormanRicciCurvature(i, j);
            }
            
            return curvatures;
        }
        
        /// <summary>
        /// Get node masses as a flat array for GPU processing.
        /// 
        /// RQ-HYPOTHESIS FIX: Now returns NodeMasses[i].TotalMass which includes
        /// ALL field contributions (fermion, scalar, gauge, correlation, vacuum),
        /// not just topological correlation mass.
        /// 
        /// This fixes the critical gravity-matter decoupling issue where gravity
        /// was ignoring scalar field (Higgs), fermion (Dirac), and gauge field
        /// energy contributions. The Einstein equations G_?? = 8?G T_??
        /// require the FULL stress-energy tensor as source.
        /// </summary>
        public float[] GetNodeMasses()
        {
            // RQ-FIX: Update unified node mass models from all fields BEFORE gravity step
            // This ensures gravity sees matter from scalar, fermion, gauge fields
            // UpdateNodeMasses is defined in RQGraph.UnifiedMass.cs
            UpdateNodeMasses();
            
            float[] masses = new float[N];
            
            // Use GetNodeTotalMass from UnifiedMass.cs which includes ALL field contributions
            for (int i = 0; i < N; i++)
            {
                masses[i] = (float)GetNodeTotalMass(i);
            }
            
            return masses;
        }
        
        /// <summary>
        /// Update graph weights from flat array (after GPU processing).
        /// Order must match FlatEdgesFrom/FlatEdgesTo arrays.
        /// </summary>
        /// <exception cref="ArgumentException">Thrown when weights array length doesn't match edge count</exception>
        public void UpdateWeightsFromFlat(float[] weights)
        {
            var edgesFrom = FlatEdgesFrom;
            var edgesTo = FlatEdgesTo;
            
            // Validate array lengths match
            if (weights.Length != edgesFrom.Length || weights.Length != edgesTo.Length)
            {
                throw new ArgumentException(
                    $"Array length mismatch: weights.Length={weights.Length}, " +
                    $"edgesFrom.Length={edgesFrom.Length}, edgesTo.Length={edgesTo.Length}. " +
                    $"All arrays must have the same length.");
            }
            
            int edgeCount = weights.Length;
            
            for (int e = 0; e < edgeCount; e++)
            {
                int i = edgesFrom[e];
                int j = edgesTo[e];
                double w = weights[e];
                Weights[i, j] = w;
                Weights[j, i] = w; // Symmetric
            }
        }
        
        // === Edge Index Lookup for GPU ===
        
        /// <summary>
        /// Cached dictionary for fast edge index lookup by node pair.
        /// Key: (min(i,j), max(i,j)), Value: edge index in flat arrays.
        /// </summary>
        private Dictionary<(int, int), int>? _edgeIndexCache;
        private int _edgeIndexCacheVersion = -1;
        
        /// <summary>
        /// Get the flat array index for an edge (i,j).
        /// Returns -1 if the edge doesn't exist.
        /// 
        /// This is required by the GPU curvature shader for CSR structure building.
        /// </summary>
        /// <param name="i">First node of edge</param>
        /// <param name="j">Second node of edge</param>
        /// <returns>Edge index in FlatEdgesFrom/FlatEdgesTo arrays, or -1 if not found</returns>
        public int GetEdgeIndex(int i, int j)
        {
            // Ensure cache is up to date
            if (_edgeIndexCache == null || _edgeIndexCacheVersion != _topologyVersion)
            {
                RebuildEdgeIndexCache();
            }
            
            // Normalize to (min, max) for undirected edge lookup
            var key = i < j ? (i, j) : (j, i);
            
            if (_edgeIndexCache!.TryGetValue(key, out int edgeIndex))
            {
                return edgeIndex;
            }
            
            return -1; // Edge not found
        }
        
        /// <summary>
        /// Rebuild the edge index cache from flat edge arrays.
        /// Called automatically when topology version changes.
        /// </summary>
        private void RebuildEdgeIndexCache()
        {
            // Ensure flat edge arrays are built
            var edgesFrom = FlatEdgesFrom;
            var edgesTo = FlatEdgesTo;
            
            _edgeIndexCache = new Dictionary<(int, int), int>(edgesFrom.Length);
            
            for (int e = 0; e < edgesFrom.Length; e++)
            {
                int i = edgesFrom[e];
                int j = edgesTo[e];
                
                // Store with normalized key (min, max)
                var key = i < j ? (i, j) : (j, i);
                _edgeIndexCache[key] = e;
            }
            
            _edgeIndexCacheVersion = _topologyVersion;
        }
    }
}
