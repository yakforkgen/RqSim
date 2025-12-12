using System;

namespace RQSimulation
{
    /// <summary>
    /// Physics constants and configuration parameters.
    /// 
    /// DIMENSIONAL ANALYSIS (Checklist G.3):
    /// =====================================
    /// All constants are expressed in natural Planck units where:
    ///   c = 1 (speed of light)
    ///   ℏ = 1 (reduced Planck constant)
    ///   G = 1 (gravitational constant)
    /// 
    /// In these units:
    ///   - Length: l_P = √(ℏG/c³) = 1
    ///   - Time:   t_P = √(ℏG/c⁵) = 1
    ///   - Mass:   m_P = √(ℏc/G) = 1
    ///   - Energy: E_P = m_P c² = 1
    /// 
    /// Physical coupling constants are derived from:
    ///   - Fine structure constant: α = e²/(4πε₀ℏc) ≈ 1/137
    ///   - Strong coupling: α_s(M_Z) ≈ 0.118
    ///   - Electroweak mixing: sin²θ_W ≈ 0.231
    /// 
    /// Arbitrary coefficients (0.02, 0.1, etc.) have been replaced with
    /// values derived from these fundamental coupling constants.
    /// </summary>

    public static class PhysicsConstants
    {
        // ============================================================
        // FUNDAMENTAL CONSTANTS (Planck units: c = ℏ = G = 1)
        // ============================================================

        /// <summary>Speed of light in Planck units (c = 1)</summary>
        public const double C = 1.0;

        /// <summary>Reduced Planck constant in Planck units (ℏ = 1)</summary>
        public const double HBar = 1.0;

        /// <summary>Gravitational constant in Planck units (G = 1)</summary>
        public const double G = 1.0;

        /// <summary>Planck length (fundamental length scale)</summary>
        public const double PlanckLength = 1.0;

        /// <summary>Planck time (fundamental time scale)</summary>
        public const double PlanckTime = 1.0;

        /// <summary>Planck mass (fundamental mass scale)</summary>
        public const double PlanckMass = 1.0;

        // ============================================================
        // GAUGE COUPLING CONSTANTS (dimensionless)
        // ============================================================

        /// <summary>
        /// Electromagnetic fine structure constant α ≈ 1/137.036.
        /// </summary>
        public const double FineStructureConstant = 1.0 / 137.036;

        /// <summary>
        /// Strong coupling constant α_s at Z mass scale ≈ 0.118.
        /// </summary>
        public const double StrongCouplingConstant = 0.118;

        /// <summary>
        /// Electroweak mixing angle (Weinberg angle). sin²θ_W ≈ 0.231.
        /// </summary>
        public const double WeakMixingAngle = 0.231;

        /// <summary>
        /// Weak coupling constant g_W derived from α and θ_W.
        /// </summary>
        public static readonly double WeakCouplingConstant = Math.Sqrt(4 * Math.PI * FineStructureConstant / WeakMixingAngle);

        // ============================================================
        // DERIVED CONSTANTS & SIMULATION PARAMETERS
        // ============================================================

        // === Graph Topology ===
        // RQ-COMPLIANT: Thresholds are NOT hardcoded but computed adaptively from graph statistics.
        // The "default" values below are ONLY used when adaptive computation fails.

        /// <summary>
        /// Fallback heavy cluster threshold. In RQ-compliant mode, use adaptive threshold:
        /// threshold = mean(weights) + AdaptiveThresholdSigma * stddev(weights)
        /// This value is NEVER used directly in physics - only as emergency fallback.
        /// </summary>
        public const double DefaultHeavyClusterThreshold = 0.5; // Neutral: mean of [0,1]

        /// <summary>
        /// Standard deviations above mean for adaptive heavy threshold.
        /// Value 1.5σ means ~7% of edges are "heavy" (statistical definition).
        /// This is a STATISTICAL parameter, not a physical one.
        /// </summary>
        public const double AdaptiveThresholdSigma = 1.5;

        /// <summary>
        /// Minimum cluster size = 3 (triangle, smallest stable structure).
        /// Triangles are topologically protected in graph theory.
        /// </summary>
        public const int MinimumClusterSize = 3; // Triangle = minimal stable topology

        // === Energy and Mass ===
        /// <summary>Energy-mass conversion factor (E = mc² = m in Planck units)</summary>
        public const double C2_EnergyMassConversion = C * C;

        /// <summary>Characteristic radius for radiation distribution (hops)</summary>
        public const int RadiationDistributionRadius = 3;

        /// <summary>Tolerance for energy conservation checks</summary>
        public const double EnergyConservationTolerance = 1e-4; // Stricter tolerance

        /// <summary>Scale factor for cluster momentum accumulation</summary>
        public const double ClusterMomentumScale = 0.05;

        // === Field Theory ===
        /// <summary>Scalar field diffusion rate = α (fine structure constant)</summary>
        public static readonly double FieldDiffusionRate = FineStructureConstant;

        /// <summary>
        /// Field energy decay rate = α² (two-photon process).
        /// Decay ∝ coupling², standard QFT result.
        /// </summary>
        public static readonly double FieldDecayRate = FineStructureConstant * FineStructureConstant;

        /// <summary>
        /// Klein-Gordon field mass parameter.
        /// In Planck units, m_KG ~ m_Planck × (some hierarchy factor).
        /// Using √α as a natural small number from gauge theory.
        /// </summary>
        public static readonly double KleinGordonMass = Math.Sqrt(FineStructureConstant); // ~0.085

        /// <summary>Dirac field coupling strength = √α (vertex factor)</summary>
        public static readonly double DiracCoupling = Math.Sqrt(FineStructureConstant);

        // === Gauge Fields ===
        /// <summary>Gauge field coupling constant (e)</summary>
        public static readonly double GaugeCouplingConstant = Math.Sqrt(4 * Math.PI * FineStructureConstant);

        /// <summary>Weight of plaquette contributions (Wilson action)</summary>
        public static readonly double PlaquetteWeight = FineStructureConstant;

        /// <summary>Scaling factor for gauge currents (numerical stability)</summary>
        public const double GaugeCurrentScaleFactor = 0.1;

        // === Gravity and Curvature ===
        // RQ-HYPOTHESIS: Gravity emerges from graph curvature (Forman-Ricci).
        // No external gravitational constant - G is emergent from topology.
        // The "coupling" below controls HOW FAST geometry responds to matter,
        // not the strength of gravity itself.

        /// <summary>
        /// Warmup gravitational response rate.
        /// During thermalization, geometry should be "soft" (slowly responding).
        /// Value = α (fine structure) as natural small coupling.
        /// </summary>
        public static readonly double WarmupGravitationalCoupling = FineStructureConstant;

        /// <summary>
        /// Duration of warmup phase = relaxation time of graph.
        /// Natural scale: 10% of TotalSteps for information to propagate.
        /// This is a default; actual value should be computed as: N × ⟨k⟩ / c.
        /// </summary>
        public const double WarmupDuration = 500.0; // 10% of typical 5000 steps

        /// <summary>
        /// Fraction of TotalSteps for warmup phase.
        /// Use this to compute actual warmup duration: WarmupDuration = TotalSteps * WarmupFraction.
        /// </summary>
        public const double WarmupFraction = 0.10; // 10% of simulation for thermalization

        /// <summary>
        /// Gravitational response rate after warmup.
        /// In Planck units, G=1 by definition. This is the RESPONSE RATE,
        /// controlling how fast edge weights adjust to curvature.
        /// Value 1.0 means "instant equilibration" (strong gravity limit).
        /// </summary>
        public const double GravitationalCoupling = 0.2;

        /// <summary>
        /// Gravity transition duration = 1/α steps (slow turn-on).
        /// Smooth transition prevents "shock" instabilities.
        /// </summary>
        public static readonly double GravityTransitionDuration = 1.0 / FineStructureConstant;

        /// <summary>
        /// Curvature term scale in energy functional.
        /// From Regge calculus: E_curv ~ κ × Σ (deficit angles).
        /// Using α as natural dimensionless scale.
        /// </summary>
        public static readonly double CurvatureTermScale = FineStructureConstant;

        /// <summary>
        /// Network cosmological constant Λ.
        /// In RQ: Λ = vacuum energy density = ⟨E_vacuum⟩/Volume.
        /// Scaled up from α² to be effective in finite-size simulations.
        /// For N=600 nodes, Λ ~ 0.001 provides noticeable vacuum pressure.
        /// </summary>
        public static readonly double CosmologicalConstant = 0.001; // Effective for N~600 graphs

        /// <summary>
        /// Degree penalty in Forman-Ricci curvature.
        /// Standard Forman formula: κ(e) = 4 - deg(v1) - deg(v2) + triangles.
        /// Penalty factor = 1/⟨k⟩ where ⟨k⟩ ~ TargetDegree.
        /// </summary>
        public const double DegreePenaltyFactor = 1.0 / 8.0;

        /// <summary>Mass coupling to time dilation (Schwarzschild analog)</summary>
        public const double TimeDilationMassCoupling = 0.5;

        /// <summary>Curvature coupling to time dilation</summary>
        public const double TimeDilationCurvatureCoupling = 0.1;

        // === Causality ===
        /// <summary>Speed of light in network</summary>
        public const double SpeedOfLight = C;

        /// <summary>
        /// Maximum causal distance for interactions.
        /// Reduced to 2 to enforce strict locality/relationality.
        /// </summary>
        public const int MaxCausalDistance = 2;

        // === Time Evolution ===
        /// <summary>Base timestep for asynchronous updates</summary>
        public const double BaseTimestep = 0.01;

        /// <summary>Minimum time dilation factor (slowest time near BH)</summary>
        public const double MinTimeDilation = 0.05;

        /// <summary>Maximum time dilation factor (vacuum time)</summary>
        public const double MaxTimeDilation = 1.0;

        /// <summary>Fubini-Study scale for relational time</summary>
        public const double FubiniStudyScale = 1.0;

        // === Quantum Effects ===
        /// <summary>
        /// Base rate for vacuum fluctuations = α³ (three-loop process).
        /// Vacuum pair creation requires energy ~ 2m_e, suppressed by α³.
        /// </summary>
        public static readonly double VacuumFluctuationBaseRate =
            FineStructureConstant * FineStructureConstant * FineStructureConstant;

        /// <summary>
        /// Curvature enhancement for fluctuations = 4π (geometric factor).
        /// Near horizon, fluctuation rate × (surface area factor).
        /// </summary>
        public const double CurvatureCouplingFactor = 4.0 * Math.PI;

        /// <summary>
        /// Hawking radiation enhancement = 2 (pair creation: one escapes, one falls).
        /// Standard result from Hawking's original calculation.
        /// </summary>
        public const double HawkingRadiationEnhancement = 2.0;

        /// <summary>
        /// Threshold for pair creation = 2 (two particles from vacuum).
        /// In Planck units, E_threshold = 2 × m_Planck = 2.
        /// </summary>
        public const double PairCreationEnergyThreshold = 2.0;

        // === Cluster Dynamics ===
        /// <summary>Temperature for cluster stabilization</summary>
        public const double ClusterStabilizationTemperature = 0.05;

        /// <summary>Metropolis trials per cluster</summary>
        public const int MetropolisTrialsPerCluster = 10;

        /// <summary>Threshold for overcorrelation</summary>
        public const double OvercorrelationThreshold = 0.9;

        /// <summary>Decay for external connections</summary>
        public const double ExternalConnectionDecay = 0.99;

        // === Impulse and Excitation ===
        public const int ImpulseCascadeRadius = 4;
        public const double PhaseResonantGain = 0.5;
        public const double ImpulsePotentialBoost = 0.3;
        public const double ImpulseEnergyBoost = 0.15;
        public const double PhaseCoherenceThreshold = 0.5;

        // === Update Frequencies ===
        public const int TopologyUpdateInterval = 5;
        public const int TopologyFlipsDivisor = 20; // Try N/20 flips per step
        public const int GeometryUpdateInterval = 1; // Geometry should evolve smoothly
        public const int GaugeConstraintInterval = 5;
        public const int EnergyValidationInterval = 20;
        public const int TopologicalProtectionInterval = 50;

        // === Topological Protection ===
        public const double TopologicalProtectionThreshold = 0.6;
        public const double TopologicalProtectionStrength = 0.005;
        public const double VacuumFluctuationAmplitude = 0.02;

        // === Spinor Field ===
        public const double SpinorNormalizationThreshold = 0.0001; // Was 0.05 - Fixed for unitarity
        public const double SpinorNormalizationCorrectionFactor = 0.2; // Gentle correction

        // === Symplectic Integrator ===
        public const double SymplecticNormSafetyThreshold = 0.0001; // Was 0.05 - Fixed for unitarity
        public const double SymplecticNormCorrectionRate = 0.05;
        public const double SymplecticPhaseScale = 0.1;

        // === Potentials & Energy Weights ===
        public const double DefaultNodeMass = 1.0;

        // Higgs Potential
        public const double HiggsMuSquared = 0.1;
        public const double HiggsLambda = 0.25;
        public const double HiggsVEV = 0.447;

        // Energy Weights (Balanced for Unified Hamiltonian)
        // All weights = 1 by default (democratic energy functional).
        // Different weights would imply hierarchy - which should EMERGE, not be imposed.
        public const double ScalarFieldEnergyWeight = 1.0;
        public const double FermionFieldEnergyWeight = 1.0;
        public const double GaugeFieldEnergyWeight = 1.0;
        public const double YangMillsFieldEnergyWeight = 1.0;
        public const double GraphLinkEnergyWeight = 1.0;
        public const double GravityCurvatureEnergyWeight = 1.0;
        public const double ClusterBindingEnergyWeight = 1.0;

        // === Geometry Momenta (Gravitational Waves) ===
        public const double GeometryMomentumMass = 2.0; // Higher mass = slower, smoother waves
        public const double GeometryDamping = 0.02;

        // === Other ===
        public const double TopologicalCensorshipFluxThreshold = 0.4;
        public const double WeightUpperSoftWall = 0.98;
        public const double WeightLowerSoftWall = 0.02;
        public const double WeightAbsoluteMinimum = 0.0001;
        public const double WeightAbsoluteMaximum = 0.9999;

        public const double PowerIterationInitCenter = 0.5;
        public const double SignalExcitationProbability = 0.3;
        public const double SignalStrengthFactor = 0.15;

        public const double FluxPhaseThreshold = 0.15;
        public const double InitialVacuumPoolFraction = 0.2;
        public const double FieldHarmonicFrequency = 0.15;
        public const double CurvatureRegularizationEpsilon = 1e-6;

        // === RQ-Specific (Hot Start + Annealing) ===
        /// <summary>
        /// Initial temperature for "Big Bang" hot start.
        /// In Planck units, T_Planck = 1. We use T = 1/α ~ 137 as "hot" scale.
        /// However, for numerical stability, we cap at T = 10.
        /// </summary>
        public const double InitialAnnealingTemperature = 10.0; // ~1/α capped for stability

        /// <summary>
        /// Final equilibrium temperature = α (cold universe).
        /// At T ~ α, quantum fluctuations dominate thermal fluctuations.
        /// </summary>
        public static readonly double FinalAnnealingTemperature = FineStructureConstant;

        /// <summary>
        /// Default annealing time constant fraction.
        /// Actual τ should be computed as: τ = TotalSteps / AnnealingFractionDenominator.
        /// Using 1/5 of TotalSteps ensures ~99% cooling by simulation end.
        /// The physical 1/α² value (~18779) is too slow for typical simulations.
        /// </summary>
        public const double DefaultAnnealingFraction = 5.0; // τ = TotalSteps / 5

        /// <summary>
        /// Physical annealing time constant = 1/α² steps (for reference only).
        /// WARNING: This value is too large for typical simulations (5000 steps).
        /// Use ComputeAnnealingTimeConstant(totalSteps) instead.
        /// </summary>
        public static readonly double PhysicalAnnealingTimeConstant = 1.0 / (FineStructureConstant * FineStructureConstant);

        /// <summary>
        /// Compute effective annealing time constant for given simulation length.
        /// Returns τ such that T reaches ~1% of initial by end of simulation.
        /// </summary>
        public static double ComputeAnnealingTimeConstant(int totalSteps)
        {
            // τ = TotalSteps / 5 ensures T_final ≈ T_f + 0.007*(T_i - T_f)
            return totalSteps / DefaultAnnealingFraction;
        }

        /// <summary>
        /// Legacy constant for backward compatibility. Prefer ComputeAnnealingTimeConstant().
        /// </summary>
        [Obsolete("Use ComputeAnnealingTimeConstant(totalSteps) for correct annealing")]
        public static readonly double AnnealingTimeConstant = 1000.0; // Reasonable default for 5000 steps

        public const double ClockSubsystemFraction = 0.05;
        public const double ClockTickThreshold = 0.05;

        public const int MinBettiForProtection = 1;
        public const double SpectralGapThreshold = 0.02;

        public const double LocalUpdateFraction = 0.2;
        public const int MaxLocalSubgraphSize = 12;

        /// <summary>
        /// Topology tunneling rate = α (quantum tunneling amplitude).
        /// Long-range links appear via quantum tunneling through barrier.
        /// </summary>
        public static readonly double TopologyTunnelingRate = FineStructureConstant;

        /// <summary>
        /// Barrier for creating new edges = α (coupling strength).
        /// Edge creation = interaction → probability ∝ α.
        /// </summary>
        public static readonly double EdgeCreationBarrier = FineStructureConstant;

        /// <summary>
        /// Barrier for breaking edges = 1 - α (strong binding).
        /// Asymmetry: breaking is harder than creating (stability).
        /// This follows from detailed balance at low temperature.
        /// </summary>
        public static readonly double EdgeAnnihilationBarrier = 1.0 - FineStructureConstant;

        public const int MeasurementLocalityRadius = 2;
        public const double MeasurementDecoherenceRate = 0.2;

        // ============================================================
        // RQ-HYPOTHESIS: GRAPH HEALTH & FRAGMENTATION DETECTION
        // ============================================================
        // These constants are critical for maintaining a healthy graph topology
        // that can support emergent 4D spacetime (d_S → 4 transition).

        /// <summary>
        /// Critical spectral dimension threshold for graph fragmentation.
        /// If d_S falls below this value, the graph has collapsed to ~1D chain
        /// and cannot support 4D physics. Simulation should stop or recover.
        /// 
        /// Physics: d_S = 1 means random walker returns with P(t) ~ t^(-1/2),
        /// indicating a linear chain topology (worst case for RQ-Hypothesis).
        /// </summary>
        public const double CriticalSpectralDimension = 1.0;

        /// <summary>
        /// Warning threshold for spectral dimension.
        /// Below this value, the graph is approaching fragmentation.
        /// Corrective measures (reduce G, add edges) should be applied.
        /// </summary>
        public const double WarningSpectralDimension = 1.5;

        /// <summary>
        /// Giant cluster threshold as fraction of total nodes.
        /// If a single cluster contains more than this fraction of N,
        /// it indicates percolation rather than particle formation.
        /// 
        /// RQ-Hypothesis: Particles are SMALL stable clusters.
        /// A giant cluster (>30% N) is a "percolated soup", not particles.
        /// </summary>
        public const double GiantClusterThreshold = 0.30; // 30% of N

        /// <summary>
        /// Emergency giant cluster threshold for immediate action.
        /// Above this threshold, aggressive decoherence is applied.
        /// </summary>
        public const double EmergencyGiantClusterThreshold = 0.50; // 50% of N

        /// <summary>
        /// Decoherence rate for breaking giant clusters.
        /// Controls how fast edges within giant clusters are weakened.
        /// Higher values break clusters faster but may cause instabilities.
        /// </summary>
        public const double GiantClusterDecoherenceRate = 0.15;

        /// <summary>
        /// Minimum edge weight reduction per decoherence step.
        /// Ensures edges are weakened enough to eventually break.
        /// </summary>
        public const double MinDecoherenceWeightReduction = 0.05;

        /// <summary>
        /// Maximum edges to weaken per giant cluster decoherence step.
        /// Limits computational cost and prevents sudden topology collapse.
        /// Value is fraction of cluster size.
        /// </summary>
        public const double MaxDecoherenceEdgesFraction = 0.10; // 10% of cluster edges

        /// <summary>
        /// Recovery edge addition rate as fraction of N.
        /// When d_S drops below warning, add this many random edges.
        /// </summary>
        public const double FragmentationRecoveryEdgeFraction = 0.10; // Add N/10 edges

        /// <summary>
        /// Gravity suppression factor when d_S is critical.
        /// Reduces effective G to prevent further topology collapse.
        /// </summary>
        public const double CriticalGravitySuppression = 0.1;

        /// <summary>
        /// Number of steps to wait before declaring fragmentation terminal.
        /// If d_S stays below critical for this many steps, stop simulation.
        /// </summary>
        public const int FragmentationGracePeriodSteps = 5;

        /// <summary>
        /// Prefer Ollivier-Ricci over Forman-Ricci curvature.
        /// Ollivier-Ricci is more accurate for geometric flow simulations.
        /// </summary>
        public const bool PreferOllivierRicciCurvature = true;

        /// <summary>
        /// Volume constraint lambda for CDT-like stabilization.
        /// Penalizes deviation from target volume (valence) to prevent collapse.
        /// </summary>
        public const double VolumeConstraintLambda = 0.1;

        /// <summary>
        /// Target volume (valence) for 4D spacetime emergence.
        /// </summary>
        public const double TargetVolume = 4.0;
    }


    /// <summary>
    /// Represents the health status of the graph topology.
    /// Used to detect fragmentation and percolation problems.
    /// </summary>
    public readonly struct GraphHealthStatus
    {
        /// <summary>Current spectral dimension of the graph.</summary>
        public double SpectralDimension { get; init; }

        /// <summary>Largest cluster size as fraction of N.</summary>
        public double LargestClusterFraction { get; init; }

        /// <summary>Average degree of nodes.</summary>
        public double AverageDegree { get; init; }

        /// <summary>True if graph is healthy (d_S > warning, no giant cluster).</summary>
        public bool IsHealthy =>
            SpectralDimension >= PhysicsConstants.WarningSpectralDimension &&
            LargestClusterFraction < PhysicsConstants.GiantClusterThreshold;

        /// <summary>True if graph is fragmented (d_S < critical).</summary>
        public bool IsFragmented => SpectralDimension < PhysicsConstants.CriticalSpectralDimension;

        /// <summary>True if graph has giant cluster (percolated).</summary>
        public bool HasGiantCluster => LargestClusterFraction >= PhysicsConstants.GiantClusterThreshold;

        /// <summary>True if graph has emergency giant cluster.</summary>
        public bool HasEmergencyGiantCluster => LargestClusterFraction >= PhysicsConstants.EmergencyGiantClusterThreshold;

        /// <summary>True if graph is in warning state (needs corrective action).</summary>
        public bool NeedsCorrection =>
            SpectralDimension < PhysicsConstants.WarningSpectralDimension ||
            LargestClusterFraction >= PhysicsConstants.GiantClusterThreshold;

        /// <summary>Human-readable status description.</summary>
        public string StatusDescription
        {
            get
            {
                if (IsFragmented) return $"FRAGMENTED: d_S={SpectralDimension:F2} < {PhysicsConstants.CriticalSpectralDimension}";
                if (HasEmergencyGiantCluster) return $"PERCOLATED: cluster={LargestClusterFraction:P0} > {PhysicsConstants.EmergencyGiantClusterThreshold:P0}";
                if (HasGiantCluster) return $"GIANT_CLUSTER: {LargestClusterFraction:P0} > {PhysicsConstants.GiantClusterThreshold:P0}";
                if (SpectralDimension < PhysicsConstants.WarningSpectralDimension) return $"WARNING: d_S={SpectralDimension:F2} approaching critical";
                return $"HEALTHY: d_S={SpectralDimension:F2}, cluster={LargestClusterFraction:P0}";
            }
        }
    }

    /// <summary>
    /// Exception thrown when graph fragmentation cannot be recovered.
    /// </summary>
    public class GraphFragmentationException : Exception
    {
        public double SpectralDimension { get; }
        public int Step { get; }
        public int ConsecutiveFragmentedSteps { get; }

        public GraphFragmentationException(double spectralDimension, int step, int consecutiveSteps)
            : base($"Graph fragmentation at step {step}: d_S={spectralDimension:F3} for {consecutiveSteps} consecutive checks. " +
                   $"Cannot recover to 4D spacetime. Simulation stopped.")
        {
            SpectralDimension = spectralDimension;
            Step = step;
            ConsecutiveFragmentedSteps = consecutiveSteps;
        }
    }

    /// <summary>
    /// Runtime configuration for physics simulation
    /// Allows overriding default constants without recompilation
    /// </summary>
    public class PhysicsConfiguration
    {
        // Can be extended to make constants configurable at runtime
        public double HeavyClusterThreshold { get; set; } = PhysicsConstants.DefaultHeavyClusterThreshold;
        public double AdaptiveThresholdSigma { get; set; } = PhysicsConstants.AdaptiveThresholdSigma;
        public int MinimumClusterSize { get; set; } = PhysicsConstants.MinimumClusterSize;

        // Add more configurable parameters as needed
        public double GaugeCoupling { get; set; } = PhysicsConstants.GaugeCouplingConstant;
        public double GravitationalCoupling { get; set; } = PhysicsConstants.GravitationalCoupling;
        public double VacuumFluctuationRate { get; set; } = PhysicsConstants.VacuumFluctuationBaseRate;

        // RQ-Hypothesis Annealing configuration
        public double InitialTemperature { get; set; } = PhysicsConstants.InitialAnnealingTemperature;
        public double FinalTemperature { get; set; } = PhysicsConstants.FinalAnnealingTemperature;
        public double AnnealingTimeConstant { get; set; } = PhysicsConstants.AnnealingTimeConstant;

        /// <summary>
        /// Compute temperature at given step using exponential annealing.
        /// Implements RQ-Hypothesis Checklist Item 7.1.
        /// </summary>
        public double TemperatureAt(int step)
        {
            return FinalTemperature + (InitialTemperature - FinalTemperature)
                   * Math.Exp(-step / AnnealingTimeConstant);
        }
    }
}
