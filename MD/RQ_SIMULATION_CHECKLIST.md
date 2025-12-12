# RQ-Hypothesis Simulation — AI Session Checklist

**Version**: 2025-12-08  
**Project**: RqSim  
**Solution**: `C:\Users\Yuriy Kurnosov\Desktop\RqSim\RqSim.csproj`  
**Target**: .NET 10, C# 14

---

## QUICK REFERENCE FOR AI SESSIONS

### What is RQ-Hypothesis?
**Relational-Quantum (RQ)** hypothesis proposes that:
1. **No background spacetime** — space emerges from quantum correlations
2. **Nodes** = local quantum states (NOT points in space)
3. **Edges** = quantum correlations with weights as correlation strength
4. **Time** = relational (Page-Wootters), derived from internal clock subsystem
5. **Mass** = topological binding energy + Yukawa coupling to scalar field (Higgs mechanism)
6. **Gravity** = curvature flow on graph (Ollivier-Ricci / Forman-Ricci)

### Key Prediction
**Spectral dimension d_S** should:
- Start at **d_S ~ 2** (Planck/UV scale)  
- Evolve to **d_S ~ 4** (macroscopic/IR scale)
- Transition smoothly via quantum correlations forming clusters

---

## GPU OPTIMIZATION STATUS (December 2025)

### Available GPU Engines (GPUOptimized/)

| Engine | File | Purpose | Status |
|--------|------|---------|--------|
| **ScalarFieldEngine** | `ScalarFieldEngine.cs` | Klein-Gordon GPU evolution | Complete |
| **SpectralWalkEngine** | `SpectralWalkEngine.cs` | Parallel random walks for d_S | Complete |
| **GpuGravityEngine** | `ImprovedNetworkGravity.cs` | Forman-Ricci + weight update | Complete |
| **StatisticsEngine** | `StatisticsEngine.cs` | GPU reductions (sum, histogram) | Complete |
| **ParallelEventEngine** | `ParallelEventEngine.cs` | Graph coloring parallelism | Complete |
| **GpuCurvatureEngine** | `OllivierRicciCurvature.cs` | Ollivier-Ricci Jaccard GPU | Complete |
| **PhaseCoherenceCorrector** | `PhaseCoherenceCorrector.cs` | Async time phase twists | Complete |
| **EventDrivenEngine** | `EventDrivenEngine.cs` | DES with local proper time | Complete |
| **GaussLawProjection** | `GaussLawProjection.cs` | Gauge constraint enforcement | Complete |

### GPU Usage Examples

```csharp
// Initialize GPU gravity
graph.InitGpuGravity();

// Use in simulation loop - batch mode for maximum performance
graph.EvolveGravityGpuBatch(batchSize: 10, dt, G, lambda);

// Or parallel event processing with graph coloring
var eventEngine = new ParallelEventEngine(graph);
eventEngine.ProcessParallelSweep(dt);

// GPU spectral dimension computation
var walkEngine = new SpectralWalkEngine();
walkEngine.Initialize(walkerCount: 10000, nodeCount: N, totalEdges: E);
walkEngine.UpdateTopology(offsets, neighbors, weights);
walkEngine.InitializeWalkersUniform();
int[] returns = walkEngine.RunSteps(100);
double d_S = walkEngine.ComputeSpectralDimension(returns);
```

---

## PHYSICS FIXES APPLIED (December 2025)

### 1. Local Dirac Operator for Event-Driven Mode
**File**: `GPUOptimized/RQGraph.EventDrivenExtensions.cs`
- Replaced empty `UpdateSpinorFieldAtNode` with proper local evolution
- Gauge-covariant parallel transport with U(1) phase
- Staggered fermion signs for correct chirality

### 2. Yukawa Mass Coupling (Higgs Mechanism)
**Files**: `Fields/RQGraph.Spinor.cs`, `GPUOptimized/RQGraph.EventDrivenExtensions.cs`
- Fermion mass from scalar field: `m = g_Y * |phi|`
- No more hardcoded `_correlationMass` in physics
- Rate-limited mass evolution (10% max change per step)

### 3. U(1) Gauge Transformation Fix
**File**: `GPUOptimized/RQGraph.EventDrivenExtensions.cs`
- `ApplyGaugeTransformation` now updates BOTH directions of edge
- Preserves antisymmetry: theta_ij = -theta_ji
- Maintains Wilson loop unitarity

### 4. Scalar-Gauge Back-Reaction
**File**: `Gauge/RQGraph.FieldTheory.cs`
- `ComputeScalarFieldCurrent` implements Higgs back-reaction
- Current: J = g * phi_i * phi_j * sin(theta_ij)
- Enables proper Higgs mechanism

### 5. Coordinate Isolation
**File**: `Spacetime/RQGraph.Spacetime.cs`
- Removed `ComputeGeodesicDeviation` using X/Y coordinates
- Replaced with `ComputeGeodesicDeviationScalar` using topological distance
- Removed `ComputeFrameDragging` (requires embedding)
- `Coordinates[]` marked Obsolete, used ONLY for rendering

### 6. Mass Unification
**File**: `Core/RQGraph.cs`
- Removed `StructuralMass` field
- All physics now uses `_correlationMass` (derived from edge weights)
- Ensures matter cannot exist without graph geometry

// ? CORRECT - Relational time from clock subsystem
double dt = graph.ComputeRelationalDtExtended(); // Fubini-Study metric
```

### ? ENERGY MUST BE CONSERVED
All energy changes must be tracked. No "magic" injections:
```csharp
// ? WRONG
LocalPotential[i] += 0.05; // Energy from nowhere!

// ? CORRECT
// Use Metropolis acceptance, energy ledger, or explicit source
```

### ? GAUGE INVARIANCE
Yang-Mills phases must satisfy Gauss law: `?·E = ?`
```csharp
// GaussLawProjection.EnforceGaussLaw(graph) after phase updates
```

### ? CAUSAL REWIRING ONLY
Edge creation must respect causality:
```csharp
// Only connect nodes within causal neighborhood
// No FTL links across disconnected graph regions
```

---

## ?? PROJECT STRUCTURE

```
RQSimulation/
??? Analysis/          (6 files) - Diagnostics, modern simulation examples
??? Core/              (14 files) - RQGraph partials, PhysicsConstants
?   ??? RQGraph.cs              - Main class, d_S computation
?   ??? RQGraph.GraphHealth.cs  - Giant cluster detection, recovery
?   ??? PhysicsConstants.cs     - All physics parameters
?   ??? SimulationEngine.cs     - Graph initialization
??? Fields/            (5 files) - Scalar, spinor fields
??? Gauge/             (9 files) - U(1), SU(2), SU(3) gauge fields
??? GPUOptimized/      (14 files) - GPU shaders, Ollivier-Ricci
??? Gravity/           (3 files) - Network gravity evolution
??? Physics/           (7 files) - Energy, clusters, black holes
??? Quantum/           (7 files) - Wavefunction, decoherence
??? Spacetime/         (8 files) - Spectral geometry, relational time
??? Topology/          (16 files) - Updates, neighbors, clusters
??? OLD/               (EMPTY) - Legacy code removed ?

Forms/
??? Form_Main.cs           - Main UI, simulation loop
??? Form_Main.Designer.cs  - WinForms designer
??? FormSimAPI.cs      - Simulation API (thread-safe metrics)
??? PartialForm.cs     - Dashboard updates
```

---

## ?? KEY FILES FOR MODIFICATIONS

| Task | Primary Files |
|------|---------------|
| Spectral dimension | `Core/RQGraph.cs` (ComputeSpectralDimension), `Spacetime/RQGraph.SpectralGeometry.cs` |
| Graph health | `Core/RQGraph.GraphHealth.cs`, `Core/PhysicsConstants.cs` |
| Gravity evolution | `GPUOptimized/ImprovedNetworkGravity.cs`, `Gravity/RQGraph.NetworkGravity.cs` |
| Quantum state | `Quantum/RQGraph.QuantumDynamics.cs`, `Topology/RQGraph.Updates.cs` |
| Gauge fields | `Gauge/RQGraph.YangMills.cs`, `GPUOptimized/GaussLawProjection.cs` |
| Cluster detection | `Physics/RQGraph.EnergyBasedClusters.cs`, `Topology/RQGraph.CoreHelpers.cs` |
| Constructor | `Topology/RQGraph.Updates.cs` (9-arg constructor) |

---

## ?? CURRENT STATUS (2025-12-07)

### Latest Snapshot Analysis
```json
{
  "step": 250,
  "spectralDimension": 2.2,      // Was 3.87 ? jumped to 2.2
  "networkTemperature": 4.67,
  "effectiveG": 0.0073,           // ? (warmup phase)
  "largestCluster": 559,          // 55.9% of N=1000 ?? GIANT CLUSTER
  "strongEdges": 5570
}
```

### Issues Identified
1. **d_S jumps from 3.87 to 2.21** — Method switching (HeatKernel?RandomWalk)
2. **Giant cluster at 55.9%** — Exceeds 50% threshold, indicates percolation
3. **d_S = 3.87 from step 0** — Initial graph too dense (edgeProb=0.03 ? E?15000)

### Fixes Applied (This Session)
1. ? **Fixed RQGraph constructor** — Was not initializing arrays
2. ? **Added d_S smoothing** — EMA filter (?=0.3) prevents method-switch jumps
3. ? **GraphHealth methods** — Already implemented in `RQGraph.GraphHealth.cs`
4. ? **Removed duplicate methods** — Cleaned up conflicts between partials

---

## ?? PHYSICS CONSTANTS (PhysicsConstants.cs)

### Fundamental
| Constant | Value | Source |
|----------|-------|--------|
| `FineStructureConstant` | 1/137.036 ? 0.0073 | QED |
| `StrongCouplingConstant` | 0.118 | QCD at M_Z |
| `WeakMixingAngle` | 0.231 | Electroweak |

### Annealing & Warmup
| Parameter | Value | Notes |
|-----------|-------|-------|
| `InitialAnnealingTemperature` | 10.0 | Hot start (Big Bang) |
| `FinalAnnealingTemperature` | ? ? 0.0073 | Cold universe |
| `WarmupDuration` | 500 steps | G stays at ? |
| `GravityTransitionDuration` | 137 steps | Smooth G increase |
| `WarmupGravitationalCoupling` | ? | Very soft gravity |
| `GravitationalCoupling` | 0.2 | Target after warmup |

### Graph Health Thresholds
| Parameter | Value | Action |
|-----------|-------|--------|
| `CriticalSpectralDimension` | 1.0 | Stop if d_S stays below |
| `WarningSpectralDimension` | 1.5 | Add edges for recovery |
| `GiantClusterThreshold` | 0.30 | Apply decoherence |
| `EmergencyGiantClusterThreshold` | 0.50 | Aggressive breakup |
| `GiantClusterDecoherenceRate` | 0.15 | Edge weakening rate |
| `FragmentationGracePeriodSteps` | 5 | Checks before exception |

---

## ?? PRIORITY TASKS

### ?? HIGH PRIORITY
1. **Verify d_S evolution** — Should go 2 ? 4, not 4 ? 2
   - Check `InitialEdgeProb` (should be ~0.005 for sparse start)
   - Verify `ComputeSpectralDimension()` returns smoothed value
   
2. **Giant cluster breakup** — Largest cluster > 50%
   - `ApplyGiantClusterDecoherence()` should run automatically
   - Check `GraphHealthLive` settings in UI

3. **Warmup phase too long** — G stays at ? for 500 steps
   - Consider reducing if d_S never reaches 4

### ?? MEDIUM PRIORITY
4. **Energy conservation validation** — Enable `ValidateEnergyConservationEnabled`
5. **Gauss law enforcement** — Enable `EnforceGaugeConstraintsEnabled`
6. **Topology tunneling** — For persistent giant clusters

### ?? LOW PRIORITY
7. **GPU shaders** — ComputeSharp descriptors need fixing
8. **Spectral coordinates** — For better visualization
9. **Ollivier-Ricci curvature** — More accurate than Forman-Ricci

---

## ?? EXPECTED SIMULATION BEHAVIOR

### Healthy Evolution
```
Step 0:     d_S ? 2.0, T = 10, largestCluster < 10%
Step 500:   d_S ? 2.5, T = 5, clusters forming
Step 2000:  d_S ? 3.5, T = 1, stable particles
Step 5000:  d_S ? 4.0, T = 0.01, frozen structure
```

### Pathological Signs
| Symptom | Cause | Fix |
|---------|-------|-----|
| d_S < 0 | Graph fragmented | Add edges, reduce G |
| d_S > 6 | Tree-like (hyperbolic) | Increase clustering |
| d_S stuck at 2 | No structure formation | Increase G, reduce T |
| d_S stuck at 4+ from start | Initial graph too dense | Reduce InitialEdgeProb |
| Giant cluster > 50% | Percolation | Apply decoherence |
| Energy exploding | Non-conservation | Track energy ledger |

---

## ?? TESTING COMMANDS

### Build
```bash
dotnet build RqSimForms.csproj
```

### Run Simulation (Programmatic)
```csharp
var config = new SimulationConfig {
    NodeCount = 1000,
    TotalSteps = 5000,
    InitialEdgeProb = 0.005, // Sparse start!
    GravitationalCoupling = 0.05,
    UseSpectralGeometry = true,
    UseNetworkGravity = true,
    UseHotStartAnnealing = true
};
var engine = new SimulationEngine(config);
// ... run loop
```

### Check Graph Health
```csharp
var status = graph.CheckGraphHealth(spectralDim, largestClusterSize);
if (status.NeedsCorrection)
{
    string actions = graph.PerformGraphRecovery(status);
    AppendConsole(actions);
}
```

---

## ?? RELATED DOCUMENTS

- `MD/DEEP_ANALYSIS_RQ_HYPOTHESIS.md` — Full mathematical formalism
- `RQSimulation/GPUOptimized/README.md` — GPU implementation details
- `RQSimulation/GPUOptimized/IMPLEMENTATION_SUMMARY.md` — Checklist status

---

## ?? COMMON PITFALLS

1. **Don't use `Coordinates[]` for physics** — Rendering only! (Enforced by Obsolete attribute)
2. **Don't add energy without tracking** — Use Metropolis or ledger
3. **Don't create FTL edges** — Respect causal structure
4. **Don't force-stabilize clusters** — Use energy-based methods
5. **Don't use global time for proper time** — Use `ComputeRelationalDtExtended()`
6. **Don't normalize spinors every step** — Use symplectic integrators
7. **Don't use `StructuralMass`** — Use `_correlationMass` or `ComputePerNodeCorrelationMass()`

---

## ?? SESSION NOTES TEMPLATE

```
### Session: [DATE]

**Goal**: [What you're trying to achieve]

**Files Modified**:
- [ ] file1.cs — reason
- [ ] file2.cs — reason

**Issues Found**:
1. Issue description
   - File: `path/to/file.cs`
   - Line: XXX
   - Fix: description

**Tests Performed**:
- [ ] Build succeeds
- [ ] d_S in range [1.0, 6.0]
- [ ] Energy conserved within tolerance
- [ ] No giant cluster (< 30%)

**Next Steps**:
1. Task 1
2. Task 2
```

---

## ?? CHECKLIST FOR CODE CHANGES

Before committing any physics-related change:

- [ ] Does NOT use external coordinates for physics
- [ ] Does NOT add energy without tracking source
- [ ] Does NOT create edges between causally disconnected nodes
- [ ] Uses `ComputeRelationalDtExtended()` for time steps (if applicable)
- [ ] Respects gauge invariance (if touching gauge fields)
- [ ] Build succeeds with 0 errors
- [ ] No duplicate method definitions across partials

---

*Last updated: 2025-12-07 by GitHub Copilot*
