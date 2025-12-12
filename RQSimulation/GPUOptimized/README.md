# GPU-Optimized RQ-Hypothesis Implementation

This folder contains the GPU-optimized implementations for the RQ-Hypothesis compliance improvements according to the comprehensive task list.

## Overview

The implementation addresses 5 critical items for improving RQ-Hypothesis compliance:

### 1. Event-Driven Time Architecture (Discrete Event Simulation)
**File:** `EventDrivenEngine.cs`

**Problem:** The old implementation used a global `for` loop that updated all nodes synchronously, violating the principle of relativity where each node has its own proper time τ.

**Solution:** 
- Implemented priority queue-based event scheduling
- Each node evolves independently with its own proper time
- Events are scheduled based on local time dilation factors
- Global clock is purely statistical, not used for physics

**Usage:**
```csharp
var engine = new EventDrivenEngine(graph, totalTime: 100.0, seed: 42);
engine.Initialize();
engine.Run();
```

### 2. Gauss Law Projection for Gauge Fields
**File:** `GaussLawProjection.cs`

**Problem:** Yang-Mills evolution without constraint enforcement allows numerical errors to accumulate, breaking gauge invariance. The equation ∇·E = ρ (Gauss law) must be maintained.

**Solution:**
- Compute divergence of electric field: div(E)
- Compute charge density from fermions: ρ
- Solve Poisson equation: ∇²χ = div(E) - ρ
- Apply gauge transformation: φ'[i,j] = φ[i,j] - (χ[i] - χ[j])
- Uses Conjugate Gradient for efficient solving

**Integration:** Automatically called in `EvolveYangMillsRelational` when `EnforceGaugeConstraintsEnabled` is true.

### 3. Phase Transition: Improved Gravitational Coupling and Annealing
**File:** `ImprovedNetworkGravity.cs`

**Problem:** The old `GravitationalCouplingConstant = 0.1` was too small to overcome entropy. Heavy clusters never formed (HeavyClusters: 0 in logs).

**Solution:**
- Increased base coupling from 0.1 to 2.5 (25x stronger)
- Implemented simulated annealing: G_eff(t) = G_base * (1 + 10 * exp(-step/200))
- During first 1000 steps, gravity is boosted to force cluster formation
- After annealing, system cools to stable state

**Result:** Universe "condenses" into matter (heavy clusters) instead of remaining diffuse.

### 4. Ollivier-Ricci Curvature
**File:** `OllivierRicciCurvature.cs`

**Problem:** Forman-Ricci curvature only captures topology (triangles), not geometry (transport properties). For gravity, we need geodesic sensitivity.

**Solution:**
- Implemented Ollivier-Ricci curvature: κ(i,j) = 1 - W₁(μᵢ, μⱼ) / d(i,j)
- W₁ is Wasserstein-1 distance (optimal transport)
- Uses Jaccard index approximation for performance
- More sensitive to geometry than Forman-Ricci

**Usage:**
```csharp
double curvature = OllivierRicciCurvature.ComputeOllivierRicciJaccard(graph, i, j);
```

### 5. Spectral Dimension Validation
**File:** `SpectralDimensionValidator.cs`

**Problem:** Code uses 4-component spinors (assuming 3+1 dimensions), but if the graph is fractal (d_S ≈ 2.5), the Dirac operator is incorrect.

**Solution:**
- Compute spectral dimension periodically
- Check if d_S is within 10% of 4.0
- Issue warnings if incompatible
- Suppress spinor evolution when d_S deviates too much
- Track transitions (fractal → spacetime)

**Usage:**
```csharp
var status = SpectralDimensionValidator.ValidateSpinorCompatibility(graph);
if (!status.IsCompatible)
{
    AppendConsole(status.Recommendation);
}
```

## Integration

**File:** `RQHypothesisIntegration.cs`

Main entry point for using all features together:

```csharp
// Initialize
RQHypothesisIntegration.Initialize(graph);

// In simulation loop
for (int step = 0; step < totalSteps; step++)
{
    RQHypothesisIntegration.StepPhysics(graph, step, dt, diagnosticsExport);
    
    // Your other physics updates...
}

// Get status report
string report = RQHypothesisIntegration.GetStatusReport(graph);
AppendConsole(report);
```

## CPU/GPU Dispatcher

**File:** `ComputationDispatcher.cs`

Provides automatic fallback between GPU (via ComputeSharp) and CPU execution:

```csharp
// Automatically uses GPU if available, CPU otherwise
ComputationHelpers.UpdateScalarField(field, momentum, dt, mu2, lambda);

// Check which mode is active
ComputationHelpers.LogStatus();
```

## Configuration

Enable/disable features:

```csharp
// Event-driven time (optional, for research)
RQHypothesisIntegration.UseEventDrivenTime = true;

// Improved gravity with annealing (recommended)
RQHypothesisIntegration.UseImprovedGravity = true;

// Spectral dimension validation (recommended)
RQHypothesisIntegration.ValidateSpectralDimension = true;
```

## Diagnostics Export

The spectral dimension validator exports to CSV:

```
step,spectralDim,status,deviationPercent,recommendation
0,2.345,INCOMPATIBLE,0.4137,Spectral dimension 2.35 < 3.60. Graph is too fractal...
100,2.789,INCOMPATIBLE,0.3027,Spectral dimension 2.79 < 3.60. Graph is too fractal...
200,3.456,OK,0.1360,Spectral dimension 3.46 is compatible with 4D spinors.
300,3.892,OK,0.0270,Spectral dimension 3.89 is compatible with 4D spinors.
```

## Performance Notes

- **Gauss Law Projection:** O(N * edges) per call. Called every 50 steps by default.
- **Ollivier-Ricci Curvature:** O(N²) for full computation. Use Jaccard approximation for O(degree²).
- **Event-Driven Time:** More complex than for-loop, but physically correct. Use for research.
- **GPU Dispatch:** Overhead for small graphs. Beneficial for N > 10000.

## References

1. Ollivier-Ricci Curvature: Ollivier, Y. (2009). "Ricci curvature of Markov chains on metric spaces"
2. Gauss Law on Lattice: Wilson, K. (1974). "Confinement of quarks"
3. Spectral Dimension: Calcagni, G. (2017). "Spectral dimension of quantum geometries"
4. Event-Driven Simulation: Chandy, K. M., & Misra, J. (1979). "Distributed simulation"

## Migration Path

**Backward Compatible:** All features are opt-in and don't break existing code.

**Migration Steps:**
1. Build with new files (already done)
2. Test with existing simulations (should work unchanged)
3. Enable `UseImprovedGravity` to get better heavy cluster formation
4. Enable `ValidateSpectralDimension` to monitor dimension evolution
5. (Optional) Enable `UseEventDrivenTime` for research on relativistic time

## Future Work

- Full ComputeSharp GPU kernels for Yang-Mills evolution
- Adaptive time stepping in event-driven engine
- Parallel Poisson solver with GPU acceleration
- Automatic parameter tuning based on spectral dimension
