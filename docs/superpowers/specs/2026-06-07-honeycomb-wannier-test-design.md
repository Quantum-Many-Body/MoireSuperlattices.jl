# Honeycomb Wannier Test & `_su2_layer_polarization!` Inline

## Context

The honeycomb Wannier constructor was added to `src/analysis.jl` (lines 109-199) but has no tests. The algorithm follows the SI of Zhou et al. (2026), *Itinerant topological magnons and spin excitons in twisted transition metal dichalcogenides*, targeting the top two moiré bands of twisted MoTe₂.

Reference: [nwaf354_supplemental_file.pdf](file:///C:/Users/waltergu/Desktop/nwaf354_supplemental_file.pdf)

## Changes

### 1. Inline `_su2_layer_polarization!`

**File**: `src/analysis.jl`

Remove the standalone function `_su2_layer_polarization!` (lines 172-199) and move its logic directly into the `MoireWannier` honeycomb constructor. The function has exactly one call site (line 138), making inlining straightforward:

- Replace lines 137-138:
  ```julia
  Ũ = zeros(ComplexF64, nband, nband, nk)
  _su2_layer_polarization!(Ũ, bloch, nG, nlayer, nk)
  ```
  with the inlined SU(2) layer-polarization computation.

- Delete lines 172-199 (the function definition and docstring).

No behavioral change — pure refactoring.

### 2. "MoireWannier-honeycomb" Test Set

**File**: `test/analysis.jl`

Add a new `@testset "MoireWannier-honeycomb"` block after the existing triangular test set (after line 67).

#### Parameters

BLTMD continuum model for twisted MoTe₂:
```
(a₀=3.52, m=0.60, θ=2.94, Vᶻ=0.0, μ=0.0, V=20.8, ψ=107.7, w=-23.80)
```
with `truncation=4`.

#### Test Steps

1. **Wannier construction** — `MoireHoneycomb(6)` lattice, `nk=24`, `bands=dimension(bltmd)-1:dimension(bltmd)` (top two bands). Assert `count(wannier) == 2`.

2. **Real-space display** — Two heatmaps (one per sublattice) of `|W(r)|²` on `RealZone([1.0 0.0; 0.0 1.0], -2=>2, -2=>2)`. Generate Plots and Makie PNGs.

3. **HoppingIntegral** — Construct `HoppingIntegral(wannier)` and extract `terms(hopping; tol=1e-6)`. Assert against recorded baseline values. Honeycomb has SOC (λ) terms with `SpinOrbitalCouplingAmplitude` — both `t` and `λ` coefficients are part of the tight-binding model.

4. **Energy bands** — Compare continuum model bands (top 2) vs TBA tight-binding bands along the `hexagon"Γ-K₁-K₂-Γ"` reciprocal path.

5. **CoulombIntegral** — Test with `BareCoulomb(10.0)`, `ImageCoulomb(10.0, 20.0)`, and `TanhCoulomb(10.0, 20.0)`. Since 2 bands produce 2×2 Coulomb matrices, the onsite Hubbard terms may differ per sublattice. Assert against recorded baseline values.

#### Numerical Assertions

Since the algorithm hasn't been validated against reference data, numerical assertions use **recorded baseline values** from the first successful run. These values will be embedded in the test and serve as regression guards. They should be updated if the algorithm is corrected or recalibrated against the reference paper.
