# Honeycomb Wannier Test & `_su2_layer_polarization!` Inline — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Inline the `_su2_layer_polarization!` helper into its sole call site and add a comprehensive "MoireWannier-honeycomb" test set following the existing triangular test pattern.

**Architecture:** Two files are modified — `src/analysis.jl` (inline refactor, no behavior change) and `test/analysis.jl` (new testset appended). The test follows the triangular pattern: Wannier construction → real-space display → HoppingIntegral + band comparison → CoulombIntegral. Numerical baseline values are captured from the first successful run and then embedded as regression assertions.

**Tech Stack:** Julia with QuantumLattices, TightBindingApproximation, StaticArrays, Plots (RecipesBase), CairoMakie

---

### Task 1: Inline `_su2_layer_polarization!` into the honeycomb constructor

**Files:**
- Modify: `src/analysis.jl:136-199`

- [ ] **Step 1: Replace the call site with inlined body**

Replace lines 136-138:
```julia
    # SU(2) rotation to maximize layer polarization
    Ũ = zeros(ComplexF64, nband, nband, nk)
    _su2_layer_polarization!(Ũ, bloch, nG, nlayer, nk)
```

with this inlined block:

```julia
    # SU(2) rotation to maximize layer polarization:
    # Diagonalize the bottom-layer projection Pb at each k.
    # Column 1 → eigenvector with max eigenvalue (maximizes bottom-layer weight)
    # Column 2 → eigenvector with min eigenvalue (maximizes top-layer weight, since Pb + Pt ≈ I)
    Ũ = zeros(ComplexF64, nband, nband, nk)
    @assert nlayer == 2 "SU(2) layer polarization requires 2 layers."
    for ik in 1:nk
        Pb = zeros(ComplexF64, 2, 2)
        for ig in 1:nG
            for μ in 1:2, ν in 1:2
                Pb[μ, ν] += bloch[1, ig, μ, ik] * conj(bloch[1, ig, ν, ik])
            end
        end
        vals, vecs = eigen(Hermitian(Pb))
        perm = sortperm(vals; rev=true)
        Ũ[:, :, ik] .= vecs[:, perm]
    end
```

Note: the previous lines 140-170 (U(1) gauge fix through end of constructor) remain unchanged but shift down by ~8 lines due to the expanded SU(2) block.

- [ ] **Step 2: Delete the standalone `_su2_layer_polarization!` function**

Delete lines 172-199 (the docstring and function definition):
```julia
"""
    _su2_layer_polarization!(Ũ, bloch, nG, nlayer, nk)

Compute the SU(2) rotation Ũ(k) at each k that maximizes layer polarization:
- Column 1: maximize bottom-layer projection ⟨P_b⟩
- Column 2: maximize top-layer projection  ⟨P_t⟩

P_b = diag(1, 0), P_t = diag(0, 1) acting on the layer index.
"""
function _su2_layer_polarization!(Ũ::Array{ComplexF64,3}, bloch::Array{ComplexF64,4}, nG::Int, nlayer::Int, nk::Int)
    @assert nlayer == 2 "SU(2) layer polarization requires 2 layers."
    for ik in 1:nk
        Pb = zeros(ComplexF64, 2, 2)
        for ig in 1:nG
            # bottom-layer component (layer=1) at G-vector ig
            for μ in 1:2, ν in 1:2
                Pb[μ, ν] += bloch[1, ig, μ, ik] * conj(bloch[1, ig, ν, ik])
            end
        end
        # diagonalize Pb (2×2 Hermitian)
        vals, vecs = eigen(Hermitian(Pb))
        # sort: eigenvector for max eigenvalue → column 1 (maximizes bottom-layer weight)
        #        eigenvector for min eigenvalue → column 2 (maximizes top-layer weight, since Pb + Pt ≈ I)
        perm = sortperm(vals; rev=true)
        Ũ[:, :, ik] .= vecs[:, perm]
    end
    return Ũ
end
```

Make sure there is exactly one blank line between the end of the constructor's `end` and the start of the `HoppingIntegral` docstring.

- [ ] **Step 3: Verify the file parses correctly**

Run:
```bash
cd f:/GitHub/Julia/MoireSuperlattices && julia --project=. -e "using MoireSuperlattices; println('Module loaded OK')"
```
Expected: prints "Module loaded OK" with no errors.

- [ ] **Step 4: Commit**

```bash
cd f:/GitHub/Julia/MoireSuperlattices && git add src/analysis.jl && git commit -m "refactor: inline _su2_layer_polarization! into honeycomb MoireWannier constructor"
```

---

### Task 2: Add the honeycomb test set structure (without numerical assertions)

**Files:**
- Modify: `test/analysis.jl` — append after line 67

- [ ] **Step 1: Append the honeycomb testset skeleton**

Add the following code after line 67 (after the closing `end` of the triangular testset, with one blank line separator):

```julia

@testset "MoireWannier-honeycomb" begin
    # Twisted MoTe₂ parameters from Zhou et al. (2026)
    parameters = (a₀=3.52, m=0.60, θ=2.94, Vᶻ=0.0, μ=0.0, V=20.8, ψ=107.7, w=-23.80)
    bltmd = Algorithm(:BLTMD, BLTMD(values(parameters)...; truncation=4), parameters)

    # Wannier function W — top two moiré bands on honeycomb effective lattice
    lattice = MoireHoneycomb(6)
    hilbert = Hilbert(Fock{:f}(1, 2), length(lattice))
    recipls = reciprocals(lattice)
    top = dimension(bltmd)
    wannier = MoireWannier(bltmd, lattice; nk=24, bands=top-1:top)
    @test count(wannier) == 2

    # (|W|^2) in the real space — two sublattices (MX at lattice[1], XM at lattice[2])
    rz = RealZone([[1.0, 0.0], [0.0, 1.0]], -2=>2, -2=>2)
    result₁ = zeros(length(rz), 2)  # sublattice 1 (MX-centered)
    result₂ = zeros(length(rz), 2)  # sublattice 2 (XM-centered)
    for (i, r) in enumerate(rz)
        result₁[i, :] = abs.(wannier(r, 1))
        result₂[i, :] = abs.(wannier(r, 2))
    end
    result₁ = reshape(result₁, map(length, reverse(shape(rz)))..., 2)
    result₂ = reshape(result₂, map(length, reverse(shape(rz)))..., 2)
    Plots.savefig(Plots.plot(rz, result₁), "Plots-MoTe₂-wannier-MX.png")
    Plots.savefig(Plots.plot(rz, result₂), "Plots-MoTe₂-wannier-XM.png")
    Makie.save("Makie-MoTe₂-wannier-MX.png", Makie.plot(rz, result₁))
    Makie.save("Makie-MoTe₂-wannier-XM.png", Makie.plot(rz, result₂))

    # HoppingIntegral and TBA
    hopping = HoppingIntegral(wannier)
    tba = Algorithm(:tba, TBA(lattice, hilbert, terms(hopping; tol=1e-6)))
    # Numerical assertions to be filled after capture run (Task 3)
    @test length(tba.parameters) > 0

    # Energy band comparison: continuum model vs TBA
    bands₁ = bltmd(:EB, EnergyBands(ReciprocalPath(recipls, hexagon"Γ-K₁-K₂-Γ", length=100)))
    bands₂ = bltmd(:EB, EnergyBands(ReciprocalPath(recipls, hexagon"Γ-K₂-K₁-Γ", length=100)))
    bands = tba(:EB, EnergyBands(ReciprocalPath(recipls, hexagon"Γ-K₁-K₂-Γ", length=100)))

    plt = Plots.plot()
    emin, emax = -40.0, 40.0
    Plots.plot!(plt, bands₁, ylims=(emin, emax), color="blue", title="")
    Plots.plot!(plt, bands₂, ylims=(emin, emax), color="green", title="")
    Plots.plot!(plt, bands, ylims=(emin, emax), ls=:dash, color="red", size=(400, 300), title="")
    Plots.savefig(plt, "Plots-MoTe₂-continuum-tba.png")
    fig = Makie.Figure()
    ax = Makie.Axis(fig[1, 1])
    Makie.plot!(ax, bands₁; ylims=(emin, emax), color=:blue, title="")
    Makie.plot!(ax, bands₂; ylims=(emin, emax), color=:green, title="")
    Makie.plot!(ax, bands; ylims=(emin, emax), linestyle=:dash, color=:red, title="")
    Makie.save("Makie-MoTe₂-continuum-tba.png", fig)

    # CoulombIntegral with three screening models
    update!(bltmd; θ=1.0, Vᶻ=0.0)
    wannier = MoireWannier(bltmd, lattice; nk=24, bands=top-1:top)
    coulomb = CoulombIntegral(wannier)
    compare(ts, vs) = all(map((t, v)->isapprox(value(t), v; atol=1e-4), ts, vs))
    # Numerical assertions to be filled after capture run (Task 3)
    @test length(terms(coulomb, BareCoulomb(10.0); order=4, tol=1e-5)) > 0
    @test length(terms(coulomb, ImageCoulomb(10.0, 20.0); order=4, tol=1e-5)) > 0
    @test length(terms(coulomb, TanhCoulomb(10.0, 20.0); order=4, tol=1e-5)) > 0
end
```

- [ ] **Step 2: Commit the skeleton**

```bash
cd f:/GitHub/Julia/MoireSuperlattices && git add test/analysis.jl && git commit -m "test: add MoireWannier-honeycomb test skeleton (assertions pending capture)"
```

---

### Task 3: Capture numerical baseline values

**Files:**
- Modify: `test/analysis.jl` — the honeycomb testset (to add println-based capture)

No permanent file changes — this task is a temporary instrumented run to capture reference values.

- [ ] **Step 1: Create a capture script**

Create a temporary script `test/capture_values.jl`:

```julia
using MoireSuperlattices
using QuantumLattices
using TightBindingApproximation
using StaticArrays: SVector

parameters = (a₀=3.52, m=0.60, θ=2.94, Vᶻ=0.0, μ=0.0, V=20.8, ψ=107.7, w=-23.80)
bltmd = Algorithm(:BLTMD, BLTMD(values(parameters)...; truncation=4), parameters)

lattice = MoireHoneycomb(6)
hilbert = Hilbert(Fock{:f}(1, 2), length(lattice))
top = dimension(bltmd)
wannier = MoireWannier(bltmd, lattice; nk=24, bands=top-1:top)
@assert count(wannier) == 2

# Capture HoppingIntegral parameters
hopping = HoppingIntegral(wannier)
tba = Algorithm(:tba, TBA(lattice, hilbert, terms(hopping; tol=1e-6)))
println("=== HOPPING PARAMETERS ===")
println("tba.parameters = ", tba.parameters)

# Capture CoulombIntegral values
update!(bltmd; θ=1.0, Vᶻ=0.0)
wannier = MoireWannier(bltmd, lattice; nk=24, bands=top-1:top)
coulomb = CoulombIntegral(wannier)

bare_terms = terms(coulomb, BareCoulomb(10.0); order=4, tol=1e-5)
bare_values = [value(t) for t in bare_terms]
println("=== BARE COULOMB ===")
println("bare_values = ", bare_values)

image_terms = terms(coulomb, ImageCoulomb(10.0, 20.0); order=4, tol=1e-5)
image_values = [value(t) for t in image_terms]
println("=== IMAGE COULOMB ===")
println("image_values = ", image_values)

tanh_terms = terms(coulomb, TanhCoulomb(10.0, 20.0); order=4, tol=1e-5)
tanh_values = [value(t) for t in tanh_terms]
println("=== TANH COULOMB ===")
println("tanh_values = ", tanh_values)
```

- [ ] **Step 2: Run the capture script**

```bash
cd f:/GitHub/Julia/MoireSuperlattices && julia --project=. test/capture_values.jl
```

Expected: script runs to completion and prints the four captured value arrays.

- [ ] **Step 3: Record the captured values**

Copy the printed values from the output. They will look something like:
```
=== HOPPING PARAMETERS ===
tba.parameters = [t₁, λ₁, t₂, λ₂, ..., μ₁, μ₂]
=== BARE COULOMB ===
bare_values = [U₁, U₂, V₁, ...]
=== IMAGE COULOMB ===
image_values = [...]
=== TANH COULOMB ===
tanh_values = [...]
```

Save these values — they will be embedded in Task 4.

- [ ] **Step 4: Clean up the capture script**

```bash
rm f:/GitHub/Julia/MoireSuperlattices/test/capture_values.jl
```

---

### Task 4: Fill in numerical assertions with captured values

**Files:**
- Modify: `test/analysis.jl` — replace placeholder assertions with captured values

Replace the placeholder lines in the honeycomb testset:

- [ ] **Step 1: Replace the HoppingIntegral placeholder assertion**

Replace:
```julia
    # Numerical assertions to be filled after capture run (Task 3)
    @test length(tba.parameters) > 0
```

with (fill in `[...]` with the captured `tba.parameters` values):

```julia
    @test all(map((x, y)->isapprox(x, y; atol=10^-4), tba.parameters,
        [/* CAPTURED: paste tba.parameters array here */]
    ))
```

- [ ] **Step 2: Replace the CoulombIntegral placeholder assertions**

Replace:
```julia
    # Numerical assertions to be filled after capture run (Task 3)
    @test length(terms(coulomb, BareCoulomb(10.0); order=4, tol=1e-5)) > 0
    @test length(terms(coulomb, ImageCoulomb(10.0, 20.0); order=4, tol=1e-5)) > 0
    @test length(terms(coulomb, TanhCoulomb(10.0, 20.0); order=4, tol=1e-5)) > 0
```

with:

```julia
    @test compare(
        terms(coulomb, BareCoulomb(10.0); order=4, tol=1e-5),
        [/* CAPTURED: paste bare_values array here */]
    )
    @test compare(
        terms(coulomb, ImageCoulomb(10.0, 20.0); order=4, tol=1e-5),
        [/* CAPTURED: paste image_values array here */]
    )
    @test compare(
        terms(coulomb, TanhCoulomb(10.0, 20.0); order=4, tol=1e-5),
        [/* CAPTURED: paste tanh_values array here */]
    )
```

- [ ] **Step 3: Commit**

```bash
cd f:/GitHub/Julia/MoireSuperlattices && git add test/analysis.jl && git commit -m "test: add numerical assertions from capture run for MoireWannier-honeycomb"
```

---

### Task 5: Verify the full test suite passes

- [ ] **Step 1: Run the honeycomb testset alone**

```bash
cd f:/GitHub/Julia/MoireSuperlattices && julia --project=. -e '
using Test
include("test/analysis.jl")
'
```

Expected: all `@test` assertions pass, no errors, PNG files generated.

- [ ] **Step 2: Run the full test suite**

```bash
cd f:/GitHub/Julia/MoireSuperlattices && julia --project=. -e '
using Test
include("test/runtests.jl")
'
```

Expected: all tests pass (`lattices`, `systems`, `analysis` with both triangular and honeycomb testsets).

- [ ] **Step 3: Commit any final adjustments**

```bash
cd f:/GitHub/Julia/MoireSuperlattices && git diff
# If clean, done. If adjustments needed:
git add -A && git commit -m "test: finalize MoireWannier-honeycomb test with verified assertions"
```

---

### Generated Artifacts

The test generates 6 PNG files in the working directory:
- `Plots-MoTe₂-wannier-MX.png` — Real-space |W|² for MX-centered Wannier (Plots)
- `Plots-MoTe₂-wannier-XM.png` — Real-space |W|² for XM-centered Wannier (Plots)
- `Makie-MoTe₂-wannier-MX.png` — Real-space |W|² for MX-centered Wannier (Makie)
- `Makie-MoTe₂-wannier-XM.png` — Real-space |W|² for XM-centered Wannier (Makie)
- `Plots-MoTe₂-continuum-tba.png` — Band structure comparison (Plots)
- `Makie-MoTe₂-continuum-tba.png` — Band structure comparison (Makie)
