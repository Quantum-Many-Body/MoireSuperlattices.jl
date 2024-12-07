```@meta
CurrentModule = MoireSuperlattices
```

# Twisted Homobilayer Transition Metal Dichalcogenides

For twisted homobilayer transition metal dichalcogenides, such as the twisted bilayer WSe₂ or MoTe₂, the continuum model [[Phys. Rev. Lett. 122, 086402 (2019)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.122.086402)] has two material parameters, i.e., `a₀` (monolayer lattice constant), and `m` (effective mass of the conduction band), one device parameter, i.e., `θ` (twist angle), and five phenomenological parameters, i.e., `Vᶻ` (perpendicular displacement field), `μ` (chemical potential), `V` (amplitude of Moire potential), `ψ` (phase of Moire potential), and `w` (interlayer hopping amplitude).

## Trivial Band Topology

The band topology of the moire bands depends on the aforementioned parameters. For a typical set of parameters [[Phys. Rev. Research 2, 033087 (2020)](https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.2.033087)] that may be related to twisted bilayer WSe₂, the topmost moire band is topologically trivial.

First, we calculate the energy bands with such a set of parameters in both valleys:
```@example WSe₂
using MoireSuperlattices
using Plots
using QuantumLattices
using TightBindingApproximation

# `a₀`: monolayer lattice constant (Å)
# `m`: effective mass of the conduction band (mₑ)
# `θ`: twist angle (°)
# `Vᶻ`: perpendicular displacement field (meV)
# `μ`: chemical potential (meV)
# `V`: amplitude of Moire potential (meV)
# `ψ`: phase of Moire potential (°)
# `w`: interlayer hopping amplitude (meV)
parameters = (a₀=3.28, m=0.45, θ=3.70, Vᶻ=0.0, μ=0.0, V=4.4, ψ=5.9, w=20.0)
WSe₂ = Algorithm(
    :WSe₂, BLTMD(values(parameters)...; truncation=4);
    parameters=parameters, map=bltmdmap
)
recipls = WSe₂.frontend.reciprocallattice.translations
bands₁ = WSe₂(:EB, EnergyBands(ReciprocalPath(recipls, hexagon"Γ-K₁-M₁-Γ", length=100)))
bands₂ = WSe₂(:EB, EnergyBands(ReciprocalPath(recipls, hexagon"Γ-K₄-M₄-Γ", length=100)))

emin, emax = -100.0, 30.0
plt = plot()
plot!(plt, bands₁, ylim=(emin, emax), color="blue", title="")
plot!(plt, bands₂, ylim=(emin, emax), color="blue", title="")
```
The bands are two-fold degenerate due to the valley degrees of freedom. Note that due to the spin-valley lock in WSe₂, such a degeneracy can also be viewed as the spin degeneracy.

Then we compute the Berry curvature and Chern number of the topmost moire band:
```@example WSe₂
berry = WSe₂(:BC, BerryCurvature(BrillouinZone(recipls, 48), [dimension(WSe₂)]))
plot(berry)
```

The perpendicular displacement field can induced a layer dependent potential, which will lead to the splitting of the two spins:
```@example WSe₂
update!(WSe₂; Vᶻ=8.0)
bands₁ = WSe₂(:EB, EnergyBands(ReciprocalPath(recipls, hexagon"Γ-K₁-M₁-Γ", length=100)))
bands₂ = WSe₂(:EB, EnergyBands(ReciprocalPath(recipls, hexagon"Γ-K₄-M₄-Γ", length=100)))

plt = plot()
plot!(plt, bands₁, ylim=(emin, emax), color="blue", title="")
plot!(plt, bands₂, ylim=(emin, emax), color="blue", title="")
```

For the topologically trivial topmost moire band, it forms a triangular lattice. The effective tight-binding model of this sole band can be obtained:
```@example WSe₂
lattice = MoireTriangular(6, reciprocals(recipls))
hilbert = Hilbert(Fock{:f}(1, 2), length(lattice))
tba = Algorithm(
    :tba,
    TBA(lattice, hilbert, terms(WSe₂, lattice, BrillouinZone(recipls, 24); tol=10^-6))
)
nothing # hide
```
Here, the hopping parameters are truncated up to the 6th order. We can compare the energy bands of this tight-binding model with those of the continuum model:

```@example WSe₂
bands₃ = tba(:EB, EnergyBands(ReciprocalPath(recipls, hexagon"Γ-K-M-Γ", length=100)))
plt = plot()
plot!(plt, bands₁, ylim=(emin, emax), color="blue", title="")
plot!(plt, bands₂, ylim=(emin, emax), color="blue", title="")
plot!(plt, bands₃, ylim=(emin, emax), color="red", lw=4, alpha=0.3, title="")
```

## Nontrivial Band Topology

With another set of parameters [[Phys. Rev. Lett. 132, 036501 (2024)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.132.036501)] that was proposed to be suitable for twisted bilayer MoTe₂, the top two moire bands have opposite Chern numbers.

First, we calculate the energy bands with such a set of parameters in a single valley:
```@example MoTe₂
using MoireSuperlattices
using Plots
using QuantumLattices
using TightBindingApproximation

parameters = (a₀=3.52, m=0.6, θ=3.89, Vᶻ=0.0, μ=30.0, V=20.8, ψ=107.7, w=-23.8)
MoTe₂ = Algorithm(
    :MoTe₂, BLTMD(values(parameters)...; truncation=4);
    parameters=parameters, map=bltmdmap
)
recipls = MoTe₂.frontend.reciprocallattice.translations
bands = MoTe₂(:EB, EnergyBands(ReciprocalPath(recipls, hexagon"Γ-K₁-M₁-Γ", length=100)))
plot(bands, ylim=(-50.0, 10.0), color="blue", title="")
```

Then we compute the Berry curvature and Chern numbers of the top two moire bands:
```@example MoTe₂
berry = MoTe₂(
    :BC, BerryCurvature(BrillouinZone(recipls, 48), [dimension(MoTe₂), dimension(MoTe₂)-1])
)
plot(berry)
```
