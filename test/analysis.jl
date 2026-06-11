using MoireSuperlattices
using QuantumLattices
using StaticArrays: SVector
using TightBindingApproximation
import CairoMakie as Makie
import Plots

@testset "MoireWannier-triangular" begin
    parameters = (a₀=3.30, m=0.45, θ=4.0, Vᶻ=10.0, μ=0.0, V=4.4, ψ=5.9, w=20.0)
    bltmd = Algorithm(:BLTMD, BLTMD(values(parameters)...; truncation=4), parameters)

    # Wannier function W
    lattice = MoireTriangular(6)
    hilbert = Hilbert(Fock{:f}(1, 2), length(lattice))
    recipls = reciprocals(lattice)
    wannier = MoireWannier(bltmd, lattice; nk=24, band=dimension(bltmd))
    @test count(wannier) == 1
    # (|W|^2) in the real space
    rz = RealZone([[1.0, 0.0], [0.0, 1.0]], -2=>2, -2=>2)
    result = zeros(length(rz), 2)
    for (i, r) in enumerate(rz)
        result[i, :] = abs.(wannier(r))
    end
    result = reshape(result, map(length, reverse(shape(rz)))..., 2)
    Plots.savefig(Plots.plot(rz, result), "Plots-WeSe₂-wannier.png")
    Makie.save("Makie-WeSe₂-wannier.png", Makie.plot(rz, result))

    # HoppingIntegral
    hopping = HoppingIntegral(wannier)
    tba = Algorithm(:tba, TBA(lattice, hilbert, terms(hopping; tol=1e-6)))
    @test all(map((x, y)->isapprox(x, y; atol=10^-4), tba.parameters, [-3.13361279, -1.35704902, -0.25580374, 0.0, -0.29205446, -0.24716631, -0.06155122, -0.01488868, -0.03976008, -0.05775271, -0.01915438, 0.0, 10.60934397]))

    bands₁ = bltmd(:EB, EnergyBands(ReciprocalPath(recipls, hexagon"Γ-K₁-K₂-Γ", length=100)))
    bands₂ = bltmd(:EB, EnergyBands(ReciprocalPath(recipls, hexagon"Γ-K₂-K₁-Γ", length=100)))
    bands = tba(:EB, EnergyBands(ReciprocalPath(recipls, hexagon"Γ-K₁-K₂-Γ", length=100)))

    plt = Plots.plot()
    emin, emax = -40.0, 40.0
    Plots.plot!(plt, bands₁, ylims=(emin, emax), color="blue", title="")
    Plots.plot!(plt, bands₂, ylims=(emin, emax), color="green", title="")
    Plots.plot!(plt, bands, ylims=(emin, emax), ls=:dash, color="red", size=(400, 300), title="")
    Plots.savefig(plt, "Plots-WeSe₂-continuum-tba.png")
    fig = Makie.Figure()
    ax = Makie.Axis(fig[1, 1])
    Makie.plot!(ax, bands₁; ylims=(emin, emax), color=:blue, title="")
    Makie.plot!(ax, bands₂; ylims=(emin, emax), color=:green, title="")
    Makie.plot!(ax, bands; ylims=(emin, emax), linestyle=:dash, color=:red, title="")
    Makie.save("Makie-WeSe₂-continuum-tba.png", fig)

    # CoulombIntegral
    update!(bltmd; θ=1.0, Vᶻ=0.0)
    wannier = MoireWannier(bltmd, lattice; nk=24, band=dimension(bltmd))
    coulomb = CoulombIntegral(wannier)
    compare(ts, vs) = all(map((t, v)->isapprox(value(t), v; atol=1e-4), ts, vs))
    @test compare(
        terms(coulomb, BareCoulomb(10.0); order=4, tol=1e-5),
        [63.46943858, 6.36851152, 3.08094620, 2.48762849, 1.55677789]
    )
    @test compare(
        terms(coulomb, ImageCoulomb(10.0, 20.0); order=4, tol=1e-5),
        [41.21364408, 3.99961333, 3.85050967, 3.83863948, 3.82626882]
    )
    @test compare(
        terms(coulomb, TanhCoulomb(10.0, 20.0); order=4, tol=1e-5),
        [31.44096775, 1.90862459, 1.90842110, 1.90842320, 1.90842098]
    )
end

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
    rz = RealZone([[1.0, 0.0], [0.0, 1.0]], -1=>1, -1=>1)
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
    tba = Algorithm(:tba, TBA(lattice, hilbert, terms(hopping; tol=1e-4)))
    @test all(map((x, y)->isapprox(x, y; atol=10^-4), tba.parameters,
        [-1.20822104, 2.09270444, -0.40709772, 0.60058736, -0.00592401, 0.01026946, -0.00181249, 0.00314168, -0.01524060, 0.00000102, -0.00689480, 0.03103737, 44.58298090, 44.58284580]
    ))

    # Energy band comparison: continuum model vs TBA
    bands₁ = bltmd(:EB, EnergyBands(ReciprocalPath(recipls, hexagon"Γ-K₁-K₂-Γ", length=100)))
    bands₂ = bltmd(:EB, EnergyBands(ReciprocalPath(recipls, hexagon"Γ-K₂-K₁-Γ", length=100)))
    bands = tba(:EB, EnergyBands(ReciprocalPath(recipls, hexagon"Γ-K₁-K₂-Γ", length=100)))

    plt = Plots.plot()
    emin, emax = -10.0, 50.0
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
    @test compare(
        terms(coulomb, BareCoulomb(10.0); order=4, tol=1e-3),
        [81.80992834, 81.20328224, 11.33693336, 5.93224071, 5.93094850, 3.40304785, 2.51544556]
    )
    @test compare(
        terms(coulomb, ImageCoulomb(10.0, 20.0); order=4, tol=1e-3),
        [56.98253589, 56.43494002, 4.35892023, 3.72565672, 3.63170873, 3.60396218]
    )
    @test compare(
        terms(coulomb, TanhCoulomb(10.0, 20.0); order=4, tol=1e-3),
        [45.90855073, 45.39858355, 1.81623487, 1.79161142, 1.78974441, 1.78975394]
    )
end
