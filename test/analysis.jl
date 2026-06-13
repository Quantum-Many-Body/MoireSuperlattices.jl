using MoireSuperlattices
using QuantumLattices
using TightBindingApproximation
import CairoMakie as Makie
import Plots

@testset "MoireWannier-triangular" begin
    parameters = (a₀=3.30, m=0.45, θ=4.0, Vᶻ=10.0, μ=0.0, V=4.4, ψ=5.9, w=20.0)
    bltmd = Algorithm(:BLTMD, BLTMD(values(parameters)...; truncation=4), parameters)

    # Wannier function W
    lattice = MoireTriangular()
    hilbert = Hilbert(Fock{:f}(1, 2), length(lattice))
    wannier = MoireWannier(bltmd, lattice; nk=18)
    @test count(wannier) == 1
    # (|W|^2) in the real space
    realzone = RealZone([[1.0, 0.0], [0.0, 1.0]], -1=>1, -1=>1)
    Plots.savefig(Plots.plot(realzone, wannier, 1), "Plots-WeSe₂-wannier.png")
    Makie.save("Makie-WeSe₂-wannier.png", Makie.plot(realzone, wannier, 1))

    # HoppingIntegral and TBA
    hopping = HoppingIntegral(wannier)
    tba = Algorithm(:tba, TBA(lattice, hilbert, terms(hopping; order=6, atol=1e-3, rtol=1e-3)))
    @test all(map(
        (x, y)->isapprox(x, y; atol=1e-3, rtol=1e-3),
        tba.parameters,
        [-3.134, -1.357, -0.256, 0.0, -0.292, -0.247, -0.062, -0.015, -0.040, -0.058, -0.019, 0.0, 10.609]
    ))
    # Energy band comparison: continuum model vs TBA
    recipls = reciprocals(lattice)
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
    coulomb = CoulombIntegral(wannier)
    compare(ts, vs) = all(map((t, v)->isapprox(value(t), v; atol=1e-3, rtol=1e-3), ts, vs))
    @test compare(
        terms(coulomb, BareCoulomb(10.0); order=4, atol=1e-3, rtol=1e-3),
        [98.816, 26.177, 11.024, 8.548, 4.634]
    )
    @test compare(
        terms(coulomb, ImageCoulomb(10.0, 20.0); order=4, atol=1e-3, rtol=1e-3),
        [101.519, 37.284, 28.939, 28.206, 27.338]
    )
    @test compare(
        terms(coulomb, TanhCoulomb(10.0, 20.0); order=4, atol=1e-3, rtol=1e-3),
        [77.818, 18.932, 13.717, 13.536, 13.369]
    )
end

@testset "MoireWannier-honeycomb" begin
    # Twisted MoTe₂ parameters from Zhou et al. (2026)
    parameters = (a₀=3.52, m=0.60, θ=2.94, Vᶻ=0.0, μ=0.0, V=20.8, ψ=107.7, w=-23.80)
    bltmd = Algorithm(:BLTMD, BLTMD(values(parameters)...; truncation=4), parameters)

    # Wannier function W — top two moire bands on honeycomb effective lattice
    lattice = MoireHoneycomb()
    hilbert = Hilbert(Fock{:f}(1, 2), length(lattice))
    wannier = MoireWannier(bltmd, lattice; nk=18)
    @test count(wannier) == 2
    # (|W|^2) in the real space — two sublattices (XM at lattice[1], MX at lattice[2])
    realzone = RealZone([[1.0, 0.0], [0.0, 1.0]], -1=>1, -1=>1)
    Plots.savefig(Plots.plot(realzone, wannier, 1), "Plots-MoTe₂-wannier-XM.png")
    Plots.savefig(Plots.plot(realzone, wannier, 2), "Plots-MoTe₂-wannier-MX.png")
    Makie.save("Makie-MoTe₂-wannier-XM.png", Makie.plot(realzone, wannier, 1))
    Makie.save("Makie-MoTe₂-wannier-MX.png", Makie.plot(realzone, wannier, 2))

    # HoppingIntegral and TBA
    hopping = HoppingIntegral(wannier)
    tba = Algorithm(:tba, TBA(lattice, hilbert, terms(hopping; order=6, atol=1e-3, rtol=1e-3)))
    @test all(map(
        (x, y)->isapprox(x, y; atol=1e-3, rtol=1e-3),
        tba.parameters,
        [-1.208, -2.093, -0.407, -0.601, 0.144, 0.25, -0.048, -0.083, -0.015, 0.0, -0.007, -0.031, 44.583]
    ))
    # Energy band comparison: continuum model vs TBA
    recipls = reciprocals(lattice)
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

    # CoulombIntegral
    coulomb = CoulombIntegral(wannier)
    compare(ts, vs) = all(map((t, v)->isapprox(value(t), v; atol=1e-3, rtol=1e-3), ts, vs))
    @test compare(
        terms(coulomb, BareCoulomb(10.0); order=4, atol=1e-3, rtol=1e-3),
        [114.733, 34.607, 16.841, 13.779, 9.045]
    )
end
