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
