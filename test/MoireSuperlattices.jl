using MoireSuperlattices
using Plots
using QuantumLattices
using TightBindingApproximation

@testset "CommensurateBilayerHoneycomb" begin
    moire = CommensurateBilayerHoneycomb((8, 1); stack=:AA, center=:carbon)
    savefig(plot(moire, :lattice; xlim=(-20, 20), ylim=(-10, 20)), "twisted-graphene-coordinate.png")
    savefig(plot(moire, :reciprocal; xlim=(-5, 5), ylim=(-4, 4)), "twisted-graphene-reciprocal.png")
end

@testset "MoireReciprocalLattice" begin
    savefig(plot(MoireReciprocalLattice(4), 1), "moire-reciprocal-lattice.png")
end

@testset "BLTMD" begin
    parameters = (a₀=3.28, m=0.45, θ=3.70, Vᶻ=38.0, μ=8.31, V=-1.28, ψ=22.7, w=-12.9)
    bltmd = Algorithm(:BLTMD, BLTMD(values(parameters)...; truncation=4); parameters=parameters, map=bltmdmap)

    emin, emax = -40.0, 40.0
    recipls = bltmd.frontend.reciprocallattice.translations
    lattice = MoireTriangular(6, reciprocals(recipls))
    hilbert = Hilbert(1=>Fock{:f}(1, 2))
    brillouinzone = BrillouinZone(recipls, 24)
    tba = Algorithm(:tba, TBA(lattice, hilbert, terms(bltmd, lattice, brillouinzone; atol=10^-6)))
    @test all(map((x, y)->isapprox(x, y, atol=10^-6), collect(tba.parameters), [-2.7598267, -4.3678292, -1.3035002, 0.0, 0.2447067, -0.6094541, -0.2700574, -0.3888003, 0.0260020, -0.0308282, -0.2999706, 0.0, 10.2102190]))
    plt = plot()
    plot!(plt, bltmd(:EB, EnergyBands(ReciprocalPath(recipls, hexagon"Γ-K₁-M₁-Γ", length=100))), ylim=(emin, emax), color="blue", title="")
    plot!(plt, bltmd(:EB, EnergyBands(ReciprocalPath(recipls, hexagon"Γ-K₄-M₄-Γ", length=100))), ylim=(emin, emax), color="green", title="")
    plot!(plt, tba(:EB, EnergyBands(ReciprocalPath(recipls, hexagon"Γ-K-M-Γ", length=100))), ylim=(emin, emax), ls=:dash, color="red", size=(400, 300), title="")
    savefig("WeSe₂-AA-stack.png")
end