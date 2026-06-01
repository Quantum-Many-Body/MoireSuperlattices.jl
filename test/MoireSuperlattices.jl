using MoireSuperlattices
import Plots
import CairoMakie as Makie
using QuantumLattices
using QuantumLattices: contentnames, getcontent
using TightBindingApproximation

@time @testset "CommensurateBilayerHoneycomb" begin
    moire = CommensurateBilayerHoneycomb((20, 3); stack=:AB, center=:carbon)
    @test rad2deg(angle(moire)) ≈ 4.613282862403654
    top, bottom = Lattice(moire, :top), Lattice(moire, :bottom)
    @test volume(top.vectors...) ≈ volume(bottom.vectors...)
    @test volume(vectors(moire)...)/volume(top.vectors...) ≈ count(moire) == 463
    Plots.savefig(Plots.plot(moire, :real; xlims=(-20, 20), ylims=(-10, 20)), "Plots-twisted-honeycomb-(20, 3)-coordinate.png")
    Makie.save("Makie-twisted-honeycomb-(20, 3)-coordinate.png", Makie.plot(moire, :real; limits=(-20, 20, -10, 20)))
    Plots.savefig(Plots.plot(moire, :reciprocal; xlims=(-5, 5), ylims=(-4, 4)), "Plots-twisted-honeycomb-(20, 3)-reciprocal.png")
    Makie.save("Makie-twisted-honeycomb-(20, 3)-reciprocal.png", Makie.plot(moire, :reciprocal; limits=(-5, 5, -4, 4)))

    moire = CommensurateBilayerHoneycomb((8, 1); stack=:AA, center=:carbon)
    @test rad2deg(angle(moire)) ≈ 3.8902381690076835
    top, bottom = Lattice(moire, :top), Lattice(moire, :bottom)
    @test volume(top.vectors...) ≈ volume(bottom.vectors...)
    @test volume(vectors(moire)...)/volume(top.vectors...) ≈ count(moire) == 217
    Plots.savefig(Plots.plot(moire, :real; xlims=(-20, 20), ylims=(-10, 20)), "Plots-twisted-honeycomb-(8, 1)-coordinate.png")
    Makie.save("Makie-twisted-honeycomb-(8, 1)-coordinate.png", Makie.plot(moire, :real; limits=(-20, 20, -10, 20)))
    Plots.savefig(Plots.plot(moire, :reciprocal; xlims=(-5, 5), ylims=(-4, 4)), "Plots-twisted-honeycomb-(8, 1)-reciprocal.png")
    Makie.save("Makie-twisted-honeycomb-(8, 1)-reciprocal.png", Makie.plot(moire, :reciprocal; limits=(-5, 5, -4, 4)))
end

@time @testset "MoireReciprocalLattice" begin
    lattice = MoireTriangularReciprocal(4)
    @test getcontent(lattice, :name) == :MoireTriangularReciprocal
    @test getcontent(lattice, :vectors) == []
    Plots.savefig(Plots.plot(lattice, 1), "Plots-moire-reciprocal-lattice.png")
    Makie.save("Makie-moire-reciprocal-lattice.png", Makie.plot(lattice, 1))
end

@time @testset "MoireTriangular" begin
    lattice = MoireTriangular(6)
    @test truncation(lattice) == truncation(typeof(lattice)) == 6
    @test lattice.coordinates == [0.0; 0.0;;]
    @test lattice.vectors ≈ reciprocals(reciprocals(C₆))
    @test length(lattice.neighbors) == 6
end

@time @testset "MoireSpinor and MoireSpace" begin
    @test MoireSpinor(1, 1, 1, 1//2, 1)' == MoireSpinor(1, 1, 1, 1//2, 2)
    @test MoireSpinor(2, 1, 1, -1//2, 2)' == MoireSpinor(2, 1, 1, -1//2, 1)
    @test string(MoireSpinor(:, 1, 1, 1//2, 1)) == "MoireSpinor(:, 1, 1, 1//2, 1)"
    @test string(MoireSpinor(1, :, 1, 0//2, 2)) == "MoireSpinor(1, :, 1, 0, 2)"
    @test statistics(MoireSpinor(:, :, :, :, 1)) == statistics(MoireSpinor) == :f
    @test !isdefinite(MoireSpinor(:, :, :, :, 1)) && isdefinite(MoireSpinor(1, 1, 1, 1//2, 2))
    @test indextype(MoireSpinor, Colon, Colon, Colon, Colon) == MoireSpinor{Colon, Colon, Colon, Colon}
    @test MoireSpinor{Colon, Colon, Colon, Colon}(1, 1, 1, 1//2, 2) == MoireSpinor(1, 1, 1, 1//2, 2)
    @test script(MoireSpinor(1, 1, 1, 1//2, 2), Val(:valley))=="1" && script(MoireSpinor(:, 1, 1, 1//2, 2), Val(:valley))==":"
    @test script(MoireSpinor(1, 1, 1, 1//2, 2), Val(:layer))=="1" && script(MoireSpinor(1, :, 1, 1//2, 2), Val(:layer))==":"
    @test script(MoireSpinor(1, 1, 1, 1//2, 2), Val(:sublattice))=="1" && script(MoireSpinor(1, 1, :, 1//2, 2), Val(:sublattice))==":"
    @test script(MoireSpinor(1, 1, 1, 1//2, 2), Val(:spin))=="↑" && script(MoireSpinor(1, 1, 1, -1//2, 2), Val(:spin))=="↓"
    @test script(MoireSpinor(1, 1, 1, 0, 2), Val(:spin))=="" && script(MoireSpinor(1, 1, 1, :, 2), Val(:spin))==":"
    @test script(MoireSpinor(1, 1, 1, 0, 2), Val(:nambu))=="\\dagger" && script(MoireSpinor(1, 1, 1, 0, 1), Val(:nambu))==""
    @test latexname(MoireSpinor{Colon, Colon, Colon, Colon}) == :MoireSpinor
    @test latexname(Index{MoireSpinor{Colon, Colon, Colon, Colon}, Int}) == Symbol("Index{MoireSpinor}")
    @test latexname(CompositeIndex{<:Index{<:MoireSpinor}}) == Symbol("CompositeIndex{Index{MoireSpinor}}")

    moire = MoireSpace(2, 2, 2, 2)
    @test shape(moire) == (1:2, 1:2, 1:2, 1:2, 1:2)
    for spinor in moire
        @test convert(MoireSpinor, convert(CartesianIndex, spinor, moire), moire) == spinor
    end
    @test shape(moire, MoireSpinor(1, 1, 1, 1//2, 1)) == (1:1, 1:1, 1:1, 2:2, 1:1)
    @test shape(moire, MoireSpinor(:, :, :, -1//2, 2)) == (1:2, 1:2, 1:2, 1:1, 2:2)
    @test shape(moire, MoireSpinor(:, :, 2, :, 2)) == (1:2, 1:2, 2:2, 1:2, 2:2)
end

@time @testset "BLTMD" begin
    parameters = (a₀=3.28, m=0.45, θ=3.70, Vᶻ=38.0, μ=0.0, V=-1.28, ψ=22.7, w=-12.9)
    bltmd = Algorithm(:BLTMD, BLTMD(values(parameters)...; truncation=4), parameters)
    @test contentnames(typeof(bltmd.frontend)) == (:parameters, :reciprocallattice, :diagonal!, :system, :quadraticization, :H)
    @test Parameters(bltmd.frontend) == (a₀=3.28, m=0.45, θ=3.7, Vᶻ=38.0, μ=0.0, potentialᵣ=-1.1808487545039954, potentialᵢ=-0.4939597341751278, interlayer₁=-12.9, interlayer₂=-12.9, interlayer₃=-12.9)
    @test dimension(bltmd.frontend) == 122
    @test count(bltmd.frontend) == 2

    update!(bltmd; μ=8.31)
    recipls = bltmd.frontend.reciprocallattice.translations
    lattice = MoireTriangular(6)
    hilbert = Hilbert(Fock{:f}(1, 2), length(lattice))
    tba = Algorithm(:tba, TBA(lattice, hilbert, terms(bltmd, lattice, BrillouinZone(recipls, 24); tol=10^-6)))
    @test coefficients(bltmd, lattice, BrillouinZone(recipls, 24)) == coefficients(bltmd.frontend, lattice, BrillouinZone(recipls, 24))
    @test all(map((x, y)->isapprox(x, y; atol=10^-6), tba.parameters, [-2.7598267, -4.3678292, -1.3035002, 0.0, 0.2447067, -0.6094541, -0.2700574, -0.3888003, 0.0260020, -0.0308282, -0.2999706, 0.0, 10.2102190]))

    plt = Plots.plot()
    emin, emax = -40.0, 40.0
    Plots.plot!(plt, bltmd(:EB, EnergyBands(ReciprocalPath(recipls, hexagon"Γ-K₁-M₁-Γ", length=100))), ylims=(emin, emax), color="blue", title="")
    Plots.plot!(plt, bltmd(:EB, EnergyBands(ReciprocalPath(recipls, hexagon"Γ-K₄-M₄-Γ", length=100))), ylims=(emin, emax), color="green", title="")
    Plots.plot!(plt, tba(:EB, EnergyBands(ReciprocalPath(recipls, hexagon"Γ-K-M-Γ", length=100))), ylims=(emin, emax), ls=:dash, color="red", size=(400, 300), title="")
    Plots.savefig("Plots-WeSe₂-AA-stack.png")

    fig = Makie.Figure()
    ax = Makie.Axis(fig[1, 1])
    Makie.plot!(ax, bltmd(:EB, EnergyBands(ReciprocalPath(recipls, hexagon"Γ-K₁-M₁-Γ", length=100))); ylims=(emin, emax), color=:blue, title="")
    Makie.plot!(ax, bltmd(:EB, EnergyBands(ReciprocalPath(recipls, hexagon"Γ-K₄-M₄-Γ", length=100))); ylims=(emin, emax), color=:green, title="")
    Makie.plot!(ax, tba(:EB, EnergyBands(ReciprocalPath(recipls, hexagon"Γ-K-M-Γ", length=100))); ylims=(emin, emax), linestyle=:dash, color=:red, title="")
    Makie.save("Makie-WeSe₂-AA-stack.png", fig)
end
