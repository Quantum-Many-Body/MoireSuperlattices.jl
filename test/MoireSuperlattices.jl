using MoireSuperlattices
using Plots
using QuantumLattices
using QuantumLattices: contentnames, getcontent, shape
using TightBindingApproximation

@time @testset "CommensurateBilayerHoneycomb" begin
    moire = CommensurateBilayerHoneycomb((20, 3); stack=:AB, center=:carbon)
    @test rad2deg(angle(moire)) ≈ 4.613282862403654
    top, bottom = Lattice(moire, :top), Lattice(moire, :bottom)
    @test volume(top.vectors...) ≈ volume(bottom.vectors...)
    @test volume(vectors(moire)...)/volume(top.vectors...) ≈ count(moire) == 463
    savefig(plot(moire, :real; xlim=(-20, 20), ylim=(-10, 20)), "twisted-honeycomb-(20, 3)-coordinate.png")
    savefig(plot(moire, :reciprocal; xlim=(-5, 5), ylim=(-4, 4)), "twisted-honeycomb-(20, 3)-reciprocal.png")

    moire = CommensurateBilayerHoneycomb((8, 1); stack=:AA, center=:carbon)
    @test rad2deg(angle(moire)) ≈ 3.8902381690076835
    top, bottom = Lattice(moire, :top), Lattice(moire, :bottom)
    @test volume(top.vectors...) ≈ volume(bottom.vectors...)
    @test volume(vectors(moire)...)/volume(top.vectors...) ≈ count(moire) == 217
    savefig(plot(moire, :real; xlim=(-20, 20), ylim=(-10, 20)), "twisted-honeycomb-(8, 1)-coordinate.png")
    savefig(plot(moire, :reciprocal; xlim=(-5, 5), ylim=(-4, 4)), "twisted-honeycomb-(8, 1)-reciprocal.png")
end

@time @testset "MoireReciprocalLattice" begin
    lattice = MoireReciprocalLattice(4)
    @test getcontent(lattice, :name) == :truncation
    @test getcontent(lattice, :vectors) == []
    savefig(plot(lattice, 1), "moire-reciprocal-lattice.png")
end

@time @testset "MoireTriangular" begin
    lattice = MoireTriangular(6, [[1.0, 0.0], [0.5, √3/2]]; origin=[0.0, 0.0])
    @test truncation(lattice) == truncation(typeof(lattice)) == 6
    @test lattice.coordinates == [0.0; 0.0;;]
    @test lattice.vectors == [[1.0, 0.0], [0.5, √3/2]]
    @test all(map(≈, lattice.neighbors, ([[0.5, √3/2]], [[0.0, √3]], [[1.0, √3]], [[0.5, 1.5√3], [-0.5, 1.5√3]], [[1.5, 1.5√3]], [[0.0, 2√3]])))
end

@time @testset "MoireSpinor and MoireSpace" begin
    @test MoireSpinor(1, 1, 1, 1//2, 1)' == MoireSpinor(1, 1, 1, 1//2, 2)
    @test MoireSpinor(2, 1, 1, -1//2, 2)' == MoireSpinor(2, 1, 1, -1//2, 1)
    @test string(MoireSpinor(:, 1, 1, 1//2, 1)) == "MoireSpinor(:, 1, 1, 1//2, 1)"
    @test string(MoireSpinor(1, :, 1, 0//2, 2)) == "MoireSpinor(1, :, 1, 0, 2)"
    @test statistics(MoireSpinor(:, :, :, :, :)) == statistics(MoireSpinor) == :f
    @test !isdefinite(MoireSpinor(:, :, :, :, :)) && isdefinite(MoireSpinor(1, 1, 1, 1//2, 2))
    @test allequalfields(MoireSpinor) == (:valley, :layer, :sublattice, :spin)
    @test indextype(MoireSpinor, Colon, Colon, Colon, Colon, Colon) == MoireSpinor{Colon, Colon, Colon, Colon, Colon}
    @test patternrule((:, :, :, :), Val(:), MoireSpinor, Val(:nambu)) == (2, 1, 2, 1)
    @test MoireSpinor{Colon, Colon, Colon, Colon, Colon}(1, 1, 1, 1//2, 2) == MoireSpinor(1, 1, 1, 1//2, 2)
    @test script(MoireSpinor(1, 1, 1, 1//2, 2), Val(:valley))=="1" && script(MoireSpinor(:, 1, 1, 1//2, 2), Val(:valley))==":"
    @test script(MoireSpinor(1, 1, 1, 1//2, 2), Val(:layer))=="1" && script(MoireSpinor(1, :, 1, 1//2, 2), Val(:layer))==":"
    @test script(MoireSpinor(1, 1, 1, 1//2, 2), Val(:sublattice))=="1" && script(MoireSpinor(1, 1, :, 1//2, 2), Val(:sublattice))==":"
    @test script(MoireSpinor(1, 1, 1, 1//2, 2), Val(:spin))=="↑" && script(MoireSpinor(1, 1, 1, -1//2, 2), Val(:spin))=="↓"
    @test script(MoireSpinor(1, 1, 1, 0, 2), Val(:spin))=="" && script(MoireSpinor(1, 1, 1, :, 2), Val(:spin))==":"
    @test script(MoireSpinor(1, 1, 1, 0, 2), Val(:nambu))=="\\dagger" && script(MoireSpinor(1, 1, 1, 0, 1), Val(:nambu))=="" && script(MoireSpinor(1, 1, 1, 0, :), Val(:nambu))==":"
    @test latexname(MoireSpinor{Colon, Colon, Colon, Colon, Colon}) == :MoireSpinor
    @test latexname(Index{MoireSpinor{Colon, Colon, Colon, Colon, Colon}, Int}) == Symbol("Index{MoireSpinor}")
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
    bltmd = Algorithm(:BLTMD, BLTMD(values(parameters)...; truncation=4); parameters=parameters, map=bltmdmap)
    @test contentnames(typeof(bltmd.frontend)) == (:parameters, :reciprocallattice, :diagonal!, :system, :quadraticization, :H)
    @test Parameters(bltmd.frontend) == (a₀=3.28, m=0.45, θ=3.7, Vᶻ=38.0, μ=0.0, potentialᵣ=-1.1808487545039954, potentialᵢ=-0.4939597341751278, interlayer₁=-12.9, interlayer₂=-12.9, interlayer₃=-12.9)
    @test dimension(bltmd.frontend) == 122
    @test count(bltmd.frontend) == 2

    update!(bltmd; μ=8.31)
    recipls = bltmd.frontend.reciprocallattice.translations
    lattice = MoireTriangular(6, reciprocals(recipls))
    hilbert = Hilbert(Fock{:f}(1, 2), length(lattice))
    tba = Algorithm(:tba, TBA(lattice, hilbert, terms(bltmd, lattice, BrillouinZone(recipls, 24); tol=10^-6)))
    @test coefficients(bltmd, lattice, BrillouinZone(recipls, 24)) == coefficients(bltmd.frontend, lattice, BrillouinZone(recipls, 24))
    @test all(map((x, y)->isapprox(x, y; atol=10^-6), tba.parameters, [-2.7598267, -4.3678292, -1.3035002, 0.0, 0.2447067, -0.6094541, -0.2700574, -0.3888003, 0.0260020, -0.0308282, -0.2999706, 0.0, 10.2102190]))

    plt = plot()
    emin, emax = -40.0, 40.0
    plot!(plt, bltmd(:EB, EnergyBands(ReciprocalPath(recipls, hexagon"Γ-K₁-M₁-Γ", length=100))), ylim=(emin, emax), color="blue", title="")
    plot!(plt, bltmd(:EB, EnergyBands(ReciprocalPath(recipls, hexagon"Γ-K₄-M₄-Γ", length=100))), ylim=(emin, emax), color="green", title="")
    plot!(plt, tba(:EB, EnergyBands(ReciprocalPath(recipls, hexagon"Γ-K-M-Γ", length=100))), ylim=(emin, emax), ls=:dash, color="red", size=(400, 300), title="")
    savefig("WeSe₂-AA-stack.png")
end
