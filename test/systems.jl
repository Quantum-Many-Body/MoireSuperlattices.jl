using MoireSuperlattices
using QuantumLattices
using QuantumLattices: contentnames
using TightBindingApproximation
import CairoMakie as Makie
import Plots

@testset "MoireSpinor and MoireSpace" begin
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

@testset "BLTMD" begin
    parameters = (a₀=3.28, m=0.45, θ=3.70, Vᶻ=38.0, μ=0.0, V=-1.28, ψ=22.7, w=-12.9)
    bltmd = Algorithm(:BLTMD, BLTMD(values(parameters)...; truncation=4), parameters)
    @test contentnames(typeof(bltmd.frontend)) == (:parameters, :reciprocallattice, :diagonal!, :system, :quadraticization, :H)
    @test Parameters(bltmd.frontend) == (a₀=3.28, m=0.45, θ=3.7, Vᶻ=38.0, μ=0.0, potentialᵣ=-1.1808487545039954, potentialᵢ=-0.4939597341751278, interlayer₁=-12.9, interlayer₂=-12.9, interlayer₃=-12.9)
    @test dimension(bltmd.frontend) == 122
    @test count(bltmd) == count(bltmd.frontend) == 2

    update!(bltmd; a₀=3.30, m=0.45, θ=4.0, Vᶻ=0.0, μ=0.0, V=4.4, ψ=5.9, w=20.0)
    recipls = reciprocals(bltmd.frontend.reciprocallattice)
    bands₁ = bltmd(:EB, EnergyBands(ReciprocalPath(recipls, hexagon"Γ-K₁-K₂-Γ", length=100)))
    bands₂ = bltmd(:EB, EnergyBands(ReciprocalPath(recipls, hexagon"Γ-K₂-K₁-Γ", length=100)))

    plt = Plots.plot()
    emin, emax = -100.0, 30.0
    Plots.plot!(plt, bands₁, ylims=(emin, emax), color="blue", title="")
    Plots.plot!(plt, bands₂, ylims=(emin, emax), color="green", title="")
    Plots.savefig(plt, "Plots-WeSe₂-continuum.png")

    fig = Makie.Figure()
    ax = Makie.Axis(fig[1, 1])
    Makie.plot!(ax, bands₁; ylims=(emin, emax), color=:blue, title="")
    Makie.plot!(ax, bands₂; ylims=(emin, emax), color=:green, title="")
    Makie.save("Makie-WeSe₂-continuum.png", fig)
end
