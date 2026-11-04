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
    parameters = (a₀=3.28, m=0.45, θ=3.70, Vᶻ=76.0, μ=0.0, V=-1.28, ψ=22.7, w=-12.9)
    bltmd = Algorithm(:BLTMD, BLTMD(values(parameters)...; truncation=4), parameters)
    @test contentnames(typeof(bltmd.frontend)) == (:parameters, :reciprocallattice, :diagonal!, :system, :quadraticization, :H)
    @test Parameters(bltmd.frontend) == (a₀=3.28, m=0.45, θ=3.7, Vᶻ=76.0, μ=0.0, potentialᵣ=-1.1808487545039954, potentialᵢ=-0.4939597341751278, interlayer₁=-12.9, interlayer₂=-12.9, interlayer₃=-12.9)
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

@testset "TLTMD" begin
    parameters = (a₀=3.317, m=0.43, θ=3.9, Vᶻ=0.0, μ=0.0, V=5.966, ψ=77.02, w=14.32)
    tltmd = Algorithm(:TLTMD, TLTMD(values(parameters)...; stacking=:αβ, truncation=4), parameters)
    @test contentnames(typeof(tltmd.frontend)) == (:parameters, :reciprocallattice, :diagonal!, :system, :quadraticization, :H)
    @test Parameters(tltmd.frontend) == (
        a₀=3.317, m=0.43, θ=3.9, Vᶻ=0.0, μ=0.0, ϕ=120.0,
        potential21minusᵣ=-5.704704954571016, potential21minusᵢ=-1.7462810144111702,
        potential21plusᵣ=4.364676197912042, potential21plusᵢ=-4.067278904547865,
        potential32minusᵣ=4.364676197912042, potential32minusᵢ=-4.067278904547865,
        potential32plusᵣ=-5.704704954571016, potential32plusᵢ=-1.7462810144111702,
        interlayer21₀=14.32, interlayer21ᵣ=-7.16, interlayer21ᵢ=12.40148378219316,
        interlayer32₀=14.32, interlayer32ᵣ=-7.16, interlayer32ᵢ=12.40148378219316,
    )
    @test dimension(tltmd.frontend) == 183
    @test count(tltmd) == count(tltmd.frontend) == 3
    @test_throws AssertionError TLTMD(values(parameters)...; stacking=:ab, truncation=1)

    recipls = reciprocals(tltmd.frontend.reciprocallattice)
    path = ReciprocalPath(recipls, hexagon"M₆-K₂-Γ-K₅-M₃", length=100)
    bands_αβ = tltmd(:EB, EnergyBands(path, dimension(tltmd)-4:dimension(tltmd)))

    update!(tltmd; Vᶻ=0.0, μ=0.0, V=5.966, ψ=77.02, w=14.32, ϕ=-120.0)
    @test tltmd.frontend.parameters.ϕ == -120.0
    bands_βα = tltmd(:EB, EnergyBands(path, dimension(tltmd)-4:dimension(tltmd)))

    emin, emax = -70.0, 25.0
    plt = Plots.plot(layout=(2, 1), size=(700, 700), link=:y)
    Plots.plot!(plt[1], bands_αβ, ylims=(emin, emax), color="blue")
    Plots.plot!(plt[2], bands_βα, ylims=(emin, emax), color="green")
    Plots.savefig(plt, "Plots-TLTMD-continuum.png")

    fig = Makie.Figure()
    ax_αβ = Makie.Axis(fig[1, 1], title="αβ")
    ax_βα = Makie.Axis(fig[2, 1], title="βα")
    Makie.plot!(ax_αβ, bands_αβ; ylims=(emin, emax), color=:blue)
    Makie.plot!(ax_βα, bands_βα; ylims=(emin, emax), color=:green)
    Makie.save("Makie-TLTMD-continuum.png", fig)
end
