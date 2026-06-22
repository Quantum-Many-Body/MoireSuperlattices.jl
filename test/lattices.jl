using MoireSuperlattices
using QuantumLattices
using QuantumLattices: getcontent
using StaticArrays: SVector
import CairoMakie as Makie
import Plots

@testset "CommensurateBilayerHoneycomb" begin
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

@testset "MoireReciprocalLattice" begin
    lattice = MoireTriangularReciprocal(4)
    @test getcontent(lattice, :name) == :MoireTriangularReciprocal
    @test getcontent(lattice, :vectors) == []
    @test truncation(lattice) == 4
    @test reciprocals(lattice) == lattice.translations
    @test isapprox(lattice.Γ, SVector(2π/√3, 0.0); atol=1e-12)
    @test isapprox(lattice.K₊, SVector(0.0, 2π/3); atol=1e-12)
    @test isapprox(lattice.K₋, SVector(0.0, -2π/3); atol=1e-12)
    @test lattice.translations ≈ reciprocals(MoireTriangularReciprocal)
    Plots.savefig(Plots.plot(lattice, 1), "Plots-moire-reciprocal-lattice.png")
    Makie.save("Makie-moire-reciprocal-lattice.png", Makie.plot(lattice, 1))
end

@testset "MoireTriangular" begin
    lattice = MoireTriangular()
    @test lattice.coordinates == [0.0; 0.0;;]
    @test lattice.vectors ≈ reciprocals(reciprocals(MoireTriangular))
end

@testset "MoireHoneycomb" begin
    lattice = MoireHoneycomb()
    vectors = reciprocals(reciprocals(MoireHoneycomb))
    v₁, v₂ = vectors[1], vectors[2]
    @test lattice.coordinates[:, 1] ≈ (2v₁ - v₂) / 3
    @test lattice.coordinates[:, 2] ≈ (v₁ + v₂) / 3
    @test lattice.vectors ≈ vectors
end
