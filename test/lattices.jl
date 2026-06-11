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

@testset "PointGroup" begin
    c₃ = C₃()
    @test angle(c₃) == angle(typeof(c₃)) == 2π/3
    @test reciprocals(c₃) == reciprocals(typeof(c₃)) == 4π/√3 * SVector(SVector(1.0, 0.0), SVector(-1/2, √3/2))

    # test sign for C₃ (uses 60° azimuth discretization to distinguish forward +1 vs reverse -1)
    ref = Bond(1, Point(1, (0.0, 0.0)), Point(2, (1.0, 0.0)))
    @test sign(C₃, ref, ref, 2) == sign(c₃, ref, ref, 2) == 1
    @test sign(C₃, ref, Bond(1, [Point(1, (0.0, 0.0)), Point(2, (-0.5, √3/2))]), 2) == 1
    @test sign(C₃, ref, Bond(1, [Point(2, (0.0, 0.0)), Point(1, (0.5, √3/2))]), 2) == -1
    @test sign(C₃, ref, Bond(1, [Point(4, (0.0, 0.0)), Point(3, (-1.0, 0.0))]), 2) == -1
    @test sign(C₃, ref, Bond(1, [Point(3, (0.0, 0.0)), Point(4, (-0.5, -√3/2))]), 2) == 1
    @test sign(C₃, ref, Bond(2, [Point(1, (0.0, 0.0)), Point(2, (1.0, 0.0))]), 2) == 0
    @test sign(C₃, ref, Bond(1, [Point(1, (0.0, 0.0)), Point(3, (1.0, 0.0))]), 2) == 0
    @test sign(C₃, ref, Bond(1, [Point(1, (0.0, 0.0)), Point(2, (√2/2, √2/2))]), 2) == 0
end

@testset "MoireReciprocalLattice" begin
    lattice = MoireTriangularReciprocal(4)
    @test getcontent(lattice, :name) == :MoireTriangularReciprocal
    @test getcontent(lattice, :vectors) == []
    @test PointGroup(lattice) == PointGroup(typeof(lattice)) == C₃()
    @test truncation(lattice) == 4
    @test reciprocals(lattice) == lattice.translations
    @test lattice.Γ ≈ SVector(2π/√3, 0.0) atol=1e-12
    @test lattice.K₊ ≈ SVector(0.0, 2π/3) atol=1e-12
    @test lattice.K₋ ≈ SVector(0.0, -2π/3) atol=1e-12
    @test lattice.translations ≈ reciprocals(C₃) atol=1e-12
    Plots.savefig(Plots.plot(lattice, 1), "Plots-moire-reciprocal-lattice.png")
    Makie.save("Makie-moire-reciprocal-lattice.png", Makie.plot(lattice, 1))
end

@testset "MoireNeighbors" begin
    # Build a Kagome lattice (3 sites per unit cell, C₃-symmetric triangular Bravais vectors)
    kagome = Lattice(
        [0.0, 0.0], [0.5, 0.0], [0.25, √3/4];
        name=:Kagome, vectors=[[1.0, 0.0], [1/2, √3/2]]
    )
    mn = MoireNeighbors{C₃}(kagome, 3)
    @test nsublattice(mn) == 3
    @test truncation(mn) == 3
    @test PointGroup(mn) == PointGroup(typeof(mn)) == C₃()
    # 11 bonds total: kind=1 (3), kind=2 (5), kind=3 (3)
    # kind=2 has 5 (not 3) because reversed-sublattice bonds with same spatial
    # direction are distinct physical bonds — see sign() fix for Kagome.
    @test length(bonds(mn)) == 11
    @test length(bonds(mn, 1)) == 3
    @test length(bonds(mn, 2)) == 5
    @test length(bonds(mn, 3)) == 3

    # in: C₃ symmetry-equivalent bonds are recognized
    # 1→2 at az=0° is equivalent to whichever 1→2 or 2→1 bond is stored for kind=1
    @test Bond(1, Point(1, (0.0, 0.0)), Point(2, (0.5, 0.0))) in mn
    # 1→2 at az=60°: with the sign() fix, reversed sublattice + even Δ is NOT
    # equivalent (it is a physically distinct bond), so this is correctly rejected
    @test !(Bond(1, Point(1, (0.0, 0.0)), Point(2, (0.5, √3/2))) in mn)
    # different kind → not equivalent
    @test !(Bond(4, Point(1, (0.0, 0.0)), Point(2, (0.5, 0.0))) in mn)
    # azimuth not at a multiple of 60°
    @test !(Bond(1, Point(1, (0.0, 0.0)), Point(2, (√2/2, √2/2))) in mn)
    # single-point bond (length≠2)
    @test !(Bond(Point(1, (0.0, 0.0))) in mn)

    # push!: deduplication against existing bonds
    n = length(bonds(mn))
    push!(mn, Bond(1, Point(1, (0.0, 0.0)), Point(2, (0.5, 0.0))))  # equivalent → skipped
    @test length(bonds(mn)) == n
    push!(mn, Bond(1, Point(1, (0.0, 0.0)), Point(2, (cosd(30), sind(30)))))  # non-equiv → added
    @test length(bonds(mn)) == n + 1
    push!(mn, Bond(Point(1, (0.0, 0.0))))  # single-point → skipped
    @test length(bonds(mn)) == n + 1
end

@testset "MoireTriangular" begin
    lattice = MoireTriangular(6)
    @test truncation(lattice) == 6
    @test lattice.coordinates == [0.0; 0.0;;]
    @test lattice.vectors ≈ reciprocals(reciprocals(C₃))
    @test truncation(lattice.neighbors) == 6

    # Neighbors bond structure (nsublattice=1, truncation=6): 7 bonds across kinds 1–6
    mn = lattice.neighbors
    @test nsublattice(mn) == 1
    @test length(bonds(mn)) == 7
    @test length(bonds(mn, 1)) == length(bonds(mn, 2)) == length(bonds(mn, 3)) == 1
    @test length(bonds(mn, 4)) == 2
    @test length(bonds(mn, 5)) == length(bonds(mn, 6)) == 1
end

@testset "MoireHoneycomb" begin
    lattice = MoireHoneycomb(6)
    vectors = reciprocals(reciprocals(C₃))
    v₁, v₂ = vectors[1], vectors[2]
    @test truncation(lattice) == 6
    @test lattice.coordinates[:, 1] ≈ (2v₁ - v₂) / 3
    @test lattice.coordinates[:, 2] ≈ (v₁ + v₂) / 3
    @test lattice.vectors ≈ vectors
    @test truncation(lattice.neighbors) == 6

    # Neighbors bond structure (nsublattice=2, truncation=6): 10 bonds across kinds 1–6
    mn = lattice.neighbors
    @test nsublattice(mn) == 2
    @test length(bonds(mn)) == 10
    # Per-kind bond counts
    @test length(bonds(mn, 1)) == 1  # inter-sublattice only
    @test length(bonds(mn, 2)) == 2  # same-sublattice
    @test length(bonds(mn, 3)) == 1  # inter-sublattice only
    @test length(bonds(mn, 4)) == 2  # inter-sublattice, 2 inequivalent directions
    @test length(bonds(mn, 5)) == 2  # same-sublattice
    @test length(bonds(mn, 6)) == 2  # same-sublattice
end
