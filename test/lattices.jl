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
    c₆ = C₆()
    @test angle(c₆) == angle(typeof(c₆)) == π/3
    @test reciprocals(c₆) == reciprocals(typeof(c₆)) == 4π/√3 * SVector(SVector(1.0, 0.0), SVector(-1/2, √3/2))

    # test sign for C₆
    ref = Bond(1, Point(1, (0.0, 0.0)), Point(2, (1.0, 0.0)))
    @test sign(C₆, ref, ref, 2) == sign(c₆, ref, ref, 2) == 1
    @test sign(C₆, ref, Bond(1, [Point(1, (0.0, 0.0)), Point(2, (-0.5, √3/2))]), 2) == 1
    @test sign(C₆, ref, Bond(1, [Point(2, (0.0, 0.0)), Point(1, (0.5, √3/2))]), 2) == -1
    @test sign(C₆, ref, Bond(1, [Point(4, (0.0, 0.0)), Point(3, (-1.0, 0.0))]), 2) == -1
    @test sign(C₆, ref, Bond(1, [Point(3, (0.0, 0.0)), Point(4, (-0.5, -√3/2))]), 2) == 1
    @test sign(C₆, ref, Bond(2, [Point(1, (0.0, 0.0)), Point(2, (1.0, 0.0))]), 2) == 0
    @test sign(C₆, ref, Bond(1, [Point(1, (0.0, 0.0)), Point(3, (1.0, 0.0))]), 2) == 0
    @test sign(C₆, ref, Bond(1, [Point(1, (0.0, 0.0)), Point(2, (√2/2, √2/2))]), 2) == 0
end

@testset "MoireReciprocalLattice" begin
    lattice = MoireTriangularReciprocal(4)
    @test getcontent(lattice, :name) == :MoireTriangularReciprocal
    @test getcontent(lattice, :vectors) == []
    @test PointGroup(lattice) == PointGroup(typeof(lattice)) == C₆()
    @test truncation(lattice) == 4
    @test reciprocals(lattice) == lattice.translations
    @test lattice.Γ ≈ SVector(2π/√3, 0.0) atol=1e-12
    @test lattice.K₊ ≈ SVector(0.0, 2π/3) atol=1e-12
    @test lattice.K₋ ≈ SVector(0.0, -2π/3) atol=1e-12
    @test lattice.translations ≈ reciprocals(C₆) atol=1e-12
    Plots.savefig(Plots.plot(lattice, 1), "Plots-moire-reciprocal-lattice.png")
    Makie.save("Makie-moire-reciprocal-lattice.png", Makie.plot(lattice, 1))
end

@testset "MoireNeighbors" begin
    # Build a Kagome lattice (3 sites per unit cell, C₆-symmetric triangular Bravais vectors)
    kagome = Lattice(
        [0.0, 0.0], [0.5, 0.0], [0.25, √3/4];
        name=:Kagome, vectors=[[1.0, 0.0], [1/2, √3/2]]
    )
    mn = MoireNeighbors{C₆}(kagome, 3)
    @test nsublattice(mn) == 3
    @test truncation(mn) == 3
    @test PointGroup(mn) == PointGroup(typeof(mn)) == C₆()
    # 9 bonds total: 3 kinds × 3 sublattice pairs (1 bond per pair, per kind)
    @test length(bonds(mn)) == 9
    # kind=1 (dist≈0.5): different-sublattice pairs, 3 bonds
    @test length(bonds(mn, 1)) == 3
    @test Set(pairs(mn, 1)) == Set([(1, 2), (1, 3), (3, 2)])
    for pair in [(1, 2), (1, 3), (2, 3)]
        @test length(bonds(mn, 1, pair)) == 1
    end
    # kind=2 (dist≈0.866): different-sublattice pairs, 3 bonds
    @test length(bonds(mn, 2)) == 3
    @test Set(pairs(mn, 2)) == Set([(1, 2), (1, 3), (3, 2)])
    for pair in [(1, 2), (1, 3), (2, 3)]
        @test length(bonds(mn, 2, pair)) == 1
    end
    # kind=3 (dist≈1.0): same-sublattice pairs, 3 bonds
    @test length(bonds(mn, 3)) == 3
    @test Set(pairs(mn, 3)) == Set([(1, 1), (2, 2), (3, 3)])
    for pair in [(1, 1), (2, 2), (3, 3)]
        @test length(bonds(mn, 3, pair)) == 1
    end

    # in: C₆ symmetry-equivalent bonds are recognized
    # kind=1 site1→site2 at az=0° is equivalent to stored site2→site1 at az=180°
    @test Bond(1, Point(1, (0.0, 0.0)), Point(2, (0.5, 0.0))) in mn
    # kind=1 site1→site2 at az=60° is equivalent (swapped, Δ=-2)
    @test Bond(1, Point(1, (0.0, 0.0)), Point(2, (0.5, √3/2))) in mn
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
    @test lattice.vectors ≈ reciprocals(reciprocals(C₆))
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
    vectors = reciprocals(reciprocals(C₆))
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
    # Per-kind bond counts and sublattice pair structure
    @test length(bonds(mn, 1)) == 1 && pairs(mn, 1) == [(1, 2)] # inter-sublattice only
    @test length(bonds(mn, 2)) == 2 && Set(pairs(mn, 2)) == Set([(1, 1), (2, 2)])  # same-sublattice
    @test length(bonds(mn, 3)) == 1 && pairs(mn, 3) == [(1, 2)] # inter-sublattice only
    @test length(bonds(mn, 4)) == 2 && Set(pairs(mn, 4)) == Set([(1, 2)]) # inter-sublattice, 2 inequivalent directions
    @test length(bonds(mn, 5)) == 2 && Set(pairs(mn, 5)) == Set([(1, 1), (2, 2)])  # same-sublattice
    @test length(bonds(mn, 6)) == 2 && Set(pairs(mn, 6)) == Set([(1, 1), (2, 2)])  # same-sublattice
end
