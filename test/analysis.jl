using MoireSuperlattices
using QuantumLattices
using StaticArrays: SVector

@testset "MoireWannier-triangular" begin
    parameters = (a₀=3.28, m=0.45, θ=3.70, Vᶻ=38.0, μ=0.0, V=-1.28, ψ=22.7, w=-12.9)
    bltmd = Algorithm(:BLTMD, BLTMD(values(parameters)...; truncation=4), parameters)
    update!(bltmd; μ=8.31)
    recipls = bltmd.frontend.reciprocallattice.translations
    lattice = MoireTriangular(6)
    bz = BrillouinZone(recipls, 12)
    w = MoireWannier(bltmd.frontend, lattice, bz; band=dimension(bltmd.frontend))
    @test size(w.energies) == (1, length(bz))
    @test size(w.bloch) == (2, dimension(bltmd.frontend)÷2, 1, length(bz))
    @test size(w.U) == (1, 1, length(bz))
    val = w([0.0, 0.0], 1)
    @test length(val) == 2
    @test eltype(val) == ComplexF64
end

@testset "HoppingIntegral" begin
    parameters = (a₀=3.28, m=0.45, θ=3.70, Vᶻ=38.0, μ=0.0, V=-1.28, ψ=22.7, w=-12.9)
    bltmd = Algorithm(:BLTMD, BLTMD(values(parameters)...; truncation=4), parameters)
    update!(bltmd; μ=8.31)
    recipls = bltmd.frontend.reciprocallattice.translations
    lattice = MoireTriangular(6)
    bz = BrillouinZone(recipls, 12)
    w = MoireWannier(bltmd.frontend, lattice, bz; band=dimension(bltmd.frontend))
    h = HoppingIntegral(w)
    # onsite
    t0 = h(SVector(0.0, 0.0))
    @test t0 isa Matrix{ComplexF64}
    @test size(t0) == (1, 1)
    @test abs(imag(t0[1,1])) < 1e-10  # onsite should be real
    # nearest-neighbor hopping
    t1 = h(icoordinate(bonds(lattice.neighbors, 1)[1]))
    @test t1 isa Matrix{ComplexF64}
    @test size(t1) == (1, 1)
end
