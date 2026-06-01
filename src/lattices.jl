"""
    CommensurateBilayerHoneycomb

Commensurate Moire superlattice composed of two layers of honeycomb lattices.
"""
struct CommensurateBilayerHoneycomb
    characters::Tuple{Int, Int}
    displacement::SVector{2, Float64}
    center::SVector{2, Float64}
    coordinates::Matrix{Float64}
    vectors::SVector{2, SVector{2, Float64}}
    function CommensurateBilayerHoneycomb(characters::Tuple{Int, Int}; stack::Symbol=:AA, center::Symbol=:carbon)
        @assert gcd(characters[1], characters[2])==1 "CommensurateBilayerHoneycomb error: not coprime integers."
        @assert stack∈(:AA, :AB) "CommensurateBilayerHoneycomb error: not supported stack."
        @assert center∈(:carbon, :hexagon) "CommensurateBilayerHoneycomb error: not supported center."
        a₁, a₂ = SVector(1.0, 0.0), SVector(-0.5, √3/2)
        coordinates = [1/3*a₁+2/3*a₂;; 2/3*a₁+1/3*a₂]
        displacement = stack==:AA ? SVector(0.0, 0.0) : -1/3*a₁+1/3*a₂
        center = center==:carbon ? SVector(coordinates[1, 1], coordinates[2, 1]) : SVector(0.0, 0.0)
        new(characters, displacement, center, coordinates, SVector(a₁, a₂))
    end
end

"""
    Lattice(moire::CommensurateBilayerHoneycomb, type::Symbol)

Get the minimum unit of the top/bottom layer of a commensurate Moire superlattice composed of two layers of honeycomb lattices.
"""
@inline function Lattice(moire::CommensurateBilayerHoneycomb, type::Symbol)
    @assert type∈(:top, :bottom) "Lattice error: incorrect type (`:$type`), which should be either `:top` or `:bottom`."
    m, r, θ = moire.characters..., angle(moire)
    if type == :top
        coordinates = rotate(moire.coordinates, +θ/2; axis=(moire.center, (0, 0)))
        vectors = map(v->SVector{2}(rotate(v, +θ/2)), moire.vectors)
        return Lattice(:top, coordinates, vectors)
    else
        coordinates = rotate(moire.coordinates.+reshape(moire.displacement, :, 1), -θ/2; axis=(moire.center, (0, 0)))
        vectors = map(v->SVector{2}(rotate(v, -θ/2)), moire.vectors)
        return Lattice(:bottom, coordinates, vectors)
    end
end

"""
    angle(moire::CommensurateBilayerHoneycomb) -> Float64

Get the twist angle of a commensurate Moire superlattice composed of two layers of honeycomb lattices.
"""
@inline function Base.angle(moire::CommensurateBilayerHoneycomb)
    m, r = moire.characters
    return acos((3m^2+3m*r+r^2/2)/(3m^2+3m*r+r^2))
end

"""
    vectors(moire::CommensurateBilayerHoneycomb) -> SVector{2, SVector{2, Float64}}

Get the translation vectors of a commensurate Moire superlattice composed of two layers of honeycomb lattices.
"""
@inline function vectors(moire::CommensurateBilayerHoneycomb)
    m, r, θ = moire.characters..., angle(moire)
    v₁, v₂ = map(v->SVector{2}(rotate(v, -θ/2)), moire.vectors)
    if r%3 == 0
        return SVector((m+r÷3)*v₁+(m+2r÷3)*v₂, -r÷3*v₁+(m+r÷3)*v₂)
    else
        return SVector(m*v₁+(2m+r)*v₂, -(m+r)*v₁+m*v₂)
    end
end

"""
    count(moire::CommensurateBilayerHoneycomb) -> Int

Count the number of honeycomb unitcells contained in the unitcell of a commensurate Moire superlattice composed of two layers of honeycomb lattices.

The total number of atoms in the unitcell of the Moire superlattice is 4 times this result because of the AB sublattice and the top/bottom layer degrees of freedom.
"""
@inline function Base.count(moire::CommensurateBilayerHoneycomb)
    m, r = moire.characters
    if r%3 == 0
        return m^2 + m*r + r^2÷3
    else
        return 3m^2 + 3m*r + r^2
    end
end

#=== point group symmetries and generic methods ===#
"""
    PointGroup

Abstract type for point group symmetries used to classify bond equivalence and parameterize lattice/reciprocal-lattice types.
"""
abstract type PointGroup end

"""
    Base.angle(sym::PointGroup) -> Float64
    Base.angle(::Type{G}) where {G<:PointGroup} -> Float64

Fundamental rotation angle of the point group in radians.
"""
@inline Base.angle(sym::PointGroup) = angle(typeof(sym))

"""
    Base.sign(sym::PointGroup, ref::Bond, bond::Bond) -> Int
    Base.sign(::Type{G}, ref::Bond, bond::Bond) where {G<:PointGroup} -> Int

Compare `bond` to `ref` under the given point group.

Returns `+1` if parallel, `-1` if antiparallel, `0` if not equivalent.
"""
@inline Base.sign(sym::PointGroup, ref::Bond, bond::Bond) = Base.sign(typeof(sym), ref, bond)
function Base.sign(::Type{G}, ref::Bond, bond::Bond) where {G<:PointGroup}
    α = angle(G)
    θ, θ₀ = azimuth(rcoordinate(bond)), azimuth(rcoordinate(ref))
    Δ = (θ - θ₀) / α
    Δ_int = round(Int, Δ)
    isapprox(Δ, Δ_int; atol=atol) || return 0
    return iseven(Δ_int) ? 1 : -1
end

"""
    reciprocals(sym::PointGroup, [T=Float64]) -> SVector{2, SVector{2, T}}
    reciprocals(::Type{G}, [T=Float64]) where {G<:PointGroup} -> SVector{2, SVector{2, T}}

Reciprocal-lattice translation vectors for the given point group.
"""
@inline reciprocals(sym::PointGroup, ::Type{T}=Float64) where {T<:Number} = reciprocals(typeof(sym), T)

#=== C₆ symmetry ===#
"""
    C₆ <: PointGroup

Six-fold rotation symmetry (60° fundamental angle, π/3 rad).
"""
struct C₆ <: PointGroup end
@inline Base.angle(::Type{C₆}) = π/3
function reciprocals(::Type{C₆}, ::Type{T}=Float64) where {T<:Number}
    b₀ = 4one(T)*π/√(3one(T))
    b₁ = SVector(one(T), zero(T)) * b₀
    b₂ = SVector(-one(T)/2, √(one(T)*3)/2) * b₀
    return SVector(b₁, b₂)
end

"""
    const C6 = C₆

ASCII alias for [`C₆`](@ref).
"""
const C6 = C₆

#=== MoireReciprocalLattice ===#
"""
    MoireReciprocalLattice{G<:PointGroup, T<:Number} <: AbstractLattice{2, T, 0}

Abstract type for Moire reciprocal lattices parameterized by point group `G`.
"""
abstract type MoireReciprocalLattice{G<:PointGroup, T<:Number} <: AbstractLattice{2, T, 0} end
@inline getcontent(moire::MoireReciprocalLattice, ::Val{:vectors}) = SVector{0, SVector{2, scalartype(moire)}}()

"""
    PointGroup(moire::MoireReciprocalLattice) -> PointGroup
    PointGroup(::Type{<:MoireReciprocalLattice{G}}) where {G<:PointGroup} -> G

Get the point group of a Moire reciprocal lattice from an instance or type.
"""
@inline PointGroup(moire::MoireReciprocalLattice) = PointGroup(typeof(moire))
@inline PointGroup(::Type{<:MoireReciprocalLattice{G}}) where {G<:PointGroup} = G()

"""
    MoireTriangularReciprocal{T<:Number} <: MoireReciprocalLattice{C₆, T}

C₆-symmetric Moire reciprocal lattice with truncation.
"""
struct MoireTriangularReciprocal{T<:Number} <: MoireReciprocalLattice{C₆, T}
    Γ::SVector{2, T}
    K₊::SVector{2, T}
    K₋::SVector{2, T}
    translations::SVector{2, SVector{2, T}}
    coordinates::Matrix{T}
    function MoireTriangularReciprocal(truncation::Int, ::Type{T}=Float64) where {T<:Number}
        b₁, b₂ = reciprocals(C₆, T)
        Γ = b₁ / 2
        K₊ = (b₁ + 2b₂) / 6
        K₋ = -K₊
        b₀ = norm(b₁)
        coordinates = SVector{2, T}[]
        for i=-2truncation:2truncation, j=-2truncation:2truncation
            coordinate = i*b₁ + j*b₂
            norm(coordinate)<=truncation*b₀+atol && push!(coordinates, coordinate)
        end
        new{T}(Γ, K₊, K₋, SVector(b₁, b₂), reduce(hcat, coordinates))
    end
end
@inline getcontent(::MoireTriangularReciprocal, ::Val{:name}) = :MoireTriangularReciprocal

#=== MoireNeighbors ===#
"""
    MoireNeighbors{G<:PointGroup, N, D<:Number}

Neighbor shells of a Moire superlattice under point group `G`, storing [`Bond`](@ref) objects grouped by neighbor order (kind).
"""
struct MoireNeighbors{G<:PointGroup, N, D<:Number}
    shells::NTuple{N, Vector{Bond{Int, Point{2, D}}}}
end

"""
    MoireNeighbors{G}(shells) where {G<:PointGroup}
    MoireNeighbors{G, N, D}() where {G<:PointGroup, N, D<:Number}
    MoireNeighbors{G, N}() where {G<:PointGroup, N}
    MoireNeighbors{G, N}(lattice::Lattice) where {G<:PointGroup, N}

Construct `MoireNeighbors` from pre-built shells, empty for coordinate type `D`,
empty (default `Float64`), or by collecting symmetry-inequivalent bonds from `lattice` up to order `N`.
"""
@inline MoireNeighbors{G}(shells::NTuple{N, Vector{Bond{Int, Point{2, D}}}}) where {G<:PointGroup, N, D<:Number} = MoireNeighbors{G, N, D}(shells)
@inline function MoireNeighbors{G, N, D}() where {G<:PointGroup, N, D<:Number}
    shells = ntuple(_ -> Bond{Int, Point{2, D}}[], Val(N))
    return MoireNeighbors{G}(shells)
end
@inline MoireNeighbors{G, N}() where {G<:PointGroup, N} = MoireNeighbors{G, N, Float64}()
@inline function MoireNeighbors{G, N}(lattice::Lattice) where {G<:PointGroup, N}
    neighbors = MoireNeighbors{G, N, scalartype(lattice)}()
    for bond in bonds(lattice, N)
        push!(neighbors, bond)
    end
    return neighbors
end

"""
    PointGroup(neighbors::MoireNeighbors) -> PointGroup
    PointGroup(::Type{<:MoireNeighbors{G}}) where {G<:PointGroup} -> G

Get the point group of MoireNeighbors from an instance or type.
"""
@inline PointGroup(neighbors::MoireNeighbors) = PointGroup(typeof(neighbors))
@inline PointGroup(::Type{<:MoireNeighbors{G}}) where {G<:PointGroup} = G()

"""
    truncation(neighbors::MoireNeighbors) -> Int
    truncation(::Type{<:MoireNeighbors{<:PointGroup, N}}) where N -> N

Get the number of neighbor shells from an instance or type.
"""
@inline truncation(neighbors::MoireNeighbors) = truncation(typeof(neighbors))
@inline truncation(::Type{<:MoireNeighbors{<:PointGroup, N}}) where N = N

"""
    getindex(neighbors::MoireNeighbors, k::Int) -> Vector{<:Bond}

Get the `k`-th neighbor shell (1-indexed).
"""
@inline Base.getindex(neighbors::MoireNeighbors, k::Int) = neighbors.shells[k]
@inline Base.firstindex(::MoireNeighbors) = 1
@inline Base.lastindex(neighbors::MoireNeighbors) = truncation(neighbors)

"""
    length(neighbors::MoireNeighbors) -> Int

Number of neighbor shells (same as [`truncation`](@ref)).
"""
@inline Base.length(neighbors::MoireNeighbors) = truncation(neighbors)

"""
    in(bond::Bond, neighbors::MoireNeighbors) -> Bool

Check whether `bond` is symmetry-equivalent to any bond in `neighbors`.
"""
@inline function Base.in(bond::Bond, neighbors::MoireNeighbors{G}) where {G<:PointGroup}
    bond.kind < 1 && return false
    return any(ref -> !iszero(sign(G, bond, ref)), neighbors[bond.kind])
end

"""
    push!(neighbors::MoireNeighbors, bond::Bond) -> MoireNeighbors

Push `bond` into `neighbors` if it is not symmetry-equivalent to any existing bond
in the same shell. Bonds with `kind ≤ 0` are ignored.
"""
@inline function Base.push!(neighbors::MoireNeighbors, bond::Bond)
    bond.kind > 0 || return neighbors
    bond ∉ neighbors && push!(neighbors.shells[bond.kind], bond)
    return neighbors
end

"""
    sign(neighbors::MoireNeighbors, bond::Bond) -> Int

Find the direction of `bond` relative to the matched reference in `neighbors[bond.kind]`. Returns `+1`/`-1`/`0` (not found).
"""
function Base.sign(neighbors::MoireNeighbors{G}, bond::Bond) where {G<:PointGroup}
    bond.kind < 1 && return 0
    for ref in neighbors[bond.kind]
        s = sign(G, ref, bond)
        s == 0 || return s
    end
    return 0
end

#=== MoireSuperlattice ===#
"""
    MoireSuperlattice{G<:PointGroup, N, D<:Number} <: AbstractLattice{2, D, 2}

Abstract type of the emergent superlattices in Moire systems, parameterized by point group `G`, truncation `N`, and coordinate type `D`.
"""
abstract type MoireSuperlattice{G<:PointGroup, N, D<:Number} <: AbstractLattice{2, D, 2} end

"""
    PointGroup(lattice::MoireSuperlattice) -> PointGroup
    PointGroup(::Type{<:MoireSuperlattice{G}}) where {G<:PointGroup} -> G

Get the point group of a Moire superlattice from an instance or type.
"""
@inline PointGroup(lattice::MoireSuperlattice) = PointGroup(typeof(lattice))
@inline PointGroup(::Type{<:MoireSuperlattice{G}}) where {G<:PointGroup} = G()

"""
    truncation(lattice::MoireSuperlattice) -> Int
    truncation(::Type{<:MoireSuperlattice{<:PointGroup, N}}) where N -> N

Get the truncation (number of shells) of a Moire superlattice from an instance or type.
"""
@inline truncation(lattice::MoireSuperlattice) = truncation(typeof(lattice))
@inline truncation(::Type{<:MoireSuperlattice{<:PointGroup, N}}) where N = N

"""
    MoireTriangular{N, D<:Number} <: MoireSuperlattice{C₆, N, D}

Emergent triangular superlattice in Moire systems.
"""
struct MoireTriangular{N, D<:Number} <: MoireSuperlattice{C₆, N, D}
    coordinates::Matrix{D}
    vectors::SVector{2, SVector{2, D}}
    neighbors::MoireNeighbors{C₆, N, D}
end
@inline getcontent(::MoireTriangular, ::Val{:name}) = :MoireTriangular

"""
    MoireTriangular(truncation::Int, [T=Float64])

Construct with truncation `N = truncation` and coordinate type `T`.

The single site per unitcell is at the origin (MM stacking site).
"""
function MoireTriangular(truncation::Int, ::Type{T}=Float64) where {T<:Number}
    vectors = reciprocals(reciprocals(C₆, T))
    coordinates = zeros(T, 2, 1)
    neighbors = MoireNeighbors{C₆, truncation}(Lattice(:MoireTriangular, coordinates, vectors))
    return MoireTriangular(coordinates, vectors, neighbors)
end

"""
    MoireHoneycomb{N, D<:Number} <: MoireSuperlattice{C₆, N, D}

Emergent honeycomb superlattice in Moire systems, composed of MX and XM stacking sites.
"""
struct MoireHoneycomb{N, D<:Number} <: MoireSuperlattice{C₆, N, D}
    coordinates::Matrix{D}
    vectors::SVector{2, SVector{2, D}}
    neighbors::MoireNeighbors{C₆, N, D}
end
@inline getcontent(::MoireHoneycomb, ::Val{:name}) = :MoireHoneycomb

"""
    MoireHoneycomb(truncation::Int, [T=Float64])

Construct with truncation `N = truncation` and coordinate type `T`.

MX at `(v₁+v₂)/3`, XM at `(2v₁-v₂)/3` (MM stacking at origin).
"""
function MoireHoneycomb(truncation::Int, ::Type{T}=Float64) where {T<:Number}
    v₁, v₂ = reciprocals(reciprocals(C₆, T))
    coordinates = zeros(T, 2, 2)
    coordinates[:, 1] = (v₁ + v₂) / 3
    coordinates[:, 2] = (2v₁ - v₂) / 3
    neighbors = MoireNeighbors{C₆, truncation}(Lattice(:MoireHoneycomb, coordinates, SVector(v₁, v₂)))
    return MoireHoneycomb(coordinates, SVector(v₁, v₂), neighbors)
end

"""
    RealZone{N, S<:SVector, V<:Number}

A rectangular zone in real space.

Alias for [`ReciprocalZone`](@ref) with the space-type parameter `K = :r`.
"""
const RealZone{N, S<:SVector, V<:Number} = ReciprocalZone{:r, N, S, V}
