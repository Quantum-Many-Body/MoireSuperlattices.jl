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

#=== sublattice ===#
"""
    Sublattice(nsublattice::Int)

Callable that maps a lattice site index to its 1-based sublattice: `(site - 1) % nsublattice + 1`.
"""
struct Sublattice
    nsublattice::Int
end
@inline (s::Sublattice)(site::Int) = (site - 1) % s.nsublattice + 1

#=== point group symmetries and generic methods ===#
"""
    PointGroup

Abstract type for point group symmetries used to classify bond equivalence and parameterize lattice/reciprocal-lattice types.
"""
abstract type PointGroup end

"""
    angle(sym::PointGroup) -> Float64
    angle(::Type{G}) where {G<:PointGroup} -> Float64

Fundamental rotation angle of the point group in radians.
"""
@inline Base.angle(sym::PointGroup) = angle(typeof(sym))

"""
    sign(::Type{G}, ref::Bond, bond::Bond, nsublattice::Int; atol::Real=atol) where {G<:PointGroup} -> Int
    sign(sym::PointGroup, ref::Bond, bond::Bond, nsublattice::Int; atol::Real=atol) -> Int

Bond equivalence check under the point group `G`.

Each concrete `PointGroup` subtype must implement `Base.sign(::Type{G}, ref, bond, nsublattice)` returning:
- `+1`: `bond` is symmetry-equivalent to `ref` in the forward direction.
- `-1`: `bond` is symmetry-equivalent to the reverse of `ref`.
- `0`: not equivalent.

The instance method `sign(sym::PointGroup, ...)` dispatches to the type-specific implementation via `sign(typeof(sym), ...)`.
"""
@inline Base.sign(sym::PointGroup, ref::Bond, bond::Bond, nsublattice::Int) = sign(typeof(sym), ref, bond, nsublattice)

"""
    reciprocals(sym::PointGroup, [T=Float64]) -> SVector{2, SVector{2, T}}
    reciprocals(::Type{G}, [T=Float64]) where {G<:PointGroup} -> SVector{2, SVector{2, T}}

Reciprocal-lattice translation vectors for the given point group.
"""
@inline reciprocals(sym::PointGroup, ::Type{T}=Float64) where {T<:Number} = reciprocals(typeof(sym), T)

#=== C₃ symmetry ===#
"""
    C₃ <: PointGroup

Three-fold rotation symmetry (120° fundamental angle, 2π/3 rad).

For bond equivalence classification, [`sign`](@ref) uses a 60° azimuth discretization to distinguish bonds that are C₃-equivalent in the forward direction (``Δ`` even → `+1`) from those that are C₃-equivalent in the reverse direction (``Δ`` odd → `-1`).
"""
struct C₃ <: PointGroup end
@inline Base.angle(::Type{C₃}) = 2π/3
function Base.sign(::Type{C₃}, ref::Bond, bond::Bond, nsublattice::Int)
    ref.kind == bond.kind || return 0
    sl = Sublattice(nsublattice)
    rᵢ, rⱼ = sl(ref[1].site), sl(ref[2].site)
    bᵢ, bⱼ = sl(bond[1].site), sl(bond[2].site)
    if rᵢ == bᵢ && rⱼ == bⱼ
        reversed = false
    elseif rᵢ == bⱼ && rⱼ == bᵢ
        reversed = true
    else
        return 0
    end
    α, θ, θ₀ = angle(C₃)/2, azimuth(rcoordinate(bond)), azimuth(rcoordinate(ref))
    Δ = (θ - θ₀) / α
    Δ_int = round(Int, Δ)
    isapprox(Δ, Δ_int; atol=atol) || return 0
    # reversed sublattice + same spatial direction → distinct bonds, not equivalent
    reversed && iseven(Δ_int) && return 0
    return iseven(Δ_int) ? 1 : -1
end
function reciprocals(::Type{C₃}, ::Type{T}=Float64) where {T<:Number}
    b₀ = 4one(T)*π/√(3one(T))
    b₁ = SVector(one(T), zero(T)) * b₀
    b₂ = SVector(-one(T)/2, √(one(T)*3)/2) * b₀
    return SVector(b₁, b₂)
end

"""
    const C3 = C₃

ASCII alias for [`C₃`](@ref).
"""
const C3 = C₃

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
    truncation(lattice::MoireReciprocalLattice) -> Int

Get the truncation (number of shells) of a Moire reciprocal lattice.
"""
@inline truncation(lattice::MoireReciprocalLattice) = lattice.truncation

"""
    reciprocals(lattice::MoireReciprocalLattice) -> SVector{2, SVector{2, scalartype(lattice)}}

Get the reciprocal translation vectors of a Moire reciprocal lattice.
"""
@inline reciprocals(lattice::MoireReciprocalLattice) = lattice.translations

"""
    MoireTriangularReciprocal{T<:Number} <: MoireReciprocalLattice{C₃, T}

C₃-symmetric Moire reciprocal lattice with truncation.
"""
struct MoireTriangularReciprocal{T<:Number} <: MoireReciprocalLattice{C₃, T}
    Γ::SVector{2, T}
    K₊::SVector{2, T}
    K₋::SVector{2, T}
    translations::SVector{2, SVector{2, T}}
    coordinates::Matrix{T}
    truncation::Int
    function MoireTriangularReciprocal(truncation::Int, ::Type{T}=Float64) where {T<:Number}
        b₁, b₂ = reciprocals(C₃, T)
        Γ = b₁ / 2
        K₊ = (b₁ + 2b₂) / 6
        K₋ = -K₊
        b₀ = norm(b₁)
        coordinates = SVector{2, T}[]
        for i=-2truncation:2truncation, j=-2truncation:2truncation
            coordinate = i*b₁ + j*b₂
            norm(coordinate)<=truncation*b₀+atol && push!(coordinates, coordinate)
        end
        new{T}(Γ, K₊, K₋, SVector(b₁, b₂), reduce(hcat, coordinates), truncation)
    end
end
@inline getcontent(::MoireTriangularReciprocal, ::Val{:name}) = :MoireTriangularReciprocal

#=== MoireNeighbors ===#
"""
    MoireNeighbors{G<:PointGroup, D<:Number}

Flat collection of symmetry-inequivalent `Bond` objects under point group `G`.

Use [`bonds`](@ref) to query by neighbor order.
"""
struct MoireNeighbors{G<:PointGroup, D<:Number}
    bonds::Vector{Bond{Int, Point{2, D}, SVector{2, Point{2, D}}}}
    nsublattice::Int
    truncation::Int
end

"""
    MoireNeighbors{G}(nsublattice, truncation) where {G<:PointGroup}
    MoireNeighbors{G, D}(nsublattice, truncation) where {G<:PointGroup, D<:Number}
    MoireNeighbors{G}(lattice, truncation) where {G<:PointGroup}

Constructors:

1. `MoireNeighbors{G}(nsublattice, truncation)` — empty, default `Float64` coordinates.
2. `MoireNeighbors{G, D}(nsublattice, truncation)` — empty, coordinate type `D`.
3. `MoireNeighbors{G}(lattice, truncation)` — collect symmetry-inequivalent bonds; `nsublattice` from `length(lattice)`.
"""
@inline MoireNeighbors{G}(nsublattice::Int, truncation::Int) where {G<:PointGroup} = MoireNeighbors{G, Float64}(nsublattice, truncation)
@inline MoireNeighbors{G, D}(nsublattice::Int, truncation::Int) where {G<:PointGroup, D<:Number} = MoireNeighbors{G, D}(Bond{Int, Point{2, D}, SVector{2, Point{2, D}}}[], nsublattice, truncation)
@inline function MoireNeighbors{G}(lattice::AbstractLattice, truncation::Int) where {G<:PointGroup}
    nsublattice = length(lattice)
    neighbors = MoireNeighbors{G, scalartype(lattice)}(nsublattice, truncation)
    for bond in bonds(lattice, truncation)
        push!(neighbors, bond)
    end
    return neighbors
end

"""
    nsublattice(neighbors::MoireNeighbors) -> Int

Number of sublattice sites per unit cell of the effective lattice.
"""
@inline nsublattice(neighbors::MoireNeighbors) = neighbors.nsublattice

"""
    truncation(neighbors::MoireNeighbors) -> Int

Maximum neighbor order (kind) stored in `neighbors`.
"""
@inline truncation(neighbors::MoireNeighbors) = neighbors.truncation

"""
    PointGroup(neighbors::MoireNeighbors) -> PointGroup
    PointGroup(::Type{<:MoireNeighbors{G}}) where {G<:PointGroup} -> G

Get the point group of MoireNeighbors from an instance or type.
"""
@inline PointGroup(neighbors::MoireNeighbors) = PointGroup(typeof(neighbors))
@inline PointGroup(::Type{<:MoireNeighbors{G}}) where {G<:PointGroup} = G()

"""
    bonds(neighbors::MoireNeighbors) -> Vector{<:Bond}
    bonds(neighbors::MoireNeighbors, k::Int) -> Vector{<:Bond}

Bond access: all stored inequivalent bonds, or filtered by neighbor order `k`.
"""
@inline bonds(neighbors::MoireNeighbors) = neighbors.bonds
@inline bonds(neighbors::MoireNeighbors, k::Int) = [bond for bond in neighbors.bonds if bond.kind == k]

"""
    in(bond::Bond, neighbors::MoireNeighbors) -> Bool

Check whether `bond` is symmetry-equivalent to any bond in `neighbors`.
"""
@inline function Base.in(bond::Bond, neighbors::MoireNeighbors{G}) where {G<:PointGroup}
    return length(bond)==2 && any(ref -> !iszero(sign(G, bond, ref, neighbors.nsublattice)), neighbors.bonds)
end

"""
    push!(neighbors::MoireNeighbors, bond::Bond) -> MoireNeighbors

Push `bond` if not symmetry-equivalent to an existing bond.
"""
@inline function Base.push!(neighbors::MoireNeighbors, bond::Bond)
    if length(bond) == 2 && bond ∉ neighbors
        push!(neighbors.bonds, Bond(bond.kind, SVector(bond[1], bond[2])))
    end
    return neighbors
end

#=== MoireSuperlattice ===#
"""
    MoireSuperlattice{G<:PointGroup, D<:Number} <: AbstractLattice{2, D, 2}

Abstract type of the emergent superlattices in Moire systems, parameterized by point group `G` and coordinate type `D`.
"""
abstract type MoireSuperlattice{G<:PointGroup, D<:Number} <: AbstractLattice{2, D, 2} end

"""
    PointGroup(lattice::MoireSuperlattice) -> PointGroup
    PointGroup(::Type{<:MoireSuperlattice{G}}) where {G<:PointGroup} -> G

Get the point group of a Moire superlattice from an instance or type.
"""
@inline PointGroup(lattice::MoireSuperlattice) = PointGroup(typeof(lattice))
@inline PointGroup(::Type{<:MoireSuperlattice{G}}) where {G<:PointGroup} = G()

"""
    truncation(lattice::MoireSuperlattice) -> Int

Get the truncation (number of shells) of a Moire superlattice.
"""
@inline truncation(lattice::MoireSuperlattice) = truncation(lattice.neighbors)

"""
    MoireTriangular{D<:Number} <: MoireSuperlattice{C₃, D}

Emergent triangular superlattice in Moire systems.
"""
struct MoireTriangular{D<:Number} <: MoireSuperlattice{C₃, D}
    coordinates::Matrix{D}
    vectors::SVector{2, SVector{2, D}}
    neighbors::MoireNeighbors{C₃, D}
end
@inline getcontent(::MoireTriangular, ::Val{:name}) = :MoireTriangular

"""
    MoireTriangular(truncation::Int, [T=Float64])

Construct with coordinate type `T`. The single site per unitcell is at the origin (MM stacking site).
"""
function MoireTriangular(truncation::Int, ::Type{T}=Float64) where {T<:Number}
    vectors = reciprocals(reciprocals(C₃, T))
    coordinates = zeros(T, 2, 1)
    neighbors = MoireNeighbors{C₃}(Lattice(:MoireTriangular, coordinates, vectors), truncation)
    return MoireTriangular(coordinates, vectors, neighbors)
end

"""
    MoireHoneycomb{D<:Number} <: MoireSuperlattice{C₃, D}

Emergent honeycomb superlattice in Moire systems, composed of MX and XM stacking sites.
"""
struct MoireHoneycomb{D<:Number} <: MoireSuperlattice{C₃, D}
    coordinates::Matrix{D}
    vectors::SVector{2, SVector{2, D}}
    neighbors::MoireNeighbors{C₃, D}
end
@inline getcontent(::MoireHoneycomb, ::Val{:name}) = :MoireHoneycomb

"""
    MoireHoneycomb(truncation::Int, [T=Float64])

Construct with coordinate type `T`. XM at `(2v₁-v₂)/3`, MX at `(v₁+v₂)/3`.
"""
function MoireHoneycomb(truncation::Int, ::Type{T}=Float64) where {T<:Number}
    v₁, v₂ = reciprocals(reciprocals(C₃, T))
    coordinates = zeros(T, 2, 2)
    coordinates[:, 1] = (2v₁ - v₂) / 3
    coordinates[:, 2] = (v₁ + v₂) / 3
    neighbors = MoireNeighbors{C₃}(Lattice(:MoireHoneycomb, coordinates, SVector(v₁, v₂)), truncation)
    return MoireHoneycomb(coordinates, SVector(v₁, v₂), neighbors)
end

"""
    RealZone{N, S<:SVector, V<:Number}

A rectangular zone in real space.

Alias for `ReciprocalZone` with the space-type parameter `K = :r`.
"""
const RealZone{N, S<:SVector, V<:Number} = ReciprocalZone{:r, N, S, V}
