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

#=== MoireReciprocalLattice ===#
"""
    MoireReciprocalLattice{T<:Number} <: AbstractLattice{2, T, 0}

Abstract type for Moire reciprocal lattices.
"""
abstract type MoireReciprocalLattice{T<:Number} <: AbstractLattice{2, T, 0} end
@inline getcontent(moire::MoireReciprocalLattice, ::Val{:vectors}) = SVector{0, SVector{2, scalartype(moire)}}()

"""
    truncation(lattice::MoireReciprocalLattice) -> Int

Get the truncation (number of shells) of a Moire reciprocal lattice.
"""
@inline truncation(lattice::MoireReciprocalLattice) = lattice.truncation

"""
    reciprocals(lattice::MoireReciprocalLattice) -> SVector{2, SVector{2, scalartype(lattice)}}

Get the reciprocal translation vectors of a Moire reciprocal lattice.
"""
@inline reciprocals(lattice::MoireReciprocalLattice) = reciprocals(typeof(lattice))

"""
    MoireTriangularReciprocal{T<:Number} <: MoireReciprocalLattice{T}

Moire reciprocal lattice with truncation.
"""
struct MoireTriangularReciprocal{T<:Number} <: MoireReciprocalLattice{T}
    Γ::SVector{2, T}
    K₊::SVector{2, T}
    K₋::SVector{2, T}
    translations::SVector{2, SVector{2, T}}
    coordinates::Matrix{T}
    truncation::Int
    function MoireTriangularReciprocal(truncation::Int, ::Type{T}=Float64) where {T<:Number}
        b₁, b₂ = reciprocals(MoireTriangularReciprocal, T)
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

"""
    reciprocals(::Type{<:MoireTriangularReciprocal}, [T=Float64]) -> SVector{2, SVector{2, T}}

Reciprocal-lattice translation vectors for the triangular Bravais lattice underlying Moire superlattices.
"""
function reciprocals(::Type{<:MoireTriangularReciprocal}, ::Type{T}=Float64) where {T<:Number}
    b₀ = 4one(T)*π/√(3one(T))
    b₁ = SVector(one(T), zero(T)) * b₀
    b₂ = SVector(-one(T)/2, √(one(T)*3)/2) * b₀
    return SVector(b₁, b₂)
end

#=== MoireSuperlattice ===#
"""
    MoireSuperlattice{D<:Number} <: AbstractLattice{2, D, 2}

Abstract type of the emergent superlattices in Moire systems, parameterized by coordinate type `D`.
"""
abstract type MoireSuperlattice{D<:Number} <: AbstractLattice{2, D, 2} end

"""
    MoireTriangular{D<:Number} <: MoireSuperlattice{D}

Emergent triangular superlattice in Moire systems.
"""
struct MoireTriangular{D<:Number} <: MoireSuperlattice{D}
    coordinates::Matrix{D}
    vectors::SVector{2, SVector{2, D}}
    function MoireTriangular(::Type{T}=Float64) where {T<:Number}
        vectors = reciprocals(reciprocals(MoireTriangularReciprocal, T))
        coordinates = zeros(T, 2, 1)
        return new{T}(coordinates, vectors)
    end
end
@inline getcontent(::MoireTriangular, ::Val{:name}) = :MoireTriangular

"""
    reciprocals(::Type{<:MoireTriangular}, [T=Float64]) -> SVector{2, SVector{2, T}}

Delegates to `MoireTriangularReciprocal` (shared triangular Bravais lattice).
"""
@inline reciprocals(::Type{<:MoireTriangular}, ::Type{T}=Float64) where {T<:Number} = reciprocals(MoireTriangularReciprocal, T)

"""
    MoireHoneycomb{D<:Number} <: MoireSuperlattice{D}

Emergent honeycomb superlattice in Moire systems, composed of XM and MX stacking sites.
"""
struct MoireHoneycomb{D<:Number} <: MoireSuperlattice{D}
    coordinates::Matrix{D}
    vectors::SVector{2, SVector{2, D}}
    function MoireHoneycomb(::Type{T}=Float64) where {T<:Number}
        v₁, v₂ = reciprocals(reciprocals(MoireTriangularReciprocal, T))
        coordinates = zeros(T, 2, 2)
        coordinates[:, 1] = (2v₁ - v₂) / 3  # XM
        coordinates[:, 2] = (v₁ + v₂) / 3   # MX
        return new{T}(coordinates, SVector(v₁, v₂))
    end
end
@inline getcontent(::MoireHoneycomb, ::Val{:name}) = :MoireHoneycomb

"""
    reciprocals(::Type{<:MoireHoneycomb}, [T=Float64]) -> SVector{2, SVector{2, T}}

Delegates to `MoireTriangularReciprocal` (shared triangular Bravais lattice).
"""
@inline reciprocals(::Type{<:MoireHoneycomb}, ::Type{T}=Float64) where {T<:Number} = reciprocals(MoireTriangularReciprocal, T)

"""
    RealZone{N, S<:SVector, V<:Number}

A rectangular zone in real space.

Alias for `ReciprocalZone` with the space-type parameter `K = :r`.
"""
const RealZone{N, S<:SVector, V<:Number} = ReciprocalZone{:r, N, S, V}
