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
    angle(moire::CommensurateBilayerHoneycomb) -> Float64

Get the twist angle of a commensurate Moire superlattice composed of two layers of honeycomb lattices.
"""
@inline function Base.angle(moire::CommensurateBilayerHoneycomb)
    m, r = moire.characters
    return acos((3m^2+3m*r+r^2/2)/(3m^2+3m*r+r^2))
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

"""
    MoireReciprocalLattice{T<:Number} <: AbstractLattice{2, T, 0}

Moire reciprocal lattice with truncation.
"""
struct MoireReciprocalLattice{T<:Number} <: AbstractLattice{2, T, 0}
    Γ::SVector{2, T}
    K₊::SVector{2, T}
    K₋::SVector{2, T}
    translations::SVector{2, SVector{2, T}}
    coordinates::Matrix{T}
    function MoireReciprocalLattice(truncation::Int, ::Type{T}=Float64) where {T<:Number}
        b₀ = 4one(T)*pi/√(3one(T))
        b₁, b₂ = SVector(one(T), zero(T))*b₀, SVector(-one(T)/2, √(one(T)*3)/2)*b₀
        Γ, K₊, K₋ = SVector(one(T)/2, zero(T))*b₀, SVector(zero(T), +√(one(T)*3)/6)*b₀, SVector(zero(T), -√(one(T)*3)/6)*b₀
        coordinates = SVector{2, T}[]
        for i=-2truncation:2truncation, j=-2truncation:2truncation
            coordinate = i*b₁ + j*b₂
            norm(coordinate)<=truncation*b₀+atol && push!(coordinates, coordinate)
        end
        new{T}(Γ, K₊, K₋, SVector(b₁, b₂), reduce(hcat, coordinates))
    end
end
@inline getcontent(moire::MoireReciprocalLattice, ::Val{:name}) = :truncation
@inline getcontent(moire::MoireReciprocalLattice, ::Val{:vectors}) = SVector{0, SVector{2, scalartype(moire)}}()

"""
    MoireSuperlattice{D<:Number} <: AbstractLattice{2, D, 2}

Abstract type of the emergent superlattices in Moire systems.
"""
abstract type MoireSuperlattice{D<:Number} <: AbstractLattice{2, D, 2} end

"""
    MoireTriangular{N, D<:Number} <: MoireSuperlattice{D}

Emergent triangular superlattice in Moire systems.
"""
struct MoireTriangular{N, D<:Number} <: MoireSuperlattice{D}
    name::Symbol
    coordinates::Matrix{D}
    vectors::SVector{2, SVector{2, D}}
    neighbors::NTuple{N, Vector{SVector{2, D}}}
end
@inline truncation(lattice::MoireTriangular) = truncation(typeof(lattice))
@inline truncation(::Type{<:MoireTriangular{N}}) where N = N

"""
    MoireTriangular(truncation::Int, vectors::AbstractVector{<:AbstractVector{<:Number}}; name=:MoireTriangular, origin=nothing)

Construct the emergent triangular superlattice in Moire systems.
"""
function MoireTriangular(truncation::Int, vectors::AbstractVector{<:AbstractVector{<:Number}}; name=:MoireTriangular, origin=nothing)
    @assert length(vectors)==2 "MoireTriangular error: 2 instead of $(length(vectors)) vectors should be input."
    @assert all(v->length(v)==2, vectors) "MoireTriangular error: the length of every vector should be 2."
    datatype = eltype(eltype(vectors))
    vectors = convert(SVector{2, SVector{2, datatype}}, vectors)
    coordinates = zeros(datatype, 2, 1)
    isnothing(origin) || begin
        @assert length(origin)==2 "MoireTriangular error: the length of the origin point should be 2."
        coordinates[1, 1] = origin[1]
        coordinates[2, 1] = origin[2]
    end
    neighbors = ntuple(i->SVector{2, datatype}[], Val(truncation))
    for bond in bonds(Lattice(name, coordinates, vectors), truncation)
        if bond.kind>0
            coordinate = rcoordinate(bond)
            if length(neighbors[bond.kind])==0
                push!(neighbors[bond.kind], coordinate)
            else
                for i = 1:length(neighbors[bond.kind])
                    another = neighbors[bond.kind][i]
                    θ = acos(dot(coordinate, another)/norm(coordinate)/norm(another))/(pi/3)
                    isapprox(round(Int, convert(Float64, θ)), convert(Float64, θ); atol=atol) && break
                    if i==length(neighbors[bond.kind])
                        push!(neighbors[bond.kind], coordinate)
                    end
                end
            end
        end
    end
    return MoireTriangular(name, coordinates, vectors, neighbors)
end

"""
    RealZone{N, S<:SVector, V<:Number}

A rectangular zone in real space.

Alias for [`ReciprocalZone`](@ref) with the space-type parameter `K = :r`.
"""
const RealZone{N, S<:SVector, V<:Number} = ReciprocalZone{:r, N, S, V}
