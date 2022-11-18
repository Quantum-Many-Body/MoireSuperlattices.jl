module MoireSuperlattices

using QuantumLattices: hexagon120°map, hexagon60°map, Bond, Lattice, Neighbors, decimaltostr, distance, reciprocals, rotate
using RecipesBase: RecipesBase, @recipe, @series
using StaticArrays: SVector

export CommensurateBilayerGraphene

struct CommensurateBilayerGraphene{T, V<:AbstractVector{T}}
    characters::Tuple{Int, Int}
    displacement::V
    center::V
    coordinates::Matrix{T}
    vectors::Vector{V}
    function CommensurateBilayerGraphene(characters::Tuple{Int, Int}, displacement::V, center::V, coordinates::Matrix{T}, vectors::Vector{V}) where {T, V<:AbstractVector{T}}
        @assert gcd(characters[1], characters[2])==1 "CommensurateBilayerGraphene error: not coprime integers."
        @assert length(displacement)==2 "CommensurateBilayerGraphene error: incorrect displacement vector between the two layers."
        @assert length(center)==2 "CommensurateBilayerGraphene error: incorrect rotation center of the two layers."
        @assert size(coordinates)==(2, 2) "CommensurateBilayerGraphene error: incorrect coordinates of the points in the unitcell of a monolayer graphene."
        @assert length(vectors)==2 && length(vectors[1])==2 && length(vectors[2])==2 "CommensurateBilayerGraphene error: incorrect translation vectors of a monolayer graphene."
        new{T, V}(characters, displacement, center, coordinates, vectors)
    end
end
function CommensurateBilayerGraphene(characters::Tuple{Int, Int}; stack::Symbol=:AA, center::Symbol=:carbon)
    @assert stack∈(:AA, :AB) "CommensurateBilayerGraphene error: not supported stack."
    @assert center∈(:carbon, :hexagon) "CommensurateBilayerGraphene error: not supported center."
    a₁, a₂ = SVector(1.0, 0.0), SVector(-0.5, √3/2)
    coordinates = [1/3*a₁+2/3*a₂;; 2/3*a₁+1/3*a₂]
    displacement = stack==:AA ? SVector(0.0, 0.0) : -1/3*a₁+1/3*a₂
    center = center==:carbon ? SVector(coordinates[1, 1], coordinates[2, 1]) : SVector(0.0, 0.0)
    return CommensurateBilayerGraphene(characters, displacement, center, coordinates, [a₁, a₂])
end

@inline function Base.angle(moire::CommensurateBilayerGraphene)
    m, r = moire.characters
    return acos((3m^2+3m*r+r^2/2)/(3m^2+3m*r+r^2))
end
@inline function translations(moire::CommensurateBilayerGraphene)
    m, r = moire.characters
    θ = angle(moire)
    a₁, a₂ = rotate(moire.vectors[1], -θ/2), rotate(moire.vectors[2], -θ/2)
    if r%3==0
        t₁ = (m+r÷3)*a₁ + (m+2r÷3)*a₂
        t₂ = -r÷3*a₁ + (m+r÷3)*a₂
    else
        t₁ = m*a₁ + (2m+r)*a₂
        t₂ = -(m+r)*a₁ + m*a₂
    end
    return [t₁, t₂]
end

@recipe function plot(moire::CommensurateBilayerGraphene, choice::Symbol)
    @assert choice∈(:lattice, :reciprocal) "plot error: incorrect choice."
    θ = angle(moire)
    n = 2*ceil(Int, √((3*moire.characters[1]^2+3*moire.characters[1]*moire.characters[2]+moire.characters[2]^2)/gcd(moire.characters[2], 3)))
    t₁, t₂ = translations(moire)
    vectors₁ = map(v->rotate(v, +θ/2), moire.vectors)
    vectors₂ = map(v->rotate(v, -θ/2), moire.vectors)
    title --> "Twisted Bilayer Graphene ($(decimaltostr(rad2deg(θ)))°)"
    aspect_ratio := :equal
    legend := false
    if choice==:lattice
        coordinates₁ = rotate(moire.coordinates, +θ/2; axis=(moire.center, (0, 0)))
        coordinates₂ = rotate(moire.coordinates.+reshape(moire.displacement, :, 1), -θ/2; axis=(moire.center, (0, 0)))
        neighbors = Neighbors(1=>distance(moire.coordinates[:, 1], moire.coordinates[:, 2]))
        @series begin
            Lattice(Lattice{2}(:top, coordinates₁, vectors₁), (2n, 2n); mode=:center), neighbors
        end
        @series begin
            Lattice(Lattice{2}(:bottom, coordinates₂, vectors₂), (2n, 2n); mode=:center), neighbors
        end
        arrow := true
        linewidth := 3
        [moire.center[1], t₁[1]+moire.center[1], NaN, moire.center[1], t₂[1]+moire.center[1]], [moire.center[2], t₁[2]+moire.center[2], NaN, moire.center[2], t₂[2]+moire.center[2]]
    else
        recipls₁ = reciprocals(vectors₁)
        recipls₂ = reciprocals(vectors₂)
        @series begin
            Lattice([mapreduce(*, +, hexagon60°map[key], recipls₁) for key in ("K₁", "K₂", "K₃", "K₄", "K₅", "K₆")]...), 1, bond::Bond->bond.kind==1
        end
        @series begin
            Lattice([mapreduce(*, +, hexagon60°map[key], recipls₂) for key in ("K₁", "K₂", "K₃", "K₄", "K₅", "K₆")]...), 1, bond::Bond->bond.kind==1
        end
        recipls = reciprocals([t₁, t₂])
        Lattice(Lattice([mapreduce(*, +, hexagon120°map[key], recipls) for key in ("K₁", "K₂")]...; vectors=recipls), (n, n); mode=:center), 1, bond::Bond->bond.kind==1
    end
end

end