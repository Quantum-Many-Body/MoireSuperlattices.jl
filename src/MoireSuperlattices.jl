module MoireSuperlattices

using LinearAlgebra: dot, eigvals, norm
using Printf: @printf
using QuantumLattices: atol, hexagon120°map, hexagon60°map, plain
using QuantumLattices: AbstractLattice, Algorithm, Bond, BrillouinZone, CompositeIID, CompositeIndex, Coupling, FID, Hilbert, Hopping, ID, IIDSpace, Index
using QuantumLattices: lazy, Image, Lattice, LaTeX, MatrixCoupling, Neighbors, Onsite, OperatorGenerator, OperatorUnitToTuple, Point, SimpleIID, SimpleInternal, Table, Term
using QuantumLattices: azimuth, azimuthd, bonds, concatenate, decimaltostr, dimension, distance, dtype, latexformat, reciprocals, rank, rcoordinate, rotate, update, @σ_str
using RecipesBase: RecipesBase, @recipe, @series
using StaticArrays: SVector
using TightBindingApproximation: AbstractTBA, Fermionic, QuadraticFormalize

import QuantumLattices: diagonalizablefields, getcontent, iidtype, isdefinite, latexname, matrix, script, shape, statistics, update!

export BLTMD, CommensurateBilayerHoneycomb, MoireEmergentSuperLattice, MoireReciprocalLattice, MoireSpace, MoireSpinor, MoireSystem, MoireTriangular, bltmd!, bltmdmap, coefficients, terms

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
@recipe function plot(moire::CommensurateBilayerHoneycomb, choice::Symbol; n=nothing, topcolor=:red, bottomcolor=:blue, vector=true, vectorcolor=:green, moirecolor=:black, anglecolor=:grey)
    @assert choice∈(:lattice, :reciprocal) "plot error: incorrect choice."
    m, r = moire.characters
    θ = acos((3m^2+3m*r+r^2/2)/(3m^2+3m*r+r^2))
    isnothing(n) && (n = 2*ceil(Int, √((3*m^2+3*m*r+r^2)/gcd(r, 3))))
    vectors₁, vectors₂ = map(v->SVector{2}(rotate(v, +θ/2)), moire.vectors), map(v->SVector{2}(rotate(v, -θ/2)), moire.vectors)
    (t₁, t₂) = r%3==0 ? ((m+r÷3)*vectors₂[1]+(m+2r÷3)*vectors₂[2], -r÷3*vectors₂[1]+(m+r÷3)*vectors₂[2]) : (m*vectors₂[1]+(2m+r)*vectors₂[2], -(m+r)*vectors₂[1]+m*vectors₂[2])
    title --> "Twisted Bilayer Honeycomb ($(decimaltostr(rad2deg(θ)))°)"
    aspect_ratio := :equal
    legend := false
    if choice==:lattice
        coordinates₁ = rotate(moire.coordinates, +θ/2; axis=(moire.center, (0, 0)))
        coordinates₂ = rotate(moire.coordinates.+reshape(moire.displacement, :, 1), -θ/2; axis=(moire.center, (0, 0)))
        neighbors = Neighbors(1=>distance(moire.coordinates[:, 1], moire.coordinates[:, 2]))
        @series begin
            color --> topcolor
            Lattice(Lattice(:top, coordinates₁, vectors₁), (2n, 2n); mode=:center), neighbors
        end
        @series begin
            color --> bottomcolor
            Lattice(Lattice(:bottom, coordinates₂, vectors₂), (2n, 2n); mode=:center), neighbors
        end
        arrow := true
        linewidth := 2
        color --> vectorcolor
        alpha := (vector ? 1.0 : 0.0)
        [moire.center[1], t₁[1]+moire.center[1], NaN, moire.center[1], t₂[1]+moire.center[1]], [moire.center[2], t₁[2]+moire.center[2], NaN, moire.center[2], t₂[2]+moire.center[2]]
    else
        recipls₁ = reciprocals(vectors₁)
        recipls₂ = reciprocals(vectors₂)
        @series begin
            color --> topcolor
            Lattice([collect(mapreduce(*, +, hexagon60°map[key], recipls₁)) for key in ("K₁", "K₂", "K₃", "K₄", "K₅", "K₆")]...), 1, bond::Bond->bond.kind==1
        end
        @series begin
            color --> bottomcolor
            Lattice([collect(mapreduce(*, +, hexagon60°map[key], recipls₂)) for key in ("K₁", "K₂", "K₃", "K₄", "K₅", "K₆")]...), 1, bond::Bond->bond.kind==1
        end
        @series begin
            color --> anglecolor
            Kt = collect(mapreduce(*, +, hexagon60°map["K₂"], recipls₁))
            Kb = collect(mapreduce(*, +, hexagon60°map["K₂"], recipls₂))
            linestyle := :dot
            [0.0, Kt[1], NaN, 0.0, Kb[1]], [0.0, Kt[2], NaN, 0.0, Kb[2]]
        end
        recipls = reciprocals([t₁, t₂])
        color --> moirecolor
        Lattice(Lattice([collect(mapreduce(*, +, hexagon120°map[key], recipls)) for key in ("K₁", "K₂")]...; vectors=recipls), (n, n); mode=:center), 1, bond::Bond->bond.kind==1
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
@inline getcontent(moire::MoireReciprocalLattice, ::Val{:vectors}) = SVector{0, SVector{2, dtype(moire)}}()

"""
    MoireEmergentSuperLattice{D<:Number} <: AbstractLattice{2, D, 2}

Abstract type of the emergent superlattices in Moire systems.
"""
abstract type MoireEmergentSuperLattice{D<:Number} <: AbstractLattice{2, D, 2} end

"""
    MoireTriangular{N, D<:Number} <: MoireEmergentSuperLattice{D}

Emergent triangular superlattice in Moire systems.
"""
struct MoireTriangular{N, D<:Number} <: MoireEmergentSuperLattice{D}
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
    MoireSpinor{V<:Union{Int, Colon}, L<:Union{Int, Colon}, S<:Union{Int, Colon}, P<:Union{Rational{Int}, Colon}, N<:Union{Int, Colon}} <: SimpleIID

The index of the internal degrees of freedom of Moire systems.
"""
struct MoireSpinor{V<:Union{Int, Colon}, L<:Union{Int, Colon}, S<:Union{Int, Colon}, P<:Union{Rational{Int}, Colon}, N<:Union{Int, Colon}} <: SimpleIID
    valley::V
    layer::L
    sublattice::S
    spin::P
    nambu::N
    function MoireSpinor(valley::Union{Int, Colon}, layer::Union{Int, Colon}, sublattice::Union{Int, Colon}, spin::Union{Rational{Int}, Colon}, nambu::Union{Int, Colon})
        @assert spin∈(-1//2, 1//2, 0, :) "MoireSpinor error: incorrect spin ($spin)."
        isa(nambu, Int) && @assert nambu∈(1, 2) "MoireSpinor error: wrong nambu ($nambu)."
        new{typeof(valley), typeof(layer), typeof(sublattice), typeof(spin), typeof(nambu)}(valley, layer, sublattice, spin, nambu)
    end
end
@inline Base.adjoint(spinor::MoireSpinor{<:Union{Int, Colon}, <:Union{Int, Colon}, <:Union{Int, Colon}, <:Union{Rational{Int}, Colon}, Int}) = MoireSpinor(spinor.valley, spinor.layer, spinor.sublattice, spinor.spin, 3-spinor.nambu)
@inline Base.show(io::IO, spinor::MoireSpinor) = @printf io "MoireSpinor(%s)" join((spinor.valley|>default, spinor.layer|>default, spinor.sublattice|>default, spinor.spin|>default, spinor.nambu|>default), ", ")
@inline statistics(::Type{<:MoireSpinor}) = :f
default(::Colon) = ":"
default(value::Int) = string(value)
default(value::Rational{Int}) = value.den==1 ? string(value.num) : string(value)

"""
    MoireSpace <: SimpleInternal{MoireSpinor{Int, Int, Int, Rational{Int}, Int}}

The internal degrees of freedom of Moire systems.
"""
struct MoireSpace <: SimpleInternal{MoireSpinor{Int, Int, Int, Rational{Int}, Int}}
    nvalley::Int
    nlayer::Int
    nsublattice::Int
    nspin::Int
end
@inline shape(moire::MoireSpace) = (1:moire.nvalley, 1:moire.nlayer, 1:moire.nsublattice, 1:moire.nspin, 1:2)
@inline Base.CartesianIndex(spinor::MoireSpinor, moire::MoireSpace) = CartesianIndex(spinor.valley, spinor.layer, spinor.sublattice, Int(spinor.spin+(moire.nspin-1)//2)+1, spinor.nambu)
@inline MoireSpinor(index::CartesianIndex{5}, moire::MoireSpace) = MoireSpinor(index[1], index[2], index[3], index[4]-1-(moire.nspin-1)//2, index[5])

# LaTeX formatted outputs
@inline script(::Val{:valley}, spinor::MoireSpinor; kwargs...) = spinor.valley==(:) ? ":" : string(spinor.valley)
@inline script(::Val{:layer}, spinor::MoireSpinor; kwargs...) = spinor.layer==(:) ? ":" : string(spinor.layer)
@inline script(::Val{:sublattice}, spinor::MoireSpinor; kwargs...) = spinor.sublattice==(:) ? ":" : string(spinor.sublattice)
@inline script(::Val{:spin}, spinor::MoireSpinor; kwargs...) = spinor.spin==(:) ? ":" : spinor.spin==0 ? "" : spinor.spin==1//2 ? "↑" : "↓"
@inline script(::Val{:nambu}, spinor::MoireSpinor; kwargs...) = spinor.nambu==(:) ? ":" : spinor.nambu==2 ? "\\dagger" : ""
@inline latexname(::Type{<:MoireSpinor}) = Symbol("MoireSpinor")
@inline latexname(::Type{<:Index{<:Union{Int, Colon}, <:MoireSpinor}}) = Symbol("Index{Union{Int, Colon}, MoireSpinor}")
@inline latexname(::Type{<:CompositeIndex{<:Index{<:Union{Int, Colon}, <:MoireSpinor}}}) = Symbol("CompositeIndex{Index{Union{Int, Colon}, MoireSpinor}}")

# Coupling related
@inline isdefinite(::Type{MoireSpinor{Int, Int, Int, Rational{Int}, Int}}) = true
@inline iidtype(::Type{MoireSpinor}, ::Type{V}, ::Type{L}, ::Type{S}, ::Type{P}, ::Type{N}) where {V<:Union{Int, Colon}, L<:Union{Int, Colon}, S<:Union{Int, Colon}, P<:Union{Rational{Int}, Colon}, N<:Union{Int, Colon}} = MoireSpinor{V, L, S, P, N}
@inline diagonalizablefields(::Type{<:Index{<:Union{Int, Colon}, <:MoireSpinor}}) = (:valley, :layer, :sublattice, :spin)
@inline function CompositeIID(coupling::Coupling{V, I}, info::Val=Val(:term)) where {V, I<:ID{<:Index{<:Union{Int, Colon}, <:MoireSpinor}}}
    return CompositeIID(map((spinor::MoireSpinor, order::Int)->MoireSpinor(spinor.valley, spinor.layer, spinor.sublattice, spinor.spin, nambu(spinor.nambu, order)), coupling.indexes.iids, ntuple(i->i, Val(rank(I)))))
end
@inline nambu(::Colon, order::Int) = order%2==1 ? 2 : 1
@inline nambu(nambu::Int, ::Int) = nambu
@inline function shape(iidspace::IIDSpace{<:MoireSpinor, MoireSpace})
    valley = moireshape(iidspace.iid.valley, iidspace.internal.nvalley)
    layer = moireshape(iidspace.iid.layer, iidspace.internal.nlayer)
    sublattice = moireshape(iidspace.iid.sublattice, iidspace.internal.nsublattice)
    spin = moireshape(iidspace.iid.spin, iidspace.internal.nspin)
    nambu = iidspace.iid.nambu:iidspace.iid.nambu
    return (valley, layer, sublattice, spin, nambu)
end
@inline moireshape(::Colon, n::Int) = 1:n
@inline moireshape(v::Int, n::Int) = (@assert(v<=n, "shape error: out of range."); v:v)
@inline moireshape(v::Rational{Int}, n::Int) = (@assert(abs(v)<=(n-1)//2, "shape error: out of range."); index=Int(v+(n-1)//2)+1; index:index)

"""
    MoireSystem{H<:OperatorGenerator, Hₘ<:Image} <: AbstractTBA{Fermionic{:TBA}, Hₘ, Nothing}

The continuum model of Moire systems.
"""
abstract type MoireSystem{H<:OperatorGenerator, Hₘ<:Image} <: AbstractTBA{Fermionic{:TBA}, Hₘ, Nothing} end
@inline getcontent(moire::MoireSystem, ::Val{:commutator}) = nothing
@inline function update!(moire::MoireSystem; parameters...)
    moire.parameters = update(moire.parameters; parameters...)
    update!(getcontent(moire, :H); parameters...)
    update!(getcontent(moire, :Hₘ); parameters...)
end
@inline function matrix(moire::MoireSystem; k, kwargs...)
    nblock = count(moire)
    reciprocallattice = getcontent(moire, :reciprocallattice)
    diagonal! = getcontent(moire, :diagonal!)
    result = zeros(valtype(moire), dimension(moire), dimension(moire))
    for i = 1:length(reciprocallattice)
        diagonal!(result, moire.parameters..., k+reciprocallattice[i]+reciprocallattice.Γ, reciprocallattice.K₊, reciprocallattice.K₋; offset=(i-1)*nblock)
    end
    for operator in getcontent(moire, :Hₘ)
        result[operator.position...] += operator.value
    end
    return result
end

"""
    BLTMD{L<:MoireReciprocalLattice, D<:Function, H<:OperatorGenerator} <: MoireSystem{H}

Twisted transition metal dichalcogenide homobilayers.
"""
mutable struct BLTMD{L<:MoireReciprocalLattice, D<:Function, H<:OperatorGenerator, Hₘ<:Image} <: MoireSystem{H, Hₘ}
    parameters::NamedTuple{(:a₀, :m, :θ, :Vᶻ, :μ), NTuple{5, Float64}}
    const reciprocallattice::L
    const diagonal!::D
    const H::H
    const Hₘ::Hₘ
end
@inline Base.count(bltmd::BLTMD) = (bltmd.H.hilbert)[1].nlayer * (bltmd.H.hilbert)[1].nsublattice

"""
    BLTMD(a₀::Number, m::Number, θ::Number, Vᶻ::Number, μ::Number, V::Number, ψ::Number, w::Number; truncation::Int=4)

[Continuum model of twisted transition metal dichalcogenide homobilayers](https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.122.086402).

Here, the parameters are as follows:
* `a₀`: monolayer lattice constant (Å)
* `m`: effective mass of the conduction band (mₑ)
* `θ`: twist angle (°)
* `Vᶻ`: perpendicular displacement field (meV)
* `μ`: chemical potential (meV)
* `V`: amplitude of Moire potential (meV)
* `ψ`: phase of Moire potential (°)
* `w`: interlayer hopping amplitude (meV)
"""
function BLTMD(a₀::Number, m::Number, θ::Number, Vᶻ::Number, μ::Number, V::Number, ψ::Number, w::Number; truncation::Int=4)
    reciprocallattice = MoireReciprocalLattice(truncation)
    terms = (
        Term{:TMD}(:potentialᵣ, V*cosd(ψ), 1, Coupling{2}(:, MoireSpinor, :, :, :, :, :), false),
        Term{:TMD}(:potentialᵢ, V*sind(ψ), 1, bond::Bond->(sign=round(Int, real(exp(3im*azimuth(rcoordinate(bond)))))::Int; (Coupling(-1im*sign, :, MoireSpinor, :, (1, 1), :, :, :), Coupling(1im*sign, :, MoireSpinor, :, (2, 2), :, :, :))), false),
        Term{:TMD}(:interlayer₁, w, 0, Coupling(:, MoireSpinor, :, (2, 1), :, :, :), false),
        Term{:TMD}(:interlayer₂, w, 1, bond::Bond->(ϕ=azimuthd(rcoordinate(bond)); ϕ≈60 ? Coupling(:, MoireSpinor, :, (2, 1), :, :, :) : ϕ≈240 ? Coupling(:, MoireSpinor, :, (1, 2), :, :, :) : Coupling(0, :, MoireSpinor, :, (0, 0), :, :, :)), false),
        Term{:TMD}(:interlayer₃, w, 1, bond::Bond->(ϕ=azimuthd(rcoordinate(bond)); ϕ≈120 ? Coupling(:, MoireSpinor, :, (2, 1), :, :, :) : ϕ≈300 ? Coupling(:, MoireSpinor, :, (1, 2), :, :, :) : Coupling(0, :, MoireSpinor, :, (0, 0), :, :, :)), false),
    )
    hilbert = Hilbert(site=>MoireSpace(1, 2, 1, 1) for site=1:length(reciprocallattice))
    H = OperatorGenerator(terms, bonds(reciprocallattice, 1), hilbert, plain, lazy; half=false)
    table = Table(hilbert, OperatorUnitToTuple(:site, :layer))
    return BLTMD((a₀=a₀, m=m, θ=θ, Vᶻ=Vᶻ, μ=μ), reciprocallattice, bltmd!, H, QuadraticFormalize{Fermionic{:TBA}}(table)(H))
end
@inline function bltmd!(dest, a₀, m, θ, Vᶻ, μ, k, K₊, K₋; offset)
    m₀ = 0.0001312169949060677
    m = m₀*m*a₀^2/(2sind(θ/2))^2
    dest[offset+1, offset+1] = - mapreduce(x->x^2, +, k-K₊)/2m + Vᶻ - μ
    dest[offset+2, offset+2] = - mapreduce(x->x^2, +, k-K₋)/2m - Vᶻ - μ
    return dest
end
@inline function bltmdmap(parameters)
    return (
        a₀=parameters[:a₀],
        m=parameters[:m],
        θ=parameters[:θ],
        Vᶻ=parameters[:Vᶻ],
        μ=parameters[:μ],
        potentialᵣ=parameters[:V]*cosd(parameters[:ψ]),
        potentialᵢ=Complex(parameters[:V]*sind(parameters[:ψ])),
        interlayer₁=parameters[:w],
        interlayer₂=parameters[:w],
        interlayer₃=parameters[:w]
    )
end

"""
    coefficients(bltmd::BLTMD, lattice::MoireTriangular, brillouinzone::BrillouinZone; band::Int=dimension(bltmd)) -> Tuple{Vector{Vector{ComplexF64}}, Float64}
    coefficients(bltmd::Algorithm{<:BLTMD}, lattice::MoireTriangular, brillouinzone::BrillouinZone; band::Int=dimension(bltmd.frontend)) -> Tuple{Vector{Vector{ComplexF64}}, Float64}

Get the coefficients of the hoppings and chemical potential of a bilayer TMD on the emergent triangular lattice.
"""
function coefficients(bltmd::BLTMD, lattice::MoireTriangular, brillouinzone::BrillouinZone; band::Int=dimension(bltmd))
    hoppings, μ = [zeros(ComplexF64, length(neighbor)) for neighbor in lattice.neighbors], 0.0
    for momentum in brillouinzone
        value = eigvals(matrix(bltmd; k=momentum))[band]/length(brillouinzone)
        for i = 1:length(lattice.neighbors)
            for j = 1:length(lattice.neighbors[i])
                hoppings[i][j] += exp(-1im*dot(momentum, lattice.neighbors[i][j]))*value
            end
        end
        μ += value
    end
    return hoppings, μ
end
@inline function coefficients(bltmd::Algorithm{<:BLTMD}, lattice::MoireTriangular, brillouinzone::BrillouinZone; band::Int=dimension(bltmd.frontend))
    return coefficients(bltmd.frontend, lattice, brillouinzone; band=band)
end 

"""
    terms(bltmd::BLTMD, lattice::MoireTriangular, brillouinzone::BrillouinZone; band::Int=dimension(bltmd), ismodulatable::Bool=true, atol=atol) -> NTuple{2*truncation(lattice)+1, Term}
    terms(bltmd::Algorithm{<:BLTMD}, lattice::MoireTriangular, brillouinzone::BrillouinZone; band::Int=dimension(bltmd.frontend), ismodulatable::Bool=true, atol=atol) -> NTuple{2*truncation(lattice)+1, Term}

Get the hopping terms and chemical potential of a bilayer TMD on the emergent triangular lattice.
"""
function terms(bltmd::BLTMD, lattice::MoireTriangular, brillouinzone::BrillouinZone; band::Int=dimension(bltmd), ismodulatable::Bool=true, atol=atol)
    tvals, μval = coefficients(bltmd, lattice, brillouinzone; band=band)
    hoppings = map(NTuple{truncation(lattice), eltype(tvals)}(tvals), lattice.neighbors, ntuple(i->i, Val(truncation(lattice)))) do values, neighbor, order
        @assert all(value->isapprox(real(value), real(values[1]); atol=atol) && isapprox(abs(imag(value)), abs(imag(values[1])); atol=atol), values) "terms error: unexpected behavior."
        θs = ntuple(i->azimuthd(neighbor[i]), length(neighbor))
        signs = ntuple(i->isapprox(imag(values[i]), 0; atol=atol) ? 1 : round(Int, imag(values[1])/imag(values[i])), length(values))
        function amplitude(bond::Bond)
            θ = azimuthd(rcoordinate(bond))
            for (sign, θ₀) in zip(signs, θs)
                Δ = (θ-θ₀)/60
                isapprox(round(Int, Δ), Δ; atol=atol) && return -1im*sign*cosd(3*(θ-θ₀))
            end
            error("amplitude error: mismatched bond.")
        end
        suffix = join('₀'+d for d in digits(order))
        return (
            Hopping(Symbol("t", suffix), real(values[1]), order; ismodulatable=ismodulatable),
            Hopping(Symbol("λ", suffix), imag(values[1]), order, MatrixCoupling(:, FID, :, σ"z", :); amplitude=amplitude, ismodulatable=ismodulatable)
        )
    end
    μ = Onsite(:μ, Complex(μval))
    return (concatenate(hoppings...)..., μ)
end
@inline function terms(bltmd::Algorithm{<:BLTMD}, lattice::MoireTriangular, brillouinzone::BrillouinZone; band::Int=dimension(bltmd.frontend), ismodulatable::Bool=true, atol=atol)
    return terms(bltmd.frontend, lattice, brillouinzone; band=band, ismodulatable=ismodulatable, atol=atol)
end

# runtime initialization
function __init__()
    latexformat(MoireSpinor, LaTeX{(:nambu,), (:layer,)}('c'))
    latexformat(Index{<:Union{Int, Colon}, <:MoireSpinor}, LaTeX{(:nambu,), (:site, :layer)}('c'))
    latexformat(CompositeIndex{<:Index{<:Union{Int, Colon}, <:MoireSpinor}}, LaTeX{(:nambu,), (:site, :layer)}('c'))
    nothing
end

end