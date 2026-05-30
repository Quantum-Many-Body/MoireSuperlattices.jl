module MoireSuperlattices

using LinearAlgebra: dot, eigvals, eigvecs, norm
using Printf: @printf
using QuantumLattices: annihilation, atol, creation, σᶻ
using QuantumLattices: AbstractLattice, Bond, BrillouinZone, CategorizedGenerator, CompositeIndex, Coupling, Hilbert, Hopping, Index, InternalIndex, LaTeX, Neighbors, Onsite, OperatorGenerator, OperatorIndexToTuple, OperatorSum, ReciprocalZone, SimpleInternal, Table, Term
using QuantumLattices: azimuth, azimuthd, bonds, concatenate, decompose, distance, latexformat, periods, reciprocals, rcoordinate, rotate, scalartype, str, update, volume, 𝕔⁺𝕔
using StaticArrays: SVector
using TightBindingApproximation: TBA, Fermionic, Quadratic, Quadraticization

import QuantumLattices: Algorithm, Lattice, Parameters, contentnames, dimension, getcontent, indextype, isdefinite, latexname, matrix, script, shape, statistics, update!

export BLTMD, CommensurateBilayerHoneycomb, BareCoulomb, CoulombIntegral, ImageCoulomb, MoireReciprocalLattice, MoireSpace, MoireSpinor, MoireSuperlattice, MoireSystem, MoireTriangular, MoireTriangularWannier, RealZone, TanhCoulomb, bltmd!, bltmdmap, coefficients, terms, truncation, vectors

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
    MoireSpinor{V<:Union{Int, Colon}, L<:Union{Int, Colon}, S<:Union{Int, Colon}, P<:Union{Rational{Int}, Colon}} <: InternalIndex

The index of the internal degrees of freedom of Moire systems.
"""
struct MoireSpinor{V<:Union{Int, Colon}, L<:Union{Int, Colon}, S<:Union{Int, Colon}, P<:Union{Rational{Int}, Colon}} <: InternalIndex
    valley::V
    layer::L
    sublattice::S
    spin::P
    nambu::Int
    function MoireSpinor(valley::Union{Int, Colon}, layer::Union{Int, Colon}, sublattice::Union{Int, Colon}, spin::Union{Rational{Int}, Int, Colon}, nambu::Int)
        @assert spin∈(-1//2, 1//2, 0, :) "MoireSpinor error: incorrect spin ($spin)."
        @assert nambu∈(1, 2) "MoireSpinor error: wrong nambu ($nambu)."
        isa(spin, Int) && (spin = convert(Rational{Int}, spin))
        new{typeof(valley), typeof(layer), typeof(sublattice), typeof(spin)}(valley, layer, sublattice, spin, nambu)
    end
end
# basic methods of concrete InternalIndex
@inline Base.adjoint(spinor::MoireSpinor) = MoireSpinor(spinor.valley, spinor.layer, spinor.sublattice, spinor.spin, 3-spinor.nambu)
@inline statistics(::Type{<:MoireSpinor}) = :f
@inline isdefinite(::Type{MoireSpinor{Int, Int, Int, Rational{Int}}}) = true
@inline Base.show(io::IO, spinor::MoireSpinor) = @printf io "MoireSpinor(%s)" join((spinor.valley|>default, spinor.layer|>default, spinor.sublattice|>default, spinor.spin|>default, spinor.nambu|>default), ", ")
@inline default(::Colon) = ":"
@inline default(value::Int) = string(value)
@inline default(value::Rational{Int}) = value.den==1 ? string(value.num) : string(value)
# requested by MatrixCoupling
@inline indextype(::Type{MoireSpinor}, ::Type{V}, ::Type{L}, ::Type{S}, ::Type{P}) where {V<:Union{Int, Colon}, L<:Union{Int, Colon}, S<:Union{Int, Colon}, P<:Union{Rational{Int}, Colon}} = MoireSpinor{V, L, S, P}
# patternrule
@inline MoireSpinor{V, L, S, P}(valley, layer, sublattice, spin, nambu) where {V, L, S, P} = MoireSpinor(valley, layer, sublattice, spin, nambu)
# LaTeX format output
@inline script(spinor::MoireSpinor, ::Val{:valley}; kwargs...) = spinor.valley==(:) ? ":" : string(spinor.valley)
@inline script(spinor::MoireSpinor, ::Val{:layer}; kwargs...) = spinor.layer==(:) ? ":" : string(spinor.layer)
@inline script(spinor::MoireSpinor, ::Val{:sublattice}; kwargs...) = spinor.sublattice==(:) ? ":" : string(spinor.sublattice)
@inline script(spinor::MoireSpinor, ::Val{:spin}; kwargs...) = spinor.spin==(:) ? ":" : spinor.spin==0 ? "" : spinor.spin==1//2 ? "↑" : "↓"
@inline script(spinor::MoireSpinor, ::Val{:nambu}; kwargs...) = spinor.nambu==2 ? "\\dagger" : ""
@inline latexname(::Type{<:MoireSpinor}) = Symbol("MoireSpinor")
@inline latexname(::Type{<:Index{<:MoireSpinor}}) = Symbol("Index{MoireSpinor}")
@inline latexname(::Type{<:CompositeIndex{<:Index{<:MoireSpinor}}}) = Symbol("CompositeIndex{Index{MoireSpinor}}")

"""
    MoireSpace <: SimpleInternal{MoireSpinor{Int, Int, Int, Rational{Int}}}

The internal degrees of freedom of Moire systems.
"""
struct MoireSpace <: SimpleInternal{MoireSpinor{Int, Int, Int, Rational{Int}}}
    nvalley::Int
    nlayer::Int
    nsublattice::Int
    nspin::Int
end
@inline shape(moire::MoireSpace) = (1:moire.nvalley, 1:moire.nlayer, 1:moire.nsublattice, 1:moire.nspin, 1:2)
@inline Base.convert(::Type{<:CartesianIndex}, spinor::MoireSpinor, moire::MoireSpace) = CartesianIndex(spinor.valley, spinor.layer, spinor.sublattice, Int(spinor.spin+(moire.nspin-1)//2)+1, spinor.nambu)
@inline Base.convert(::Type{<:MoireSpinor}, index::CartesianIndex{5}, moire::MoireSpace) = MoireSpinor(index[1], index[2], index[3], index[4]-1-(moire.nspin-1)//2, index[5])
# requested by Coupling expansion
@inline function shape(moire::MoireSpace, spinor::MoireSpinor)
    valley = moireshape(spinor.valley, moire.nvalley)
    layer = moireshape(spinor.layer, moire.nlayer)
    sublattice = moireshape(spinor.sublattice, moire.nsublattice)
    spin = moireshape(spinor.spin, moire.nspin)
    nambu = spinor.nambu:spinor.nambu
    return (valley, layer, sublattice, spin, nambu)
end
@inline moireshape(::Colon, n::Int) = 1:n
@inline moireshape(v::Int, n::Int) = (@assert(v<=n, "shape error: out of range."); v:v)
@inline moireshape(v::Rational{Int}, n::Int) = (@assert(abs(v)<=(n-1)//2, "shape error: out of range."); index=Int(v+(n-1)//2)+1; index:index)

"""
    MoireSystem{P<:Parameters, L<:MoireReciprocalLattice, D<:Function, S<:OperatorGenerator, Q<:Quadraticization, H<:CategorizedGenerator{<:OperatorSum{<:Quadratic}}} <: TBA{Fermionic{:TBA}, H, Nothing}

The continuum model of Moire systems.
"""
abstract type MoireSystem{P<:Parameters, L<:MoireReciprocalLattice, D<:Function, S<:OperatorGenerator, Q<:Quadraticization, H<:CategorizedGenerator{<:OperatorSum{<:Quadratic}}} <: TBA{Fermionic{:TBA}, H, Nothing} end
@inline contentnames(::Type{<:MoireSystem}) = (:parameters, :reciprocallattice, :diagonal!, :system, :quadraticization, :H)
@inline Parameters(moire::MoireSystem) = (; moire.parameters..., Parameters(getcontent(moire, :system))...)
@inline dimension(moire::MoireSystem) = length(getcontent(moire, :quadraticization).table)
@inline function update!(moire::MoireSystem; parameters...)
    moire.parameters = update(moire.parameters; parameters...)
    update!(getcontent(moire, :system); parameters...)
    update!(getcontent(moire, :H); parameters...)
end
@inline function matrix(moire::MoireSystem, k::AbstractVector{<:Number}; kwargs...)
    nblock = count(moire)
    reciprocallattice = getcontent(moire, :reciprocallattice)
    diagonal! = getcontent(moire, :diagonal!)
    result = zeros(scalartype(moire), dimension(moire), dimension(moire))
    for i = 1:length(reciprocallattice)
        diagonal!(result, moire.parameters..., k+reciprocallattice[i]+reciprocallattice.Γ, reciprocallattice.K₊, reciprocallattice.K₋; offset=(i-1)*nblock)
    end
    for operator in getcontent(moire, :H)
        result[operator.position...] += operator.value
    end
    return result
end

"""
    BLTMD{
        L<:MoireReciprocalLattice,
        D<:Function,
        S<:OperatorGenerator,
        Q<:Quadraticization,
        H<:CategorizedGenerator{<:OperatorSum{<:Quadratic}}
    } <: MoireSystem{NamedTuple{(:a₀, :m, :θ, :Vᶻ, :μ), NTuple{5, Float64}}, L, D, S, Q, H}

Twisted transition metal dichalcogenide homobilayers.
"""
mutable struct BLTMD{
    L<:MoireReciprocalLattice,
    D<:Function,
    S<:OperatorGenerator,
    Q<:Quadraticization,
    H<:CategorizedGenerator{<:OperatorSum{<:Quadratic}}
} <: MoireSystem{NamedTuple{(:a₀, :m, :θ, :Vᶻ, :μ), NTuple{5, Float64}}, L, D, S, Q, H}
    parameters::NamedTuple{(:a₀, :m, :θ, :Vᶻ, :μ), NTuple{5, Float64}}
    const reciprocallattice::L
    const diagonal!::D
    const system::S
    const quadraticization::Q
    const H::H
end
@inline Base.count(bltmd::BLTMD) = (bltmd.system.hilbert)[1].nlayer * (bltmd.system.hilbert)[1].nsublattice

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
    nambus = (creation, annihilation)
    coupling = Coupling{MoireSpinor}(:, :, :, :, :, nambus)
    coupling₁₁ = Coupling{MoireSpinor}(:, :, (1, 1), :, :, nambus)
    coupling₁₂ = Coupling{MoireSpinor}(:, :, (1, 2), :, :, nambus)
    coupling₂₁ = Coupling{MoireSpinor}(:, :, (2, 1), :, :, nambus)
    coupling₂₂ = Coupling{MoireSpinor}(:, :, (2, 2), :, :, nambus)
    coupling₀ = Coupling{MoireSpinor}(0, :, :, (0, 0), :, :, nambus)
    terms = (
        Term{:TMD}(:potentialᵣ, V*cosd(ψ), 1, coupling, false),
        Term{:TMD}(:potentialᵢ, V*sind(ψ), 1, bond::Bond->(sign=round(Int, real(exp(3im*azimuth(rcoordinate(bond)))))::Int; (-1im*sign*coupling₁₁, 1im*sign*coupling₂₂)), false),
        Term{:TMD}(:interlayer₁, w, 0, coupling₂₁, false),
        Term{:TMD}(:interlayer₂, w, 1, bond::Bond->(ϕ=azimuthd(rcoordinate(bond)); ϕ≈60 ? coupling₂₁ : ϕ≈240 ? coupling₁₂ : coupling₀), false),
        Term{:TMD}(:interlayer₃, w, 1, bond::Bond->(ϕ=azimuthd(rcoordinate(bond)); ϕ≈120 ? coupling₂₁ : ϕ≈300 ? coupling₁₂ : coupling₀), false),
    )
    reciprocallattice = MoireReciprocalLattice(truncation)
    hilbert = Hilbert(site=>MoireSpace(1, 2, 1, 1) for site=1:length(reciprocallattice))
    system = OperatorGenerator(bonds(reciprocallattice, 1), hilbert, terms; half=false)
    table = Table(hilbert, OperatorIndexToTuple(:site, :layer))
    quadraticization = Quadraticization{Fermionic{:TBA}}(table)
    return BLTMD((a₀=a₀, m=m, θ=θ, Vᶻ=Vᶻ, μ=μ), reciprocallattice, bltmd!, system, quadraticization, quadraticization(system))
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
    Algorithm(name::Symbol, bltmd::BLTMD, parameters::Parameters; kwargs...)

Construct an `Algorithm` with a `BLTMD` as the frontend.
"""
@inline function Algorithm(name::Symbol, bltmd::BLTMD, parameters::Parameters; kwargs...)
    return Algorithm(name, bltmd, parameters, bltmdmap; kwargs...)
end

"""
    coefficients(bltmd::BLTMD, lattice::MoireTriangular, brillouinzone::BrillouinZone; band::Int=dimension(bltmd)) -> Tuple{Vector{Vector{ComplexF64}}, Float64}
    coefficients(bltmd::Algorithm{<:BLTMD}, lattice::MoireTriangular, brillouinzone::BrillouinZone; band::Int=dimension(bltmd.frontend)) -> Tuple{Vector{Vector{ComplexF64}}, Float64}

Get the coefficients of the hoppings and chemical potential of a bilayer TMD on the emergent triangular lattice.
"""
function coefficients(bltmd::BLTMD, lattice::MoireTriangular, brillouinzone::BrillouinZone; band::Int=dimension(bltmd))
    hoppings, μ = [zeros(ComplexF64, length(neighbor)) for neighbor in lattice.neighbors], 0.0
    for momentum in brillouinzone
        value = eigvals(bltmd, momentum)[band]/length(brillouinzone)
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
    terms(bltmd::BLTMD, lattice::MoireTriangular, brillouinzone::BrillouinZone; band::Int=dimension(bltmd), ismodulatable::Bool=true, tol=atol) -> NTuple{2*truncation(lattice)+1, Term}
    terms(bltmd::Algorithm{<:BLTMD}, lattice::MoireTriangular, brillouinzone::BrillouinZone; band::Int=dimension(bltmd.frontend), ismodulatable::Bool=true, tol=atol) -> NTuple{2*truncation(lattice)+1, Term}

Get the hopping terms and chemical potential of a bilayer TMD on the emergent triangular lattice.
"""
function terms(bltmd::BLTMD, lattice::MoireTriangular, brillouinzone::BrillouinZone; band::Int=dimension(bltmd), ismodulatable::Bool=true, tol=atol)
    tvals, μval = coefficients(bltmd, lattice, brillouinzone; band=band)
    hoppings = map(NTuple{truncation(lattice), eltype(tvals)}(tvals), lattice.neighbors, ntuple(i->i, Val(truncation(lattice)))) do values, neighbor, order
        @assert all(value->isapprox(real(value), real(values[1]); atol=tol) && isapprox(abs(imag(value)), abs(imag(values[1])); atol=tol), values) "terms error: unexpected behavior."
        θs = ntuple(i->azimuthd(neighbor[i]), length(neighbor))
        signs = ntuple(i->isapprox(imag(values[i]), 0; atol=tol) ? 1 : round(Int, imag(values[1])/imag(values[i])), length(values))
        function amplitude(bond::Bond)
            θ = azimuthd(rcoordinate(bond))
            for (sign, θ₀) in zip(signs, θs)
                Δ = (θ-θ₀)/60
                isapprox(round(Int, Δ), Δ; atol=tol) && return -1im*sign*cosd(3*(θ-θ₀))
            end
            error("amplitude error: mismatched bond.")
        end
        suffix = join('₀'+d for d in digits(order))
        return (
            Hopping(Symbol("t", suffix), real(values[1]), order; ismodulatable=ismodulatable),
            Hopping(Symbol("λ", suffix), imag(values[1]), order, 𝕔⁺𝕔(:, :, σᶻ); amplitude=amplitude, ismodulatable=ismodulatable)
        )
    end
    μ = Onsite(:μ, Complex(μval))
    return (concatenate(hoppings...)..., μ)
end
@inline function terms(bltmd::Algorithm{<:BLTMD}, lattice::MoireTriangular, brillouinzone::BrillouinZone; band::Int=dimension(bltmd.frontend), ismodulatable::Bool=true, tol=atol)
    return terms(bltmd.frontend, lattice, brillouinzone; band=band, ismodulatable=ismodulatable, tol=tol)
end

"""
    RealZone{N, S<:SVector, V<:Number}

A rectangular zone in real space.

Alias for [`ReciprocalZone`](@ref) with the space-type parameter `K = :r`.
"""
const RealZone{N, S<:SVector, V<:Number} = ReciprocalZone{:r, N, S, V}

"""
    MoireTriangularWannier{L<:MoireTriangular, G<:MoireReciprocalLattice, B<:BrillouinZone}

Wannier function constructed on the emergent triangular lattice of a Moire system.

Fields:
- `aₘ::Float64` — lattice constant of the moire superlattice
- `lattice::L` — emergent triangular lattice (site positions and neighbor shells)
- `reciprocallattice::G` — truncated plane-wave basis (G-vectors) from the continuum model
- `brillouinzone::B` — uniform k-point mesh over the moire Brillouin zone
- `bloch::Matrix{ComplexF64}` — gauge-fixed Bloch eigenvectors, (D × N_k)
"""
struct MoireTriangularWannier{L<:MoireTriangular, G<:MoireReciprocalLattice, B<:BrillouinZone}
    aₘ::Float64
    lattice::L
    reciprocallattice::G
    brillouinzone::B
    bloch::Matrix{ComplexF64}
end

"""
    MoireTriangularWannier(moiresystem::MoireSystem, lattice::MoireTriangular, brillouinzone::BrillouinZone; band::Int=dimension(moiresystem))

Construct the Wannier function for band `band` of `moiresystem`, localized on the sites of the emergent `lattice`.

Steps:
1. Diagonalize the Hamiltonian at each k-point in `brillouinzone`.
2. Gauge-fix: choose the phase so the bottom-layer component at r=0 is real and positive (Appendix A of PRR 2, 033087).
"""
function MoireTriangularWannier(moiresystem::MoireSystem, lattice::MoireTriangular, brillouinzone::BrillouinZone; band::Int=dimension(moiresystem))
    dim = dimension(moiresystem)
    @assert 1 <= band <= dim "MoireTriangularWannier error: band index $band out of range [1, $dim]."
    bloch = zeros(ComplexF64, dim, length(brillouinzone))
    for (i, k) in enumerate(brillouinzone)
        psi = eigvecs(moiresystem, k)[:, band]
        bottom = zero(ComplexF64)
        for j in 1:2:dim
            bottom += psi[j]
        end
        @assert abs(bottom) > atol "MoireTriangularWannier error: bottom-layer component at r=0 is zero at $k; gauge fixing failed."
        phase = conj(bottom) / abs(bottom)
        @views bloch[:, i] .= psi .* phase
    end
    aₘ = moiresystem.parameters.a₀ / (2sind(moiresystem.parameters.θ/2))
    reciprocallattice = getcontent(moiresystem, :reciprocallattice)
    return MoireTriangularWannier(aₘ, lattice, reciprocallattice, brillouinzone, bloch)
end

"""
    (wannier::MoireTriangularWannier)(r::AbstractVector{<:Number}) -> SVector{2, ComplexF64}

Evaluate the Wannier function at real-space position `r`: W(r) = (1/√N) Σ_k ψ_k(r),  ψ_k(r) = (1/√NΩ) Σ_G c_{k,G} e^{i(k+G)·r}

Returns a 2-component layer-pseudospin spinor [W_b(r), W_t(r)].
"""
function (wannier::MoireTriangularWannier)(r::AbstractVector{<:Number})
    dim, nₖ = size(wannier.bloch)
    nblock = dim ÷ length(wannier.reciprocallattice)
    @assert nblock == 2 "MoireTriangularWannier error: only 2-layer systems supported (got nblock=$nblock)."
    result = SVector(ComplexF64(0.0), ComplexF64(0.0))
    for (i, momentum) in enumerate(wannier.brillouinzone)
        k = SVector(momentum[1], momentum[2])
        for (j, G) in enumerate(wannier.reciprocallattice)
            phase = cis(dot(k + G, r))
            bot = nblock * (j - 1) + 1  # bottom-layer index within block
            top = nblock * (j - 1) + 2  # top-layer index within block
            result += SVector(wannier.bloch[bot, i], wannier.bloch[top, i]) * phase
        end
    end
    Ω = volume(wannier.lattice.vectors)
    return result / (sqrt(Ω) * nₖ)
end

"""
    CoulombIntegral{W<:MoireTriangularWannier}

Precomputed Coulomb form factor |M(q)|² for a triangular-lattice Wannier function.

Fields:
- `wannier::W` — reference to the Wannier function
- `qs::Vector{SVector{2,Float64}}` — unique q-vectors from the pairwise convolution
- `formfactor::Vector{Float64}` — |M(q)|² at each q-point
"""
struct CoulombIntegral{W<:MoireTriangularWannier}
    wannier::W
    qs::Vector{SVector{2,Float64}}
    formfactor::Vector{Float64}
end

"""
    CoulombIntegral(wannier::MoireTriangularWannier)

Construct by computing the form factor M(q) from Bloch coefficients via pairwise convolution.

Algorithm: collect all extended momenta p = k+G with their Bloch coefficients,
then for each pair (p, p'), accumulate dot(c(p), c(p')) into M(q = p-p'). Normalized by N_k at extraction.
The q-mesh emerges naturally from the G/G' truncation.
"""
function CoulombIntegral(wannier::MoireTriangularWannier)
    dim, nk = size(wannier.bloch)
    nG = length(wannier.reciprocallattice)
    nblock = dim ÷ nG
    @assert nblock == 2 "CoulombIntegral error: only 2-layer systems supported (got nblock=$nblock)."
    b₁, b₂ = wannier.reciprocallattice.translations
    N₁, N₂ = periods(wannier.brillouinzone)
    # Decompose each k-point and G-vector into integer coordinates in the (b₁, b₂) basis
    ks = Vector{Tuple{Int,Int}}(undef, nk)
    for (i, k) in enumerate(wannier.brillouinzone)
        f₁, f₂ = decompose(k, b₁, b₂)
        ks[i] = (round(Int, f₁*N₁), round(Int, f₂*N₂))
    end
    Gs = Vector{Tuple{Int,Int}}(undef, nG)
    for (i, G) in enumerate(wannier.reciprocallattice)
        g₁, g₂ = decompose(G, b₁, b₂)
        Gs[i] = (round(Int, g₁), round(Int, g₂))
    end
    # Pairwise convolution, accumulated into integer-keyed dict
    data = Vector{Tuple{Int, Int, SVector{2, ComplexF64}}}(undef, nk * nG)
    idx = 1
    for (ik, (k₁, k₂)) in enumerate(ks)
        for (ig, (g₁, g₂)) in enumerate(Gs)
            p₁ = k₁ + N₁*g₁
            p₂ = k₂ + N₂*g₂
            c₁ = wannier.bloch[nblock*(ig-1)+1, ik]
            c₂ = wannier.bloch[nblock*(ig-1)+2, ik]
            data[idx] = (p₁, p₂, SVector(c₁, c₂))
            idx += 1
        end
    end
    Ms = Dict{Tuple{Int,Int}, ComplexF64}()
    for (p₁, p₂, coeff) in data, (p₁′, p₂′, coeff′) in data
        q = (p₁′ - p₁, p₂′ - p₂)
        Ms[q] = get(Ms, q, zero(ComplexF64)) + dot(coeff, coeff′)
    end
    # Extract qs and |M(q)|²
    qs = SVector{2, Float64}[]
    formfactor = Float64[]
    for ((q₁, q₂), M) in Ms
        push!(qs, (q₁/N₁)*b₁ + (q₂/N₂)*b₂)
        push!(formfactor, abs2(M/nk))
    end
    return CoulombIntegral(wannier, qs, formfactor)
end

const e²_meV_Å = 14400.0  # e² = 1.44 eV·nm in meV·Å
"""
    BareCoulomb(ϵ::Real)

Bare (unscreened) 2D Coulomb potential: V(q) = 2π e² / (ϵ aₘ |q|).
q=0 skipped (returns 0).
"""
struct BareCoulomb
    ϵ::Float64
end
function (v::BareCoulomb)(q::Real, aₘ::Real)
    q < 1e-14 && return 0.0
    return 2π * e²_meV_Å / (v.ϵ * aₘ * q)
end

"""
    ImageCoulomb(ϵ::Real, d::Real)

Image-charge screened Coulomb: V(q) = 2π e² (1 - e^{-2dₘ|q|}) / (ϵ aₘ |q|) where dₘ = d/aₘ.
"""
struct ImageCoulomb
    ϵ::Float64
    d::Float64
end
function (v::ImageCoulomb)(q::Real, aₘ::Real)
    dₘ = v.d / aₘ
    q < 1e-14 && return 2π * e²_meV_Å * 2dₘ / v.ϵ
    return 2π * e²_meV_Å * (1 - exp(-2dₘ * q)) / (v.ϵ * aₘ * q)
end

"""
    TanhCoulomb(ϵ::Real, d::Real)

Gate-screened Coulomb: V(q) = 2π e² tanh(dₘ|q|) / (ϵ aₘ |q|) where dₘ = d/aₘ.
"""
struct TanhCoulomb
    ϵ::Float64
    d::Float64
end
function (v::TanhCoulomb)(q::Real, aₘ::Real)
    dₘ = v.d / aₘ
    q < 1e-14 && return 2π * e²_meV_Å * dₘ / v.ϵ
    return 2π * e²_meV_Å * tanh(dₘ * q) / (v.ϵ * aₘ * q)
end

"""
    (c::CoulombIntegral)(R::AbstractVector{<:Number}, V=BareCoulomb(1.0)) -> Float64

Compute U(R) = (1/(N_k Ω)) Σ_q V(|q|, aₘ) |M(q)|² e^{iq·R} in meV.

`V` is a callable `V(q::Real, aₘ::Real) -> Real`, e.g. `BareCoulomb(ϵ)`, `ImageCoulomb(ϵ, d)`, `TanhCoulomb(ϵ, d)`, or a user-defined function. Defaults to `BareCoulomb(1.0)` (unscreened, ε=1).
"""
function (c::CoulombIntegral)(R::AbstractVector{<:Number}, V=BareCoulomb(1.0))
    nk = length(c.wannier.brillouinzone)
    Ω = volume(c.wannier.lattice.vectors)
    aₘ = c.wannier.aₘ
    U = 0.0
    for (q, m²) in zip(c.qs, c.formfactor)
        Vq = V(norm(q), aₘ)
        U += Vq * m² * cos(dot(q, R))
    end
    return U / (nk * Ω)
end

# runtime initialization
function __init__()
    latexformat(MoireSpinor, LaTeX{(:nambu,), (:layer, :spin)}('c'))
    latexformat(Index{<:MoireSpinor}, LaTeX{(:nambu,), (:site, :layer, :spin)}('c'))
    latexformat(CompositeIndex{<:Index{<:MoireSpinor}}, LaTeX{(:nambu,), (:site, :layer, :spin)}('c'))
    nothing
end

end