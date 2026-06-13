#=== internal degrees of freedom and continuum model for Moire systems ===#

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
@inline Base.count(moire::Algorithm{<:MoireSystem}) = count(moire.frontend)

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
    reciprocallattice = MoireTriangularReciprocal(truncation)
    hilbert = Hilbert(site=>MoireSpace(1, 2, 1, 1) for site=1:length(reciprocallattice))
    system = OperatorGenerator(bonds(reciprocallattice, 1), hilbert, terms; half=false)
    table = Table(hilbert, OperatorIndexToTuple(:site, :layer))
    quadraticization = Quadraticization{Fermionic{:TBA}}(table)
    return BLTMD((a₀=a₀, m=m, θ=θ, Vᶻ=Vᶻ, μ=μ), reciprocallattice, bltmd!, system, quadraticization, quadraticization(system))
end
@inline function bltmd!(dest, a₀, m, θ, Vᶻ, μ, k, K₊, K₋; offset)
    m₀ = 0.0001312169949060677
    m = m₀*m*a₀^2/(2sind(θ/2))^2
    dest[offset+1, offset+1] = - mapreduce(x->x^2, +, k-K₊)/2m + Vᶻ/2 - μ
    dest[offset+2, offset+2] = - mapreduce(x->x^2, +, k-K₋)/2m - Vᶻ/2 - μ
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
