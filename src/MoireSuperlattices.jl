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

include("lattices.jl")
include("systems.jl")

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

end