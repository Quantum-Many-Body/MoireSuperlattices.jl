"""
    MoireWannier{L<:MoireSuperlattice, G<:MoireReciprocalLattice, B<:BrillouinZone}

Wannier function constructed on an emergent Moire superlattice.

Fields:
- `aₘ::Float64` — lattice constant of the moire superlattice
- `lattice::L` — emergent superlattice (MoireTriangular or MoireHoneycomb)
- `reciprocallattice::G` — truncated plane-wave basis (G-vectors) from the continuum model
- `brillouinzone::B` — uniform k-point mesh over the moire Brillouin zone
- `energies::Matrix{Float64}` — raw band energies, (nband, nk)
- `bloch::Array{ComplexF64, 4}` — raw Bloch eigenvectors, (nlayer, nG, nband, nk), pre-gauge
- `U::Array{ComplexF64, 3}` — gauge transformation matrices, (nband, nband, nk)
"""
struct MoireWannier{L<:MoireSuperlattice, G<:MoireReciprocalLattice, B<:BrillouinZone}
    aₘ::Float64
    lattice::L
    reciprocallattice::G
    brillouinzone::B
    energies::Matrix{Float64}
    bloch::Array{ComplexF64, 4}
    U::Array{ComplexF64, 3}
end

"""
    count(wannier::MoireWannier) -> Int

Return the number of Wannier orbitals (``n_{\\rm band}``).
"""
@inline Base.count(wannier::MoireWannier) = size(wannier.energies, 1)

"""
    (wannier::MoireWannier)(r::AbstractVector{<:Number}, sublattice::Int=1) -> Vector{ComplexF64}

Evaluate the Wannier function at real-space position `r` for a given sublattice.

Formula: Wₙˡ(r) = (1/N√Ω) Σ_{k, G, ν} bloch_{G, l, ν}(k) · U_{ν, n}(k) · e^{i(k+G)·r},
where N is the number of k-points and Ω is the volume of the unit cell in the real space.
"""
function (wannier::MoireWannier)(r::AbstractVector{<:Number}, sublattice::Int=1)
    nlayer, _, nband, nk = size(wannier.bloch)
    @assert 1 <= sublattice <= nband "MoireWannier error: sublattice $sublattice out of range [1, $nband]."
    result = zeros(ComplexF64, nlayer)
    for (ik, k) in enumerate(wannier.brillouinzone), (ig, G) in enumerate(wannier.reciprocallattice)
        phase = cis(dot(k + G, r))
        for ib in 1:nband, il in 1:nlayer
            result[il] += wannier.bloch[il, ig, ib, ik] * wannier.U[ib, sublattice, ik] * phase
        end
    end
    Ω = volume(wannier.lattice.vectors)
    return broadcast!(/, result, result, sqrt(Ω)*nk)
end

"""
    MoireWannier(moiresystem::Algorithm{<:MoireSystem}, args...; kwargs...)

Forward to the wrapped [`MoireSystem`](@ref) frontend.
"""
@inline MoireWannier(moiresystem::Algorithm{<:MoireSystem}, args...; kwargs...) = MoireWannier(moiresystem.frontend, args...; kwargs...)

"""
    MoireWannier(moiresystem::MoireSystem, lattice::MoireSuperlattice; nk=24, kwargs...)

Convenience constructor that auto-generates a `BrillouinZone` with `nk` k-points per dimension from `lattice` and delegates to the full constructor.
"""
@inline function MoireWannier(moiresystem::MoireSystem, lattice::MoireSuperlattice; nk=24, kwargs...)
    recipls = reciprocals(lattice)
    @assert recipls ≈ reciprocals(moiresystem.reciprocallattice) atol=1e-12 "MoireWannier error: mismatched reciprocals between input `moiresystem` and `lattice`."
    brillouinzone = BrillouinZone(reciprocals(lattice), nk)
    return MoireWannier(moiresystem, lattice, brillouinzone; kwargs...)
end

#=== triangular constructor (nband=1, U(1) gauge fix) ===#
"""
    MoireWannier(moiresystem::MoireSystem, lattice::MoireTriangular, brillouinzone::BrillouinZone; band::Int=dimension(moiresystem), tol::Real=atol)

Construct the Wannier function for a single band on a triangular lattice.

Gauge fixing: U(1) phase such that the bottom-layer component at r=0 (MM site) is real and positive.
"""
function MoireWannier(moiresystem::MoireSystem, lattice::MoireTriangular, brillouinzone::BrillouinZone; band::Int=dimension(moiresystem), tol::Real=atol)
    dim = dimension(moiresystem)
    @assert 1 <= band <= dim "MoireWannier error: band index $band out of range [1, $dim]."
    nlayer, nband, nk = 2, 1, length(brillouinzone)
    nG = dim ÷ nlayer
    # extract band energies and Bloch states
    energies = zeros(Float64, nband, nk)
    bloch = zeros(ComplexF64, nlayer, nG, nband, nk)
    for (ik, k) in enumerate(brillouinzone)
        eigensystem = eigen(moiresystem, k)
        energies[1, ik] = eigensystem.values[band]
        bloch[:, :, 1, ik] = reshape(eigensystem.vectors[:, band], nlayer, nG)
    end
    # U(1) gauge fix: ψ_k(r_MM) real positive where r_MM = origin (0,0)
    U = zeros(ComplexF64, nband, nband, nk)
    for (ik, k) in enumerate(brillouinzone)
        ψ_MM = zero(ComplexF64)
        for ig in 1:nG
            ψ_MM += bloch[1, ig, 1, ik]
        end
        @assert abs(ψ_MM) > tol "MoireWannier error: wavefunction at r=0 is zero at k=$k; U(1) gauge fixing failed."
        U[1, 1, ik] = conj(ψ_MM) / abs(ψ_MM)
    end
    aₘ = moiresystem.parameters.a₀ / (2sind(moiresystem.parameters.θ/2))
    reciprocallattice = getcontent(moiresystem, :reciprocallattice)
    return MoireWannier(aₘ, lattice, reciprocallattice, brillouinzone, energies, bloch, U)
end

# === honeycomb constructor (nband=2, SU(2) + U(1) gauge fix) ===#
"""
    MoireWannier(moiresystem::MoireSystem, lattice::MoireHoneycomb, brillouinzone::BrillouinZone; bands::UnitRange{Int}, tol::Real=atol)

Construct Wannier functions for a 2-band subspace on a honeycomb lattice.

Steps:
1. Extract raw Bloch states for the 2-band subspace
2. SU(2) rotation: maximize layer polarization via diagonalizing layer projectors
3. U(1) gauge fix: ψ₁(r_XM) real positive (Wannier 1 at XM, bottom-layer), ψ₂(r_MX) real positive (Wannier 2 at MX, top-layer)
"""
function MoireWannier(moiresystem::MoireSystem, lattice::MoireHoneycomb, brillouinzone::BrillouinZone; bands::UnitRange{Int}, tol::Real=atol)
    dim, nband = dimension(moiresystem), length(bands)
    @assert nband == 2 "MoireWannier error: honeycomb requires exactly 2 bands, got $nband."
    @assert all(band -> 1 <= band <= dim, bands) "MoireWannier error: band indices out of range [1, $dim]."
    nlayer, nk = 2, length(brillouinzone)
    nG = dim ÷ nlayer
    # extract raw Bloch states and energies
    energies = zeros(Float64, nband, nk)
    bloch = zeros(ComplexF64, nlayer, nG, nband, nk)
    for (ik, k) in enumerate(brillouinzone)
        eigensystem = eigen(moiresystem, k)
        for (ib, band) in enumerate(bands)
            energies[ib, ik] = eigensystem.values[band]
            bloch[:, :, ib, ik] .= reshape(eigensystem.vectors[:, band], nlayer, nG)
        end
    end
    # SU(2) rotation to maximize layer polarization:
    # Diagonalize top-layer projector P at each k, where P[μ, ν] = ⟨ψ_μ|[0 0; 1 0]|ψ_ν⟩ = Σ_G (ψ_μ(G))' * [0 0; 1 0] * ψ_ν(G)
    Ũ = zeros(ComplexF64, nband, nband, nk)
    for ik in 1:nk
        P₁₁, P₂₂, P₁₂ = zero(ComplexF64), zero(ComplexF64), zero(ComplexF64)
        for ig in 1:nG
            P₁₁ += conj(bloch[2, ig, 1, ik]) * bloch[2, ig, 1, ik]
            P₁₂ += conj(bloch[2, ig, 1, ik]) * bloch[2, ig, 2, ik]
            P₂₂ += conj(bloch[2, ig, 2, ik]) * bloch[2, ig, 2, ik]
        end
        Ũ[:, :, ik] = eigvecs(Hermitian(SMatrix{2, 2}(P₁₁, conj(P₁₂), P₁₂, P₂₂)))
    end
    # U(1) gauge fix:
    # W₁ (bottom-layer-dominant) → real positive at XM (lattice[1])
    # W₂ (top-layer-dominant) → real positive at MX (lattice[2])
    r_XM, r_MX= lattice[1], lattice[2]
    U = zeros(ComplexF64, nband, nband, nk)
    for (ik, k) in enumerate(brillouinzone)
        ψ₁ = zero(ComplexF64)
        ψ₂ = zero(ComplexF64)
        for (ig, G) in enumerate(moiresystem.reciprocallattice)
            phase_XM = cis(dot(k + G, r_XM))
            phase_MX = cis(dot(k + G, r_MX))
            # W₁: bottom-layer (il=1) component at XM
            # W₂: top-layer (il=2) component at MX
            for ν in 1:nband
                ψ₁ += bloch[1, ig, ν, ik] * Ũ[ν, 1, ik] * phase_XM
                ψ₂ += bloch[2, ig, ν, ik] * Ũ[ν, 2, ik] * phase_MX
            end
        end
        @assert abs(ψ₁) > tol "MoireWannier error: ψ₁(r_XM) is zero at k=$k; U(1) gauge fixing failed."
        @assert abs(ψ₂) > tol "MoireWannier error: ψ₂(r_MX) is zero at k=$k; U(1) gauge fixing failed."
        φ₁ = conj(ψ₁) / abs(ψ₁)
        φ₂ = conj(ψ₂) / abs(ψ₂)
        # U(k) = Ũ(k) × diag(e^{-iφ₁}, e^{-iφ₂})
        for ib in 1:nband
            U[ib, 1, ik] = Ũ[ib, 1, ik] * φ₁
            U[ib, 2, ik] = Ũ[ib, 2, ik] * φ₂
        end
    end
    aₘ = moiresystem.parameters.a₀ / (2sind(moiresystem.parameters.θ/2))
    reciprocallattice = getcontent(moiresystem, :reciprocallattice)
    return MoireWannier(aₘ, lattice, reciprocallattice, brillouinzone, energies, bloch, U)
end

"""
    HoppingIntegral{W<:MoireWannier}

Hopping amplitude calculator for a Wannier function.

Fields:
- `wannier::W` — reference to the MoireWannier

Callable as `(hopping::HoppingIntegral)(R::AbstractVector{<:Number}) -> Matrix{ComplexF64}`:

```math
t_{mn}(R) = (1/N) Σ_k exp(-ik·R) [U(k) diag(ε(k)) U†(k)]_{mn}
```
"""
struct HoppingIntegral{W<:MoireWannier}
    wannier::W
end
function (hopping::HoppingIntegral)(R::AbstractVector{<:Number})
    nband, nk = size(hopping.wannier.energies)
    result = zeros(ComplexF64, nband, nband)
    for (ik, k) in enumerate(hopping.wannier.brillouinzone)
        phase = exp(-1im * dot(k, R))
        U = hopping.wannier.U[:, :, ik]
        ε = hopping.wannier.energies[:, ik]
        for m in 1:nband, n in 1:nband
            acc = zero(ComplexF64)
            for i in 1:nband
                acc += conj(U[i, m]) * ε[i] * U[i, n]
            end
            result[m, n] += phase * acc
        end
    end
    return broadcast!(/, result, result, nk)
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
    CoulombIntegral{W<:MoireWannier}

Precomputed Coulomb form factor for a Wannier function.

Fields:
- `wannier::W` — reference to the Wannier function
- `qs::Vector{SVector{2,Float64}}` — unique q-vectors from the pairwise convolution
- `formfactor::Vector{Matrix{ComplexF64}}` — M(q)†M(q) / Nₖ² matrices at each q, (nband×nband)
"""
struct CoulombIntegral{W<:MoireWannier}
    wannier::W
    qs::Vector{SVector{2,Float64}}
    formfactor::Vector{Matrix{ComplexF64}}
end

"""
    CoulombIntegral(wannier::MoireWannier)

Construct by computing the form factor M(q) from Bloch coefficients via pairwise convolution.

Algorithm: compute gauge-transformed Bloch coefficients c_n(p) for each Wannier function n at
extended momenta p=k+G, then for each pair (n,m) accumulate M_{m,n}(q) = Σ_p dot(c_n(p+q), c_m(p)).
The q-mesh emerges naturally from the G/G' truncation.
"""
function CoulombIntegral(wannier::MoireWannier)
    nlayer, nG, nband, nk = size(wannier.bloch)
    b₁, b₂ = wannier.reciprocallattice.translations
    N₁, N₂ = periods(wannier.brillouinzone)
    # integer coordinates for k-points
    ks = Vector{Tuple{Int, Int}}(undef, nk)
    for (ik, k) in enumerate(wannier.brillouinzone)
        f₁, f₂ = decompose(k, b₁, b₂)
        ks[ik] = (round(Int, f₁*N₁), round(Int, f₂*N₂))
    end
    # integer coordinates for G-vectors
    Gs = Vector{Tuple{Int, Int}}(undef, nG)
    for (ig, G) in enumerate(wannier.reciprocallattice)
        g₁, g₂ = decompose(G, b₁, b₂)
        Gs[ig] = (round(Int, g₁), round(Int, g₂))
    end
    # extended momentum integer coordinates: p = k + G
    ps = Matrix{Tuple{Int,Int}}(undef, nk, nG)
    for (ik, (k₁, k₂)) in enumerate(ks), (ig, (g₁, g₂)) in enumerate(Gs)
        ps[ik, ig] = (k₁ + N₁*g₁, k₂ + N₂*g₂)
    end
    # gauge-transformed coefficients: coeff[il, ig, iw, ik], same shape as bloch
    coeff = zeros(ComplexF64, nlayer, nG, nband, nk)
    for ik in 1:nk, iw in 1:nband, ib in 1:nband, ig in 1:nG, il in 1:nlayer
        coeff[il, ig, iw, ik] += wannier.bloch[il, ig, ib, ik] * wannier.U[ib, iw, ik]
    end
    # pairwise convolution: M_n(q) = Σ_p dot(coeff_n(p+q), coeff_n(p))
    Ms = Dict{Tuple{Int,Int}, Vector{ComplexF64}}()
    for ik₁ in 1:nk, ig₁ in 1:nG
        p₁ = ps[ik₁, ig₁]
        for ik₂ in 1:nk, ig₂ in 1:nG
            p₂ = ps[ik₂, ig₂]
            M = get!(Ms, (p₂[1]-p₁[1], p₂[2]-p₁[2])) do
                zeros(ComplexF64, nband)
            end
            for n in 1:nband, il in 1:nlayer
                M[n] += conj(coeff[il, ig₂, n, ik₂]) * coeff[il, ig₁, n, ik₁]
            end
        end
    end
    # formfactor: M*_m(q) M_n(q) / nk² at each q (outer product of M(q) with its conjugate)
    qs = Vector{SVector{2, Float64}}(undef, length(Ms))
    formfactor = Vector{Matrix{ComplexF64}}(undef, length(Ms))
    for (i, ((q₁, q₂), M)) in enumerate(Ms)
        qs[i] = (q₁/N₁)*b₁ + (q₂/N₂)*b₂
        m′m = zeros(ComplexF64, nband, nband)
        for m in 1:nband, n in 1:nband
            m′m[m, n] = conj(M[m]) * M[n] / nk^2
        end
        formfactor[i] = m′m
    end
    return CoulombIntegral(wannier, qs, formfactor)
end

"""
    (coulomb::CoulombIntegral)(R::AbstractVector{<:Number}, potential=BareCoulomb(1.0))

Compute U(R) = (1/(N_k Ω)) Σ_q V(|q|, aₘ) |M(q)|² e^{iq·R} in meV.

`potential` is a callable `potential(q::Real, aₘ::Real) -> Real`, e.g. `BareCoulomb(ϵ)`, `ImageCoulomb(ϵ, d)`, `TanhCoulomb(ϵ, d)`, or a user-defined function. Defaults to `BareCoulomb(1.0)` (unscreened, ε=1).

Returns an nband×nband `Matrix{Float64}` Coulomb interaction matrix.
"""
function (coulomb::CoulombIntegral)(R::AbstractVector{<:Number}, potential=BareCoulomb(1.0))
    nk = length(coulomb.wannier.brillouinzone)
    Ω = volume(coulomb.wannier.lattice.vectors)
    nband = size(coulomb.wannier.energies, 1)
    result = zeros(Float64, nband, nband)
    for (q, m′m) in zip(coulomb.qs, coulomb.formfactor)
        Vq = potential(norm(q), coulomb.wannier.aₘ) * cos(dot(q, R))
        for i in eachindex(result, m′m)
            result[i] += Vq * real(m′m[i])
        end
    end
    return broadcast!(/, result, result, nk*Ω)
end

#=== Hopping and Coulomb terms from Wannier integrals ===#
"""
    SublatticeAmplitude{G<:PointGroup, N, D<:Number} <: Function

Spin-independent hopping amplitude. Matches a bond to the registered set of reference bonds via [`sign`](@ref) and returns `1` (match) or `0` (no match).
"""
struct SublatticeAmplitude{G<:PointGroup, N, D<:Number} <: Function
    refs::NTuple{N, Bond{Int, Point{2, D}, SVector{2, Point{2, D}}}}
    nsublattice::Int
    function SublatticeAmplitude{G}(refs::NTuple{N, Bond{Int, Point{2, D}, SVector{2, Point{2, D}}}}, nsublattice::Int) where {G<:PointGroup, N, D<:Number}
        return new{G, N, D}(refs, nsublattice)
    end
end
function (amp::SublatticeAmplitude{G})(bond::Bond) where {G<:PointGroup}
    for ref in amp.refs
        !iszero(sign(G, ref, bond, amp.nsublattice)) && return 1
    end
    return 0
end

"""
    SpinOrbitalCouplingAmplitude{G<:PointGroup, N, D<:Number} <: Function

Spin-orbital-coupling hopping amplitude under point group `G`. Matches a bond to the registered reference bonds via [`sign`](@ref), which yields a spatial parity ``r = \\pm 1``. The final SOC factor is ``-i \\cdot s \\cdot r``, where ``s`` is the relative sign of the SOC coefficient ``\\lambda`` between bonds in the same ``(t, |\\lambda|)`` group (precomputed by [`terms`](@ref) and stored in `signs`).
"""
struct SpinOrbitalCouplingAmplitude{G<:PointGroup, N, D<:Number} <: Function
    signs::NTuple{N, Int}
    refs::NTuple{N, Bond{Int, Point{2, D}, SVector{2, Point{2, D}}}}
    nsublattice::Int
    function SpinOrbitalCouplingAmplitude{G}(signs::NTuple{N, Int}, refs::NTuple{N, Bond{Int, Point{2, D}, SVector{2, Point{2, D}}}}, nsublattice::Int) where {G<:PointGroup, N, D<:Number}
        return new{G, N, D}(signs, refs, nsublattice)
    end
end
function (amp::SpinOrbitalCouplingAmplitude{G})(bond::Bond) where G<:PointGroup
    for (s, ref) in zip(amp.signs, amp.refs)
        result = sign(G, ref, bond, amp.nsublattice)
        iszero(result) && continue
        return -1im * s * result
    end
    return 0im
end

"""
    OnsiteAmplitude <: Function

Onsite amplitude. Matches a self-bond (bondkind=0) to a specific sublattice and returns `1` (match) or `0` (no match).
"""
struct OnsiteAmplitude <: Function
    sublattice::Int
    nsublattice::Int
end
function (amp::OnsiteAmplitude)(bond::Bond)
    bond.kind == 0 || return 0
    return Sublattice(amp.nsublattice)(bond[1].site) == amp.sublattice ? 1 : 0
end

"""
    terms(hopping::HoppingIntegral; order::Int=truncation(hopping.wannier.lattice), ismodulatable::Bool=true, tol=atol) -> Tuple{Vararg{Term}}

Generate spin-independent and spin-orbital-coupling `Hopping` terms from a [`HoppingIntegral`](@ref).

## Algorithm

1. **Per-bond extraction** — For each neighbor order ``k``, the coefficient ``t + iλ = hopping(R)[i,j]`` is extracted for each symmetry-inequivalent bond independently. Each bond stored in [`MoireNeighbors`](@ref) is already inequivalent under the prescribed point group, so no within-group consistency assertion is needed.

2. **Grouping** — Entries are grouped by ``(t, |λ|)``: ``t`` must match exactly (including sign), ``λ`` by magnitude. Relative signs between grouped ``λ`` values are precomputed.

3. **Term construction** — Each group yields two terms:
   - `t` → `Hopping` with [`SublatticeAmplitude`](@ref) for sublattice matching.
   - `λ` → `Hopping` with ``σᶻ`` coupling and [`SpinOrbitalCouplingAmplitude`](@ref) for spatial + sublattice matching.

4. **Onsite** — Diagonal ``hopping(0)[i, i]`` → `Onsite` chemical potentials.

## Naming

- `nband = 1` and `nband > 1` with single group per shell: `t₁`, `λ₁`, `t₂`, `λ₂`, …
- `nband > 1` with multiple groups per shell: all pairs listed, e.g. `t₁₍₁₋₂₊₂₋₃₎`.
"""
function terms(hopping::HoppingIntegral; order::Int=truncation(hopping.wannier.lattice), ismodulatable::Bool=true, tol=atol)
    lattice = hopping.wannier.lattice
    B = eltype(bonds(lattice.neighbors))
    DataEntry = @NamedTuple{t::Float64, λ::Float64, bond::B}
    GroupEntry = @NamedTuple{t::Float64, λ::Float64, signs::Vector{Int}, bonds::Vector{B}}
    # P1: extract per-bond coefficients (each bond is symmetry-inequivalent by construction)
    shells = Vector{Vector{GroupEntry}}(undef, order)
    sl = Sublattice(nsublattice(lattice.neighbors))
    for k in 1:order
        data = DataEntry[]
        for bond in bonds(lattice.neighbors, k)
            i, j = sl(bond[1].site), sl(bond[2].site)
            coeff = hopping(icoordinate(bond))[i, j]
            t, λ = real(coeff), imag(coeff)
            isapprox(t, 0; atol=tol) && isapprox(λ, 0; atol=tol) && continue
            push!(data, (t=t, λ=λ, bond=bond))
        end
        # P2: group by (t, |λ|)
        groups = GroupEntry[]
        for (t, λ, bond) in data
            found = false
            for group in groups
                isapprox(t, group.t; atol=tol) && isapprox(abs(λ), abs(group.λ); atol=tol) || continue
                push!(group.signs, sign(λ)*sign(group.λ))
                push!(group.bonds, bond)
                found = true
                break
            end
            found && continue
            push!(groups, (t=t, λ=λ, signs=[abs(sign(λ))], bonds=[bond]))
        end
        shells[k] = groups
    end
    # P3: generate Hopping terms
    G = typeof(PointGroup(lattice))
    nband = size(hopping.wannier.energies, 1)
    hoppings = map(enumerate(shells)) do (k, groups)
        result = Term[]
        for group in groups
            suffix = join('₀'+d for d in reverse(digits(k)))
            refs = Tuple(group.bonds)
            if length(groups) > 1
                pairs = String[]
                for bond in refs
                    i = join('₀'+d for d in reverse(digits(sl(bond[1].site))))
                    j = join('₀'+d for d in reverse(digits(sl(bond[2].site))))
                    push!(pairs, string(i, '₋', j))
                end
                suffix *= string('₍', join(unique!(pairs), '₊'), '₎')
            end
            push!(result, Hopping(
                Symbol("t", suffix), group.t, k;
                amplitude=SublatticeAmplitude{G}(refs, sl.nsublattice),
                ismodulatable=ismodulatable
            ))
            push!(result, Hopping(
                Symbol("λ", suffix), group.λ, k, 𝕔⁺𝕔(:, :, σᶻ);
                amplitude=SpinOrbitalCouplingAmplitude{G}(Tuple(group.signs), refs, sl.nsublattice),
                ismodulatable=ismodulatable
            ))
        end
        return Tuple(result)
    end
    # P4: onsite — group by μ value
    h₀ = hopping(SVector(0.0, 0.0))
    μs = [real(h₀[i, i]) for i in 1:nband]
    onsites = Term[]
    if all(μ->isapprox(μ, first(μs); atol=tol), μs)
        push!(onsites, Onsite(:μ, first(μs); ismodulatable=ismodulatable))
    else
        for (i, μ) in enumerate(μs)
            name = Symbol("μ", join('₀'+d for d in reverse(digits(i))))
            push!(onsites, Onsite(name, μ; amplitude=OnsiteAmplitude(i, sl.nsublattice), ismodulatable=ismodulatable))
        end
    end
    return (concatenate(hoppings...)..., onsites...)
end

"""
    terms(coulomb::CoulombIntegral, potential=BareCoulomb(1.0); order::Int=truncation(coulomb.wannier.lattice), ismodulatable::Bool=true, tol=atol) -> Tuple{Vararg{Term}}

Generate Coulomb interaction terms from a [`CoulombIntegral`](@ref).

`potential` is a callable `potential(q::Real, aₘ::Real) -> Real`, e.g. `BareCoulomb(ϵ)`, `ImageCoulomb(ϵ, d)`, `TanhCoulomb(ϵ, d)`, or a user-defined function.

## Algorithm

1. **Per-bond extraction** — For each neighbor order ``k``, the coefficient ``V = coulomb(R, potential)[i, j]`` is extracted for each symmetry-inequivalent bond independently. Each bond stored in [`MoireNeighbors`](@ref) is already inequivalent under the prescribed point group, so no within-group consistency assertion is needed.

2. **Grouping** — Entries are grouped by ``V`` value across all pairs in the same shell.

3. **Term construction** — Each group yields a `Coulomb` term with [`SublatticeAmplitude`](@ref) for sublattice matching. Coulomb interaction is purely real, so there is no spin-orbital-coupling counterpart. Parameter names use ``V`` prefix (e.g., `V₁`, `V₂`).

4. **Onsite** — Diagonal ``coulomb(0, potential)[i, i]`` → `Hubbard` onsite repulsion with ``U`` prefix.

## Naming

- R > 0, `nband = 1` or single group per shell: `V₁`, `V₂`, …
- R > 0, `nband > 1` with multiple groups per shell: all pairs listed, e.g. `V₁₍₁₋₂₊₂₋₃₎`.
- R = 0 (onsite): `U` (all equal) or `U₁`, `U₂`, … (per-sublattice).
"""
function terms(coulomb::CoulombIntegral, potential=BareCoulomb(1.0); order::Int=truncation(coulomb.wannier.lattice), ismodulatable::Bool=true, tol=atol)
    lattice = coulomb.wannier.lattice
    B = eltype(bonds(lattice.neighbors))
    DataEntry = @NamedTuple{V::Float64, bond::B}
    GroupEntry = @NamedTuple{V::Float64, bonds::Vector{B}}
    # P1: extract per-bond coefficients (R > 0)
    shells = Vector{Vector{GroupEntry}}(undef, order)
    sl = Sublattice(nsublattice(lattice.neighbors))
    for k in 1:order
        data = DataEntry[]
        for bond in bonds(lattice.neighbors, k)
            i, j = sl(bond[1].site), sl(bond[2].site)
            V = coulomb(icoordinate(bond), potential)[i, j]
            isapprox(V, 0; atol=tol) && continue
            push!(data, (V=V, bond=bond))
        end
        # P2: group by V value
        groups = GroupEntry[]
        for (V, bond) in data
            found = false
            for group in groups
                isapprox(V, group.V; atol=tol) || continue
                push!(group.bonds, bond)
                found = true
                break
            end
            found && continue
            push!(groups, (V=V, bonds=[bond]))
        end
        shells[k] = groups
    end
    # P3: generate Coulomb terms (R > 0)
    G = typeof(PointGroup(lattice))
    nband = size(coulomb.wannier.energies, 1)
    coulombs = map(enumerate(shells)) do (k, groups)
        result = Term[]
        for group in groups
            suffix = join('₀'+d for d in reverse(digits(k)))
            refs = Tuple(group.bonds)
            if length(groups) > 1
                pairs = String[]
                for bond in refs
                    i = join('₀'+d for d in reverse(digits(sl(bond[1].site))))
                    j = join('₀'+d for d in reverse(digits(sl(bond[2].site))))
                    push!(pairs, string(i, '₋', j))
                end
                suffix *= string('₍', join(unique!(pairs), '₊'), '₎')
            end
            push!(result, Coulomb(
                Symbol("V", suffix), group.V, k;
                amplitude=SublatticeAmplitude{G}(refs, sl.nsublattice),
                ismodulatable=ismodulatable
            ))
        end
        return Tuple(result)
    end
    # P4: onsite (R=0) — Hubbard terms
    U₀ = coulomb(SVector(0.0, 0.0), potential)
    Us = [U₀[i, i] for i in 1:nband]
    hubbards = Term[]
    if all(U->isapprox(U, first(Us); atol=tol), Us)
        push!(hubbards, Hubbard(:U, first(Us); ismodulatable=ismodulatable))
    else
        for (i, U) in enumerate(Us)
            name = Symbol("U", join('₀'+d for d in reverse(digits(i))))
            push!(hubbards, Hubbard(name, U; amplitude=OnsiteAmplitude(i, sl.nsublattice), ismodulatable=ismodulatable))
        end
    end
    return (hubbards..., concatenate(coulombs...)...)
end
