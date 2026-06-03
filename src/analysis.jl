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

#=== triangular constructor (nband=1, U(1) gauge fix) ===#
"""
    MoireWannier(moiresystem::MoireSystem, lattice::MoireTriangular, brillouinzone::BrillouinZone; band::Int=dimension(moiresystem))

Construct the Wannier function for a single band on a triangular lattice.

Gauge fixing: U(1) phase such that the bottom-layer component at r=0 (MM site) is real and positive.
"""
function MoireWannier(moiresystem::MoireSystem, lattice::MoireTriangular, brillouinzone::BrillouinZone; band::Int=dimension(moiresystem))
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
        @assert abs(ψ_MM) > atol "MoireWannier error: wavefunction at r=0 is zero at k=$k; gauge fixing failed."
        U[1, 1, ik] = conj(ψ_MM) / abs(ψ_MM)
    end
    aₘ = moiresystem.parameters.a₀ / (2sind(moiresystem.parameters.θ/2))
    reciprocallattice = getcontent(moiresystem, :reciprocallattice)
    return MoireWannier(aₘ, lattice, reciprocallattice, brillouinzone, energies, bloch, U)
end

# === honeycomb constructor (nband=2, SU(2) + U(1) gauge fix) ===#
# """
#     MoireWannier(moiresystem::MoireSystem, lattice::MoireHoneycomb, brillouinzone::BrillouinZone; bands::UnitRange{Int})

# Construct Wannier functions for a 2-band subspace on a honeycomb lattice.

# Steps:
# 1. Extract raw Bloch states for the 2-band subspace
# 2. SU(2) rotation: maximize layer polarization via diagonalizing layer projectors
# 3. U(1) gauge fix: ψ̃₁(r_MX) real positive, ψ̃₂(r_XM) real positive
# """
# function MoireWannier(moiresystem::MoireSystem, lattice::MoireHoneycomb, brillouinzone::BrillouinZone; bands::UnitRange{Int})
#     dim = dimension(moiresystem)
#     nband = length(bands)
#     @assert nband == 2 "MoireWannier error: honeycomb requires exactly 2 bands, got $nband."
#     @assert all(b -> 1 <= b <= dim, bands) "MoireWannier error: band indices out of range [1, $dim]."
#     nk = length(brillouinzone)
#     nlayer = 2
#     nG = dim ÷ nlayer
#     # extract raw Bloch states and energies
#     bloch = zeros(ComplexF64, nlayer, nG, nband, nk)
#     energies = zeros(Float64, nband, nk)
#     band_indices = collect(bands)
#     for (ik, k) in enumerate(brillouinzone)
#         eigensystem = eigen(moiresystem, k)
#         for (ib, b) in enumerate(band_indices)
#             energies[ib, ik] = eigensystem.values[b]
#             bloch[:, :, ib, ik] .= reshape(eigensystem.vectors[:, b], nlayer, nG)
#         end
#     end
#     # SU(2) rotation to maximize layer polarization
#     U_tilde = zeros(ComplexF64, nband, nband, nk)
#     _su2_layer_polarization!(U_tilde, bloch, nG, nlayer, nk)
#     # U(1) gauge fix:
#     # ψ̃₁ at r_MX (MX position, sublattice 1) → real positive
#     # ψ̃₂ at r_XM (XM position, sublattice 2) → real positive
#     # MX = coordinates[:, 1], XM = coordinates[:, 2] in the honeycomb lattice
#     r_MX = SVector(lattice.coordinates[1, 1], lattice.coordinates[2, 1])
#     r_XM = SVector(lattice.coordinates[1, 2], lattice.coordinates[2, 2])
#     U = zeros(ComplexF64, nband, nband, nk)
#     for ik in 1:nk
#         k = SVector(brillouinzone[ik][1], brillouinzone[ik][2])
#         psi1_MX = zero(ComplexF64)
#         psi2_XM = zero(ComplexF64)
#         for ig in 1:nG
#             Gvec = SVector(moiresystem.reciprocallattice[ig][1], moiresystem.reciprocallattice[ig][2])
#             phase_MX = cis(dot(k + Gvec, r_MX))
#             phase_XM = cis(dot(k + Gvec, r_XM))
#             for il in 1:nlayer
#                 for ν in 1:nband
#                     psi1_MX += bloch[il, ig, ν, ik] * U_tilde[ν, 1, ik] * phase_MX
#                     psi2_XM += bloch[il, ig, ν, ik] * U_tilde[ν, 2, ik] * phase_XM
#                 end
#             end
#         end
#         @assert abs(psi1_MX) > atol "MoireWannier error: ψ̃₁(r_MX) is zero at k=$k; gauge fixing failed."
#         @assert abs(psi2_XM) > atol "MoireWannier error: ψ̃₂(r_XM) is zero at k=$k; gauge fixing failed."
#         phi1 = conj(psi1_MX) / abs(psi1_MX)
#         phi2 = conj(psi2_XM) / abs(psi2_XM)
#         # U(k) = Ũ(k) × diag(e^{iφ₁k}, e^{iφ₂k})
#         # U[:, n, ik] = U_tilde[:, n, ik] * exp(-i*phi_n)  for each column n
#         U[:, 1, ik] .= U_tilde[:, 1, ik] .* phi1
#         U[:, 2, ik] .= U_tilde[:, 2, ik] .* phi2
#     end
#     aₘ = moiresystem.parameters.a₀ / (2sind(moiresystem.parameters.θ/2))
#     reciprocallattice = getcontent(moiresystem, :reciprocallattice)
#     return MoireWannier(aₘ, lattice, reciprocallattice, brillouinzone, energies, bloch, U)
# end

# """
#     _su2_layer_polarization!(U_tilde, bloch, nG, nlayer, nk)

# Compute the SU(2) rotation Ũ(k) at each k that maximizes layer polarization:
# - Column 1: maximize bottom-layer projection ⟨P_b⟩
# - Column 2: maximize top-layer projection  ⟨P_t⟩

# P_b = diag(1, 0), P_t = diag(0, 1) acting on the layer index.
# """
# function _su2_layer_polarization!(U_tilde::Array{ComplexF64,3}, bloch::Array{ComplexF64,4}, nG::Int, nlayer::Int, nk::Int)
#     @assert nlayer == 2 "SU(2) layer polarization requires 2 layers."
#     for ik in 1:nk
#         Pb = zeros(ComplexF64, 2, 2)
#         for ig in 1:nG
#             # bottom-layer component (layer=1) at G-vector ig
#             for μ in 1:2, ν in 1:2
#                 Pb[μ, ν] += bloch[1, ig, μ, ik] * conj(bloch[1, ig, ν, ik])
#             end
#         end
#         # diagonalize Pb (2×2 Hermitian)
#         vals, vecs = eigen(Hermitian(Pb))
#         # sort: eigenvector for max eigenvalue → column 1 (maximizes bottom-layer weight)
#         #        eigenvector for min eigenvalue → column 2 (maximizes top-layer weight, since Pb + Pt ≈ I)
#         perm = sortperm(vals; rev=true)
#         U_tilde[:, :, ik] .= vecs[:, perm]
#     end
#     return U_tilde
# end

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
                acc += U[m, i] * ε[i] * conj(U[n, i])
            end
            result[m, n] += phase * acc
        end
    end
    return broadcast!(/, result, result, nk)
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
    (c::CoulombIntegral)(R::AbstractVector{<:Number}, V=BareCoulomb(1.0))

Compute U(R) = (1/(N_k Ω)) Σ_q V(|q|, aₘ) |M(q)|² e^{iq·R} in meV.

`V` is a callable `V(q::Real, aₘ::Real) -> Real`, e.g. `BareCoulomb(ϵ)`, `ImageCoulomb(ϵ, d)`, `TanhCoulomb(ϵ, d)`, or a user-defined function. Defaults to `BareCoulomb(1.0)` (unscreened, ε=1).

Returns an nband×nband `Matrix{Float64}` Coulomb interaction matrix.
"""
function (c::CoulombIntegral)(R::AbstractVector{<:Number}, V=BareCoulomb(1.0))
    nk = length(c.wannier.brillouinzone)
    Ω = volume(c.wannier.lattice.vectors)
    nband = size(c.wannier.energies, 1)
    result = zeros(Float64, nband, nband)
    for (q, m′m) in zip(c.qs, c.formfactor)
        Vq = V(norm(q), c.wannier.aₘ) * cos(dot(q, R))
        for i in eachindex(result, m′m)
            result[i] += Vq * real(m′m[i])
        end
    end
    return broadcast!(/, result, result, nk*Ω)
end

#=== Hopping and Coulomb terms from Wannier integrals ===#

"""
    MoireAmplitude{N, G<:PointGroup}

SOC hopping amplitude for Moiré superlattices under point group `G`.

Fields:
- `signs::NTuple{N, Int}` — relative polarities between symmetry-inequivalent stars within a shell
- `θs::NTuple{N, Float64}` — reference azimuthal angles of each star (°)
- `ℓ::Int` — angular momentum channel (3 for Moiré C₆ systems)

When called with a [`Bond`](@ref), matches its azimuth to the correct star and returns
`-1im * sign * cosd(ℓ * Δθ)` where Δθ is the angular deviation from the reference.
"""
struct MoireAmplitude{G<:PointGroup, N} <: Function
    signs::NTuple{N, Int}
    θs::NTuple{N, Float64}
    ℓ::Int
    function MoireAmplitude{G}(λs::AbstractVector{<:Real}, shell::AbstractVector{<:Bond}; ℓ::Int=3, atol::Real=atol) where G<:PointGroup
        N = length(λs)
        signs = ntuple(b -> isapprox(λs[b], 0; atol=atol) ? 1 : round(Int, λs[1]/λs[b]), N)
        θs = ntuple(i -> azimuthd(rcoordinate(shell[i])), N)
        return new{G, N}(signs, θs, ℓ)
    end
end

function (amp::MoireAmplitude{G})(bond::Bond) where G<:PointGroup
    αd = rad2deg(angle(G))
    θ = azimuthd(rcoordinate(bond))
    for (sign, θ₀) in zip(amp.signs, amp.θs)
        Δ = (θ - θ₀) / αd
        isapprox(round(Int, Δ), Δ; atol=atol) || continue
        return -1im * sign * cosd(amp.ℓ * (θ - θ₀))
    end
    error("amplitude error: mismatched bond.")
end

"""
    terms(h::HoppingIntegral; order::Int=truncation(h.wannier.lattice), ismodulatable::Bool=true, tol=atol) -> NTuple{...}, Term}

Generate hopping terms from a [`HoppingIntegral`](@ref) for any Moiré superlattice.

For each neighbor shell, collects the nband×nband hopping matrices for symmetry-inequivalent
bonds, decomposes each matrix element ``(i,j)`` into spin-independent (real part) and SOC
(imaginary part) components, and generates corresponding [`Hopping`](@ref) terms.

For triangular lattices (nband=1), each shell produces 2 terms (spin-independent + SOC).
For multi-band lattices (nband>1), terms are generated per sublattice pair with subscripted
names (e.g., `t₁₁₂` for shell 1, sublattice pair (1,2)).
"""
function terms(h::HoppingIntegral; order::Int=truncation(h.wannier.lattice), ismodulatable::Bool=true, tol=atol)
    lattice = h.wannier.lattice
    nband = size(h.wannier.energies, 1)
    G = typeof(PointGroup(lattice))
    # collect hopping matrices and validate symmetry per shell
    shells = Vector{Tuple{Int, Vector{Bond}, Vector{Any}}}(undef, order)
    for k in 1:order
        shell = lattice.neighbors[k]
        nstar = length(shell)
        hmats = [h(icoordinate(bond)) for bond in shell]
        elements = Vector{Any}()
        for i in 1:nband, j in 1:nband
            ts = [real(hmats[b][i,j]) for b in 1:nstar]
            λs = [imag(hmats[b][i,j]) for b in 1:nstar]
            # skip zero element
            all(t -> isapprox(t, 0; atol=tol), ts) && all(λ -> isapprox(λ, 0; atol=tol), λs) && continue
            # assert |t| and |λ| are the same magnitude across all stars
            tref, λref = abs(ts[1]), abs(λs[1])
            @assert all(b -> isapprox(abs(ts[b]), tref; atol=tol), 1:nstar) "terms error: |t| mismatch for ($i,$j) in shell $k."
            @assert all(b -> isapprox(abs(λs[b]), λref; atol=tol), 1:nstar) "terms error: |λ| mismatch for ($i,$j) in shell $k."
            # name suffix
            suffix_str = join('₀'+d for d in digits(k))
            if nband > 1
                suffix_str *= string(Char(0x2080+i), Char(0x2080+j))
            end
            push!(elements, (ts=ts, λs=λs, suffix=suffix_str))
        end
        shells[k] = (k, shell, elements)
    end
    # generate Hopping terms from validated data
    hoppings = map(shells) do (k, shell, elements)
        terms_list = Term[]
        for elem in elements
            push!(terms_list, Hopping(Symbol("t", elem.suffix), elem.ts[1], k; ismodulatable=ismodulatable))
            push!(terms_list, Hopping(
                Symbol("λ", elem.suffix), elem.λs[1], k, 𝕔⁺𝕔(:, :, σᶻ);
                amplitude=MoireAmplitude{G}(elem.λs, shell; ℓ=3, atol=tol),
                ismodulatable=ismodulatable
            ))
        end
        return (terms_list...,)
    end
    # onsite chemical potentials
    h0 = h(SVector(0.0, 0.0))
    μ_terms = Term[]
    for i in 1:nband
        isapprox(real(h0[i,i]), 0; atol=tol) && continue
        name = nband == 1 ? :μ : Symbol("μ", Char(0x2080+i))
        push!(μ_terms, Onsite(name, real(h0[i,i])))
    end
    return (concatenate(hoppings...)..., μ_terms...)
end

"""
    terms(c::CoulombIntegral; order::Int=truncation(c.wannier.lattice), tol=atol) -> NTuple{...}, Term}

Generate Coulomb interaction terms from a CoulombIntegral.

Not yet implemented.
"""
function terms(c::CoulombIntegral; order::Int=truncation(c.wannier.lattice), tol=atol)
    error("Coulomb terms not yet implemented.")
end
