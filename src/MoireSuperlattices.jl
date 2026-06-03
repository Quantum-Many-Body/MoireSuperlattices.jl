module MoireSuperlattices

using LinearAlgebra: Hermitian, dot, eigen, eigvals, eigvecs, norm
using Printf: @printf
using QuantumLattices: annihilation, atol, creation, σᶻ
using QuantumLattices: AbstractLattice, Bond, BrillouinZone, CategorizedGenerator, CompositeIndex, Coupling, Hilbert, Hopping, icoordinate, Index, InternalIndex, LaTeX, Neighbors, Onsite, OperatorGenerator, OperatorIndexToTuple, OperatorSum, Point, ReciprocalZone, SimpleInternal, Table, Term
using QuantumLattices: azimuth, azimuthd, bonds, concatenate, decompose, distance, latexformat, periods, rcoordinate, rotate, scalartype, str, update, volume, 𝕔⁺𝕔

using StaticArrays: SVector
using TightBindingApproximation: TBA, Fermionic, Quadratic, Quadraticization

import QuantumLattices: Algorithm, Lattice, Parameters, contentnames, dimension, getcontent, indextype, isdefinite, latexname, matrix, reciprocals, script, shape, statistics, update!

export BLTMD, C6, C₆, CommensurateBilayerHoneycomb, BareCoulomb, CoulombIntegral, HoppingIntegral, ImageCoulomb, MoireAmplitude, MoireHoneycomb, MoireNeighbors, MoireReciprocalLattice, MoireSpace, MoireSpinor, MoireSuperlattice, MoireSystem, MoireTriangular, MoireTriangularReciprocal, MoireWannier, PointGroup, RealZone, TanhCoulomb, bltmd!, bltmdmap, terms, truncation, vectors

include("lattices.jl")
include("systems.jl")
include("analysis.jl")

# runtime initialization
function __init__()
    latexformat(MoireSpinor, LaTeX{(:nambu,), (:layer, :spin)}('c'))
    latexformat(Index{<:MoireSpinor}, LaTeX{(:nambu,), (:site, :layer, :spin)}('c'))
    latexformat(CompositeIndex{<:Index{<:MoireSpinor}}, LaTeX{(:nambu,), (:site, :layer, :spin)}('c'))
    nothing
end

end