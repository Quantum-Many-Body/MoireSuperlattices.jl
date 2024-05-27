var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = MoireSuperlattices","category":"page"},{"location":"#MoireSuperlattices","page":"Home","title":"MoireSuperlattices","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for MoireSuperlattices.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [MoireSuperlattices]","category":"page"},{"location":"#MoireSuperlattices.BLTMD","page":"Home","title":"MoireSuperlattices.BLTMD","text":"BLTMD{L<:MoireReciprocalLattice, D<:Function, H<:OperatorGenerator} <: MoireSystem{H}\n\nTwisted transition metal dichalcogenide homobilayers.\n\n\n\n\n\n","category":"type"},{"location":"#MoireSuperlattices.BLTMD-NTuple{8, Number}","page":"Home","title":"MoireSuperlattices.BLTMD","text":"BLTMD(a₀::Number, m::Number, θ::Number, Vᶻ::Number, μ::Number, V::Number, ψ::Number, w::Number; truncation::Int=4)\n\nContinuum model of twisted transition metal dichalcogenide homobilayers.\n\nHere, the parameters are as follows:\n\na₀: monolayer lattice constant (Å)\nm: effective mass of the conduction band (mₑ)\nθ: twist angle (°)\nVᶻ: perpendicular displacement field (meV)\nμ: chemical potential (meV)\nV: amplitude of Moire potential (meV)\nψ: phase of Moire potential (°)\nw: interlayer hopping amplitude (meV)\n\n\n\n\n\n","category":"method"},{"location":"#MoireSuperlattices.CommensurateBilayerHoneycomb","page":"Home","title":"MoireSuperlattices.CommensurateBilayerHoneycomb","text":"CommensurateBilayerHoneycomb\n\nCommensurate Moire superlattice composed of two layers of honeycomb lattices.\n\n\n\n\n\n","category":"type"},{"location":"#MoireSuperlattices.MoireEmergentSuperLattice","page":"Home","title":"MoireSuperlattices.MoireEmergentSuperLattice","text":"MoireEmergentSuperLattice{D<:Number} <: AbstractLattice{2, D, 2}\n\nAbstract type of the emergent superlattices in Moire systems.\n\n\n\n\n\n","category":"type"},{"location":"#MoireSuperlattices.MoireReciprocalLattice","page":"Home","title":"MoireSuperlattices.MoireReciprocalLattice","text":"MoireReciprocalLattice{T<:Number} <: AbstractLattice{2, T, 0}\n\nMoire reciprocal lattice with truncation.\n\n\n\n\n\n","category":"type"},{"location":"#MoireSuperlattices.MoireSpace","page":"Home","title":"MoireSuperlattices.MoireSpace","text":"MoireSpace <: SimpleInternal{MoireSpinor{Int, Int, Int, Rational{Int}, Int}}\n\nThe internal degrees of freedom of Moire systems.\n\n\n\n\n\n","category":"type"},{"location":"#MoireSuperlattices.MoireSpinor","page":"Home","title":"MoireSuperlattices.MoireSpinor","text":"MoireSpinor{V<:Union{Int, Colon}, L<:Union{Int, Colon}, S<:Union{Int, Colon}, P<:Union{Rational{Int}, Colon}, N<:Union{Int, Colon}} <: SimpleIID\n\nThe index of the internal degrees of freedom of Moire systems.\n\n\n\n\n\n","category":"type"},{"location":"#MoireSuperlattices.MoireSystem","page":"Home","title":"MoireSuperlattices.MoireSystem","text":"MoireSystem{H<:OperatorGenerator, Hₘ<:Image} <: AbstractTBA{Fermionic{:TBA}, Hₘ, Nothing}\n\nThe continuum model of Moire systems.\n\n\n\n\n\n","category":"type"},{"location":"#MoireSuperlattices.MoireTriangular","page":"Home","title":"MoireSuperlattices.MoireTriangular","text":"MoireTriangular{N, D<:Number} <: MoireEmergentSuperLattice{D}\n\nEmergent triangular superlattice in Moire systems.\n\n\n\n\n\n","category":"type"},{"location":"#MoireSuperlattices.MoireTriangular-Tuple{Int64, AbstractVector{<:AbstractVector{<:Number}}}","page":"Home","title":"MoireSuperlattices.MoireTriangular","text":"MoireTriangular(truncation::Int, vectors::AbstractVector{<:AbstractVector{<:Number}}; name=:MoireTriangular, origin=nothing)\n\nConstruct the emergent triangular superlattice in Moire systems.\n\n\n\n\n\n","category":"method"},{"location":"#MoireSuperlattices.coefficients-Tuple{BLTMD, MoireTriangular, QuantumLattices.Spatials.BrillouinZone}","page":"Home","title":"MoireSuperlattices.coefficients","text":"coefficients(bltmd::BLTMD, lattice::MoireTriangular, brillouinzone::BrillouinZone; band::Int=dimension(bltmd)) -> Tuple{Vector{Vector{ComplexF64}}, Float64}\ncoefficients(bltmd::Algorithm{<:BLTMD}, lattice::MoireTriangular, brillouinzone::BrillouinZone; band::Int=dimension(bltmd.frontend)) -> Tuple{Vector{Vector{ComplexF64}}, Float64}\n\nGet the coefficients of the hoppings and chemical potential of a bilayer TMD on the emergent triangular lattice.\n\n\n\n\n\n","category":"method"},{"location":"#MoireSuperlattices.terms-Tuple{BLTMD, MoireTriangular, QuantumLattices.Spatials.BrillouinZone}","page":"Home","title":"MoireSuperlattices.terms","text":"terms(bltmd::BLTMD, lattice::MoireTriangular, brillouinzone::BrillouinZone; band::Int=dimension(bltmd), ismodulatable::Bool=true, atol=atol) -> NTuple{2*truncation(lattice)+1, Term}\nterms(bltmd::Algorithm{<:BLTMD}, lattice::MoireTriangular, brillouinzone::BrillouinZone; band::Int=dimension(bltmd.frontend), ismodulatable::Bool=true, atol=atol) -> NTuple{2*truncation(lattice)+1, Term}\n\nGet the hopping terms and chemical potential of a bilayer TMD on the emergent triangular lattice.\n\n\n\n\n\n","category":"method"}]
}