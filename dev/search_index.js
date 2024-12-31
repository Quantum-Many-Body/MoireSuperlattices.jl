var documenterSearchIndex = {"docs":
[{"location":"examples/HomobilayerTMD/","page":"Twisted Homobilayer Transition Metal Dichalcogenides","title":"Twisted Homobilayer Transition Metal Dichalcogenides","text":"CurrentModule = MoireSuperlattices","category":"page"},{"location":"examples/HomobilayerTMD/#Twisted-Homobilayer-Transition-Metal-Dichalcogenides","page":"Twisted Homobilayer Transition Metal Dichalcogenides","title":"Twisted Homobilayer Transition Metal Dichalcogenides","text":"","category":"section"},{"location":"examples/HomobilayerTMD/","page":"Twisted Homobilayer Transition Metal Dichalcogenides","title":"Twisted Homobilayer Transition Metal Dichalcogenides","text":"For twisted homobilayer transition metal dichalcogenides, such as the twisted bilayer WSe₂ or MoTe₂, the continuum model [Phys. Rev. Lett. 122, 086402 (2019)] has two material parameters, i.e., a₀ (monolayer lattice constant), and m (effective mass of the conduction band), one device parameter, i.e., θ (twist angle), and five phenomenological parameters, i.e., Vᶻ (perpendicular displacement field), μ (chemical potential), V (amplitude of Moire potential), ψ (phase of Moire potential), and w (interlayer hopping amplitude).","category":"page"},{"location":"examples/HomobilayerTMD/#Trivial-Band-Topology","page":"Twisted Homobilayer Transition Metal Dichalcogenides","title":"Trivial Band Topology","text":"","category":"section"},{"location":"examples/HomobilayerTMD/","page":"Twisted Homobilayer Transition Metal Dichalcogenides","title":"Twisted Homobilayer Transition Metal Dichalcogenides","text":"The band topology of the moire bands depends on the aforementioned parameters. For a typical set of parameters [Phys. Rev. Research 2, 033087 (2020)] that may be related to twisted bilayer WSe₂, the topmost moire band is topologically trivial.","category":"page"},{"location":"examples/HomobilayerTMD/","page":"Twisted Homobilayer Transition Metal Dichalcogenides","title":"Twisted Homobilayer Transition Metal Dichalcogenides","text":"First, we calculate the energy bands with such a set of parameters in both valleys:","category":"page"},{"location":"examples/HomobilayerTMD/","page":"Twisted Homobilayer Transition Metal Dichalcogenides","title":"Twisted Homobilayer Transition Metal Dichalcogenides","text":"using MoireSuperlattices\nusing Plots\nusing QuantumLattices\nusing TightBindingApproximation\n\n# `a₀`: monolayer lattice constant (Å)\n# `m`: effective mass of the conduction band (mₑ)\n# `θ`: twist angle (°)\n# `Vᶻ`: perpendicular displacement field (meV)\n# `μ`: chemical potential (meV)\n# `V`: amplitude of Moire potential (meV)\n# `ψ`: phase of Moire potential (°)\n# `w`: interlayer hopping amplitude (meV)\nparameters = (a₀=3.28, m=0.45, θ=3.70, Vᶻ=0.0, μ=0.0, V=4.4, ψ=5.9, w=20.0)\nWSe₂ = Algorithm(\n    :WSe₂, BLTMD(values(parameters)...; truncation=4);\n    parameters=parameters, map=bltmdmap\n)\nrecipls = WSe₂.frontend.reciprocallattice.translations\nbands₁ = WSe₂(:EB, EnergyBands(ReciprocalPath(recipls, hexagon\"Γ-K₁-M₁-Γ\", length=100)))\nbands₂ = WSe₂(:EB, EnergyBands(ReciprocalPath(recipls, hexagon\"Γ-K₄-M₄-Γ\", length=100)))\n\nemin, emax = -100.0, 30.0\nplt = plot()\nplot!(plt, bands₁, ylim=(emin, emax), color=\"blue\", title=\"\")\nplot!(plt, bands₂, ylim=(emin, emax), color=\"blue\", title=\"\")","category":"page"},{"location":"examples/HomobilayerTMD/","page":"Twisted Homobilayer Transition Metal Dichalcogenides","title":"Twisted Homobilayer Transition Metal Dichalcogenides","text":"The bands are two-fold degenerate due to the valley degrees of freedom. Note that due to the spin-valley lock in WSe₂, such a degeneracy can also be viewed as the spin degeneracy.","category":"page"},{"location":"examples/HomobilayerTMD/","page":"Twisted Homobilayer Transition Metal Dichalcogenides","title":"Twisted Homobilayer Transition Metal Dichalcogenides","text":"Then we compute the Berry curvature and Chern number of the topmost moire band:","category":"page"},{"location":"examples/HomobilayerTMD/","page":"Twisted Homobilayer Transition Metal Dichalcogenides","title":"Twisted Homobilayer Transition Metal Dichalcogenides","text":"berry = WSe₂(:BC, BerryCurvature(BrillouinZone(recipls, 48), [dimension(WSe₂)]))\nplot(berry)","category":"page"},{"location":"examples/HomobilayerTMD/","page":"Twisted Homobilayer Transition Metal Dichalcogenides","title":"Twisted Homobilayer Transition Metal Dichalcogenides","text":"The perpendicular displacement field can induced a layer dependent potential, which will lead to the splitting of the two spins:","category":"page"},{"location":"examples/HomobilayerTMD/","page":"Twisted Homobilayer Transition Metal Dichalcogenides","title":"Twisted Homobilayer Transition Metal Dichalcogenides","text":"update!(WSe₂; Vᶻ=8.0)\nbands₁ = WSe₂(:EB, EnergyBands(ReciprocalPath(recipls, hexagon\"Γ-K₁-M₁-Γ\", length=100)))\nbands₂ = WSe₂(:EB, EnergyBands(ReciprocalPath(recipls, hexagon\"Γ-K₄-M₄-Γ\", length=100)))\n\nplt = plot()\nplot!(plt, bands₁, ylim=(emin, emax), color=\"blue\", title=\"\")\nplot!(plt, bands₂, ylim=(emin, emax), color=\"blue\", title=\"\")","category":"page"},{"location":"examples/HomobilayerTMD/","page":"Twisted Homobilayer Transition Metal Dichalcogenides","title":"Twisted Homobilayer Transition Metal Dichalcogenides","text":"For the topologically trivial topmost moire band, it forms a triangular lattice. The effective tight-binding model of this sole band can be obtained:","category":"page"},{"location":"examples/HomobilayerTMD/","page":"Twisted Homobilayer Transition Metal Dichalcogenides","title":"Twisted Homobilayer Transition Metal Dichalcogenides","text":"lattice = MoireTriangular(6, reciprocals(recipls))\nhilbert = Hilbert(Fock{:f}(1, 2), length(lattice))\ntba = Algorithm(\n    :tba,\n    TBA(lattice, hilbert, terms(WSe₂, lattice, BrillouinZone(recipls, 24); tol=10^-6))\n)\nnothing # hide","category":"page"},{"location":"examples/HomobilayerTMD/","page":"Twisted Homobilayer Transition Metal Dichalcogenides","title":"Twisted Homobilayer Transition Metal Dichalcogenides","text":"Here, the hopping parameters are truncated up to the 6th order. We can compare the energy bands of this tight-binding model with those of the continuum model:","category":"page"},{"location":"examples/HomobilayerTMD/","page":"Twisted Homobilayer Transition Metal Dichalcogenides","title":"Twisted Homobilayer Transition Metal Dichalcogenides","text":"bands₃ = tba(:EB, EnergyBands(ReciprocalPath(recipls, hexagon\"Γ-K-M-Γ\", length=100)))\nplt = plot()\nplot!(plt, bands₁, ylim=(emin, emax), color=\"blue\", title=\"\")\nplot!(plt, bands₂, ylim=(emin, emax), color=\"blue\", title=\"\")\nplot!(plt, bands₃, ylim=(emin, emax), color=\"red\", lw=4, alpha=0.3, title=\"\")","category":"page"},{"location":"examples/HomobilayerTMD/#Nontrivial-Band-Topology","page":"Twisted Homobilayer Transition Metal Dichalcogenides","title":"Nontrivial Band Topology","text":"","category":"section"},{"location":"examples/HomobilayerTMD/","page":"Twisted Homobilayer Transition Metal Dichalcogenides","title":"Twisted Homobilayer Transition Metal Dichalcogenides","text":"With another set of parameters [Phys. Rev. Lett. 132, 036501 (2024)] that was proposed to be suitable for twisted bilayer MoTe₂, the top two moire bands have opposite Chern numbers.","category":"page"},{"location":"examples/HomobilayerTMD/","page":"Twisted Homobilayer Transition Metal Dichalcogenides","title":"Twisted Homobilayer Transition Metal Dichalcogenides","text":"First, we calculate the energy bands with such a set of parameters in a single valley:","category":"page"},{"location":"examples/HomobilayerTMD/","page":"Twisted Homobilayer Transition Metal Dichalcogenides","title":"Twisted Homobilayer Transition Metal Dichalcogenides","text":"using MoireSuperlattices\nusing Plots\nusing QuantumLattices\nusing TightBindingApproximation\n\nparameters = (a₀=3.52, m=0.6, θ=3.89, Vᶻ=0.0, μ=30.0, V=20.8, ψ=107.7, w=-23.8)\nMoTe₂ = Algorithm(\n    :MoTe₂, BLTMD(values(parameters)...; truncation=4);\n    parameters=parameters, map=bltmdmap\n)\nrecipls = MoTe₂.frontend.reciprocallattice.translations\nbands = MoTe₂(:EB, EnergyBands(ReciprocalPath(recipls, hexagon\"Γ-K₁-M₁-Γ\", length=100)))\nplot(bands, ylim=(-50.0, 10.0), color=\"blue\", title=\"\")","category":"page"},{"location":"examples/HomobilayerTMD/","page":"Twisted Homobilayer Transition Metal Dichalcogenides","title":"Twisted Homobilayer Transition Metal Dichalcogenides","text":"Then we compute the Berry curvature and Chern numbers of the top two moire bands:","category":"page"},{"location":"examples/HomobilayerTMD/","page":"Twisted Homobilayer Transition Metal Dichalcogenides","title":"Twisted Homobilayer Transition Metal Dichalcogenides","text":"berry = MoTe₂(\n    :BC, BerryCurvature(BrillouinZone(recipls, 48), [dimension(MoTe₂), dimension(MoTe₂)-1])\n)\nplot(berry)","category":"page"},{"location":"examples/Visualization/","page":"Visualization of Moire Superlattices by Twist","title":"Visualization of Moire Superlattices by Twist","text":"CurrentModule = MoireSuperlattices","category":"page"},{"location":"examples/Visualization/#Visualization-of-Moire-Superlattices-by-Twist","page":"Visualization of Moire Superlattices by Twist","title":"Visualization of Moire Superlattices by Twist","text":"","category":"section"},{"location":"examples/Visualization/","page":"Visualization of Moire Superlattices by Twist","title":"Visualization of Moire Superlattices by Twist","text":"Here follows a quick view of the commensurate moire superlattices composed of twisted homobilayer graphene:","category":"page"},{"location":"examples/Visualization/","page":"Visualization of Moire Superlattices by Twist","title":"Visualization of Moire Superlattices by Twist","text":"using MoireSuperlattices\nusing Plots\n\nanim = @animate for i ∈ 2:10\n    plot(\n        CommensurateBilayerHoneycomb((i, 1)), :real, 30;\n        xlim=(-50, 50), ylim=(-30, 30), axis=false, size=(600, 400), vector=false\n    )\nend\ngif(anim; fps=1.5)","category":"page"},{"location":"examples/Visualization/","page":"Visualization of Moire Superlattices by Twist","title":"Visualization of Moire Superlattices by Twist","text":"The crystal structure of a commensurate moire superlattice composed of twisted homobilayer graphene is characterized by two coprime positive integers (m, r) [Phys. Rev. B 86, 155449 (2012)]. In general, we need to distinguish two cases, i.e., whether gcd(r, 3) is 1 or 3, where gcd denotes the greatest common divisor.","category":"page"},{"location":"examples/Visualization/#When-gcd(r,-3)-is-1","page":"Visualization of Moire Superlattices by Twist","title":"When gcd(r, 3) is 1","text":"","category":"section"},{"location":"examples/Visualization/","page":"Visualization of Moire Superlattices by Twist","title":"Visualization of Moire Superlattices by Twist","text":"For example, m=8, r=1:","category":"page"},{"location":"examples/Visualization/","page":"Visualization of Moire Superlattices by Twist","title":"Visualization of Moire Superlattices by Twist","text":"plot(CommensurateBilayerHoneycomb((8, 1)), :real; xlim=(-20, 20), ylim=(-10, 20))","category":"page"},{"location":"examples/Visualization/","page":"Visualization of Moire Superlattices by Twist","title":"Visualization of Moire Superlattices by Twist","text":"In the reciprocal space, the K points of the top and bottom layers correspond to two inequivalent K points of the Moire Brillouin zone:","category":"page"},{"location":"examples/Visualization/","page":"Visualization of Moire Superlattices by Twist","title":"Visualization of Moire Superlattices by Twist","text":"plot(CommensurateBilayerHoneycomb((8, 1)), :reciprocal; xlim=(-5, 5), ylim=(-4, 4))","category":"page"},{"location":"examples/Visualization/#When-gcd(r,-3)-is-3","page":"Visualization of Moire Superlattices by Twist","title":"When gcd(r, 3) is 3","text":"","category":"section"},{"location":"examples/Visualization/","page":"Visualization of Moire Superlattices by Twist","title":"Visualization of Moire Superlattices by Twist","text":"For example, m=20, r=3:","category":"page"},{"location":"examples/Visualization/","page":"Visualization of Moire Superlattices by Twist","title":"Visualization of Moire Superlattices by Twist","text":"plot(CommensurateBilayerHoneycomb((20, 3)), :real; xlim=(-20, 20), ylim=(-10, 20))","category":"page"},{"location":"examples/Visualization/","page":"Visualization of Moire Superlattices by Twist","title":"Visualization of Moire Superlattices by Twist","text":"In the reciprocal space, the K points of the top and bottom layers correspond to two equivalent K points of the Moire Brillouin zone:","category":"page"},{"location":"examples/Visualization/","page":"Visualization of Moire Superlattices by Twist","title":"Visualization of Moire Superlattices by Twist","text":"plot(CommensurateBilayerHoneycomb((20, 3)), :reciprocal; xlim=(-5, 5), ylim=(-4, 4))","category":"page"},{"location":"examples/Introduction/","page":"Introduction","title":"Introduction","text":"CurrentModule = MoireSuperlattices","category":"page"},{"location":"examples/Introduction/#examples","page":"Introduction","title":"Introduction","text":"","category":"section"},{"location":"examples/Introduction/","page":"Introduction","title":"Introduction","text":"Here are some examples to illustrate how this package could be used.","category":"page"},{"location":"examples/Introduction/","page":"Introduction","title":"Introduction","text":"Pages = [\n        \"Visualization.md\",\n        \"HomobilayerTMD.md\",\n        ]\nDepth = 2","category":"page"},{"location":"manual/","page":"Manual","title":"Manual","text":"CurrentModule = MoireSuperlattices","category":"page"},{"location":"manual/","page":"Manual","title":"Manual","text":"Modules = [MoireSuperlattices]","category":"page"},{"location":"manual/#MoireSuperlattices.BLTMD","page":"Manual","title":"MoireSuperlattices.BLTMD","text":"BLTMD{\n    L<:MoireReciprocalLattice,\n    D<:Function,\n    S<:OperatorGenerator,\n    Q<:Quadraticization,\n    H<:CategorizedGenerator{<:OperatorSum{<:Quadratic}}\n} <: MoireSystem{NamedTuple{(:a₀, :m, :θ, :Vᶻ, :μ), NTuple{5, Float64}}, L, D, S, Q, H}\n\nTwisted transition metal dichalcogenide homobilayers.\n\n\n\n\n\n","category":"type"},{"location":"manual/#MoireSuperlattices.BLTMD-NTuple{8, Number}","page":"Manual","title":"MoireSuperlattices.BLTMD","text":"BLTMD(a₀::Number, m::Number, θ::Number, Vᶻ::Number, μ::Number, V::Number, ψ::Number, w::Number; truncation::Int=4)\n\nContinuum model of twisted transition metal dichalcogenide homobilayers.\n\nHere, the parameters are as follows:\n\na₀: monolayer lattice constant (Å)\nm: effective mass of the conduction band (mₑ)\nθ: twist angle (°)\nVᶻ: perpendicular displacement field (meV)\nμ: chemical potential (meV)\nV: amplitude of Moire potential (meV)\nψ: phase of Moire potential (°)\nw: interlayer hopping amplitude (meV)\n\n\n\n\n\n","category":"method"},{"location":"manual/#MoireSuperlattices.CommensurateBilayerHoneycomb","page":"Manual","title":"MoireSuperlattices.CommensurateBilayerHoneycomb","text":"CommensurateBilayerHoneycomb\n\nCommensurate Moire superlattice composed of two layers of honeycomb lattices.\n\n\n\n\n\n","category":"type"},{"location":"manual/#MoireSuperlattices.MoireReciprocalLattice","page":"Manual","title":"MoireSuperlattices.MoireReciprocalLattice","text":"MoireReciprocalLattice{T<:Number} <: AbstractLattice{2, T, 0}\n\nMoire reciprocal lattice with truncation.\n\n\n\n\n\n","category":"type"},{"location":"manual/#MoireSuperlattices.MoireSpace","page":"Manual","title":"MoireSuperlattices.MoireSpace","text":"MoireSpace <: SimpleInternal{MoireSpinor{Int, Int, Int, Rational{Int}, Int}}\n\nThe internal degrees of freedom of Moire systems.\n\n\n\n\n\n","category":"type"},{"location":"manual/#MoireSuperlattices.MoireSpinor","page":"Manual","title":"MoireSuperlattices.MoireSpinor","text":"MoireSpinor{V<:Union{Int, Colon}, L<:Union{Int, Colon}, S<:Union{Int, Colon}, P<:Union{Rational{Int}, Colon}, N<:Union{Int, Colon}} <: SimpleInternalIndex\n\nThe index of the internal degrees of freedom of Moire systems.\n\n\n\n\n\n","category":"type"},{"location":"manual/#MoireSuperlattices.MoireSuperlattice","page":"Manual","title":"MoireSuperlattices.MoireSuperlattice","text":"MoireSuperlattice{D<:Number} <: AbstractLattice{2, D, 2}\n\nAbstract type of the emergent superlattices in Moire systems.\n\n\n\n\n\n","category":"type"},{"location":"manual/#MoireSuperlattices.MoireSystem","page":"Manual","title":"MoireSuperlattices.MoireSystem","text":"MoireSystem{P<:Parameters, L<:MoireReciprocalLattice, D<:Function, S<:OperatorGenerator, Q<:Quadraticization, H<:CategorizedGenerator{<:OperatorSum{<:Quadratic}}} <: TBA{Fermionic{:TBA}, H, Nothing}\n\nThe continuum model of Moire systems.\n\n\n\n\n\n","category":"type"},{"location":"manual/#MoireSuperlattices.MoireTriangular","page":"Manual","title":"MoireSuperlattices.MoireTriangular","text":"MoireTriangular{N, D<:Number} <: MoireSuperlattice{D}\n\nEmergent triangular superlattice in Moire systems.\n\n\n\n\n\n","category":"type"},{"location":"manual/#MoireSuperlattices.MoireTriangular-Tuple{Int64, AbstractVector{<:AbstractVector{<:Number}}}","page":"Manual","title":"MoireSuperlattices.MoireTriangular","text":"MoireTriangular(truncation::Int, vectors::AbstractVector{<:AbstractVector{<:Number}}; name=:MoireTriangular, origin=nothing)\n\nConstruct the emergent triangular superlattice in Moire systems.\n\n\n\n\n\n","category":"method"},{"location":"manual/#QuantumLattices.Spatials.Lattice-Tuple{CommensurateBilayerHoneycomb, Symbol}","page":"Manual","title":"QuantumLattices.Spatials.Lattice","text":"Lattice(moire::CommensurateBilayerHoneycomb, type::Symbol)\n\nGet the minimum unit of the top/bottom layer of a commensurate Moire superlattice composed of two layers of honeycomb lattices.\n\n\n\n\n\n","category":"method"},{"location":"manual/#Base.angle-Tuple{CommensurateBilayerHoneycomb}","page":"Manual","title":"Base.angle","text":"angle(moire::CommensurateBilayerHoneycomb) -> Float64\n\nGet the twist angle of a commensurate Moire superlattice composed of two layers of honeycomb lattices.\n\n\n\n\n\n","category":"method"},{"location":"manual/#Base.count-Tuple{CommensurateBilayerHoneycomb}","page":"Manual","title":"Base.count","text":"count(moire::CommensurateBilayerHoneycomb) -> Int\n\nCount the number of honeycomb unitcells contained in the unitcell of a commensurate Moire superlattice composed of two layers of honeycomb lattices.\n\nThe total number of atoms in the unitcell of the Moire superlattice is 4 times this result because of the AB sublattice and the top/bottom layer degrees of freedom.\n\n\n\n\n\n","category":"method"},{"location":"manual/#MoireSuperlattices.coefficients-Tuple{BLTMD, MoireTriangular, QuantumLattices.Spatials.BrillouinZone}","page":"Manual","title":"MoireSuperlattices.coefficients","text":"coefficients(bltmd::BLTMD, lattice::MoireTriangular, brillouinzone::BrillouinZone; band::Int=dimension(bltmd)) -> Tuple{Vector{Vector{ComplexF64}}, Float64}\ncoefficients(bltmd::Algorithm{<:BLTMD}, lattice::MoireTriangular, brillouinzone::BrillouinZone; band::Int=dimension(bltmd.frontend)) -> Tuple{Vector{Vector{ComplexF64}}, Float64}\n\nGet the coefficients of the hoppings and chemical potential of a bilayer TMD on the emergent triangular lattice.\n\n\n\n\n\n","category":"method"},{"location":"manual/#MoireSuperlattices.terms-Tuple{BLTMD, MoireTriangular, QuantumLattices.Spatials.BrillouinZone}","page":"Manual","title":"MoireSuperlattices.terms","text":"terms(bltmd::BLTMD, lattice::MoireTriangular, brillouinzone::BrillouinZone; band::Int=dimension(bltmd), ismodulatable::Bool=true, tol=tol) -> NTuple{2*truncation(lattice)+1, Term}\nterms(bltmd::Algorithm{<:BLTMD}, lattice::MoireTriangular, brillouinzone::BrillouinZone; band::Int=dimension(bltmd.frontend), ismodulatable::Bool=true, tol=tol) -> NTuple{2*truncation(lattice)+1, Term}\n\nGet the hopping terms and chemical potential of a bilayer TMD on the emergent triangular lattice.\n\n\n\n\n\n","category":"method"},{"location":"manual/#MoireSuperlattices.vectors-Tuple{CommensurateBilayerHoneycomb}","page":"Manual","title":"MoireSuperlattices.vectors","text":"vectors(moire::CommensurateBilayerHoneycomb) -> SVector{2, SVector{2, Float64}}\n\nGet the translation vectors of a commensurate Moire superlattice composed of two layers of honeycomb lattices.\n\n\n\n\n\n","category":"method"},{"location":"manual/#RecipesBase.apply_recipe","page":"Manual","title":"RecipesBase.apply_recipe","text":"plot(moire::CommensurateBilayerHoneycomb, choice::Symbol, n=2*ceil(Int, √count(moire)); topcolor=:red, bottomcolor=:blue, vector=true, vectorcolor=:green, moirecolor=:black, anglecolor=:grey)\n\nPlot a Moire superlattice composed of two layers of honeycomb lattices in the real space or in the reciprocal space.\n\n\n\n\n\n","category":"function"},{"location":"#MoireSuperlattices","page":"Home","title":"MoireSuperlattices","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"(Image: Stable) (Image: Dev) (Image: Build Status) (Image: Coverage) (Image: 996.icu) (Image: LICENSE) (Image: LICENSE) (Image: Code Style: Blue) (Image: ColPrac: Contributor's Guide on Collaborative Practices for Community Packages)","category":"page"},{"location":"#Introduction","page":"Home","title":"Introduction","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package aims at the construction of Moire Superlattices, acquisition of their effective models, and investigation of their low-energy physics, based on the QuantumLattices pack and the TightBindingApproximation pack.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"In Julia v1.10+, please type ] in the REPL to use the package mode, then type this command:","category":"page"},{"location":"","page":"Home","title":"Home","text":"pkg> add MoireSuperlattices","category":"page"},{"location":"#Getting-Started","page":"Home","title":"Getting Started","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Examples of Moire Superlattice","category":"page"},{"location":"#Note","page":"Home","title":"Note","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Due to the fast development of this package, releases with different minor version numbers are not guaranteed to be compatible with previous ones before the release of v1.0.0. Comments are welcomed in the issues.","category":"page"},{"location":"#Contact","page":"Home","title":"Contact","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"waltergu1989@gmail.com","category":"page"}]
}