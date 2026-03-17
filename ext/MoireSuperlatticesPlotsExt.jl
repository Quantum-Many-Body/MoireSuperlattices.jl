module MoireSuperlatticesPlotsExt

using RecipesBase: @recipe, @series
using MoireSuperlattices: CommensurateBilayerHoneycomb, Lattice, Neighbors, Bond, vectors, angle, distance, reciprocals, hexagon120°map, hexagon60°map, rcoordinate, str
using QuantumLattices: bonds

"""
    plot(moire::CommensurateBilayerHoneycomb, choice::Symbol, n=2*ceil(Int, √count(moire)); topcolor=:red, bottomcolor=:blue, vector=true, vectorcolor=:green, moirecolor=:black, anglecolor=:grey)

Plot a Moire superlattice composed of two layers of honeycomb lattices in the real space or in the reciprocal space.
"""
@recipe function plot(moire::CommensurateBilayerHoneycomb, choice::Symbol, n=2*ceil(Int, √count(moire)); topcolor=:red, bottomcolor=:blue, vector=true, vectorcolor=:green, moirecolor=:black, anglecolor=:grey)
    @assert choice∈(:real, :reciprocal) "plot error: incorrect choice (`:$choice`), which should be either `:real` or `:reciprocal`."
    θ = angle(moire)
    top, bottom, (t₁, t₂) = Lattice(moire, :top), Lattice(moire, :bottom), vectors(moire)
    title --> "Twisted Bilayer Honeycomb ($(str(rad2deg(θ)))°)"
    aspect_ratio := :equal
    legend := false
    if choice == :real
        neighbors = Neighbors(1=>distance(moire.coordinates[:, 1], moire.coordinates[:, 2]))
        @series begin
            color --> topcolor
            Lattice(top, (2n, 2n); mode=:center), neighbors
        end
        @series begin
            color --> bottomcolor
            Lattice(bottom, (2n, 2n); mode=:center), neighbors
        end
        arrow := true
        linewidth := 2
        color --> vectorcolor
        alpha := (vector ? 1.0 : 0.0)
        [moire.center[1], t₁[1]+moire.center[1], NaN, moire.center[1], t₂[1]+moire.center[1]], [moire.center[2], t₁[2]+moire.center[2], NaN, moire.center[2], t₂[2]+moire.center[2]]
    else
        recipls₁ = reciprocals(top)
        recipls₂ = reciprocals(bottom)
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

end # module
