module MoireSuperlatticesMakieExt

using QuantumLattices: Lattice, Neighbors, distance, hexagon120°map, hexagon60°map, reciprocals, str
using MoireSuperlattices: CommensurateBilayerHoneycomb, vectors
import Makie

# Plot a Moire superlattice composed of two layers of honeycomb lattices in the real space or in the reciprocal space.
@inline Makie.plot(moire::CommensurateBilayerHoneycomb, choice::Symbol, n=2*ceil(Int, √count(moire)); kwargs...) = Makie.plot!(Makie.Figure(), moire, choice, n; kwargs...)
function Makie.plot!(fig::Makie.Figure, moire::CommensurateBilayerHoneycomb, choice::Symbol, n=2*ceil(Int, √count(moire)); kwargs...)
    Makie.plot!(Makie.Axis(fig[1, 1]), moire, choice, n; kwargs...)
    return fig
end
function Makie.plot!(
    ax::Makie.AbstractAxis, moire::CommensurateBilayerHoneycomb, choice::Symbol, n=2*ceil(Int, √count(moire));
    topcolor=:red, bottomcolor=:blue, vector=true, vectorcolor=:green, moirecolor=:black, anglecolor=:grey, kwargs...
    )
    @assert choice in (:real, :reciprocal) "plot error: incorrect choice (`:$choice`), which should be either `:real` or `:reciprocal`."
    θ = angle(moire)
    title = "Twisted Bilayer Honeycomb ($(str(rad2deg(θ)))°)"
    top, bottom, (t₁, t₂) = Lattice(moire, :top), Lattice(moire, :bottom), vectors(moire)
    if choice == :real
        neighbors = Neighbors(1=>distance(moire.coordinates[:, 1], moire.coordinates[:, 2]))
        Makie.plot!(ax, Lattice(top, (2n, 2n); mode=:center), neighbors; title, color=topcolor, kwargs...)
        Makie.plot!(ax, Lattice(bottom, (2n, 2n); mode=:center), neighbors; title, color=bottomcolor, kwargs...)
        vector && Makie.arrows2d!(ax, [moire.center[1], moire.center[1]], [moire.center[2], moire.center[2]], [t₁[1], t₂[1]], [t₁[2], t₂[2]]; color=vectorcolor)
    else
        recipls₁ = reciprocals(top)
        recipls₂ = reciprocals(bottom)
        filter = bond->bond.kind==1
        Makie.plot!(ax, Lattice([collect(mapreduce(*, +, hexagon60°map[key], recipls₁)) for key in ("K₁", "K₂", "K₃", "K₄", "K₅", "K₆")]...), 1, filter; title, color=topcolor, kwargs...)
        Makie.plot!(ax, Lattice([collect(mapreduce(*, +, hexagon60°map[key], recipls₂)) for key in ("K₁", "K₂", "K₃", "K₄", "K₅", "K₆")]...), 1, filter; title, color=bottomcolor, kwargs...)
        Kt = collect(mapreduce(*, +, hexagon60°map["K₂"], recipls₁))
        Kb = collect(mapreduce(*, +, hexagon60°map["K₂"], recipls₂))
        Makie.lines!(ax, [0.0, Kt[1]], [0.0, Kt[2]]; color=anglecolor, linestyle=:dot, linewidth=2)
        Makie.lines!(ax, [0.0, Kb[1]], [0.0, Kb[2]]; color=anglecolor, linestyle=:dot, linewidth=2)
        recipls = reciprocals([t₁, t₂])
        Makie.plot!(ax, Lattice(Lattice([collect(mapreduce(*, +, hexagon120°map[key], recipls)) for key in ("K₁", "K₂")]...; vectors=recipls), (n, n); mode=:center), 1, filter; title, color=moirecolor, kwargs...)
    end
    return ax
end

end # module
