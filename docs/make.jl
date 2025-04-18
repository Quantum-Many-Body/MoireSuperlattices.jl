using MoireSuperlattices
using Documenter

DocMeta.setdocmeta!(MoireSuperlattices, :DocTestSetup, :(using MoireSuperlattices); recursive=true)

makedocs(;
    modules=[MoireSuperlattices],
    authors="waltergu <waltergu1989@gmail.com> and contributors",
    repo="https://github.com/Quantum-Many-Body/MoireSuperlattices.jl/blob/{commit}{path}#{line}",
    sitename="MoireSuperlattices.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Quantum-Many-Body.github.io/MoireSuperlattices.jl",
        edit_link="main",
        assets=["assets/favicon.ico"],
    ),
    pages=[
        "Home" => "index.md",
        "Examples" => [
            "examples/Introduction.md",
            "examples/Visualization.md",
            "examples/HomobilayerTMD.md",
        ],
        "Manual" => "manual.md",
    ],
)

deploydocs(;
    repo="github.com/Quantum-Many-Body/MoireSuperlattices.jl",
    devbranch="main",
)
