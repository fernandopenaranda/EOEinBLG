using EOEinBLG
using Documenter

DocMeta.setdocmeta!(EOEinBLG, :DocTestSetup, :(using EOEinBLG); recursive=true)

makedocs(;
    modules=[EOEinBLG],
    authors="Fernando Peñaranda",
    repo="https://github.com/fernandopenaranda/EOEinBLG.jl/blob/{commit}{path}#{line}",
    sitename="EOEinBLG.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://fernandopenaranda.github.io/EOEinBLG.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/fernandopenaranda/EOEinBLG.jl",
    devbranch="main",
)
