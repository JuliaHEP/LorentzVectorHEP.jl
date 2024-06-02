using Documenter, LorentzVectorHEP

makedocs(;
    modules=[LorentzVectorHEP],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        assets=String[],
    ),
    pages=[
        "Introduction" => "index.md",
    ],
    repo="https://github.com/JuliaHEP/LorentzVectorHEP.jl/blob/{commit}{path}#L{line}",
    sitename="LorentzVectorHEP.jl",
    authors="Jerry Ling and contributors",
)

deploydocs(;
    repo="github.com/JulieHEP/LorentzVectorHEP.jl",
)
