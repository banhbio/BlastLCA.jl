using BlastLCA
using Documenter

DocMeta.setdocmeta!(BlastLCA, :DocTestSetup, :(using BlastLCA); recursive=true)

makedocs(;
    modules=[BlastLCA],
    authors="banhbio <ban@kuicr.kyoto-u.ac.jp>",
    repo="https://github.com/banhbio/BlastLCA.jl/blob/{commit}{path}#{line}",
    sitename="BlastLCA.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://banhbio.github.io/BlastLCA.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/banhbio/BlastLCA.jl",
)
