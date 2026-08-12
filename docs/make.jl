using Documenter
using Kneading

DocMeta.setdocmeta!(
    Kneading,
    :DocTestSetup,
    :(using Kneading);
    recursive = true,
)

makedocs(
    modules = [
        Kneading,
        Kneading.OneDimensionalMaps,
        Kneading.Diagrams,
    ],
    authors = "Carter Hinsley and contributors",
    sitename = "Kneading.jl",
    format = Documenter.HTML(
        prettyurls = true,
        canonical = "https://hinsley.github.io/Kneading.jl/dev/",
    ),
    checkdocs = :exports,
    pages = [
        "Home" => "index.md",
        "Tutorials" => [
            "Chebyshev cubic kneading diagram" =>
                "tutorials/chebyshev-cubic-kneading-diagram.md",
        ],
        "One-dimensional maps" => "one-dimensional-maps.md",
        "Kneading diagrams" => "diagrams.md",
    ],
)

deploydocs(
    repo = "github.com/hinsley/Kneading.jl.git",
    devbranch = "main",
)
