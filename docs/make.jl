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
        canonical = "https://hinsley.github.io/Kneading.jl/stable/",
    ),
    checkdocs = :exports,
    pages = [
        "Home" => "index.md",
        "One-dimensional maps" => "one-dimensional-maps.md",
        "Kneading diagrams" => "diagrams.md",
    ],
)

deploydocs(
    repo = "github.com/hinsley/Kneading.jl.git",
    devbranch = "main",
)
