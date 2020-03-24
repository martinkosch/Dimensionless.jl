using Documenter, Dimensionless

makedocs(sitename="Dimensionless.jl",
pages = ["Home" => "index.md",
        "example.md",
        "library.md"
        ],
format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true")
)

deploydocs(
    repo = "github.com/martinkosch/Dimensionless.jl.git",
)
