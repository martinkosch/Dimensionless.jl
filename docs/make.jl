using Documenter, Dimensionless

makedocs(sitename="Dimensionless.jl",
pages = ["Home" => "index.md",
        "example.md",
        "library.md"
        ],
format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true")
)
