using Documenter
using Fluxion


makedocs(
    modules = [Fluxion],
    format = :html,
    sitename = "Fluxion.jl",
    authors = "Bogdan Opanchuk",
    pages = [
        "Home" => "index.md",
    ],
)

deploydocs(
    repo = "github.com/fluxion/Fluxion.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
)
