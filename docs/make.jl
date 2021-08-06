using Documenter
using EnergySystemModeling

makedocs(
    sitename = "EnergySystemModeling.jl",
    format = Documenter.HTML(
        assets = ["assets/favicon.ico"]
    ),
    modules = [EnergySystemModeling],
    authors = "Lucas Condeixa, Fabricio Oliveira, Jaan Tollander de Balsch",
    pages = [
        "Home" => "index.md",
        "plotting.md",
        "api.md",
        "datageneration.md"
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/gamma-opt/EnergySystemModeling.jl.git"
)
