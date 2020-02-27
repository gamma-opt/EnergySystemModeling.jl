using Documenter

push!(LOAD_PATH, dirname(@__DIR__))
using EnergySystemModel

makedocs(
    sitename = "EnergySystemModel",
    format = Documenter.HTML(
        assets = ["assets/favicon.ico"]
    ),
    modules = [EnergySystemModel],
    authors = "Lucas Condeixa, Fabricio Oliveira, Jaan Tollander de Balsch",
    pages = [
        "Home" => "index.md",
        "api.md",
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
