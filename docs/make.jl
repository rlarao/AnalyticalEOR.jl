using Documenter
using AnalyticalEOR

makedocs(
    sitename = "AnalyticalEOR",
    format = Documenter.HTML(),
    modules = [AnalyticalEOR],
    sitename="AnalyticalEOR.jl",
    format=Documenter.HTML(),
    pages =[
        "Home" => "index.md",
        "Fractional Flow" => [
            "Water Flooding" => "fractionalflow/waterflooding.md"
            "Chemical Flooding" => "fractionalflow/chemicalflooding.md"
        ],
        "Single Phase Reactive Transport" => [
            "Theory" => "reactivetransport/theory.md"
            "Solution with Ion Exchange Reactions" => "reactivetransport/ionexchange.md"
            "Plotting Recipes" => "reactivetransport/plotting.md"
        ]
        
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = //github.com/rlarao/AnalyticalEOR.jl"
)
