using Documenter
using AnalyticalEOR

makedocs(;
    sitename = "AnalyticalEOR",
    format = Documenter.HTML(),
    modules = [AnalyticalEOR],
    sitename = "AnalyticalEOR.jl",
    format = Documenter.HTML(),
    pages = [
        "Home" => "index.md",
        "Fractional Flow" => [
            "Water Flooding" => "fractional_flow/waterflooding.md",
            "Chemical Flooding" => "fractional_flow/chemicalflooding.md"
                                ],
        "Single Phase Reactive Transport" => [
            "Theory" => "reactive_transport/theory.md",
            "Solution with Ion Exchange Reactions" => "reactive_transport/ionexchange.md",
            "Plotting Recipes" => "reactive_transport/plotting.md",
                                                ],
            ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = //github.com/rlarao/AnalyticalEOR.jl"
)
