using Documenter
using AnalyticalEOR

makedocs(;
    modules = [AnalyticalEOR],
    sitename = "AnalyticalEOR.jl",
    format = Documenter.HTML(sidebar_sitename=false),
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
        "Types and Methods" => "types_and_methods.md"
            ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/rlarao/AnalyticalEOR.jl"
)
