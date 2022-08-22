using Documenter
using AnalyticalEOR

makedocs(
    sitename = "AnalyticalEOR",
    format = Documenter.HTML(),
    modules = [AnalyticalEOR]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/rlarao/AnalyticalEOR.jl"
)
