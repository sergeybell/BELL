# =============================================================================
# docs/make.jl - build the HTML Beam Physics Manual with Documenter.jl
# =============================================================================
# WHAT : renders docs/src/*.md (Markdown + LaTeX math, MathJax) to docs/build/.
# HOW  : julia --project=docs docs/make.jl ; open docs/build/index.html
# NOTE : physics-only manual; no docstrings are pulled from the code.
# =============================================================================
using Documenter
const PAGES = [
    "Home" => "index.md",
    "Coordinates" => "coordinates.md",
    "Coordinate transformation" => "transform.md",
    "Elements and Lie maps" => "elements.md",
    "Linear optics" => "linear.md",
    "Second order" => "second_order.md",
    "Chromaticity" => "chromaticity.md",
    "RF cavity" => "rf.md",
    "Acceleration" => "acceleration.md",
    "Longitudinal motion" => "longitudinal.md",
    "Symplecticity and Taylor maps" => "symplectic.md",
    "Sextupoles and optimiser" => "sextupoles.md",
    "Assumptions and hypotheses" => "assumptions.md",
    "Code map" => "codemap.md",
    "Programmer manual" => "programmer.md",
]
makedocs(sitename = "LieRing Beam Physics Manual", root = @__DIR__, source = "src", build = "build",
         format = Documenter.HTML(prettyurls = false, mathengine = Documenter.MathJax3()),
         pages = PAGES, remotes = nothing, checkdocs = :none, doctest = false,
         warnonly = true)
