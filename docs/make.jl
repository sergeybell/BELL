using Documenter
using Bell

# Site sources live in docs/site. docs/src keeps the Markdown of the Beam Physics
# Manual, which is published as the PDF (docs/make_pdf.jl), not as HTML pages.
const PAGES = [
    "Home" => "index.md",
    "Physics Reference" => "physics.md",
    "User's Guide" => [
        "Loading a lattice" => "guide/lattice.md",
        "Fringe-field models" => "guide/fringe.md",
        "Optics and chromaticity" => "guide/optics.md",
        "Tracking and acceleration" => "guide/tracking.md",
        "Normal forms" => "guide/normal_form.md",
    ],
    "Developer's Guide" => [
        "Architecture & conventions" => "dev/architecture.md",
        "Lattice & structure files" => "dev/lattice.md",
        "Elements: explicit maps" => "dev/elements.md",
        "Hamiltonians, maps & Lie algebra" => "dev/maps.md",
        "Tracking & integrators" => "dev/tracking.md",
        "Optics" => "dev/optics.md",
        "RF & longitudinal" => "dev/rf.md",
        "Normal form" => "dev/normal_form.md",
        "Fringe fields" => "dev/fringe.md",
        "Differential algebra" => "dev/da.md",
        "ODE integrators" => "dev/integrators.md",
        "Extending Bell" => "dev/extending.md",
    ],
    "API Reference" => [
        "Core" => "api/core.md",
        "Maps & Lie algebra" => "api/maps.md",
        "Optics" => "api/optics.md",
        "Tracking, RF, acceleration" => "api/tracking.md",
        "Normal form" => "api/normal_form.md",
        "Fringe" => "api/fringe.md",
        "DA & integrators" => "api/da.md",
        "Plotting & optimisation" => "api/utilities.md",
        "Index" => "api/index.md",
    ],
]
# All submodules are listed so that `checkdocs = :exports` verifies that every
# exported name has a docstring; the API pages collect them with @autodocs.
const MODULES = [Bell, Bell.Coordinates, Bell.Elements, Bell.Structure, Bell.LieAlgebra,
                 Bell.Hamiltonians, Bell.OptimConverter, Bell.Linear, Bell.Symplectic,
                 Bell.Tracking, Bell.HigherOrder, Bell.DifferentialAlgebra, Bell.Longitudinal,
                 Bell.ExactIntegrator, Bell.Acceleration, Bell.NormalForm, Bell.Fringe,
                 Bell.Plotting, Bell.Optimizers, Bell.DiffAlgebra, Bell.Integrators]

makedocs(sitename = "Bell",
         modules = MODULES,
         root = @__DIR__,
         source = "site",
         build = "build",
         format = Documenter.HTML(prettyurls = false, mathengine = Documenter.MathJax3()),
         pages = PAGES,
         remotes = nothing,
         checkdocs = :exports,
         doctest = false,
         warnonly = true)

# The Physics Reference page links to the PDF next to it in the build.
let pdf = joinpath(@__DIR__, "pdf", "Bell_Beam_Physics_Manual.pdf")
    isfile(pdf) && cp(pdf, joinpath(@__DIR__, "build", "Bell_Beam_Physics_Manual.pdf"); force = true)
end
