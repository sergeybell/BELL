using Documenter
using Bell

const PAGES = [
    "Architecture" => "architecture.md",
    "Beam Physics" => [
        "Home" => "index.md",
        "Canonical phase space" => "canonical.md",
        "Element Hamiltonians" => "hamiltonians.md",
        "Lie product of the ring" => "lie.md",
        "Linear optics" => "linear.md",
        "Second order" => "second_order.md",
        "Chromaticity" => "chromaticity.md",
        "Fringe fields" => "fringe.md",
        "RF cavity" => "rf.md",
        "Acceleration" => "acceleration.md",
        "Longitudinal motion" => "longitudinal.md",
        "Sextupoles" => "sextupoles.md",
        "Equilibrium momentum shift" => "equilibrium.md",
        "Assumptions" => "assumptions.md",
    ],
    "Programmer Manual" => [
        "Overview" => "programmer/index.md",
        "Architecture & conventions" => "programmer/architecture.md",
        "Lattice, elements & I/O" => "programmer/lattice.md",
        "Hamiltonians, maps & Lie algebra" => "programmer/maps.md",
        "Tracking & integrators" => "programmer/tracking.md",
        "Optics" => "programmer/optics.md",
        "RF & longitudinal" => "programmer/rf.md",
        "Normal form" => "programmer/normal_form.md",
        "Fringe fields" => "programmer/fringe.md",
        "Coordinates & utilities" => "programmer/utilities.md",
        "API index" => "programmer/api.md",
    ],
]
# All submodules are listed so that `checkdocs = :exports` verifies that every
# exported name of every module appears in a @docs block of the Programmer Manual.
const MODULES = [Bell, Bell.Coordinates, Bell.Elements, Bell.Structure, Bell.LieAlgebra,
                 Bell.Hamiltonians, Bell.OptimConverter, Bell.Linear, Bell.Symplectic,
                 Bell.Tracking, Bell.HigherOrder, Bell.DifferentialAlgebra, Bell.Longitudinal,
                 Bell.ExactIntegrator, Bell.Acceleration, Bell.NormalForm, Bell.Fringe,
                 Bell.Plotting, Bell.Optimizers]

makedocs(sitename = "Bell Manuals",
         modules = MODULES,
         root = @__DIR__,
         source = "src",
         build = "build",
         format = Documenter.HTML(prettyurls = false, mathengine = Documenter.MathJax3()),
         pages = PAGES,
         remotes = nothing,
         checkdocs = :exports,
         doctest = false,
         warnonly = true)
