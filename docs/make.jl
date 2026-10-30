using Documenter

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
    ],
]
makedocs(sitename = "Bell Manuals",
         root = @__DIR__,
         source = "src",
         build = "build",
         format = Documenter.HTML(prettyurls = false, mathengine = Documenter.MathJax3()),
         pages = PAGES,
         remotes = nothing,
         checkdocs = :none,
         doctest = false,
         warnonly = true)
