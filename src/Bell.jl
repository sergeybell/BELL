# =============================================================================
# Bell.jl - top-level module
# BELL - Beam Evolution in Accelerator Lattice
# =============================================================================
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# Loads every block in dependency order and re-exports public names.
# A script only needs `using Bell`. NormalForm is included as research and
# is not re-exported; load it with `using Bell.NormalForm`.
# =============================================================================
"""
    Bell

BELL - Beam Evolution in Accelerator Lattice.

Beam-dynamics package for a storage ring. The package works in 6D canonical
coordinates (x, px, y, py, l, δ). One formulation of the element Hamiltonian
is truncated as a polynomial or integrated exactly; Lie maps, Taylor maps and
the exact tracker are three uses of that formulation, not three independent
trackers. The averaged turn-by-turn model in (φ, δ) is Smooth.

`Bell` includes the submodules in dependency order and re-exports every
exported name except the research module `NormalForm`. See the Programmer
Manual for the module map.
"""
module Bell

# pure-mathematics submodules (no accelerator concepts); not re-exported into Bell
# because their names (compose, integrate, ...) are generic: `using Bell.DiffAlgebra`,
# `using Bell.Integrators`.
include("da/DiffAlgebra.jl")
include("integrators/Integrators.jl")

include("coordinates.jl")
include("elements.jl")
include("structure.jl")
include("lie_algebra.jl")
include("hamiltonians.jl")
include("optim_converter.jl")
include("linear.jl")
include("symplectic.jl")
include("exact_integrator.jl")
include("tracking.jl")
include("lie_flow.jl")
include("map_readings.jl")
include("smooth.jl")
include("acceleration.jl")
include("fringe.jl")
include("sextupole_families.jl")
include("plotting/Plotting.jl")
include("optimizers/Optimizers.jl")
include("research/normal_form.jl")

using .Coordinates, .Elements, .Structure, .LieAlgebra, .Hamiltonians, .OptimConverter,
      .Linear, .Symplectic, .ExactIntegrator, .Tracking, .LieFlow, .MapReadings, .Smooth,
      .Acceleration, .Fringe, .SextupoleFamilies, .Plotting, .Optimizers

for mod in (Coordinates, Elements, Structure, LieAlgebra, Hamiltonians, OptimConverter,
            Linear, Symplectic, ExactIntegrator, Tracking, LieFlow, MapReadings, Smooth,
            Acceleration, Fringe, SextupoleFamilies, Plotting, Optimizers)
    for n in names(mod)
        n === nameof(mod) && continue
        @eval export $n
    end
end

end # module
