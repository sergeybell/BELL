# =============================================================================
# Bell.jl - top-level module
# =============================================================================
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# WHAT
#   Loads every block in dependency order and re-exports public names.
#   A script only needs `using Bell`.
#
# =============================================================================
"""
    Bell

Beam-dynamics package for a storage ring. The package works in 6D canonical
coordinates (x, px, y, py, l, δ) with Lie-algebraic transfer maps of
order 3, exact symplectic integration, RF/longitudinal motion, normal forms
and fringe-field models.

`Bell` only includes the submodules in dependency order and re-exports every
exported name, so `using Bell` is all a script needs. See the Programmer
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
include("tracking.jl")
include("higher_order.jl")
include("differential_algebra.jl")
include("longitudinal.jl")
include("exact_integrator.jl")
include("acceleration.jl")
include("normal_form.jl")
include("fringe.jl")
include("lie_flow.jl")
include("plotting.jl")
include("optimizers.jl")

using .Coordinates, .Elements, .Structure, .LieAlgebra, .Hamiltonians, .OptimConverter,
      .Linear, .Symplectic, .Tracking, .HigherOrder, .DifferentialAlgebra, .NormalForm,
      .Longitudinal, .ExactIntegrator, .Acceleration, .Plotting, .Optimizers, .Fringe, .LieFlow

for mod in (Coordinates, Elements, Structure, LieAlgebra, Hamiltonians, OptimConverter,
            Linear, Symplectic, Tracking, HigherOrder, DifferentialAlgebra, NormalForm,
            Longitudinal, ExactIntegrator, Acceleration, Plotting, Optimizers, Fringe, LieFlow)
    for n in names(mod)
        n === nameof(mod) && continue
        @eval export $n
    end
end

end # module
