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
module Bell

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
include("plotting.jl")
include("optimizers.jl")

using .Coordinates, .Elements, .Structure, .LieAlgebra, .Hamiltonians, .OptimConverter,
      .Linear, .Symplectic, .Tracking, .HigherOrder, .DifferentialAlgebra, .NormalForm,
      .Longitudinal, .ExactIntegrator, .Acceleration, .Plotting, .Optimizers

for mod in (Coordinates, Elements, Structure, LieAlgebra, Hamiltonians, OptimConverter,
            Linear, Symplectic, Tracking, HigherOrder, DifferentialAlgebra, NormalForm,
            Longitudinal, ExactIntegrator, Acceleration, Plotting, Optimizers)
    for n in names(mod)
        n === nameof(mod) && continue
        @eval export $n
    end
end

end # module
