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
#
# Layers: math → core → io → physics → maps → analysis →
# longitudinal_models → optimizers → plotting.
# A new `using ..X` follows that direction. Include order is the existing
# load order: Hamiltonians and Symplectic stay after LieAlgebra, Fringe
# after Tracking, SextupoleFamilies after MapReadings.
# =============================================================================
"""
    Bell

BELL - Beam Evolution in Accelerator Lattice.

Two independent physics stacks share elements, coordinates and integrators:

- **Lie path** — polynomial `hamiltonian`, Lie algebra / Dragt–Finn, `Tracking`.
- **Field / DA path** — equations of motion in `Fields`, `FieldTrack` (Float64)
  and `DATracking` (any order N). Thin `RFCavity` uses `rf_kick` in both stacks.

Smooth `(φ, δ)` is a separate averaged model, not 6D RF tracking.
`NormalForm` is research-only (`using Bell.NormalForm`).
"""
module Bell

# pure-mathematics DiffAlgebra / Integrators are not re-exported into Bell
# (generic names): `using Bell.DiffAlgebra`, `using Bell.Integrators`.
include("math/DiffAlgebra.jl")
include("math/integrators/Integrators.jl")

include("core/coordinates.jl")
include("core/elements.jl")
include("core/structure.jl")
include("maps/lie/lie_algebra.jl")
include("physics/hamiltonians.jl")
include("physics/fields.jl")
include("io/optim_converter.jl")
include("maps/lie/linear.jl")
include("core/symplectic.jl")
include("maps/da/field_track.jl")
include("maps/lie/tracking.jl")
include("maps/da/tracking.jl")
include("analysis/map_readings.jl")
include("longitudinal_models/smooth.jl")
include("longitudinal_models/acceleration.jl")
include("physics/fringe.jl")
include("core/sextupole_families.jl")
include("plotting/Plotting.jl")
include("optimizers/Optimizers.jl")
include("maps/lie/normal_form.jl")

using .Coordinates, .Elements, .Structure, .LieAlgebra, .Hamiltonians, .Fields,
      .OptimConverter, .Linear, .Symplectic, .FieldTrack, .Tracking, .DATracking,
      .MapReadings, .Smooth, .Acceleration, .Fringe, .SextupoleFamilies, .Plotting, .Optimizers

for mod in (Coordinates, Elements, Structure, LieAlgebra, Hamiltonians, Fields,
            OptimConverter, Linear, Symplectic, FieldTrack, Tracking, DATracking,
            MapReadings, Smooth, Acceleration, Fringe, SextupoleFamilies, Plotting, Optimizers)
    for n in names(mod)
        n === nameof(mod) && continue
        @eval export $n
    end
end

end # module
