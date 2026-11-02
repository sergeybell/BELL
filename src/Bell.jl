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
include("da/DiffAlgebra.jl")
include("integrators/Integrators.jl")

include("common/coordinates.jl")
include("common/elements.jl")
include("common/structure.jl")
include("lie/lie_algebra.jl")
include("lie/hamiltonians.jl")
include("da/fields.jl")
include("common/optim_converter.jl")
include("lie/linear.jl")
include("common/symplectic.jl")
include("da/field_track.jl")
include("lie/tracking.jl")
include("da/tracking.jl")
include("common/map_readings.jl")
include("longitudinal_models/smooth.jl")
include("longitudinal_models/acceleration.jl")
include("lie/fringe.jl")
include("common/sextupole_families.jl")
include("plotting/Plotting.jl")
include("optimizers/Optimizers.jl")
include("lie/normal_form.jl")

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
