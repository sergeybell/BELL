# da

Steps of the field equations of motion, the comparison formulation in `src/physics/fields.jl`. `field_track.jl` integrates them in `Float64` and applies the thin `rf_kick`. `tracking.jl` (`DATracking`) integrates the same right-hand side on the `DiffAlgebra` algebra to a chosen order.

`FieldTrack` depends on `Coordinates`, `Elements`, `Structure`, `Fields`, and `Integrators`. `DATracking` adds `DiffAlgebra` and `FieldTrack`. `Trackers` and `Acceleration` call `FieldTrack`; `MapReadings` calls `DATracking`. Edges: [src/README.md](../../README.md).

A new step of this formulation goes here. The s-Hamiltonian stays in `src/physics/hamiltonians.jl`.
