# da

Point tracking and DA maps. Magnets and a thick RF cavity are integrated from the s-Hamiltonian by `ElementFlow` (`../element_flow.jl`). `field_track.jl` is the point tracker: a magnet or a thick cavity calls that flow, a marker is the identity, and a thin cavity is a kick. `tracking.jl` (`DATracking`) builds the DA map from the same flow. The field equations in `src/physics/fields.jl` remain for comparison and for tests.

`FieldTrack` depends on `Coordinates`, `Elements`, `Structure`, `Integrators`, and `ElementFlow`. `DATracking` depends on `DiffAlgebra`, `Integrators`, `Coordinates`, `Elements`, `Structure`, `ElementFlow`, and `FieldTrack`. `Trackers` and `Acceleration` call `FieldTrack`; `MapReadings` calls `DATracking`. Edges: [src/README.md](../../README.md).

A new magnetic integrator belongs in `element_flow.jl`. An RF kick stays in this folder. The s-Hamiltonian stays in `src/physics/hamiltonians.jl`.
