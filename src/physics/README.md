# physics

The model lives here. `hamiltonians.jl` is the element s-Hamiltonian. `fields.jl` is the field equations of motion, an alternative formulation kept for comparison. `fringe.jl` is the edge-field model chosen when a ring is built (`FringeModel()` leaves the body Hamiltonian unchanged).

`Fields` depends only on `Coordinates` and `Elements`. `Hamiltonians` also depends on `LieAlgebra`. `Fringe` also depends on `Structure`, `Symplectic`, `Integrators`, and on map modules listed in [src/README.md](../README.md). Lie `Linear` and `Tracking` take the Hamiltonian; `FieldTrack` and `DATracking` take the field equations. `DiffAlgebra` and Taylor series are algebras. The Hamiltonian is the model.

A new element formula goes in `hamiltonians.jl`. A matching comparison right-hand side goes in `fields.jl`.
