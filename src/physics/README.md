# physics

The model lives here. `hamiltonians.jl` is the element s-Hamiltonian. `fields.jl` is the field equations of motion, an alternative formulation kept for comparison and for tests. `fringe.jl` is the edge-field model: types, profiles, and the s-dependent fringe Hamiltonian. `FringeModel()` leaves the body Hamiltonian unchanged. The edge maps are built in `maps/lie/fringe_maps.jl`.

`Fields` depends only on `Coordinates` and `Elements`. `Hamiltonians` also depends on `LieAlgebra`. `Fringe` depends on `Coordinates`, `Elements`, and `Structure`. `ElementFlow` integrates a magnet, a thick RF cavity, or an explicit `H(z, s)`. `FieldTrack` and `DATracking` call that flow. A thin cavity stays a kick. `DiffAlgebra` and Taylor series are algebras. The Hamiltonian is the model.

A new element formula goes in `hamiltonians.jl`. A matching comparison right-hand side goes in `fields.jl`.
