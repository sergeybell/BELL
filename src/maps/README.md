# maps

Map construction. `abstract.jl` (`Maps`) holds `AbstractMap`, `LinearMap`, and `DAMap`, and depends only on `DiffAlgebra`. `lie/` applies the Taylor-series algebra to the Hamiltonian. `element_flow.jl` integrates a magnet from that Hamiltonian, `dz/ds = S ∇H`, on numbers, on Taylor series, or on `DiffAlgebra`. `da/` tracks a point and builds the DA map from that flow. The field equations stay for comparison and for tests. A thick RF cavity still uses them for its drift halves.

`analysis`, the trackers in `src/tracking/`, and the `Bell` facade use these modules. The edge list is [src/README.md](../README.md).

A new wrapper type goes in `abstract.jl`. Taylor-series factor code goes in `lie/`. The magnetic flow is `element_flow.jl`. An RF kick stays in `da/`.
