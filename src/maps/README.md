# maps

Map construction. `abstract.jl` (`Maps`) holds `AbstractMap`, `LinearMap`, and `DAMap`, and depends only on `DiffAlgebra`. `lie/` applies the Taylor-series algebra to the Hamiltonian. `lie/fringe_maps.jl` builds edge maps from the fringe Hamiltonian. `element_flow.jl` integrates a magnet, a thick RF cavity, or an explicit `H(z, s)`, on numbers, on Taylor series, or on `DiffAlgebra`. `da/` tracks a point and builds the DA map from that flow. A thin cavity stays a kick. The field equations stay for comparison and for tests.

`analysis`, the trackers in `src/tracking/`, and the `Bell` facade use these modules. The edge list is [src/README.md](../README.md).

A new wrapper type goes in `abstract.jl`. Taylor-series factor code goes in `lie/`. The magnetic flow is `element_flow.jl`. An RF kick stays in `da/`.
