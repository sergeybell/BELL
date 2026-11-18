# maps

Map construction. `abstract.jl` (`Maps`) holds `AbstractMap`, `LinearMap`, and `DAMap`, and depends only on `DiffAlgebra`. `lie/` applies the Taylor-series algebra to the Hamiltonian. `da/` integrates the field equations, the comparison formulation, on numbers or on `DiffAlgebra`.

`analysis`, the trackers in `src/tracking/`, and the `Bell` facade use these modules. The edge list is [src/README.md](../README.md).

A new wrapper type goes in `abstract.jl`. Taylor-series factor code goes in `lie/`. A step of the field equations goes in `da/`.
