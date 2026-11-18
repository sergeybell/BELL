# math

Algebras and generic ODE methods. `DiffAlgebra.jl` is the truncated power-series algebra. `lie_series.jl` is included into that module (internal canonical pairs and `da_symplectic_error`). ODE methods live in `integrators/`.

`DiffAlgebra` depends only on `Printf`. `Integrators` depends on `DiffAlgebra` and `LinearAlgebra`. Map code, `Symplectic`, `Fringe`, and the trackers use them. The layer diagram is [src/README.md](../README.md).

A new series operation goes in `DiffAlgebra.jl` or a file it includes. A new ODE method goes in `integrators/`.
