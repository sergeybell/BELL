# lie

The Taylor-series algebra for maps. `lie_algebra.jl` is that algebra: Poisson bracket and Dragt–Finn factors on `TaylorSeries`. `linear.jl` builds the linear matrix from the Hamiltonian. `tracking.jl` builds and applies the factors. `fringe_maps.jl` builds hard-edge generators and the soft-edge flow from the fringe Hamiltonian. `normal_form.jl` (`NormalForm`) is research: included last and not re-exported (`using Bell.NormalForm`).

`LieAlgebra` depends on `Maps`, and it imports `Coordinates` (`SYMPLECTIC_S`) and `Symplectic` so it can add the `TaylorN` method of `symplectic_error`. `Linear` and `Tracking` also depend on `Coordinates`, `Elements`, `Structure`, `Hamiltonians`, and `LieAlgebra`; `Tracking` also on `Integrators` and `ElementFlow`. `FringeMaps` depends on `Fringe`, `Hamiltonians`, `LieAlgebra`, `Symplectic`, `Tracking`, `MapReadings`, `Integrators`, and `ElementFlow`. `NormalForm` also uses `MapReadings`, `Smooth`, and `Linear`. Callers include `Hamiltonians`, `FringeMaps`, `MapReadings`, and `Trackers`. See [src/README.md](../../README.md).

New factor or linear-matrix code goes in this folder and is included from `Bell.jl` in the current load order. `NormalForm` stays out of the facade re-export.
