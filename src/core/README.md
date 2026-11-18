# core

Lattice data and checks that do not build maps: coordinates and the reference particle (`SYMPLECTIC_S` is the coordinate pairing), element types, sequences and `Lattice`, the numeric and DA methods of `symplectic_error`, the `clear_caches!` registry, and plot-name stubs. `caches.jl` and `plotting.jl` are included into `Bell`. The map modules own the cache dictionaries; `Bell.__init__` registers their emptiers. Sextupole families live in `src/analysis/`.

`Coordinates` and `Elements` import no Bell module. `Structure` uses those two. `Symplectic` uses `Coordinates` and `DiffAlgebra`. The `TaylorN` method of `symplectic_error` is added by `LieAlgebra`. Later layers take these types; see [src/README.md](../README.md).

A new element or lattice query goes here. The Hamiltonian belongs in `src/physics/`. A map belongs in `src/maps/`. A reading of a map that is already built belongs in `src/analysis/`.
