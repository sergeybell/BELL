# core

Lattice data and checks that do not build maps: coordinates and the reference particle, element types, sequences and `Lattice`, `symplectic_error`, sextupole families, map caches, and plot-name stubs. `caches.jl` and `plotting.jl` are included into `Bell`.

`Coordinates` and `Elements` import no Bell module. `Structure` uses those two. `Symplectic` also uses `DiffAlgebra` and `LieAlgebra`. `SextupoleFamilies` also uses `Structure`, `LieAlgebra`, and `MapReadings`. Later layers take these types; see [src/README.md](../README.md).

A new element or lattice query goes here. The Hamiltonian belongs in `src/physics/`. A map belongs in `src/maps/`.
