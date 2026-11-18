# analysis

`map_readings.jl` (`MapReadings`) reads coefficients of a map that is already built: closed orbit, dispersion and chromaticity orders, momentum compaction. It does not integrate the lattice.

`sextupole_families.jl` (`SextupoleFamilies`) places and fits sextupole families from an optics table that is already built. It checks that `opt.bounds` lines up with `lattice.elements`. It depends on `Coordinates`, `Elements`, `Structure`, and `MapReadings`.

`MapReadings` depends on `Coordinates`, `Elements`, `Structure`, `LieAlgebra`, `Linear`, Lie `Tracking`, `DATracking`, and `DiffAlgebra`. `FringeMaps` and `NormalForm` call it. The layer order is in [src/README.md](../README.md).

A new reading of an existing map goes in this folder. Building a map stays under `src/maps/`.
