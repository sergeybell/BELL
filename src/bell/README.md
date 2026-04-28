# `bell` package index

This file is the **folder-level index**: it describes what each file in this package is for.

## Files

- `__init__.py`: public package exports, version exposure.
- `_version.py`: package version.
- `elements.py`: element definitions (drift, quadrupole, bend, …) and their current **linear** transfer maps (`map6()`).
- `lattice.py`: lattice container + matrix composition + single-particle tracking (current MVP).
- `optics.py`: basic linear optics helpers (e.g. Twiss from 2×2).
- `particle.py`: `Particle` representation (6D phase space vector).

## Planned additions (next)

- `maps/`: transfer maps (linear part + higher-order Lie map machinery)
- `models/`: element physics models (parameter -> Map), traced to `docs/ref/`
- `beam.py`: beam/distribution helpers (multi-particle tracking)

