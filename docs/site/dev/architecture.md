# Architecture and conventions

```@meta
CurrentModule = Bell
```

This page explains how the package is organised. Two independent physics stacks
share elements, coordinates and ODE integrators; they do not call each other.

## Package layout

```
src/
  Bell.jl                 top module: include order + re-export
  common/                 shared by both stacks
    coordinates.jl        Coordinates        particle, reference, laboratory maps
    elements.jl           Elements           element types (RFCavity thin or thick)
    structure.jl          Structure          Lattice, .bell parser/writer, magnetic_lattice
    optim_converter.jl    OptimConverter     OptiM -> .bell
    map_readings.jl       MapReadings        closed orbit, chromaticity, compaction
    symplectic.jl         Symplectic         symplecticity checks
    sextupole_families.jl SextupoleFamilies
  lie/                    Hamiltonian / Lie-algebra path
    lie_algebra.jl        LieAlgebra         polynomial space, Poisson, Lie factors
    hamiltonians.jl       Hamiltonians       polynomial H of every element (RF: rf_potential)
    linear.jl             Linear             linear matrices / optics
    tracking.jl           Tracking           Lie-factor tracking, Taylor map of H
    fringe.jl             Fringe             HardEdge / SoftEdge
    normal_form.jl        NormalForm         research (using Bell.NormalForm)
  da/                     field / DA path
    DiffAlgebra.jl        DiffAlgebra        truncated series (pure maths)
    lie_series.jl         (in DiffAlgebra)   Poisson / Lie series on DA
    fields.jl             Fields             dz/ds from multipole fields
    field_track.jl        FieldTrack         Float64 field ODEs, rf_kick, renormalize
    tracking.jl           DATracking         field ODEs on DA (any order N)
  longitudinal_models/    approximate models kept for comparison
    smooth.jl             Smooth             averaged (phi, delta)
    acceleration.jl       Acceleration       ramp: rf_kick + renormalize
  integrators/            Integrators        one file per ODE method
  optimizers/, plotting/
  ...
```

## Dependency graph

```
Coordinates ─► Elements ─► Structure ─► OptimConverter
     │             │           │
     │             │           ├─► Hamiltonians ─► Linear ─► Tracking ─► MapReadings
     │             │           └─► Fields ─► FieldTrack ─► Acceleration
     │             │                     └─► DATracking ─► MapReadings
     └──► LieAlgebra ◄─────────┘
DiffAlgebra, Integrators (pure mathematics)
research/NormalForm (not re-exported)
```

## Conventions

### Thin RF

Only [`RFCavity`](@ref) (zero length). [`rf_kick`](@ref) uses closed-form
relativistic formulas in canonical `z` (numbers or DA).

### Polynomials (Lie path)

Hamiltonians are expanded to degree [`JET_ORDER`](@ref) = 4; maps are cut at
[`MAP_ORDER`](@ref) = 3. Arbitrary order lives on the field / DA path.

## Data flow

```
 .opt ──convert_optim──► .bell ──load_structure──► Lattice
                                                    │
     Lie path: hamiltonian ─► lie_map ─► Dragt-Finn ─► track / readings
     Field path: element_field ─► FieldTrack (Float64) / DATracking (order N)
     RF: rf_kick on both paths; accelerate = FieldTrack + renormalize
```
