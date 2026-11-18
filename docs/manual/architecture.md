# Bell architecture

**Author:** Kolokolchikov Sergey Dmitrievich  
**Email:** sergey.bell13@gmail.com

Bell is a Julia program for symplectic modelling of a magnetic ring. Two
independent physics stacks share elements, coordinates and ODE integrators.
They do not call each other; comparison is done from the outside.

## Blocks

```mermaid
flowchart TD
  coords[coordinates]
  elements[elements]
  structure[structure_file_and_Lattice]
  H[hamiltonians]
  fields[da_fields]
  da[DiffAlgebra]
  integrators[integrators_one_file_each]
  lie[lie_algebra]
  datrack[da_tracking]
  linear[linear_matrix]
  readings[map_readings]
  track[lie_factor_tracking]
  ftrack[field_track]
  long[smooth]
  fringe[fringe]
  coords --> elements --> structure
  structure --> H
  structure --> fields
  H --> lie --> track
  fields --> ftrack
  fields --> datrack
  da --> datrack
  integrators --> ftrack
  integrators --> track
  integrators --> datrack
  integrators --> fringe
  lie --> linear
  lie --> readings
  datrack --> readings
  ftrack --> accel[acceleration]
  readings --> long
  research[research_normal_form]
```

| File | Block | Role |
|------|--------|------|
| `src/core/coordinates.jl` | Coordinates | Particle, reference, canonical and laboratory coordinates |
| `src/core/elements.jl` | Elements | Data types: drift, bend, quad, sext, thin RF, marker |
| `src/core/structure.jl` | Structure | Lattice and structure file |
| `src/physics/hamiltonians.jl` | Hamiltonians | s-Hamiltonian of every element |
| `src/physics/fields.jl` | Fields | Comparison equations of motion `dz/ds` |
| `src/maps/element_flow.jl` | ElementFlow | Integrate a magnet from `H`: `dz/ds = S ∇H` |
| `src/maps/lie/lie_algebra.jl` | Lie algebra | Poisson bracket, Lie series, Dragt–Finn |
| `src/maps/lie/tracking.jl` | Tracking | Lie factors of truncated `H`; thin RF via `rf_kick` |
| `src/maps/da/field_track.jl` | FieldTrack | Float64 tracking of magnets from `H`, plus `rf_kick` |
| `src/maps/da/tracking.jl` | DATracking | Magnet maps from `H` on DA (any order N) |
| `src/math/DiffAlgebra.jl` | DiffAlgebra | Truncated power series (pure maths) |
| `src/longitudinal_models/smooth.jl` | Smooth | Averaged `(φ, δ)` — not 6D RF |
| `src/longitudinal_models/acceleration.jl` | Acceleration | Ramp with thin `RFCavity` + `renormalize` |

## Two stacks

**Lie path.** Polynomial Hamiltonian → Dragt–Finn → `Tracking`.

**Magnet flow.** `dz/ds = S ∇H` from the element Hamiltonian (`ElementFlow`).
`FieldTrack` evaluates it on `Float64`. `DATracking` evaluates it on DA to any
order `N`. The field equations remain for comparison. Thin RF is the same
`rf_kick` on both paths. A thick cavity is drift/kick/drift.

## Thin RF in 6D

Only zero-length [`RFCavity`](@ref) and [`rf_kick`](@ref). On a ramp,
[`renormalize`](@ref) updates `p0` after the kick of the synchronous particle.

## Smooth model

`smooth.jl` does not track the ring in 6D. Acceleration does not use it.
