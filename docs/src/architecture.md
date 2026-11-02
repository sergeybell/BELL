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
| `src/coordinates.jl` | Coordinates | Particle, reference, canonical and laboratory coordinates |
| `src/elements.jl` | Elements | Data types: drift, bend, quad, sext, thin RF, marker |
| `src/structure.jl` | Structure | Lattice and structure file |
| `src/hamiltonians.jl` | Hamiltonians | Polynomial `H` for the Lie path only |
| `src/da/fields.jl` | Fields | Equations of motion `dz/ds` from multipole fields |
| `src/lie_algebra.jl` | Lie algebra | Poisson bracket, Lie series, Dragt–Finn |
| `src/tracking.jl` | Tracking | Lie factors of truncated `H`; thin RF via `rf_kick` |
| `src/field_track.jl` | FieldTrack | Float64 tracking of field ODEs + `rf_kick` |
| `src/da/tracking.jl` | DATracking | Same field ODEs on DA (any order N) |
| `src/da/DiffAlgebra.jl` | DiffAlgebra | Truncated power series (pure maths) |
| `src/smooth.jl` | Smooth | Averaged `(φ, δ)` — not 6D RF |
| `src/acceleration.jl` | Acceleration | Ramp with thin `RFCavity` + `renormalize` |

## Two stacks

**Lie path.** Polynomial Hamiltonian → Dragt–Finn → `Tracking`.

**Field / DA path.** Write `dz/ds = f(z)` from hard-edge fields
(`da/fields.jl`), not from `H`. Evaluate on `Float64` (`FieldTrack`) or on DA
(`DATracking`) to any order `N`. Same thin `rf_kick`.

## Thin RF in 6D

Only zero-length [`RFCavity`](@ref) and [`rf_kick`](@ref). On a ramp,
[`renormalize`](@ref) updates `p0` after the kick of the synchronous particle.

## Smooth model

`smooth.jl` does not track the ring in 6D. Acceleration does not use it.
