# Coordinates, constants and utilities

```@meta
CurrentModule = Bell
```

Modules: `Coordinates` (`src/coordinates.jl`), `Plotting` (`src/plotting.jl`),
`Optimizers` (`src/optimizers.jl`).

## Particles, reference and coordinate transformations

The reference particle fixes p0, which every normalised quantity refers to.
The energy term ``e(\delta)`` and the velocity ratio ``v/v_0`` are the only
places where the particle's energy enters the Hamiltonians.

```@docs
C_LIGHT
Particle
DEUTERON
RefParticle
reference_from_kinetic
describe_reference
energy_term
velocity_ratio
revolution_frequency
canonical_to_user
user_to_canonical
IX
IPX
IY
IPY
IL
IDELTA
VAR_NAMES
USER_NAMES
```

## Plotting

`Plotting` uses Plots.jl with the GR backend in headless mode
(`GKSwstype = 100`). Every function writes a file and returns its path. The
module is a leaf: no model code depends on it.

```@docs
plot_twiss
plot_dispersion2
plot_chromatic
plot_acceleration
separatrix!
gif_ramp
gif_elements
```

## Optional optimisation

`Optimizers` is never called by the model modules. It provides response
matrices of the ring metrics with respect to the sextupole families, and two
example fitting targets. The targets are illustrations, not physics
requirements.

```@docs
ring_metrics
response
solve_A
solve_B
```
