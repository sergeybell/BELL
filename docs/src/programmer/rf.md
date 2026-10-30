# RF and longitudinal dynamics

```@meta
CurrentModule = Bell
```

Modules: `Longitudinal` (`src/longitudinal.jl`), `Acceleration`
(`src/acceleration.jl`). The RF elements and integrators are documented under
[Lattice](lattice.md) ([`RFCavity`](@ref), [`Cavity`](@ref)) and
[Tracking](tracking.md) ([`rf_kick`](@ref), [`exact_track_element`](@ref)).

## Smooth longitudinal model

With turn number n, ``\phi = -h\omega_0\Delta t`` and slip factors η0, η1 read
from the one-turn map ([`path_expansion`](@ref)):

```math
\frac{d\phi}{dn} = -2\pi h(\eta_0\delta + \eta_1\delta^2),\qquad
\frac{d\delta}{dn} = A\,(\sin(\phi_s+\phi) - \sin\phi_s),\quad A = \frac{ZVT}{\beta^2E}.
```

Its invariant is [`long_H`](@ref). The module also provides the fixed points,
the separatrix, the bucket height and the linear synchrotron tune. The η1 term
creates the second (α-bucket) fixed point and makes the bucket asymmetric.

```@docs
LongitudinalParams
long_H
synchrotron_tune
fixed_points
separatrix
bucket_height
longitudinal_hamiltonian_text
```

## Acceleration

[`accelerate`](@ref) tracks an ensemble with the exact integrator. Particle 1
is the synchronous particle. After every passage through the cavity, the
reference is moved onto it by [`renormalize`](@ref). This ramps p0 and produces
adiabatic damping of px, py and δ automatically. At the target energy the phase
switches to a flat top. [`measure_tune`](@ref) extracts tunes from turn-by-turn
data without needing the exact fixed point.

```@docs
accelerate
turn_snapshots
measure_tune
```
