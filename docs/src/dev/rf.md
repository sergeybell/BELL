# RF, acceleration and longitudinal dynamics

```@meta
CurrentModule = Bell
```

Modules: `FieldTrack`, `Acceleration`, `Smooth`. The only RF element is thin
[`RFCavity`](@ref). Smooth `(φ, δ)` is analysis only — not 6D tracking.

## Is the longitudinal motion truncated in 6D field tracking? — No.

| where | longitudinal treatment | truncated? |
|:--|:--|:--|
| [`exact_track_element`](@ref) (magnets) | field ODEs with exact roots | **no** |
| [`rf_kick`](@ref) | exact energy kick; ``E \leftrightarrow \delta``, ``l \leftrightarrow \Delta t`` | **no** |
| [`accelerate`](@ref), [`turn_snapshots`](@ref) | field track + [`renormalize`](@ref) after each thin kick | **no** |
| [`element_da_map`](@ref) / DA `rf_kick` | same formulas on DA → jet to order N | series cut is the DA order only |
| [`track`](@ref) with Lie factors | magnets truncated at `MAP_ORDER`; RF still exact `rf_kick` | Lie magnets yes |
| `Smooth` | ``\eta_0\delta + \eta_1\delta^2`` | yes (averaged model) |

Assumptions in tracking: magnets ramp with ``p_0`` (fixed normalised strengths);
RF frequency follows ``h f_0(p_0)``.

## Smooth model

See [`longitudinal_hamiltonian_text`](@ref). Not used by `accelerate`.
