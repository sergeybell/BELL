# RF, acceleration and longitudinal dynamics

```@meta
CurrentModule = Bell
```

Modules: `FieldTrack` (`src/maps/da/field_track.jl`), `Acceleration`
(`src/longitudinal_models/acceleration.jl`), `Smooth`
(`src/longitudinal_models/smooth.jl`). The RF element is [`RFCavity`](@ref),
thin (`L = 0`) or thick (`L > 0`). `Smooth` in `(φ, δ)` is a comparison. It
is not 6D tracking. The Hamiltonian is in
[RF Hamiltonian and acceleration](rf_acceleration.md).

## Is the longitudinal motion truncated?

| where | longitudinal treatment | series cut? |
|:--|:--|:--|
| [`exact_track_element`](@ref), magnet or thick cavity | [`element_flow`](@ref) of ``H``: square roots in [`energy_term`](@ref) and ``p_s`` | no; the cut is the ODE step |
| thin [`rf_kick`](@ref) | exact energy kick; ``E \leftrightarrow \delta``, ``l \leftrightarrow \Delta t`` | no |
| thick cavity, Lie path | ``\exp(-L{:}H{:})`` | ``H`` to degree [`JET_ORDER`](@ref), map to [`MAP_ORDER`](@ref) |
| [`accelerate`](@ref), [`turn_snapshots`](@ref) | the same point step, then [`renormalize`](@ref) after each thin kick | no |
| [`element_da_map`](@ref) | flow of ``H`` on a DA jet; a thin cavity is `rf_kick` on DA | the DA order ``N`` |
| [`track`](@ref) with Lie factors | generators ``f_3``, ``f_4``; a thin cavity is ``\exp(-{:}G{:})`` | yes, at [`MAP_ORDER`](@ref) |
| `Smooth` | ``\eta_0\delta + \eta_1\delta^2`` | yes (averaged model) |

[`element_flow`](@ref) rejects a thin cavity. The field equations in
`src/physics/fields.jl` are the comparison formulation. [`element_rhs`](@ref)
throws for an [`RFCavity`](@ref). A thick cavity is the flow of ``H``.

On a ramp the magnets keep their normalised strengths (they are assumed to
follow ``p_0``). The RF frequency follows ``h f_0(p_0)`` because the next
element is tracked with the reference rebuilt by [`renormalize`](@ref).

## Smooth model

See [`longitudinal_hamiltonian_text`](@ref). [`accelerate`](@ref) does not call it.
