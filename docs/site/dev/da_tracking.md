# Field ODEs on DA (`Bell.Fields`, `Bell.DATracking`)

```@meta
CurrentModule = Bell
```

Element transfer maps of any order ``N`` by integrating the field equations of
motion ([`element_field`](@ref)) on DA initial conditions. This path does **not**
build a Hamiltonian in DA and does **not** call the Lie / Dragt–Finn stack.

```julia
using Bell, Bell.DiffAlgebra
ref = reference_from_kinetic(deuteron, 270.0)
M  = element_da_map(Quadrupole("QF", 0.3, 1.5), ref, 7)   # Vector{DA}
Mr = one_turn_map_da(lat, 5)
```

## Physics (`src/da/fields.jl`)

[`element_field`](@ref) writes ``dz/ds`` from hard-edge multipole strengths and
relativistic kinematics in canonical ``z``. The same function evaluates on
`Float64` or on DA. Thin RF is [`rf_kick`](@ref) (shared with `FieldTrack`).

## Tracking

| call | role |
|:--|:--|
| [`map_space`](@ref)`(N)` | cached 6-variable DA space of order N |
| [`element_da_map`](@ref)`(e, ref, N)` | integrate field ODEs on the identity jet |
| [`one_turn_map_da`](@ref)`(lat, N)` | compose element maps (drops thin RF) |
| [`exact_track_element`](@ref) | same ODEs on a Float64 point |

## Relation to the Lie path

At low order the DA field map of a short magnet agrees with
`lie_map(hamiltonian(e, ref), L)` to working precision; the codes do not share
`H`. Symplecticity of a DA map is checked with [`da_symplectic_error`](@ref).

## Chromaticity

[`chromaticity_series`](@ref)`(lattice, N)` builds [`one_turn_map_da`](@ref) and
reads ``Q(\delta)`` as a series in δ.
