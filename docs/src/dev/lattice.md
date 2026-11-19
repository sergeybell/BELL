# Lattice and structure files (internals)

```@meta
CurrentModule = Bell
```

Modules: `Elements` (`src/core/elements.jl`), `Structure`
(`src/core/structure.jl`), `OptimConverter` (`src/io/optim_converter.jl`).

## Elements are data

An element is an immutable struct holding parameters only. It has no methods
for physics. Each physics module adds its own method:

| question | function | module |
|:--|:--|:--|
| Hamiltonian | [`hamiltonian`](@ref) | Hamiltonians |
| linear matrix | [`linear_matrix`](@ref) | Linear |
| Lie map / factor | [`element_lie_map`](@ref), [`build_factors`](@ref) | Tracking |
| flow of ``H`` | [`element_flow`](@ref) | ElementFlow |
| Taylor check | [`taylor_map_element`](@ref) | Tracking (`src/maps/lie/tracking.jl`) |
| point tracking | [`exact_track_element`](@ref) | FieldTrack (`src/maps/da/field_track.jl`) |
| DA map | [`element_da_map`](@ref) | DATracking (`src/maps/da/tracking.jl`) |
| edge maps | [`edge_map`](@ref) | FringeMaps (`src/maps/lie/fringe_maps.jl`) |

How to add a new element type: [Extending Bell](extending.md).

The explicit Hamiltonian, matrix, polynomial map and Lie generators of every
element type are on the page [Elements: explicit maps](elements.md).

## The `Lattice`

A [`Lattice`](@ref) is the expanded ring: a flat vector of elements in tracking
order, the reference particle and the harmonic number. Repeated cells are the
same element values, so caches compute each distinct element once.
Editing helpers never mutate a lattice; they return a new one.

## Structure files (`.bell`)

[`load_structure`](@ref) is a small hand-written parser. It is line-based and
uses four sections: `reference`, `parameters`, `elements`, `sequence`. It has a
recursive-descent expression evaluator (`_eval_sum` → `_eval_prod` →
`_eval_atom`) that supports `+ - * /`, parentheses, numbers, parameters, and a
trailing `deg`. Sequences may be named and repeated (`8 * cell`); expansion is
recursive with cycle detection.
[`write_structure`](@ref) writes numbers with `%.16g`, so `write → load` is
an exact round trip.

## Lattice editing

The editing helpers are used by the sextupole and acceleration examples. They
change the sequence (a named drift is split around a thin cavity), and
[`insert_cavity`](@ref) keeps [`path_length`](@ref) when the cavity length is
zero.

## `.opt` converter

[`convert_optim`](@ref) reads a `.opt` file once and writes a `.bell` file.
Tracking does not read `.opt`. The element kind comes from the first letter
of the name. Lengths and fields are converted with the reference magnetic
rigidity ``B\rho``. Anything ignored or guessed is listed in `warnings`.
