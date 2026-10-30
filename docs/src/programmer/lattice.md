# Lattice, structure files and I/O

```@meta
CurrentModule = Bell
```

Modules: `Elements` (`src/elements.jl`), `Structure` (`src/structure.jl`),
`OptimConverter` (`src/optim_converter.jl`).

## Elements are data

An element is an immutable struct holding parameters only. It has no methods
for physics. Each physics module adds its own method:

| question | function | module |
|:--|:--|:--|
| Hamiltonian | [`hamiltonian`](@ref) | Hamiltonians |
| linear matrix | [`linear_matrix`](@ref) | Linear |
| Lie map / factor | [`element_lie_map`](@ref), [`build_factors`](@ref) | Tracking |
| RK4 Taylor map | [`taylor_map_element`](@ref) | DifferentialAlgebra |
| exact tracking | [`exact_track_element`](@ref) | ExactIntegrator |
| edge maps | [`edge_map`](@ref) | Fringe |

To add a new element type:

1. Define `struct NewElem <: AbstractElement` with `name::String` and, if it has
   length, `L::Float64`.
2. Add `kind_string`.
3. Add `hamiltonian(::NewElem, ref)`. Linear, Lie and DA maps then follow
   automatically.
4. Add an `_params` method in ExactIntegrator, or a dedicated field.
5. Add the card in `_make_element` and `_card_line` of `structure.jl`.

The explicit Hamiltonian, matrix, polynomial map and Lie generators of every
element type are on the page [Elements: explicit maps](elements.md); the
element types are documented there as well.


## The `Lattice`

A [`Lattice`](@ref) is the expanded ring: a flat vector of elements in tracking
order, the reference particle and the harmonic number. Repeated cells are the
same element values, so caches compute each distinct element once.
Editing helpers never mutate a lattice; they return a new one.

```@docs
Lattice
circumference
total_bend_angle
unique_elements
with_reference
```

## Structure files (`.bell`)

[`load_structure`](@ref) is a small hand-written parser. It is line-based and
uses four sections: `reference`, `parameters`, `elements`, `sequence`. It has a
recursive-descent expression evaluator (`_eval_sum` → `_eval_prod` →
`_eval_atom`) that supports `+ - * /`, parentheses, numbers, parameters, and a
trailing `deg`. Sequences may be named and repeated (`8 * cell`); expansion is
recursive with cycle detection.
[`write_structure`](@ref) writes numbers with `%.16g`, so `write → load` is
an exact round trip.

```@docs
load_structure
write_structure
StructureFile
```

## Lattice editing

The editing helpers are used by the sextupole and acceleration examples. They
change the sequence (drift splitting), never the element types, and always
keep the circumference.

```@docs
SEXT_FAMILIES
find_optics_maxima
place_sextupole_families
set_family_strengths
insert_cavity
set_cavity_phase
```

## OptiM converter

[`convert_optim`](@ref) is a one-way converter that runs once to produce a
`.bell` file. Bell never reads `.opt` files during a computation. The element
kind comes from the first letter of the name, and units are converted with Bρ
of the reference particle (details in [`OptiMFile`](@ref)).
Anything ignored or guessed is listed in `warnings`.

```@docs
convert_optim
OptiMFile
```
