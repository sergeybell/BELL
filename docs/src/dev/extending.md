# Extending Bell

```@meta
CurrentModule = Bell
```

## A new element type

1. `src/core/elements.jl`: `struct NewElem <: AbstractElement` with `name::String`
   and, if it has length, `L::Float64`; add `kind_string(::NewElem)`.
2. `src/physics/hamiltonians.jl`: `hamiltonian(e::NewElem, ref, z)` (and
   `hamiltonian(e, ref)` for the Lie polynomial). That scalar is the model.
   [`linear_matrix`](@ref), [`element_lie_map`](@ref), and
   [`element_flow`](@ref) follow from it. [`exact_track_element`](@ref),
   [`element_da_map`](@ref), and [`taylor_map_element`](@ref) call that flow.
   See [Elements](elements.md).
3. Optional check — `src/physics/fields.jl`: `magnetic_parameters` and a
   method of [`element_field`](@ref) / [`element_rhs`](@ref), so tests can
   compare that right-hand side with ``S\nabla H``. Those functions are not
   the step `FieldTrack` or `DATracking` takes.
4. `src/core/structure.jl`: the card in `_make_element` and `_card_line`.
5. A docstring; it appears in the API Reference via `@autodocs`.

## A new fringe model

| new … | implement |
|:--|:--|
| profile | `struct P <: FringeProfile` in `src/physics/fringe.jl`, `_raw(p::P, z)` (rising, 0 → 1) and `_halfwidth(p::P)`; `profile_function` centres it on the effective field boundary |
| integration method | `struct M <: FringeMethod` in `src/physics/fringe.jl`, and `_soft` in `src/maps/lie/fringe_maps.jl` returning the degree-3 map |
| edge model | `struct X <: EdgeModel` and a method of `edge_map` in `src/maps/lie/fringe_maps.jl` |
| element family | an entry in `FringeModel` and in `edge_model` (`src/physics/fringe.jl`) |

The body Hamiltonians are not touched; edges are separate thin items
(see [Fringe fields](fringe.md)). Check new models with `fringe_observables`
(symplecticity of every edge map) and the ``g \to 0`` limit in `examples/15_fringe_fields`.
