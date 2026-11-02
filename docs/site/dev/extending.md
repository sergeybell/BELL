# Extending Bell

```@meta
CurrentModule = Bell
```

## A new element type

1. `src/elements.jl`: `struct NewElem <: AbstractElement` with `name::String`
   and, if it has length, `L::Float64`; add `kind_string(::NewElem)`.
2. Lie path — `src/hamiltonians.jl`: `hamiltonian(e::NewElem, ref)`. Linear
   matrix, Lie map and Dragt–Finn then follow (see [Elements](elements.md)).
3. Field path — `src/da/fields.jl`: `magnetic_parameters` and/or a dedicated
   branch of [`element_field`](@ref) / [`element_rhs`](@ref). `FieldTrack` and
   `DATracking` then use the same ODEs.
4. `src/structure.jl`: the card in `_make_element` and `_card_line`.
5. A docstring; it appears in the API Reference via `@autodocs`.

## A new fringe model

| new … | implement |
|:--|:--|
| profile | `struct P <: FringeProfile`, `_raw(p::P, z)` (rising, 0 → 1) and `_halfwidth(p::P)`; `profile_function` centres it on the EFB |
| integration method | `struct M <: FringeMethod` and `_soft(E::_Edge, m::M, ref)` returning the degree-3 map over ``[-W, W]`` |
| edge model | `struct X <: EdgeModel` and a branch in `edge_map(E, ref)` returning a polynomial map |
| element family | an entry in `FringeModel` and in `edge_model` |

The body Hamiltonians are not touched; edges are separate thin items
(see [Fringe fields](fringe.md)). Check new models with `fringe_observables`
(symplecticity of every edge map) and the ``g \to 0`` limit of example 19.
