# Extending Bell

```@meta
CurrentModule = Bell
```

## A new element type

1. `src/elements.jl`: `struct NewElem <: AbstractElement` with `name::String`
   and, if it has length, `L::Float64`; add `kind_string(::NewElem)`.
2. `src/hamiltonians.jl`: `hamiltonian(e::NewElem, ref)` built from `Z`,
   `energy_term` and `p_s` exactly like the existing methods. Linear matrix, Lie
   map, Dragt–Finn factor and RK4 map then follow automatically
   (see [Elements](elements.md)).
3. `src/exact_integrator.jl`: an `_params(e::NewElem)` method if the element is a
   special case of ``(h, k_1, k_2)``, otherwise a dedicated vector field and a
   branch in `exact_track_element`.
4. `src/structure.jl`: the card in `_make_element` and `_card_line`, so the
   element can be read and written.
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
