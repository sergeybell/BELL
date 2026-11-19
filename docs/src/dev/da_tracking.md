# DA maps from the s-Hamiltonian (`Bell.DATracking`)

```@meta
CurrentModule = Bell
```

[`element_da_map`](@ref) builds a transfer map of any order ``N`` by integrating
``dz/ds = S\nabla H`` on the DA identity jet. The integrator is
[`element_flow`](@ref) (`src/maps/element_flow.jl`). This path does not
integrate [`element_field`](@ref), and it does not call the Lie / Dragt–Finn
stack. Both computations read the same [`hamiltonian`](@ref). The field
equations in `src/physics/fields.jl` stay for comparison and for tests.

```julia
using Bell, Bell.DiffAlgebra
ref = reference_from_kinetic(deuteron, 270.0)
M  = element_da_map(Quadrupole("QF", 0.3, 1.5), ref, 7)   # Vector{DA}
Mr = one_turn_map_da(lat, 5)                              # RF included
```

## How the jet is differentiated

On a DA state, [`element_flow`](@ref) builds the polynomial vector field once
with [`vector_field_on_identity`](@ref Bell.DiffAlgebra.vector_field_on_identity), in a space one order above the map.
[`deriv`](@ref Bell.DiffAlgebra.deriv) drops the top degree of a jet that is already truncated, so
the extra order is what keeps degree ``N`` of ``S\nabla H``. The integrator
evaluates that polynomial on the current jet.

On numbers and on `TaylorN` the same call evaluates
[`symplectic_gradient`](@ref Bell.DiffAlgebra.symplectic_gradient) at each stage. That gradient is the forward-mode
[`Dual`](@ref Bell.DiffAlgebra.Dual) partial of the scalar `H`. The dual partial and `deriv` agree
through order ``N-1``, and at order ``N`` when `H` is a polynomial of degree
at most ``N``.

A [`Marker`](@ref) is the identity. A thin [`RFCavity`](@ref) is
[`rf_kick`](@ref) on the DA variables (`path_length` is required).
[`element_flow`](@ref) rejects a thin cavity. A thick cavity is the flow of
[`hamiltonian`](@ref)`(e, ref, path_length, z)` and needs `path_length > 0`.

## Calls

| call | role |
|:--|:--|
| [`map_space`](@ref)`(N)` | cached 6-variable DA space of order N (`src/maps/da/tracking.jl`) |
| [`element_da_map`](@ref)`(e, ref, N)` | flow of ``H`` on the identity jet; thin RF is `rf_kick` |
| [`element_da_maps`](@ref)`(lat, N)` | one map per element; an `RFCavity` uses `path_length(lat)` and is not cached |
| [`one_turn_map_da`](@ref)`(lat, N)` | composition of those maps, including RF |
| [`exact_track_element`](@ref) | the same flow on a `Float64` point (`src/maps/da/field_track.jl`) |
| [`magnetic_lattice`](@ref) | drop a thin cavity; replace a thick cavity by a drift of the same length |

The default method is [`DormandPrince853`](@ref Bell.Integrators.DormandPrince853) with `tol = 1e-14`. Pass
`nsteps` for a fixed step. An [`SBend`](@ref) map is the body only.
[`element_da_maps`](@ref) warns when an edge angle is nonzero.

[`element_field`](@ref) (`src/physics/fields.jl`) writes ``dz/ds`` from
hard-edge multipole strengths. `test/test_physics_consistency.jl` checks that
this right-hand side matches ``S\nabla H``. Production maps do not call it.
[`element_rhs`](@ref) throws for a [`Marker`](@ref) and for an [`RFCavity`](@ref).

## Relation to the Lie path

At low order the DA map of a short magnet agrees with
`lie_map(hamiltonian(e, ref), L)` to working precision. A DA map is symplectic
only up to its truncation. Symplecticity of a Jacobian, a `TaylorN` map, a DA
map, or a tracking step is [`symplectic_error`](@ref). The DA method of that
check is [`da_symplectic_error`](@ref Bell.DiffAlgebra.da_symplectic_error).

## Chromaticity

[`chromaticity_series`](@ref)`(lattice, N)` builds [`one_turn_map_da`](@ref)
and reads ``Q(\delta)`` as a series in ``\delta``.
