# Architecture and conventions

```@meta
CurrentModule = Bell
```

The model of an element is its s-Hamiltonian `H`. The equations are
`dz/ds = S ∇H`. [`element_flow`](@ref) (`src/maps/element_flow.jl`) integrates
that equation for a magnet, a thick [`RFCavity`](@ref), or an explicit
`H(z, s)`. A thin cavity stays [`rf_kick`](@ref). The field equations in
`src/physics/fields.jl` are a second formulation, kept for comparison and for
tests. They are not the map `element_da_map` or `exact_track_element` builds.

Two computations stay in the tree and do not call each other. The Lie path
(`src/maps/lie/`) applies Taylor series to `H`. The DA path
(`src/maps/da/`) integrates the same `H` through `element_flow`. Comparison
is done from the outside.

## Layers

`src/Bell.jl` includes files in this order, except where a `using` forces a
file out of sequence (`Hamiltonians` after `LieAlgebra`, `FringeMaps` after
`MapReadings`). `NormalForm` is included last and is not re-exported
(`using Bell.NormalForm`).

```
math → core → io → physics → maps → analysis → tracking → longitudinal_models → optimizers → research
```

| layer | where | role |
|:--|:--|:--|
| math | `src/math` | `DiffAlgebra` (truncated series and `Dual`), one file per ODE method in `src/math/integrators`. No elements and no Hamiltonian. |
| core | `src/core` | coordinates, `RefParticle`, element data, `Lattice`, `symplectic_error`, caches, plot names |
| io | `src/io` | `OptimConverter`: a `.opt` file → a `.bell` structure file |
| physics | `src/physics` | one s-Hamiltonian per element, field equations kept for comparison, fringe model (no maps) |
| maps | `src/maps` | `element_flow`; `lie/` (Taylor factors); `da/` (point track and DA map of that flow); `lie/fringe_maps.jl` (edge maps) |
| analysis | `src/analysis` | readings of a map that is already built; sextupole families |
| tracking | `src/tracking` | `ParticleState`, `Beam`, and trackers that call an existing engine |
| longitudinal_models | `src/longitudinal_models` | averaged `(φ, δ)` model and the energy ramp, kept beside the `RFCavity` Hamiltonian |
| optimizers | `src/optimizers` | derivative-free fit of a scalar `f(x) -> Real` |
| research | `src/maps/lie/normal_form.jl` | `NormalForm`. Not part of `using Bell`. |

[`hamiltonian`](@ref)`(e, ref, z)` evaluates the same s-Hamiltonian on numbers,
`TaylorN`, and DA. [`hamiltonian`](@ref)`(e, ref)` returns the Lie-path
polynomial (degree [`JET_ORDER`](@ref)). A thick [`RFCavity`](@ref) also takes
the ring path length. A thin cavity has no continuous Hamiltonian: the map is
[`rf_kick`](@ref), and the Lie generator is [`rf_potential`](@ref).

### `element_flow`, `Dual`, and `deriv`

[`element_flow`](@ref)`(e, ref, z)` is the autonomous element. On numbers and
on `TaylorN`, each stage evaluates [`symplectic_gradient`](@ref): the
forward-mode [`Dual`](@ref) partial of the scalar `H`, then the pairing `S`.
On a DA jet that same method builds the polynomial vector field once with
[`vector_field_on_identity`](@ref), in a space one order above the map,
because [`deriv`](@ref) drops the top degree of a jet that is already
truncated. The integrator then evaluates the polynomial. The dual partial and
`deriv` agree through order `N - 1`, and at order `N` when `H` is a polynomial
of degree at most `N`.

[`element_flow`](@ref)`(H, z; s0, s1)` is the s-dependent scalar, used by a
soft edge. Every stage evaluates `symplectic_gradient` of `w -> H(w, s)`,
including on DA.

`src/maps/lie/` is the Lie computation of `H`. `src/maps/da/` is the point
tracker and the DA map of the same flow. `src/physics/fields.jl` is the check
that [`element_field`](@ref) matches `S ∇H` (`test/test_physics_consistency.jl`).

Plots is a weak dependency. The names `plot_twiss`, `plot_dispersion2`,
`plot_chromatic`, `plot_acceleration`, `separatrix!`, `gif_ramp`, and
`gif_elements` are declared in `src/core/plotting.jl`. The methods live in
`ext/BellPlotsExt.jl` and appear after `using Plots`.

People follow `CONTRIBUTING.md`: one issue, one branch, one pull request.
Agents follow `AGENTS.md`. A user-visible change is recorded in the root
`CHANGELOG.md` under `[Unreleased]`, with the pull request number. `version`
in `Project.toml` is bumped, and that section is renamed, only in a release
pull request.

## Package layout

```
src/
  Bell.jl                 facade: include order and re-export
  math/
    DiffAlgebra.jl        truncated series
    dual.jl               Dual, symplectic_gradient (included into DiffAlgebra)
    lie_series.jl         canonical pairs, Hamiltonian field (included into DiffAlgebra)
    integrators/          one file per ODE method; Integrators.jl is the driver
  core/
    coordinates.jl        particle, reference, laboratory ↔ canonical
    elements.jl           element data (RFCavity thin or thick)
    structure.jl          Lattice, .bell reader/writer, path_length, magnetic_lattice
    symplectic.jl         symplectic_error
    caches.jl             clear_caches! registry (included into Bell)
    plotting.jl           plot_* names (included into Bell)
  io/
    optim_converter.jl    .opt file → .bell file
  physics/
    hamiltonians.jl       hamiltonian(e, ref, z); rf_potential
    fields.jl             element_field, the comparison right-hand side
    fringe.jl             edge types, profiles, fringe_hamiltonian, FringeModel
  maps/
    element_flow.jl       dz/ds = S ∇H
    abstract.jl           AbstractMap, LinearMap, DAMap
    lie/
      lie_algebra.jl      Poisson bracket, Lie factors
      linear.jl           R = exp(L S Hess H)
      tracking.jl         Lie-factor tracking; taylor_map_element
      fringe_maps.jl      edge maps, cache, fringe_observables
      normal_form.jl      research (using Bell.NormalForm)
    da/
      field_track.jl      exact_track_element, rf_kick, renormalize
      tracking.jl         element_da_map from element_flow
  analysis/
    map_readings.jl       closed orbit, dispersion, chromaticity, compaction
    sextupole_families.jl family placement and fit
  tracking/
    state.jl              ParticleState, Beam
    trackers.jl           FieldTracker, LieTracker, MapTracker
  longitudinal_models/
    smooth.jl             averaged (φ, δ)
    acceleration.jl       ramp: exact_track_element, rf_kick, renormalize
  optimizers/
    Optimizers.jl         derivative-free fitting
ext/
  BellPlotsExt.jl         plot_*; loaded with Plots
```

## Dependency graph

```
math: DiffAlgebra, Integrators
        │
core: Coordinates ─► Elements ─► Structure ─► io: OptimConverter
                      │              │
                      │              └─► physics: Hamiltonians, Fields, Fringe
                      │                        │
                      │                        ▼
                      │                 maps: ElementFlow
                      │                   ├─► lie: Linear, Tracking, FringeMaps
                      │                   └─► da: FieldTrack, DATracking
                      └──► maps/lie: LieAlgebra
analysis: MapReadings, SextupoleFamilies
tracking: BeamState, Trackers   (call FieldTrack and Lie Tracking)
longitudinal_models: Smooth, Acceleration
optimizers
research: NormalForm (`using Bell.NormalForm`)
ext/BellPlotsExt (Plots; weak dependency)
```

The full `using` graph, including the two rightward edges, is `src/README.md`.

## Conventions

### Polynomials (Lie path)

Hamiltonians on the Lie path are expanded to degree [`JET_ORDER`](@ref) = 4;
maps are cut at [`MAP_ORDER`](@ref) = 3. The current factors carry generators
`f3` and `f4`. A Lie factor applied by its flow stays symplectic at any
truncation. A Taylor map of the same `H` is symplectic only up to its
truncation. Arbitrary order lives on the DiffAlgebra path.

### RF

A zero-length [`RFCavity`](@ref) is [`rf_kick`](@ref) (closed form, numbers or
DA). [`element_flow`](@ref) rejects it. A cavity with `L > 0` has an
s-Hamiltonian. The Lie path is `exp(-L:H:)`. The point tracker, the DA map,
and [`taylor_map_element`](@ref) integrate that `H` with [`element_flow`](@ref).
On a ramp, [`renormalize`](@ref) updates `p0` after the thin synchronous kick.
The reference is not changed inside the cavity.

## Data flow

```
 .opt ──convert_optim──► .bell ──load_structure──► Lattice
                                                    │
     H = hamiltonian(e, ref[, path_length], z)
     Lie:     H ─► lie_map ─► Dragt–Finn ─► track
     numbers / DA: H ─► element_flow ─► exact_track_element / element_da_map
     thin RF: rf_kick (element_flow does not take it)
     check:   element_field  vs  S ∇H     (physics/fields.jl, tests only)
     ramp:    exact_track_element + renormalize after each thin RFCavity
```
