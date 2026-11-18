# Architecture and conventions

```@meta
CurrentModule = Bell
```

Dependencies point only toward the layers on the left. Two computations stay
in the tree: the Lie path on TaylorSeries, and the magnet flow on DiffAlgebra.
Both take the s-Hamiltonian. Comparison with the field equations is done from
the outside.

## Layers

```
math → core → io → physics → maps/{lie,da} → analysis → tracking interfaces → effects (future)
```

| layer | where | role |
|:--|:--|:--|
| math | `src/math` | `DiffAlgebra`, canonical Lie-series helpers, one file per ODE method in `Integrators`. Pure mathematics. |
| core | `src/core` | particle and reference, element types, `Lattice` and the `.bell` file, symplecticity, sextupole families, caches, plot names |
| io | `src/io` | `OptimConverter`: OptiM `.opt` → `.bell` |
| physics | `src/physics` | one s-Hamiltonian per element, comparison field ODEs `dz/ds`, fringe models |
| maps/lie | `src/maps/lie` | TaylorSeries Lie algebra, linear matrices, Lie-factor tracking, research `NormalForm` |
| maps | `src/maps/element_flow.jl` | magnet flow `dz/ds = S ∇H` on numbers, Taylor series, or DA |
| maps/da | `src/maps/da` | `FieldTrack` (Float64) and `DATracking` (DiffAlgebra, any order N), both from that flow for magnets |
| analysis | `src/analysis` | readings of a one-turn object: closed orbit, chromaticity, momentum compaction |
| tracking interfaces | public map and track names | `element_lie_map`, `taylor_map_element` (`maps/lie/tracking.jl`), `element_da_map` (`maps/da/tracking.jl`), `exact_track_element` (`maps/da/field_track.jl`), `track`. Callers use this boundary over both engines. |
| effects | — | future. No directory yet. New physics that is not an element Hamiltonian lands here, and calls the tracking interfaces. |

`src/longitudinal_models` (averaged Smooth model and the acceleration ramp) and
`src/optimizers` are included after `analysis`. They are present code beside
the chain above.

[`hamiltonian`](@ref)`(e, ref, z)` evaluates the same s-Hamiltonian on numbers,
`TaylorN` and DA. [`hamiltonian`](@ref)`(e, ref)` returns the Lie-path
polynomial (degree [`JET_ORDER`](@ref)). A thick [`RFCavity`](@ref) also takes
the ring path length; a thin cavity is the kick [`rf_kick`](@ref).

Plots is a weak dependency. The names `plot_twiss`, `plot_dispersion2`,
`plot_chromatic`, `plot_acceleration`, `separatrix!`, `gif_ramp` and
`gif_elements` are declared in `src/core/plotting.jl`. The methods live in the
package extension `ext/BellPlotsExt.jl` and appear after `using Plots`.

People follow `CONTRIBUTING.md`: one Issue, one branch, one pull request.
Agents follow `AGENTS.md`. A user-visible change is recorded in the root
`CHANGELOG.md` under `[Unreleased]`, with the pull request number. `version`
in `Project.toml` is bumped, and that section is renamed, only in a release
pull request.

## Package layout

```
src/
  Bell.jl                 top module: include order + re-export
  math/                   pure mathematics
    DiffAlgebra.jl        DiffAlgebra        truncated series (pure maths)
    lie_series.jl         (in DiffAlgebra)   canonical pairs, Hamiltonian field, DA symplectic defect
    integrators/          Integrators        one file per ODE method
  core/                   shared by both stacks
    coordinates.jl        Coordinates        particle, reference, laboratory maps
    elements.jl           Elements           element types (RFCavity thin or thick)
    structure.jl          Structure          Lattice, .bell parser/writer, magnetic_lattice
    symplectic.jl         Symplectic         symplecticity checks
    sextupole_families.jl SextupoleFamilies
    caches.jl             (in Bell)          clear_caches!
    plotting.jl           (in Bell)          plot_* names; methods in BellPlotsExt
  io/
    optim_converter.jl    OptimConverter     OptiM -> .bell
  physics/
    hamiltonians.jl       Hamiltonians       hamiltonian(e, ref, z); RF: rf_potential
    fields.jl             Fields             dz/ds from multipole fields
    fringe.jl             Fringe             HardEdge / SoftEdge
  maps/
    element_flow.jl       ElementFlow        dz/ds = S ∇H for a magnet
    lie/                  TaylorSeries Lie path
      lie_algebra.jl      LieAlgebra         polynomial space, Poisson, Lie factors
      linear.jl           Linear             linear matrices / optics
      tracking.jl         Tracking           Lie-factor tracking, Taylor map of H
      normal_form.jl      NormalForm         research (using Bell.NormalForm)
    da/                   point track and DA maps from H
      field_track.jl      FieldTrack         Float64 magnet flow, rf_kick, renormalize
      tracking.jl         DATracking         magnet maps on DA (any order N)
  analysis/
    map_readings.jl       MapReadings        closed orbit, chromaticity, compaction
  longitudinal_models/    approximate models kept for comparison
    smooth.jl             Smooth             averaged (phi, delta)
    acceleration.jl       Acceleration       ramp: rf_kick + renormalize
  optimizers/
    Optimizers.jl         Optimizers         derivative-free fitting
ext/
  BellPlotsExt.jl         BellPlotsExt       plot_*; loaded with Plots
```

## Dependency graph

```
math: DiffAlgebra, Integrators
        │
core: Coordinates ─► Elements ─► Structure ─► io: OptimConverter
                      │              │
                      │              ├─► physics: Hamiltonians ─► maps: ElementFlow ─► FieldTrack, DATracking
                      │              │                         └─► maps/lie: Linear ─► Tracking
                      │              └─► physics: Fields (comparison; thick-RF drift halves)
                      └──► maps/lie: LieAlgebra
analysis: MapReadings  (reads either one-turn object)
tracking interfaces: element_lie_map, element_da_map, exact_track_element, track
effects: future
research/NormalForm (`using Bell.NormalForm`)
ext/BellPlotsExt (Plots; weak dependency)
```

## Conventions

### Polynomials (Lie path)

Hamiltonians on the Lie path are expanded to degree [`JET_ORDER`](@ref) = 4;
maps are cut at [`MAP_ORDER`](@ref) = 3. Arbitrary order lives on the
DiffAlgebra path.

### RF

A zero-length [`RFCavity`](@ref) is [`rf_kick`](@ref) (closed form, numbers or
DA). A cavity with `L > 0` has an s-Hamiltonian: `exp(-L:H:)` on the Lie path
and drift/kick/drift on the field / DA path. On a ramp, [`renormalize`](@ref)
updates `p0` after the synchronous kick.

## Data flow

```
 .opt ──convert_optim──► .bell ──load_structure──► Lattice
                                                    │
     Lie / TaylorSeries: hamiltonian(e, ref[, z]) ─► lie_map ─► Dragt–Finn ─► track
     Magnet flow:        hamiltonian(e, ref, z) ─► element_flow ─► FieldTrack / DATracking
     tracking interfaces: element_*_map, exact_track_element, track
     RF: rf_kick (thin) or drift/kick/drift; accelerate = FieldTrack + renormalize
     effects: future, calls the tracking interfaces
```
