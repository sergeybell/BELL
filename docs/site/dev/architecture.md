# Architecture and conventions

```@meta
CurrentModule = Bell
```

This page explains how the package is organised, which conventions every
module uses, and how data flows from a lattice file to the ring observables.
The other Developer's Guide pages describe the algorithms of each area; the
function reference is in the API Reference.

## Package layout

```
src/
  Bell.jl                 top module: include order + re-export
  coordinates.jl          Coordinates        particle, reference, laboratory maps
  elements.jl             Elements           element types (data only)
  structure.jl            Structure          Lattice, .bell parser/writer
  lie_algebra.jl          LieAlgebra         polynomial space, Poisson, Lie factors
  hamiltonians.jl         Hamiltonians       polynomial H and exact fields
  optim_converter.jl      OptimConverter     OptiM .opt -> .bell
  linear.jl               Linear             R = exp(L S Hess H), periodic Twiss
  symplectic.jl           Symplectic         one check max |J^T S J - S|
  exact_integrator.jl     ExactIntegrator    unexpanded H, rf_kick, renormalize
  tracking.jl             Tracking           Lie-factor tracking, Taylor map of H
  lie_flow.jl             LieFlow            flow to the order of the call
  map_readings.jl         MapReadings        closed orbit, chromaticity, compaction
  smooth.jl               Smooth             averaged (phi, delta) through eta1
  acceleration.jl         Acceleration       energy ramp
  fringe.jl               Fringe             HardEdge thin lenses, SoftEdge Body(-W)
  sextupole_families.jl   SextupoleFamilies  optional family placement
  plotting/               Plotting           PNG/GIF figures (Plots, GR)
  optimizers/             Optimizers         GridSearch, RandomStart stubs
  research/normal_form.jl NormalForm         research; not re-exported
  da/                     DiffAlgebra        truncated series of any order
  integrators/            Integrators        one file per ODE method
```

The include order in `Bell.jl` is the dependency order. A module only uses
modules included before it (via `using ..Name`). Exported names except those of
`NormalForm` are re-exported by `Bell`. Load the research module with
`using Bell.NormalForm`.

## Dependency graph

```
Coordinates ─► Elements ─► Structure ─► OptimConverter
     │             │           │
     └──► LieAlgebra ◄─────────┘
              │
         Hamiltonians ─► Linear ─► ExactIntegrator ─► Tracking ─► LieFlow ─► MapReadings
              │                                                                   │
              └─► Smooth, Acceleration, Fringe, SextupoleFamilies, Plotting, Optimizers
DiffAlgebra, Integrators (pure mathematics)
research/NormalForm (not re-exported)
```

## Conventions

### Coordinates

| index | canonical `z` | meaning | laboratory `u` |
|:--|:--|:--|:--|
| `index_x = 1` | x | horizontal offset [m] | x [m] |
| `index_px = 2` | px = Px/p0 | horizontal momentum | x' = px/ps |
| `index_y = 3` | y | vertical offset [m] | y [m] |
| `index_py = 4` | py = Py/p0 | vertical momentum | y' |
| `index_l = 5` | l = −v·Δt | path-length / time-of-flight variable [m] | φ = −h ω0 Δt [rad] |
| `index_delta = 6` | δ = Δp/p0 | momentum deviation | δ |

Here ``p_s = \sqrt{(1+\delta)^2 - p_x^2 - p_y^2}``. The pairs (x, px), (y, py),
(l, δ) are canonical with the independent variable s (arc length along the
reference orbit). See [`canonical_to_laboratory`](@ref) and [`laboratory_to_canonical`](@ref).

### Units

| quantity | unit |
|:--|:--|
| lengths, β, D | m |
| angles, phases | rad |
| energies, masses, p0·c | MeV |
| voltages | MV |
| k1 | 1/m² |
| k2 | 1/m³ |
| curvature h | 1/m |
| Bρ | T·m |
| tunes | 1 (fractional, unless stated) |

### Polynomials and truncation

All maps are vectors of six `TaylorSeries.TaylorN` polynomials in [`Z`](@ref).
Hamiltonians are expanded to degree [`JET_ORDER`](@ref) = 4, and maps are cut
at [`MAP_ORDER`](@ref) = 3. A map `M` is applied to a point with
`evaluate(M[i], z)`; composition is [`compose_maps`](@ref).

### Generic code

The tracking functions ([`apply_factor`](@ref), [`rf_kick`](@ref),
[`exact_track_element`](@ref), [`canonical_to_laboratory`](@ref), ...) are written
without type restrictions on `z`. With `Float64` input they track particles.
With `z0 .+ Z` they return the exact Taylor expansion of themselves. This gives
Jacobians, through [`jacobian_of_function`](@ref), and polynomial maps without
any extra code.

### Caches

Element maps are cached in module-level dictionaries keyed by the element value
and `ref.p0_MeV`:

- `Linear._HCACHE`: S·Hessian;
- `Tracking._MCACHE`: Lie maps;
- `Tracking._TMCACHE`: Taylor maps of `H`;
- `Fringe._ECACHE`: edge maps.

Elements are immutable, so equal elements share entries, and a lattice at a new
energy gets new maps automatically. The caches are never invalidated (elements
cannot change).

## Data flow

```
 .opt ──convert_optim──► .bell ──load_structure──► Lattice
                                                    │
     hamiltonian(e, ref) ◄──────────────────────────┤ per element
        │                                           │
        ├─► linear_matrix = exp(L S A) ─► ring_matrix ─► periodic_optics, optics_along
        ├─► lie_map = exp(-L:H:) ─► dragt_finn_factor ─► LieFactor per element
        │        └─► build_factors ─► ring_lie_factor (concatenate, BCH)
        │                 └─► factor_to_map ─► observables: dispersion_orders,
        │                                       path_expansion, chromaticity_orders
        ├─► taylor_map_element (integrator of hamiltonian_field) ─► ring_taylor_map
        └─► magnetic_field / cavity_field ─► exact_track_element ─► exact_track_turn, accelerate
```

The **one-turn Lie factor** ([`ring_lie_factor`](@ref)) is the primary ring
object. Linear optics, slip factors and chromaticities are readings of it
(see [`observables`](@ref)). A Taylor map of the same `H` is a check. Smooth
`(φ, δ)` is a separate averaged equation.

## Building this manual

```
julia --project=docs docs/make.jl          # HTML -> docs/build/index.html
```

Docstrings live next to the code in `src/*.jl`. The API Reference pages collect
them with `@autodocs` per module, so a new docstring appears there
automatically; `checkdocs = :exports` reports exported names without a docstring.
