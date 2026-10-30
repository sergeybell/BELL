# Architecture and conventions

```@meta
CurrentModule = Bell
```

This page explains how the package is organised, which conventions every
module uses, and how data flows from a lattice file to the ring observables.
The rest of the Programmer Manual has one page per package area. Each page
gives a hand-written explanation of the algorithms, then the reference of every
type and function.

## Package layout

```
src/
  Bell.jl                 top module: include order + re-export of all names
  coordinates.jl          Coordinates        particle, reference, coordinate maps
  elements.jl             Elements           element types (data only)
  structure.jl            Structure          Lattice, .bell parser/writer, editing
  lie_algebra.jl          LieAlgebra         polynomial space, Poisson, Lie factors
  hamiltonians.jl         Hamiltonians       H of every element (single physics source)
  optim_converter.jl      OptimConverter     OptiM .opt -> .bell
  linear.jl               Linear             R = exp(L S A), Twiss, dispersion
  symplectic.jl           Symplectic         Jacobians, symplecticity error
  tracking.jl             Tracking           Lie-factor tracking, thin RF kick
  higher_order.jl         HigherOrder        readings of the one-turn map
  differential_algebra.jl DifferentialAlgebra RK4-in-TaylorN maps (cross-check)
  longitudinal.jl         Longitudinal       smooth (phi, delta) motion, bucket
  exact_integrator.jl     ExactIntegrator    exact-H symplectic integrator, cavity
  acceleration.jl         Acceleration       energy ramp, adiabatic damping
  normal_form.jl          NormalForm         6D normal forms with RF
  fringe.jl               Fringe             hard/soft edge models
  plotting.jl             Plotting           PNG/GIF figures (Plots, GR)
  optimizers.jl           Optimizers         optional sextupole fitting
```

The include order in `Bell.jl` is the dependency order. A module only uses
modules included before it (via `using ..Name`). Everything a module exports is
re-exported by `Bell`, so all names on these pages can be used directly after
`using Bell`.

## Dependency graph

```
Coordinates ─► Elements ─► Structure ─► OptimConverter
     │             │           │
     └──► LieAlgebra ◄─────────┘
              │
         Hamiltonians ─► Linear ─► Tracking ─► HigherOrder ─► DifferentialAlgebra
              │                                   │                │
              │                       Longitudinal│                │
              └─► ExactIntegrator ─► Acceleration  └──► NormalForm ◄┘
                                                   Fringe (Linear, Tracking, Hamiltonians)
                                    Plotting, Optimizers (leaf modules)
```

## Conventions

### Coordinates

| index | canonical `z` | meaning | user `u` |
|:--|:--|:--|:--|
| `IX = 1` | x | horizontal offset [m] | x [m] |
| `IPX = 2` | px = Px/p0 | horizontal momentum | x' = px/ps |
| `IY = 3` | y | vertical offset [m] | y [m] |
| `IPY = 4` | py = Py/p0 | vertical momentum | y' |
| `IL = 5` | l = −v·Δt | path-length / time-of-flight variable [m] | φ = −h ω0 Δt [rad] |
| `IDELTA = 6` | δ = Δp/p0 | momentum deviation | δ |

Here ``p_s = \sqrt{(1+\delta)^2 - p_x^2 - p_y^2}``. The pairs (x, px), (y, py),
(l, δ) are canonical with the independent variable s (arc length along the
reference orbit). See [`canonical_to_user`](@ref) and [`user_to_canonical`](@ref).

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
[`exact_track_element`](@ref), [`canonical_to_user`](@ref), ...) are written
without type restrictions on `z`. With `Float64` input they track particles.
With `z0 .+ Z` they return the exact Taylor expansion of themselves. This gives
Jacobians, through [`jacobian_of_function`](@ref), and polynomial maps without
any extra code.

### Caches

Element maps are cached in module-level dictionaries keyed by the element value
and `ref.p0_MeV`:

- `Linear._HCACHE`: S·Hessian;
- `Tracking._MCACHE`: Lie maps;
- `DifferentialAlgebra._TMCACHE`: RK4 maps;
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
        ├─► taylor_map_element (RK4, check) ─► ring_taylor_map
        └─► exact field ─► exact_track_element ─► exact_track_turn, accelerate
```

The **one-turn Lie factor** ([`ring_lie_factor`](@ref)) is the primary ring
object. Linear optics, slip factors and chromaticities are readings of it
(see [`observables`](@ref)). The RK4 Taylor map and the exact integrator are
independent cross-checks.

## Building this manual

```
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'   # once
julia --project=docs docs/make.jl          # HTML -> docs/build/index.html
julia --project=docs docs/make_pdf.jl programmer   # old PDF programmer overview
```

Docstrings live next to the code in `src/*.jl`. A new exported name needs a
docstring and an entry in a `@docs` block on the page of its area. `make.jl`
runs with `checkdocs = :exports` and reports any exported name that is missing.

## Module docstring

```@docs
Bell
```
