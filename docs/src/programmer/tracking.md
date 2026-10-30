# Tracking and integrators

```@meta
CurrentModule = Bell
```

Modules: `Tracking` (`src/tracking.jl`), `ExactIntegrator`
(`src/exact_integrator.jl`).

Bell has four integrators. All of them work on `Float64` and on `TaylorN` input.

| integrator | where | order | symplectic | used for |
|:--|:--|:--|:--|:--|
| Lie series ``e^{-L{:}H{:}}`` | [`lie_map`](@ref) | exact to MAP_ORDER | to MAP_ORDER (polynomial) | element maps |
| implicit midpoint flow of f3, f4 | `Tracking._flow` | 2 (in t ∈ [0,1], 2 steps) | yes | particle tracking with Lie factors |
| RK4 in TaylorN | [`taylor_map_element`](@ref) | 4 | no | cross-check maps |
| Yoshida-4 implicit midpoint, exact H | [`exact_track_element`](@ref) | 4 | yes | reference tracking, RF cavity, acceleration |

## Lie-factor tracking

[`build_factors`](@ref) turns each element into one [`LieFactor`](@ref).
[`track`](@ref) then applies them in order. A factor is applied to a point by
[`apply_factor`](@ref), in three steps:

1. **P4 = exp(:f4:)**: the time-1 flow of the vector field ``-S\nabla f_4``
   (`F.field4`).
2. **N3 = exp(:f3:)**: the time-1 flow of ``-S\nabla f_3`` (`F.field`).
3. **R**: multiplication by the 6×6 matrix.

Each flow is integrated by the implicit midpoint rule (`_flow`):

```math
w_{n+1} = w_n + h\,F\!\left(\tfrac{w_n + w_{n+1}}{2}\right),\qquad h = 1/n_\text{steps},
```

The implicit equation is solved by fixed-point iteration from an explicit Euler
guess, to 1e-17 for numbers or with a fixed iteration count for polynomials. The
midpoint rule is symplectic for every Hamiltonian vector field, so the tracked
map is symplectic to round-off. It is *not* equal to the truncated polynomial
`factor_to_map(F)`; the two differ at degree > MAP_ORDER.

```@docs
build_factors
ring_lie_factor
element_lie_map
apply_factor
Tracking.apply_nonlinear
Tracking._flow
track
```

## Thin RF kick

The thin cavity is an energy kick at fixed arrival time, and is symplectic by
construction (see [`rf_kick`](@ref)). Its polynomial version
[`rf_kick_taylor`](@ref) is used to build one-turn maps with RF.

```@docs
rf_kick
```

## Exact-Hamiltonian integrator

`ExactIntegrator` integrates Hamilton's equations of the **non-expanded**
Hamiltonians:

```math
H_\text{mag} = e(\delta) - (1+hx)p_s + hx + \tfrac{h^2x^2}{2} + \tfrac{k_1}{2}(x^2-y^2) + \tfrac{k_2}{6}(x^3-3xy^2),
```

```math
H_\text{cav} = e(\delta) - p_s + A\cos\Phi,\quad
\Phi = \frac{\omega(s-L/2)}{v_0} - \frac{\omega\, l}{v(\delta)} - \phi_s,\quad
A = \frac{ZV}{p_0c}\frac{c}{\omega L}.
```

The cavity Hamiltonian follows from the on-axis field ``E_s = -\partial_t A_s``.
It gives the energy gain ``ZVT\sin(\phi_s+\phi)`` with transit factor
``T = \sin(\theta/2)/(\theta/2)``, ``\theta = \omega L/v_0``.

**Step by step** (`_integrate`):

1. Split the element into n steps of length h, where n = ⌈L/ds⌉ for magnets,
   one step for a drift, and ⌈L/0.02 m⌉ (at least 4) for a cavity.
2. Each step is Yoshida's triple jump of implicit midpoint substeps
   with lengths W1·h, W0·h, W1·h, where ``W_1 = 1/(2-2^{1/3})`` and
   ``W_0 = 1-2W_1``. The s value is advanced inside each substep, so
   s-dependent fields (the cavity) are handled correctly.
3. Each implicit midpoint substep `_midpoint` solves
   ``z' = z + h F((z+z')/2, s+h/2)`` by fixed-point iteration (≤ 40 iterations,
   tolerance 1e-16).

The result is a 4th-order symplectic integrator. It is the reference for the
truncated maps, and the only integrator that handles the time-dependent cavity
and energy ramps ([`renormalize`](@ref)).

```@docs
exact_track_element
exact_track_turn
cavity_parameters
transit_factor
renormalize
ExactIntegrator.W1
ExactIntegrator._field_mag
ExactIntegrator._field_cav
ExactIntegrator._midpoint
ExactIntegrator._integrate
```
