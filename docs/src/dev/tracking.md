# Tracking and integrators

```@meta
CurrentModule = Bell
```

Production motion is ``dz/ds = S\nabla H`` of the element s-Hamiltonian.
[`element_flow`](@ref) (`src/maps/element_flow.jl`) is that integral.
[`exact_track_element`](@ref), [`element_da_map`](@ref), and
[`taylor_map_element`](@ref) call it for a magnet and for a thick
[`RFCavity`](@ref). A thin cavity is [`rf_kick`](@ref). The field equations
in `src/physics/fields.jl` are the comparison formulation. They are not the
step these functions take.

| module | file | role |
|:--|:--|:--|
| `ElementFlow` | `src/maps/element_flow.jl` | integrate ``dz/ds = S\nabla H`` |
| `FieldTrack` | `src/maps/da/field_track.jl` | one point: `exact_track_*`, `rf_kick`, `renormalize` |
| `Tracking` | `src/maps/lie/tracking.jl` | Lie factors, `track`, Taylor check `taylor_map_element` |
| `DATracking` | `src/maps/da/tracking.jl` | DA map from the same flow |
| `Acceleration` | `src/longitudinal_models/acceleration.jl` | ramp: track, then `renormalize` after each cavity |
| `Trackers` | `src/tracking/trackers.jl` | `FieldTracker`, `LieTracker`, `MapTracker` |
| `BeamState` | `src/tracking/state.jl` | `ParticleState`, `Beam` |

[`symplectic_error`](@ref) is one check: a numeric Jacobian, a `TaylorN` map,
a DA map ([`da_symplectic_error`](@ref)), or a tracking step at a point.

| route | function | what moves the state | cut | symplectic |
|:--|:--|:--|:--|:--|
| point tracking | [`exact_track_element`](@ref), [`exact_track_turn`](@ref) | [`element_flow`](@ref) of ``H``; thin RF is [`rf_kick`](@ref) | step of the ODE method | of the solver |
| Lie factors | [`track`](@ref), [`apply_factor`](@ref) | ``R``, ``f_3``, ``f_4`` | generators through degree 4 | yes for the applied flow |
| Lie series | [`element_lie_map`](@ref) | ``e^{-L{:}H{:}}z`` | map degree [`MAP_ORDER`](@ref) | up to that degree |
| Taylor check | [`taylor_map_element`](@ref) | same ``H`` through [`element_flow`](@ref) (default RK4) | degree 3 | up to truncation |
| DA map | [`element_da_map`](@ref) | same flow on a DA jet | order ``N`` | up to truncation |

## 1. Point tracking (`FieldTrack`)

[`exact_track_element`](@ref)`(z, e, ref, path_length; ds = 0.1, method = ImplicitYoshida(4))`:

- a [`Marker`](@ref), or a magnet of length 0, returns a copy of `z`;
- a thin [`RFCavity`](@ref) (`L = 0`) is [`rf_kick`](@ref). [`element_flow`](@ref) rejects a thin cavity;
- a thick cavity is [`element_flow`](@ref) of [`hamiltonian`](@ref)`(e, ref, path_length, z)` over `L`, with `nsteps = max(1, ceil(L / ds))`. `path_length` must be positive;
- a [`Drift`](@ref) is one step of `method`;
- any other magnet is the same flow of [`hamiltonian`](@ref)`(e, ref, z)`, with `nsteps = max(1, ceil(L / ds))`.

[`element_flow`](@ref) itself defaults to [`DormandPrince853`](@ref). When
`nsteps` and `ds` are both omitted, an adaptive method uses `tol`. The point
tracker always passes a step count, so its step is fixed. On numbers and on
`TaylorN`, each stage is the dual [`symplectic_gradient`](@ref) of ``H``. On a
DA jet the polynomial field is built once; see
[DA maps](da_tracking.md).

[`exact_track_turn`](@ref) walks `lattice.elements` in order and keeps the
reference it was given. Nothing on this route truncates ``H`` as a polynomial.
The cut is the step in ``s``.

The comparison right-hand side for a hard-edge magnet is [`element_field`](@ref)
in `src/physics/fields.jl`. [`element_rhs`](@ref) throws for a marker and for
an [`RFCavity`](@ref): a thick cavity is the flow of ``H``, not that field.

## 2. Acceleration ([`accelerate`](@ref))

The ramp uses a thin [`RFCavity`](@ref). Magnets are still the flow of ``H``.
The reference is rebuilt after the cavity, not inside the Hamiltonian.

```
C = path_length(lattice);  h = lattice.harmonic;  ref = injection reference
z_k = laboratory_to_canonical(u_k, ref, C, h)     # u_1 is the synchronous particle
for each turn:
    for e in lattice.elements:
        z_k = exact_track_element(z_k, e, ref, C; ds)
        (z_k, ref) = after_element!(z_k, ref, e)  # RFCavity → renormalize(..., 1)
    record ref and canonical_to_laboratory of every particle
    when ref.Ekin first reaches Ekin_target:
        set_cavity_phase(lattice, phi_flat)       # default π
    stop flat_turns after that
```

`after_element!` calls [`renormalize`](@ref) for an [`RFCavity`](@ref)
and returns the same state for every other element.

With ``(x_s, p_{xs}, y_s, p_{ys}, l_s, \delta_s)`` the synchronous particle
before the rebuild, and ``v(\delta)`` the velocity at that same ``\delta``
(not at ``\delta'``):

```math
\begin{aligned}
p_0' &= p_0(1+\delta_s),\qquad
E_\text{kin}' = \sqrt{p_0'^2 + m^2} - m,\\
\delta' &= \frac{1+\delta}{1+\delta_s} - 1,\qquad
p_x' = \frac{p_x - p_{xs}}{1+\delta_s},\qquad
p_y' = \frac{p_y - p_{ys}}{1+\delta_s},\\
x' &= x - x_s,\qquad y' = y - y_s,\\
\Delta t' &= -\frac{l}{v(\delta)} + \frac{l_s}{v(\delta_s)},\qquad
l' = -v(\delta)\,\Delta t' .
\end{aligned}
```

``p_x/(1+\delta_s)`` is the adiabatic damping. Normalised magnet strengths are
kept: the magnets are assumed to ramp with ``p_0``. The next element is
tracked with the new reference, so ``\beta_0``, ``\gamma_0``, the revolution
frequency, and the cavity scale ``ZV/p_0`` follow the beam.

[`turn_snapshots`](@ref) runs one such turn and records laboratory coordinates
after every element.

Lie-factor [`track`](@ref) does not ramp. A thin kick changes ``\delta`` of a
particle and leaves ``p_0`` fixed.

## 3. Lie-factor tracking (`src/maps/lie/tracking.jl`)

[`build_factors`](@ref) stores one [`LieFactor`](@ref) ``(R, f_3, f_4)`` per
element, from [`element_lie_map`](@ref). A thin cavity is the Lie map of
[`rf_potential`](@ref). A thick cavity is ``\exp(-L{:}H{:})``. [`track`](@ref)
only applies those factors:

```
apply_factor(F, z):
    w = z
    if f4 ≠ 0:  w = flow of −S∇f4 over t ∈ [0, 1]    # ImplicitMidpoint, 2 steps
    w = flow of −S∇f3 over t ∈ [0, 1]                # ImplicitMidpoint, 2 steps
    return R · w
```

Each flow is Hamiltonian, so the applied factor is symplectic. The generators
are truncated: ``f_3`` and ``f_4`` represent the element through degree 3 of
the map. The polynomial `factor_to_map(F)` is symplectic only up to
[`MAP_ORDER`](@ref).

### Thin kick and thick cavity

[`rf_kick`](@ref) accepts numbers or DA. With path length ``C`` and
``v_0 = \beta_0 c``:

```math
\Delta t = -\frac{l}{v(\delta)},\qquad
\varphi = -h\omega_0\Delta t = \frac{2\pi h\, v_0}{C}\,\frac{l}{v(\delta)},
```

```math
E' = \sqrt{p_0^2(1+\delta)^2 + m^2} + ZV\sin(\varphi_s + \varphi),\qquad
\delta' = \frac{\sqrt{E'^2 - m^2}}{p_0} - 1,\qquad
l' = -v(\delta')\,\Delta t .
```

``x, p_x, y, p_y`` are unchanged, and so is the reference momentum. An energy
kick at fixed arrival time is canonical in ``(t, -E)``, and the change of
variables to ``(l, \delta)`` is canonical, so the kick is symplectic.

A thick cavity is the flow of its s-Hamiltonian: ``\exp(-L{:}H{:})`` on the
Lie path, and [`element_flow`](@ref) of the same ``H`` on the point tracker
and the DA map. The Hamiltonian, and the ramp that still rebuilds the
reference with [`renormalize`](@ref) after each thin kick, are written out in
[RF Hamiltonian and acceleration](rf_acceleration.md).

## 4. Taylor check

[`taylor_map_element`](@ref) integrates ``dz/ds = S\nabla H`` in `TaylorN`,
default method RK4, and truncates to [`MAP_ORDER`](@ref). A magnet and a thick
cavity go through [`element_flow`](@ref). A marker is the identity. A thin
cavity integrates the vector field of [`rf_potential`](@ref) over a unit
interval. [`ring_taylor_map`](@ref) composes the element maps. This checks the
Lie product. It is not the tracking route.

## 5. What each route keeps exact

| quantity | point tracking | Lie-factor tracking | DA map |
|:--|:--|:--|:--|
| ``p_s``, ``v(\delta)``, [`energy_term`](@ref) | square roots in ``H`` | ``H`` to degree 4 | series to order ``N`` |
| RF | thin [`rf_kick`](@ref); thick flow of ``H`` | thin ``\exp(-{:}G{:})``; thick ``\exp(-L{:}H{:})`` | same split as the point tracker |
| reference energy | [`accelerate`](@ref) updates it after each thin cavity | fixed | fixed |
| cut | ODE step (`ImplicitYoshida(4)` by default on `exact_track_*`) | generators ``f_3``, ``f_4`` | order ``N`` |
