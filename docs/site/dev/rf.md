# RF, acceleration and longitudinal dynamics

```@meta
CurrentModule = Bell
```

Modules: `ExactIntegrator`, `Acceleration` (`src/acceleration.jl`),
`Longitudinal` (`src/longitudinal.jl`). The cavity Hamiltonian, the exact kick
and the turn-by-turn ramp algorithm are on the [Tracking](tracking.md) page; this
page states precisely **what is and what is not truncated** in the longitudinal
plane, and documents the analysis-only longitudinal model.

## Is the longitudinal motion truncated in 6D tracking? — No.

| where | longitudinal treatment | truncated? |
|:--|:--|:--|
| [`exact_track_element`](@ref) / [`exact_track_turn`](@ref) (magnets) | ``l' = v(\delta)/v_0 - (1+hx)(1+\delta)/p_s``, exact roots | **no** |
| [`exact_track_element`](@ref) (`Cavity`) | ``H_\text{cav} = e(\delta) - p_s + A\cos\Phi``, exact ``\sin\Phi``, exact ``v(\delta)``, s- and t-dependent | **no** |
| [`accelerate`](@ref), [`turn_snapshots`](@ref) | exact integrator + [`renormalize`](@ref) after each cavity | **no** |
| [`rf_kick`](@ref) (thin cavity, numbers) | exact energy kick and exact ``E \leftrightarrow \delta``, ``l \leftrightarrow \Delta t`` | **no** |
| [`track`](@ref) with Lie factors | magnets: generators of degree ≤ 4 (map exact to degree 3) — the **same** truncation in all six variables, the longitudinal plane is not treated differently | yes, by design of the map route |
| [`rf_kick_taylor`](@ref), [`one_turn_map_rf`](@ref), normal forms | polynomial of degree 3 | yes (analysis) |
| `Longitudinal` module | smooth model, ``\eta_0\delta + \eta_1\delta^2`` | yes, 2nd order (analysis only) |

So the 6D tracking used for the RF and acceleration studies (examples 14 and 16)
is **full 6D, non-truncated**: the element-by-element exact integrator with the
finite cavity. The smooth "longitudinal model" is **not** a tracking engine: no
tracking routine calls it. It is an analysis tool fed with coefficients read from
the one-turn map, used to compare tracking with the equations (bucket height,
separatrix, ``Q_s``), and to pick initial particles inside the bucket in
example 14. ``\beta_s`` of the map-based normal form [`normal_form_rf`](@ref) uses
the small-amplitude [`synchrotron_tune`](@ref).

Assumptions that *are* made in tracking (stated in the code): the cavity field is
purely longitudinal on axis (no transverse RF fields); magnets ramp with ``p_0``
(fixed ``k_1, k_2, h``); the RF frequency follows ``h f_0(p_0)``.

## The smooth longitudinal model (analysis)

With turn number ``n``, ``\varphi = -h\omega_0\Delta t`` and the slip factors
``\eta_0, \eta_1`` read from the one-turn map ([`path_expansion`](@ref)):

```math
\frac{d\varphi}{dn} = -2\pi h\,(\eta_0\delta + \eta_1\delta^2),\qquad
\frac{d\delta}{dn} = A\,\big(\sin(\varphi_s+\varphi) - \sin\varphi_s\big),\qquad A = \frac{ZVT}{\beta_0^2E_0},
```

```math
H(\varphi,\delta) = -2\pi h\Big(\eta_0\frac{\delta^2}{2} + \eta_1\frac{\delta^3}{3}\Big)
 + A\big(\cos(\varphi_s+\varphi) + \varphi\sin\varphi_s\big)\qquad(\texttt{long\_H}).
```

* ``\eta_0 = \alpha_0 - 1/\gamma_0^2``,
  ``\eta_1 = \alpha_1 - \alpha_0/\gamma_0^2 + \tfrac{3\beta_0^2}{2\gamma_0^2} + \tfrac{1}{\gamma_0^4}``
  ([`LongitudinalParams`](@ref)).
* Small-amplitude tune ``Q_s = \sqrt{h\eta_0 ZVT\cos\varphi_s/(2\pi\beta_0^2E_0)}``,
  stable iff ``\eta_0\cos\varphi_s > 0``.
* Fixed points: stable ``\varphi = 0``, unstable ``\varphi_u = \pi - 2\varphi_s``,
  second (α-bucket) momentum ``\delta_2 = -\eta_0/\eta_1``.
* Separatrix ``H = H(\varphi_u, 0)``; bucket half-heights = roots of
  ``H(0,\delta) = H_\text{sep}`` by 200 bisection steps in ``|\delta| < |\eta_0/\eta_1|/2``.

This model is truncated at second order in ``\delta`` *deliberately*, as an
analysis object to compare with the exact tracking.

## Acceleration

The full algorithm (order of operations per turn, where the kick is applied, the
kick formula, and which reference quantities are recomputed) is in
[Tracking → One turn with acceleration](tracking.md). [`measure_tune`](@ref)
extracts tunes from turn-by-turn data by fitting
``x_{n+1}+x_{n-1} = 2\cos(2\pi Q)\,x_n + b``, so an unknown fixed point does not bias ``Q``.
