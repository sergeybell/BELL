# ODE integrators (`Bell.Integrators`)

The methods live under `src/math/integrators/`, one file per family, with the
driver in `Integrators.jl`. They are generic over `Float64` and
[`DA`](@ref Bell.DiffAlgebra.DA). The method is an argument. The module holds
no element and no Hamiltonian.

| file | methods |
|:--|:--|
| `src/math/integrators/rk4.jl` | `RK4` |
| `src/math/integrators/dormand_prince.jl` | `DormandPrince5`, `DormandPrince853` |
| `src/math/integrators/gragg_bulirsch_stoer.jl` | `GraggBulirschStoer` |
| `src/math/integrators/implicit_midpoint.jl` | `ImplicitMidpoint`, `ImplicitYoshida` |
| `src/math/integrators/yoshida.jl` | `Leapfrog`, `Yoshida`, `ForestRuth`, `SplitFlows` |

```julia
using Bell.DiffAlgebra, Bell.Integrators
r = integrate(f, z0, s0, s1; method = DormandPrince853(), tol = 1e-13)   # adaptive
r = integrate(f, z0, s0, s1; method = RK4(), nsteps = 200)                # fixed step
r = integrate(SplitFlows(drift, kick), z0, s0, s1; method = Yoshida(6), nsteps = 100)
r.z, r.nsteps, r.nrejected, r.nfev
```

| method | order | type |
|:--|:--|:--|
| `RK4()` | 4 | fixed step |
| `DormandPrince5()` | 5(4) | embedded, adaptive |
| `DormandPrince853()` | 8(5,3) | embedded, adaptive |
| `GraggBulirschStoer(k)` | ``2k`` | extrapolated midpoint, adaptive, any even order |
| `ImplicitMidpoint()` | 2 | fixed step; symplectic for a Hamiltonian field |
| `ImplicitYoshida(4)` (any even order ``\ge 2``) | 4, 6, … | Yoshida composition of the implicit midpoint |
| `Leapfrog()` | 2 | symplectic splitting ``A(h/2)\,B(h)\,A(h/2)`` |
| `Yoshida(4)` = `ForestRuth()` | 4 | triple jump on two exact sub-flows |
| `Yoshida(6)`, `Yoshida(8)` | 6, 8 | recursive triple jump ``w_1 = 1/(2-2^{1/(2m+1)}),\ w_0 = 1-2w_1`` |

`is_adaptive` is true for `DormandPrince5`, `DormandPrince853`, and
`GraggBulirschStoer`. Every other method needs `nsteps` or `h`. Splitting
methods need exact sub-flows `flowA(z, h)`, `flowB(z, h)` and are always
fixed step.

**Step control for DA states.** [`errnorm`](@ref) of a DA vector is the
maximum coefficient norm over all orders, so the step is controlled on the
whole series, not only its constant part. High orders have larger coefficients
and dominate the estimate. Use an absolute `tol` scaled to those magnitudes,
or a fixed step. Implicit midpoint stops when `errnorm` of the update is
within tolerance, including a DA state and any other series `errnorm` knows.

## Where the element step calls them

`Integrators` does not know a lattice. Callers pass a right-hand side.

- `src/maps/element_flow.jl` integrates ``dz/ds = S\nabla H``. Called with no
  step count, the default method is `DormandPrince853()`.
- `exact_track_element` (`src/maps/da/field_track.jl`) passes
  `ImplicitYoshida(4)` and a fixed `nsteps` from `ds`.
- Lie-factor tracking (`src/maps/lie/tracking.jl`) flows ``f_4`` and then
  ``f_3`` with `ImplicitMidpoint`, two steps on ``t \in [0, 1]``.
- `TaylorIntegrate` (`src/maps/lie/fringe_maps.jl`) calls `element_flow` on
  the s-dependent fringe Hamiltonian with `RK4()`.
- `Leapfrog` and `Yoshida` stay available for a user-supplied `SplitFlows`.
  They are not the production element step.

## Comparison study (Hénon–Heiles, ``E \approx 0.062``)

`benchmarks/integrator_study/integrator_study.jl` writes
`benchmarks/integrator_study/output/integrator_study.md`. It is run by hand
and is not part of `Pkg.test()`. The table below is that recorded run: error
at ``T = 10`` against a GBS(16) reference, order from ``h = 0.4 \to 0.2``,
symplecticity ``\|J^\mathsf{T} S J - S\|\) of the ``T = 10`` map, maximum
energy error over ``T = 10^4`` at ``h = 0.1``, force evaluations per step,
and run time.

| method | err h=0.1 | err h=0.05 | order (h 0.4→0.2) | ‖JᵀSJ−S‖ | max ΔE (T=1e4) | f/step | time, s |
|:--|:--|:--|:--|:--|:--|:--|:--|
| Leapfrog | 1.73e-03 | 4.33e-04 | 2.04 | 1.84e-15 | 9.19e-05 | 1.0 | 0.04 |
| Yoshida4 (Forest-Ruth) | 4.16e-05 | 2.59e-06 | 4.07 | 1.67e-15 | 7.07e-07 | 3.0 | 0.03 |
| Yoshida6 | 3.04e-07 | 4.69e-09 | 6.19 | 4.14e-15 | 1.13e-08 | 9.0 | 0.05 |
| Yoshida8 | 4.91e-09 | 1.93e-11 | 7.37 | 5.85e-15 | 1.03e-10 | 27.0 | 0.14 |
| RK4 | 4.89e-06 | 3.02e-07 | 4.03 | 4.41e-06 | 1.73e-04 | 4.0 | 0.03 |
| DP5 | 2.00e-08 | 6.29e-10 | 4.86 | 1.75e-07 | 6.73e-06 | 7.0 | 1.59 |
| DOP853 | 3.30e-14 | 1.31e-14 | 8.05 | 6.21e-14 | 1.47e-12 | 13.0 | 3.12 |
| GBS(4) ord 8 | 2.76e-13 | 1.25e-14 | 8.04 | 3.56e-13 | 1.04e-11 | 36.0 | 0.93 |

Symplectic splittings have ``J^\mathsf{T} S J - S`` at roundoff and a bounded
energy error. RK methods drift secularly. At this step DOP853 and GBS keep
that drift at ``10^{-12}``–``10^{-11}``.

## Tests

`test/test_integrators.jl`, included from `test/runtests.jl`, checks
fixed-step order on the harmonic oscillator (RK4, Dormand–Prince 5(4),
DOP853, Gragg–Bulirsch–Stoer, Yoshida 2/4/6) and one adaptive Kepler period
at eccentricity 0.6. The DA initial condition ``x' = x^2`` runs only when
`BELL_FULL_TESTS=1`.
