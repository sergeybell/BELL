# ODE integrators (`Bell.Integrators`)

`src/integrators/Integrators.jl` is a standalone collection of integrators, generic
over `Float64` and [`DA`](@ref Bell.DiffAlgebra.DA) state vectors. The method is a parameter:

```julia
using Bell.DiffAlgebra, Bell.Integrators
r = integrate(f, z0, s0, s1; method = DormandPrince853(), tol = 1e-13)   # adaptive
r = integrate(f, z0, s0, s1; method = RK4(), nsteps = 200)              # fixed step
r = integrate(SplitFlows(drift, kick), z0, s0, s1; method = Yoshida(6), nsteps = 100)
r.z, r.nsteps, r.nrejected, r.nfev
```

| method | order | type |
|---|---|---|
| `RK4()` | 4 | fixed step |
| `DormandPrince5()` | 5(4) | embedded, adaptive |
| `DormandPrince853()` | 8(5,3) | embedded, adaptive (Hairer DOP853 error estimate) |
| `GraggBulirschStoer(k)` | ``2k`` | extrapolated midpoint, adaptive, any even order |
| `Leapfrog()` | 2 | symplectic splitting ``A(h/2)B(h)A(h/2)`` |
| `Yoshida(4)` = `ForestRuth()` | 4 | triple jump |
| `Yoshida(6)`, `Yoshida(8)` | 6, 8 | recursive triple jump ``w_1=1/(2-2^{1/(2m+1)}),\ w_0=1-2w_1`` |

Splitting methods need exact sub-flows `flowA(z, h)`, `flowB(z, h)` (e.g. drift and
kick of ``H=T(p)+V(q)``); they are always fixed step. Any `integrate` passing
`h` or `nsteps` is fixed step; RK4 is always fixed step.

**Step control for DA states.** The error norm of a DA vector is the maximum
coefficient norm over all orders (`errnorm`), so the step is controlled on the whole
Taylor map, not only its constant part. High orders have larger coefficients and
dominate the estimate; use an absolute `tol` scaled to the coefficient magnitudes,
or a fixed step.

## Bell's existing symplectic schemes
The element integrators of Bell (Yoshida splitting inside `fringe.jl`,
`exact_integrator.jl`, the Lie factorisation) work on accelerator-specific
objects and were left untouched. They are exposed only in the sense that any
exactly solvable drift/kick pair can be passed as `SplitFlows`; a direct wrapper
was not added, to avoid changing those code paths.

## Comparison study (Hénon–Heiles, ``E\approx0.062``)
`test/integrator_study.jl`: error at ``T=10`` vs. a GBS(16) reference, order from
``h=0.4\to0.2``, symplecticity ``\|J^TSJ-S\|`` of the ``T=10`` map (Jacobian from a
first-order DA run), max energy error over ``T=10^4`` at ``h=0.1``, force evaluations
per step and run time.
| method | err h=0.1 | err h=0.05 | order (h 0.4→0.2) | ‖JᵀSJ−S‖ | max ΔE (T=1e4) | f/step | time, s |
|---|---|---|---|---|---|---|---|
| Leapfrog | 1.73e-03 | 4.33e-04 | 2.04 | 1.84e-15 | 9.19e-05 | 1.0 | 0.04 |
| Yoshida4 (Forest-Ruth) | 4.16e-05 | 2.59e-06 | 4.07 | 1.67e-15 | 7.07e-07 | 3.0 | 0.03 |
| Yoshida6 | 3.04e-07 | 4.69e-09 | 6.19 | 4.14e-15 | 1.13e-08 | 9.0 | 0.05 |
| Yoshida8 | 4.91e-09 | 1.93e-11 | 7.37 | 5.85e-15 | 1.03e-10 | 27.0 | 0.14 |
| RK4 | 4.89e-06 | 3.02e-07 | 4.03 | 4.41e-06 | 1.73e-04 | 4.0 | 0.03 |
| DP5 | 2.00e-08 | 6.29e-10 | 4.86 | 1.75e-07 | 6.73e-06 | 7.0 | 1.59 |
| DOP853 | 3.30e-14 | 1.31e-14 | 8.05 | 6.21e-14 | 1.47e-12 | 13.0 | 3.12 |
| GBS(4) ord 8 | 2.76e-13 | 1.25e-14 | 8.04 | 3.56e-13 | 1.04e-11 | 36.0 | 0.93 |

Symplectic splittings have ``J^TSJ-S`` at roundoff and bounded energy error;
RK methods drift secularly but DOP853/GBS are so accurate that the drift is at
1e-12–1e-11 for this step. DP5/DOP853 timings include allocation overhead of the
generic implementation (not yet optimised).

## Tests
`test/runtests.jl`: orders measured on the harmonic oscillator (RK4 4.00, DP5 5.04,
DOP853 8.05, GBS(4) 8.01, GBS(5) 10.01, Yoshida 2/4/6 → 2.01/4.00/6.02); adaptive Kepler
``e=0.6``, one period (DP5 1.4e-10, DOP853 1.7e-11, GBS(8) 4.5e-12); DA initial condition
for ``x'=x^2``: Taylor coefficients to relative 4e-15 (DOP853) and 6e-14 (GBS).