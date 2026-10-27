# Programmer Manual

Implementation notes for the LieRing Julia package. Physics is in the
[Beam Physics Manual](../index.md).

## Layout

| path | role |
|---|---|
| `src/Constants.jl` | CODATA deuteron, reference particle, ``v(\delta)``, ``e(\delta)``, exact ``z\leftrightarrow u`` |
| `src/Jets.jl` | TaylorN space, `JET_ORDER=4`, `MAP_ORDER=3` |
| `src/Elements.jl` | all element types + lattice-editing helpers |
| `src/Hamiltonians.jl` | element ``H`` as polynomials |
| `src/LieOperators.jl` | Poisson bracket, ``\exp(:f:)``, Dragt–Finn |
| `src/OptiMParser.jl` | `.opt` reader (mass from particle, not file header) |
| `src/LinearOptics.jl` | Twiss, dispersion, tunes, ``\alpha_0`` |
| `src/Symplectic.jl` | ``J^\mathsf{T}SJ-S`` |
| `src/Tracking.jl` | factorised Lie tracking, thin RF kick |
| `src/TaylorMap.jl` | DA (RK4) one-turn map |
| `src/NormalForm.jl` | linear normalisation (nonlinear: planned) |
| `src/SecondOrder.jl` | ``D_2``, ``\alpha_{0,1}``, `blond_alphas`, BLonD drifts, separatrix |
| `src/ExactIntegrator.jl` | exact-``H`` integrator, finite cavity, renormalise |
| `src/Acceleration.jl` | ramp, tune fit, turn snapshots |
| `src/Plotting.jl` | figures and GIFs |
| `src/Optimizers.jl` | optional sextupole response / example solves |

Entry: `include("src/LieRing.jl"); using .LieRing`.

## Scripts

```
julia --project=. scripts/01_linear_optics.jl
julia --project=. scripts/02_lie_maps.jl
julia --project=. scripts/03_second_order.jl
julia --project=. scripts/04_optional_sextupole_optimization.jl
julia --project=. scripts/05_rf_acceleration.jl
julia --project=. scripts/06_coordinates.jl
julia --project=. scripts/07_blond_longitudinal.jl
```

## Accepted design choices

- Injection kinetic energy 10 MeV; final 270 MeV from the lattice header.
- RF: ``L=1`` m, ``V=0.2`` MV, ``h=5``, ``\phi_s=150^\circ`` (ramp) / ``180^\circ`` (flat).
- Mass = CODATA deuteron ``1875.61294257`` MeV (OptiM `Mass=2*938.272` ignored).
- Lattice sequence used once (file `Number of periods=2` is an error; bend is 360 deg).
- Closed ``D(0)=0.011`` m (OptiM initial 0 is not the closed solution).
- Test amplitudes: ``x_0=y_0=2`` mm, ``\delta_0=10^{-3}``; ensemble as in script 05.
- Optimiser targets A/B are examples only; default sextupole strengths are zero.

## Docs build

```
julia --project=docs docs/make.jl          # HTML -> docs/build/
julia docs/make_pdf.jl                     # Beam Physics PDF
julia docs/make_pdf_programmer.jl          # Programmer PDF
```