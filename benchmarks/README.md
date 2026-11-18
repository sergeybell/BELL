# Benchmarks

`benchmarks/` holds physical comparisons and performance studies. Run the scripts by hand. They are not part of `Pkg.test()`.

Numbers in this tree come from real runs. They are never fitted.

- `twiss/` — recorded optics runs for `magnetic_2`
- `pendulum/` — integrator error against the number of steps
- `integrator_study/` — Hénon–Heiles comparison of the integrators. `integrator_study.jl` writes `output/integrator_study.md`
