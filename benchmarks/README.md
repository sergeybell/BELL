# Benchmarks

`benchmarks/` holds physical comparisons and performance studies. Run the scripts by hand. They are not part of `Pkg.test()`. Numbers in this tree come from real runs and are never fitted.

- `twiss/` — recorded optics runs for `magnetic_2`; start at [twiss/README.md](twiss/README.md)
- `pendulum/` — integrator step error on a mathematical pendulum
- `integrator_study/` — Hénon–Heiles comparison of the integrators; `integrator_study.jl` writes `output/integrator_study.md`

Nothing under `src/` imports this tree. A new benchmark is a new subfolder with its own `run.jl` (or study script) and a short README.
