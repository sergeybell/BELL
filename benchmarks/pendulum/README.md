# pendulum

Step error of every ODE method in `Bell.Integrators` on a mathematical pendulum. Run `julia --project=. benchmarks/pendulum/run.jl`.

The script depends on `Bell` and `Bell.Integrators`. Nothing under `src/` calls this script.

A new check of this kind goes in this folder. A new integrator goes in `src/math/integrators/`.
