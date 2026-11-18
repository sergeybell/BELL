# benchmarks

Standalone scripts. They are not loaded with the package. `pendulum/` measures integrator step error on a mathematical pendulum. `twiss/` writes and compares optics tables; start at [twiss/README.md](twiss/README.md).

Each script depends on the packages named in its header. Nothing under `src/` imports this tree.

A new benchmark is a new subfolder with its own `run.jl` and a short README.
