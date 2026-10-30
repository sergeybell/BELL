"""Build examples/Bell_examples.ipynb (unexecuted)."""
from pathlib import Path
import nbformat as nbf

nb = nbf.v4.new_notebook()
nb.metadata["kernelspec"] = {
    "display_name": "Julia 1.13",
    "language": "julia",
    "name": "julia-1.13",
}
nb.metadata["language_info"] = {
    "file_extension": ".jl",
    "mimetype": "application/julia",
    "name": "julia",
    "version": "1.13.2",
}

cells = []


def md(src):
    cells.append(nbf.v4.new_markdown_cell(src.strip() + "\n"))


def code(src):
    cells.append(nbf.v4.new_code_cell(src.strip() + "\n"))


md(r"""
# Bell examples (v0.10)

**Author:** Kolokolchikov Sergey Dmitrievich (`sergey.bell13@gmail.com`)

This notebook runs every numbered example of the Julia package Bell.
It is not a product release (`Project.toml` version `0.10.0`).

Open this file in Jupyter with the working directory `examples/`.
The first cell activates the examples environment (`examples/Project.toml`),
where Bell is a dev dependency on the repository root.

## What is an example, and what is not

| File | Role |
|------|------|
| `10`–`18` | Program examples. They load Bell and a structure file. |
| `common/print_util.jl` | Helper: `banner` / `section` for `10`–`18`. |
| `common/reference_values.jl` | Helper: OptiM path and COSY numbers for comparison tables. |

Each example is `run_example("NN_name/run.jl")` in a fresh module, so
constants of different examples do not clash. The printed output is the
same as `julia --project=examples examples/NN_name/run.jl`.

## Blocks of the program vs examples

| `src/` block | Example that uses it |
|--------------|----------------------|
| `coordinates.jl` | `15` (main), also `10`, `14` |
| `elements.jl` | `10`, `17` |
| `structure.jl` | `10`–`16`, `18`; FODO cell below (no OptiM) |
| `hamiltonians.jl` | `11`, `17` |
| `lie_algebra.jl` | `11`, `17` |
| `differential_algebra.jl` | `11`, `12` (check only) |
| `linear.jl` | `10` |
| `higher_order.jl` | `12` |
| `normal_form.jl` | `18` |
| `tracking.jl` | `11` |
| `exact_integrator.jl` | `12`, `14` |
| `longitudinal.jl` | `12`, `14` |
| `acceleration.jl` | `14` |
| `optimizers.jl` | `13` |
| `plotting.jl` | `10`, `12`, `13`, `14` |
| `symplectic.jl` | `10`, `11`, `15` |
| `optim_converter.jl` | `10` and every script that calls `convert_optim` |

`16` compares an external BLonD `(dt, dE)` drift with Bell tracking.
Those BLonD formulae are not a `src/` block.
""")

code(r"""
import Pkg
# Jupyter cwd is the notebook folder (examples/).
Pkg.activate(pwd())
Pkg.instantiate()
ENV["GKSwstype"] = "100"   # headless Plots
using Dates

# Each example file is a script for Main.  A fresh module avoids const
# clashes between examples but must still have include().
function run_example(path)
    m = Module()
    Core.eval(m, :(include(x) = Base.include(@__MODULE__, x)))
    Base.include(m, path)
    return nothing
end

println("cwd      = ", pwd())
println("project  = ", Base.active_project())
println("started  = ", now())
""")

md(r"""
## FODO through Bell (`lattices/fodo.bell`)

This cell loads a structure file directly:
`load_structure` → `Lattice` → one-turn `R` → closed Twiss.
It is the structure-file block without OptiM.
""")
code(r"""
using Bell, Printf
lat = load_structure(joinpath(pwd(), "lattices", "fodo.bell"))
C = circumference(lat)
po = periodic_optics(ring_matrix(lat), lat.ref, C)
@printf("  name      = %s\n", lat.name)
@printf("  elements  = %d\n", length(lat.elements))
@printf("  C         = %.6f m\n", C)
@printf("  Qx        = %.6f\n", po.x.mu / 2pi)
@printf("  Qy        = %.6f\n", po.y.mu / 2pi)
@printf("  alpha_c   = %.8f\n", po.alpha_c)
@printf("  CLOSED check uses linear.jl + structure.jl + hamiltonians.jl\n")
""")

md(r"""
## 10 — Linear optics of `magnetic_2`

OptiM `.opt` is converted to `lattices/magnetic_2.bell`. Bell loads the
structure file. Closed Twiss, dispersion, tunes, comparison with OptiM
and COSY numbers.

Blocks: converter, structure, coordinates, elements, Hamiltonians,
linear, symplectic, plotting.
""")
code('run_example(joinpath("10_linear_optics", "run.jl"))')

md(r"""
## 11 — Lie maps and tracking

Element `H`, Dragt–Finn `f3`, ring Lie product versus the independent
RK4 Taylor map, chromaticity, symplectic factor-by-factor tracking.

Blocks: Hamiltonians, Lie algebra, differential algebra (check),
tracking, symplectic.
""")
code('run_example(joinpath("11_lie_maps", "run.jl"))')

md(r"""
## 12 — Higher order (`D2`, `η`, `ξ`)

Second-order dispersion, momentum compaction, slip, chromaticity.
Cross-check with the exact-H integrator.

Blocks: higher order, differential algebra, exact integrator,
longitudinal (tune formula), plotting.
""")
code('run_example(joinpath("12_higher_order", "run.jl"))')

md(r"""
## 13 — Optional sextupole families

Places `SXD` / `SXB` / `SYB` and solves two example targets. Other
examples leave sextupoles at zero.

Block: optimizers (optional; not part of the model).
""")
code('run_example(joinpath("13_sextupole_optimization", "run.jl"))')

md(r"""
## 14 — Finite RF cavity and acceleration

Cavity checks, synchrotron tune, ramp 10 → 270 MeV, adiabatic damping.
Injection energy 10 MeV is assumed (not in the OptiM file). Takes a
few minutes.

Blocks: exact integrator, acceleration, longitudinal, plotting.
""")
code('run_example(joinpath("14_rf_acceleration", "run.jl"))')

md(r"""
## 15 — User ↔ canonical coordinates

Exact nonlinear conversion `(x, x', y, y', φ, δ) ↔ (x, px, y, py, l, δ)`.
The user Jacobian is not symplectic.

Block: coordinates.
""")
code('run_example(joinpath("15_coordinates", "run.jl"))')

md(r"""
## 16 — BLonD-like longitudinal drift (external formulae)

`(dt, dE)` drift of BLonD compared with element tracking. The formulae
live in this example, not in `src/`.
""")
code('run_example(joinpath("16_blond_longitudinal", "run.jl"))')

md(r"""
## 17 — Lie product of two factors

`concatenate(Fa, Fb)` versus the Lie series of a pair of elements.
Fixes the sign of the BCH term.

Block: Lie algebra (product of Dragt–Finn factors).
""")
code('run_example(joinpath("17_lie_product", "run.jl"))')

md(r"""
## 18 — Normal form on generators

One-turn Lie factor with thin RF, then the homological equation on
`f3`, `f4`. Torus average of laboratory `δ`.

Block: normal form.
""")
code('run_example(joinpath("18_lie_normal_form", "run.jl"))')

md(r"""
## Done

Every numbered example has been included. Plots and GIFs are under
`examples/NN_name/output/`.
""")
code(r"""
println("finished = ", now())
""")

nb.cells = cells
out = Path(__file__).with_name("Bell_examples.ipynb")
nbf.write(nb, out)
print("wrote", out)
