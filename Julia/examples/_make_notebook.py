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

Open from the `Julia` folder, or open this file in Jupyter with the
working directory `Julia/examples`. The first cell activates the Bell
project one directory up.

## What is an example, and what is not

| File | Role |
|------|------|
| `01`–`06` | Julia language tutorials. They do not `using Bell`. |
| `10`–`18` | Program examples. They load Bell and a structure file. |
| `07`–`09` | Do not exist. Old `Julia/ring` scripts `07`–`09` became `16`–`18`. |
| `print_util.jl` | Helper: `banner` / `section` for `10`–`18`. |
| `reference_values.jl` | Helper: OptiM path and COSY numbers for comparison tables. |
| `linear_optics.jl` | Leftover 2×2 transfer-matrix helpers. `03` is self-contained and does not include it. |

Each numbered example is `run_example("NN_....jl")` in a fresh module, so
constants such as `N` in `03` and `05` do not clash. The printed output
is the same as `julia --project=. examples/NN_....jl`.

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
Pkg.activate(joinpath(pwd(), ".."))
Pkg.instantiate()
ENV["GKSwstype"] = "100"   # headless Plots
using Dates

# Each example file is a script for Main.  A fresh module avoids const
# clashes (N in 03 vs 05) but must still have include().
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
## 01 — Julia basics

Printing, types, arithmetic (`/` vs `div`), functions, `if`, `for`, arrays.
No physics. No Bell.
""")
code('run_example("01_hello.jl")')

md(r"""
## 02 — Vectors and matrices

`LinearAlgebra`: product vs element-wise `.*`, `det`, `inv`, `A \\ b`,
eigenvalues, a 2-D rotation. Needed later for transfer matrices.
""")
code('run_example("02_matrices.jl")')

md(r"""
## 03 — FODO cell by hand (not Bell)

Standard linear optics with 2×2 / 3×3 transfer matrices for one FODO
cell: periodic Twiss, dispersion, tune, chromaticity, plots. This is
the language tutorial. The same lattice as a Bell structure file is
`lattices/fodo.bell` (cell after example 06).
""")
code('run_example("03_fodo_linear.jl")')

md(r"""
## 04 — Dual numbers

First-order automatic differentiation: feed `x0 + ε` into an unchanged
program, read `f'(x0)` from the `ε` coefficient. Prototype of DA.
""")
code('run_example("04_dual_numbers.jl")')

md(r"""
## 05 — Taylor maps (`TaylorSeries`)

A map as a vector of truncated multivariate polynomials: coefficients
`(final | monomial)`, composition, inverse, Jacobian. This is the DA
object. In Bell it is the *check*, not the ring object.
""")
code('run_example("05_taylor_map.jl")')

md(r"""
## 06 — Poisson bracket and Lie operators

`[f,g]`, Lie operator `:f:`, flow `exp(:-t H:)`. Language of
`src/lie_algebra.jl`.
""")
code('run_example("06_poisson_bracket.jl")')

md(r"""
## FODO through Bell (`lattices/fodo.bell`)

Examples `01`–`06` do not load the package. This cell does:
`load_structure` → `Lattice` → one-turn `R` → closed Twiss.
It is the structure-file block without OptiM.
""")
code(r"""
using Bell, Printf
lat = load_structure(joinpath(pwd(), "..", "lattices", "fodo.bell"))
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
code('run_example("10_linear_optics.jl")')

md(r"""
## 11 — Lie maps and tracking

Element `H`, Dragt–Finn `f3`, ring Lie product versus the independent
RK4 Taylor map, chromaticity, symplectic factor-by-factor tracking.

Blocks: Hamiltonians, Lie algebra, differential algebra (check),
tracking, symplectic.
""")
code('run_example("11_lie_maps.jl")')

md(r"""
## 12 — Higher order (`D2`, `η`, `ξ`)

Second-order dispersion, momentum compaction, slip, chromaticity.
Cross-check with the exact-H integrator.

Blocks: higher order, differential algebra, exact integrator,
longitudinal (tune formula), plotting.
""")
code('run_example("12_higher_order.jl")')

md(r"""
## 13 — Optional sextupole families

Places `SXD` / `SXB` / `SYB` and solves two example targets. Other
examples leave sextupoles at zero.

Block: optimizers (optional; not part of the model).
""")
code('run_example("13_optional_sextupole_optimization.jl")')

md(r"""
## 14 — Finite RF cavity and acceleration

Cavity checks, synchrotron tune, ramp 10 → 270 MeV, adiabatic damping.
Injection energy 10 MeV is assumed (not in the OptiM file). Takes a
few minutes.

Blocks: exact integrator, acceleration, longitudinal, plotting.
""")
code('run_example("14_rf_acceleration.jl")')

md(r"""
## 15 — User ↔ canonical coordinates

Exact nonlinear conversion `(x, x', y, y', φ, δ) ↔ (x, px, y, py, l, δ)`.
The user Jacobian is not symplectic.

Block: coordinates.
""")
code('run_example("15_coordinates.jl")')

md(r"""
## 16 — BLonD-like longitudinal drift (external formulae)

`(dt, dE)` drift of BLonD compared with element tracking. The formulae
live in this example, not in `src/`.
""")
code('run_example("16_blond_longitudinal.jl")')

md(r"""
## 17 — Lie product of two factors

`concatenate(Fa, Fb)` versus the Lie series of a pair of elements.
Fixes the sign of the BCH term.

Block: Lie algebra (product of Dragt–Finn factors).
""")
code('run_example("17_lie_product.jl")')

md(r"""
## 18 — Normal form on generators

One-turn Lie factor with thin RF, then the homological equation on
`f3`, `f4`. Torus average of laboratory `δ`.

Block: normal form.
""")
code('run_example("18_lie_normal_form.jl")')

md(r"""
## Done

Every numbered example has been included. Plots and GIFs are under
`Julia/output/`.
""")
code(r"""
println("finished = ", now())
""")

nb.cells = cells
out = Path(__file__).with_name("Bell_examples.ipynb")
nbf.write(nb, out)
print("wrote", out)
