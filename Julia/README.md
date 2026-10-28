# Bell

**Author:** Kolokolchikov Sergey Dmitrievich (`sergey.bell13@gmail.com`)

Julia package for symplectic modelling of a magnetic ring. The object of
the ring is a product of Lie transformations. Linear optics and higher-order
numbers are readings of that product. A truncated Taylor map is an independent
check. Architecture: [`docs/architecture.md`](docs/architecture.md). Physics
manual: [`docs/src/index.md`](docs/src/index.md).

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
using Bell
lat = load_structure("lattices/fodo.bell")
```

A concrete ring is a structure file under `lattices/`. Convert an OptiM
`.opt` with `convert_optim(opt_path, structure_path)`, then
`load_structure`. The equilibrium-momentum study lives in
`C:\GIT_REPOS\SCT_2\Julia\equilibrium_energy\` and should `Pkg.develop`
this folder.

Package version `0.10.0` (not a product release). Tutorial examples `01`–`06`
teach the language and do not load the package. Program examples start at
`10_linear_optics.jl`. There are no `07`–`09`. Notebook:
`examples/Bell_examples.ipynb`.

## Install

Windows (recommended — installs `juliaup`, the Julia version manager):

```
winget install --id 9NJNWW8PVKMN -s msstore     # "Julia" from the Microsoft Store
# alternative: winget install Julialang.Juliaup
juliaup status              # installed versions
juliaup update              # update to the latest release
```

Linux/macOS: `curl -fsSL https://install.julialang.org | sh`

Open a **new** terminal afterwards so that `julia` is on the PATH, then check `julia --version`.

## Run

Always run from this folder (`C:\GIT_REPOS\BELL\Julia`) with `--project=.`, so that the
packages listed in `Project.toml` are found:

```
julia --project=. -e "using Pkg; Pkg.instantiate()"   # once: installs packages from Project.toml
julia --project=. examples/01_hello.jl                 # run one example
julia --project=.                                      # interactive REPL in this project
```

**Windows console tip:** run `chcp 65001` once in a `cmd` window (switches the console to
UTF-8) if you see garbled characters, e.g. in Julia's own messages. The examples themselves
print only ASCII text, so their output is readable in any console.

In the REPL: `include("examples/03_fodo_linear.jl")` runs a file, `?det` shows help,
`]` enters package mode, `;` shell mode, Backspace returns. Unicode names are typed
LaTeX-style: `\mu<TAB>` → `μ`, `\pi<TAB>` → `π`.

**VS Code**: install the *Julia* extension (julialang.language-julia). Open this folder,
then `Shift+Enter` runs the current line/selection in the Julia REPL, `Ctrl+F5` runs the file.
The extension also gives a plot pane, variable workspace and debugger.

## Packages

```
julia --project=.
julia> ]                      # package mode
(Julia) pkg> add TaylorSeries
(Julia) pkg> status
```

or in code: `using Pkg; Pkg.activate("."); Pkg.add("TaylorSeries")`.
`Project.toml` lists the dependencies, `Manifest.toml` pins exact versions.
`LinearAlgebra` and `Printf` are standard libraries (always available).

## Examples

Every example is self-documenting:

* **inside the file** a header block explains what the example teaches, the math/physics
  background with the key formulas, the sections of the file, how to run it and how to read
  the output; each code block has a short comment saying what it does and why;
* **when run** it first prints a banner with *What this example does*, then for every section
  a heading with 1-3 lines explaining what is computed and what to look for, then the results,
  and at the end a *Summary / what you learned*.

Program examples `10`–`18` take `banner` / `section` from `print_util.jl` and
the OptiM/COSY numbers from `reference_values.jl`. Those two files are not
examples. `linear_optics.jl` is leftover 2×2-matrix helpers; example `03` does
not include it. Printed text is plain ASCII.

| File | Content |
|------|---------|
| `examples/01_hello.jl` | Beginner basics: Hello World, variables and types, arithmetic, integer vs. float division, powers, `sqrt`, `sin`/`cos`, `pi`, functions, `if`, `for` loops (sum of squares, factorial), arrays, string interpolation |
| `examples/02_matrices.jl` | Vectors and matrices: creation, indexing, matrix vs. element-wise product (`*` vs `.*`), transpose, `det`, `tr`, `inv`, identity `I`, solving `A x = b` with `\`, eigenvalues, 2D rotation matrix |
| `examples/03_fodo_linear.jl` | Standard linear-optics procedure for one FODO cell (thick quads + sector dipoles): build the cell, periodic Twiss and dispersion, propagation, closure check, tunes / beta / D / chromaticity, plots (saved to `output/`) |
| `examples/04_dual_numbers.jl` | Hand-written dual numbers (first-order automatic differentiation = 1st-order DA); transfer matrix of a beamline by tracking `x0 + eps` |
| `examples/05_taylor_map.jl` | Taylor maps with `TaylorSeries` (`JetSpace`, `TaylorN`): truncated polynomial algebra, expansions, coefficient table (final\|monomial), composition, inverse, Jacobian |
| `examples/06_poisson_bracket.jl` | Poisson brackets and Lie operators: bracket rules, conserved quantities, Lie transformation `exp(:-tH:)` as the time-t flow, BCH, symplecticity |
| `examples/linear_optics.jl` | Not an example: leftover 2×2 transfer-matrix helpers (`03` is self-contained) |
| `examples/print_util.jl` | Not an example: `banner` / `section` used by `10`–`18` |
| `examples/reference_values.jl` | Not an example: OptiM path and COSY numbers used by `10`–`16`, `18` |
| `examples/Bell_examples.ipynb` | Notebook: descriptions + `include` of every numbered example |
| `examples/10_linear_optics.jl` | Load the Nuclotron-based structure, closed Twiss, dispersion, tunes |
| `examples/11_lie_maps.jl` | Lie factors, tracking, Taylor-map check |
| `examples/12_higher_order.jl` | D2, slip, chromaticity from the one-turn factor |
| `examples/13_optional_sextupole_optimization.jl` | Optional sextupole families |
| `examples/14_rf_acceleration.jl` | Finite cavity and ramp |
| `examples/15_coordinates.jl` | Canonical ↔ laboratory transform |
| `examples/16_blond_longitudinal.jl` | External (dt, dE) drift formulae (not the ring product) |
| `examples/17_lie_product.jl` | BCH product of two factors |
| `examples/18_lie_normal_form.jl` | One-turn factor and torus average |

Tutorials 01, 02 and 04 use only the standard library; 03 needs `Plots`, 05-06 need `TaylorSeries`.
Program examples 10–18 need `using Bell`.
All of them: `julia --project=. examples/XX_name.jl`.

### 01 - Julia basics

```
julia --project=. examples/01_hello.jl
```

The first steps in the language: `println`, variables and their types (`typeof`), arithmetic
(`/` vs `div` vs `%`), math functions and `pi`, defining functions, `if`/`else`, `for` loops,
arrays (indices start at 1, `push!`) and string interpolation (`"$x"`). Every printed line has
the form *expression = value*, often with a check (e.g. sum of squares vs. `n(n+1)(2n+1)/6`).

### 02 - Vectors and matrices

```
julia --project=. examples/02_matrices.jl
```

Linear algebra with the standard library `LinearAlgebra`: building vectors and matrices,
indexing, `A * B` vs `B * A` vs `A .* B`, transpose, `det`, `tr`, `inv`, the identity `I`,
solving `A x = b` with `A \ b`, eigenvalues/eigenvectors and a 2D rotation matrix
(det = 1, orthogonal, eigenvalues `exp(+-i theta)`). These are the tools used for transfer
matrices in example 03.

### 03 - Linear optics of a FODO cell

```
julia --project=. examples/03_fodo_linear.jl
```

The standard procedure of linear beam optics for one FODO cell
(`QF/2 - O - B - O - QD - O - B - O - QF/2`, 16 cells form a ring):

1. **build the cell** from a list of thick elements and multiply their matrices
   (3x3 horizontal with dispersion, 2x2 vertical);
2. **find the periodic solution**: stability `|tr M| < 2`, periodic Twiss parameters and
   dispersion read off the cell matrix;
3. **propagate** beta, alpha, D and the phase advance along s in fine slices
   (`B1 = M B0 M^T`);
4. **check closure**: values at the cell end equal the start values (`max |end-start| ~ 2e-13`);
5. **main parameters**: tunes (`Qx = 2.55418, Qy = 2.16017`), beta and D extrema,
   natural chromaticity;
6. **plots** of beta_x, beta_y, D_x and the cell layout, saved to `output/fodo_twiss.png`
   (one cell) and `output/fodo_twiss_3cells.png` (three cells, showing the periodicity).

*Twiss parameters*: in a periodic lattice a particle performs betatron oscillations
`x(s) = sqrt(2J beta(s)) cos(mu(s) + phi0)`. `beta(s)` is the envelope function (beam size
`sigma = sqrt(emittance * beta)`), `alpha = -beta'/2`, `gamma = (1 + alpha^2)/beta`, and
`mu(s) = int ds/beta` is the phase advance; the tune is the phase advance per turn / 2pi.

`Plots` (GR backend) is already listed in `Project.toml`; if it is missing, install it once with
`julia --project=.` and then `] add Plots` in the REPL. The **first run takes about a minute or more**
because Plots is precompiled; later runs take a few seconds.

### 04 - Dual numbers (automatic differentiation)

```
julia --project=. examples/04_dual_numbers.jl
```

Problem: a function is given by a program, and we want its value and its *exact* derivative
without deriving formulas. Solution: dual numbers `a + b eps` with `eps^2 = 0`; feeding
`x0 + eps` into the unchanged program returns `f(x0) + f'(x0) eps`, so the derivative is the eps
coefficient (e.g. `x^2` at 3 gives `9 + 6eps`, `x^3 + 2x` gives `33 + 29eps`). Then the same
trick is applied to a beamline (drift - thin lens - drift): pushing `x0 + eps` (or `x0' + eps`)
through the element functions gives the transfer matrix, which is compared with the matrix
product. This is first-order differential algebra (DA).

### 05 - Taylor maps

```
julia --project=. examples/05_taylor_map.jl
```

What a Taylor map is - a vector of multivariate polynomials truncated at order N, stored as a
coefficient table (final\|monomial) - and how to work with it in `TaylorSeries` (`JetSpace`,
`TaylorN`): polynomial arithmetic and truncation, `getcoeff`, `evaluate`, Taylor expansions of
`exp`, `sin`, `1/(1+x)`, `sqrt(1+x)`, `differentiate`/`integrate`, a nonlinear 2D map,
composition (not commutative), inverse map (`inverse_map` and fixed-point iteration), Jacobian
and its determinant.

### 06 - Poisson brackets and Lie operators

```
julia --project=. examples/06_poisson_bracket.jl
```

Poisson brackets on `TaylorN` polynomials: bracket properties (antisymmetry, linearity,
Leibniz, Jacobi), conserved quantities, the Lie transformation `exp(:-tH:)` as the time-t flow
(free particle, harmonic oscillator vs. cos/sin, anharmonic oscillator vs. a hand-written RK4),
non-commutativity and the BCH term, symplecticity `[Q,P] = 1`.
