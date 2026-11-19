# Bell

**BELL — Beam Evolution in Accelerator Lattice**

**Author:** S. D. Kolokolchikov (Kolokolchikov Sergey Dmitrievich, `sergey.bell13@gmail.com`)

Bell is the Hamiltonian dynamics of a charged particle in an accelerator lattice. The model of an element is its s-Hamiltonian `H`. The equations of motion `dz/ds = S ∇H` and the maps of the element are derived from `H`. The variables `z = (x, p_x, y, p_y, l, δ)` are normalised to the reference particle.

The model is the Hamiltonian. The field equations of motion are an alternative formulation, kept for comparison and for tests. Computation holds no physics: differential algebra and Taylor series are algebras, and the integrators are generic. A Taylor map is symplectic only up to its truncation order. A Lie factor applied exactly stays symplectic at any truncation; the current factors carry generators `f3` and `f4`. Magnets and a thick RF cavity are integrated from `H` (`dz/ds = S ∇H`). A thin cavity stays the kick `rf_kick`.

The beam-physics manual has the same chapters in English (`docs/manual/en/`) and Russian (`docs/manual/ru/`). Figures are shared in `docs/manual/assets/`. Build one language with `julia --project=docs docs/make_pdf.jl en` or `ru`. The programmer manual stays in `docs/manual/programmer/`.

## Installation

Develop the clone in the Julia environment that should see Bell. `path` is the directory of the clone that contains `Project.toml`.

```julia
using Pkg
Pkg.develop(path = "path/to/Bell")
```

After that, `using Bell` loads this clone.

## Tests

Fast tests, from the repository root:

```text
julia --project=. -e "using Pkg; Pkg.test()"
```

Full tests are the same command with `BELL_FULL_TESTS=1`.

How to propose a change: [`CONTRIBUTING.md`](CONTRIBUTING.md). Rules for agents: [`AGENTS.md`](AGENTS.md).
