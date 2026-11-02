# Bell

**BELL - Beam Evolution in Accelerator Lattice**

**Author:** S. D. Kolokolchikov (Kolokolchikov Sergey Dmitrievich, `sergey.bell13@gmail.com`)

Bell is a Julia code for symplectic modelling of beam and particle evolution
in an accelerator lattice. Element Hamiltonians are turned into Lie
transformations; their product is the one-turn object from which optics and
higher-order quantities are read, with exact integration used as a check.

## Install and use

```
julia --project=. -e "using Pkg; Pkg.instantiate()"
julia --project=examples -e "using Pkg; Pkg.instantiate()"
julia --project=examples examples/10_linear_optics/run.jl
```

```julia
using Bell
lat = load_structure("examples/lattices/fodo.bell")
```

## Documentation

* Physics manual: [`docs/src/index.md`](docs/src/index.md)
* Docs site: [`docs/site/index.md`](docs/site/index.md); build HTML with `julia --project=docs docs/make.jl`
* Version history: [`docs/CHANGELOG.md`](docs/CHANGELOG.md), details in [`docs/HISTORY.md`](docs/HISTORY.md)
