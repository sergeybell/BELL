# AGENTS.md

Rules for sessions that change the Bell repository. The human workflow is [`CONTRIBUTING.md`](CONTRIBUTING.md): one issue, one branch, one pull request.

## What Bell is

Bell (Beam Evolution in Accelerator Lattice) is a Julia package for the Hamiltonian dynamics of a charged particle in an accelerator lattice.

The physics of an element is its s-Hamiltonian `H`, not a hand-written equation of motion. The equations `dz/ds = S ∇H` and the maps of the element are derived from `H`.

The field equations of motion (`src/physics/fields.jl`) are an alternative formulation of that motion, kept for comparison. The DA path integrates those equations. The Lie path uses `H`.

The computation layer holds no physics. Differential algebra (`DiffAlgebra`) and Taylor series are algebras. The integrators in `src/math/integrators/` are generic ODE methods. The repository keeps both engines: the Lie path on TaylorSeries, and the DA path on DiffAlgebra. A Taylor map is symplectic only up to its truncation order. A Lie factor applied exactly stays symplectic at any truncation; the current factors carry generators `f3` and `f4`.

The module map and the known layering violations are in [`src/README.md`](src/README.md). Read that file before adding a `using` or an `import`.

## Code map

Dependency order used by [`src/README.md`](src/README.md):

```text
math → core → io → physics → maps → analysis → tracking → longitudinal_models → optimizers → research
```

`research` is not a directory. It is `NormalForm` in `src/maps/lie/`, included last and left out of the `Bell` facade (`using Bell.NormalForm`). The header of `src/Bell.jl` names a shorter chain and does not list `tracking` or `research`. The README is the map.

| Layer | Folder | What lives there |
| --- | --- | --- |
| math | [`src/math`](src/math/README.md) | Truncated power series (`DiffAlgebra`) and generic ODE methods ([`integrators`](src/math/integrators/README.md)). No elements and no Hamiltonian. |
| core | [`src/core`](src/core/README.md) | Coordinates, `RefParticle`, element data, sequences and `Lattice`, `symplectic_error`. |
| io | [`src/io`](src/io/README.md) | A lattice-file reader that writes a Bell structure file. |
| physics | [`src/physics`](src/physics/README.md) | The model: one s-Hamiltonian per element, the field equations kept for comparison, and the fringe model. |
| maps | [`src/maps`](src/maps/README.md) | Map construction. [`lie`](src/maps/lie/README.md) applies Taylor series to `H`. [`da`](src/maps/da/README.md) integrates the field equations on numbers or on `DiffAlgebra`. |
| analysis | [`src/analysis`](src/analysis/README.md) | Readings of a map that is already built (closed orbit, dispersion, chromaticity, momentum compaction). It does not integrate the lattice. |
| tracking | [`src/tracking`](src/tracking/README.md) | `ParticleState`, `Beam`, and trackers that call an existing engine. |
| longitudinal_models | [`src/longitudinal_models`](src/longitudinal_models/README.md) | Averaged longitudinal models kept for comparison with the `RFCavity` Hamiltonian. |
| optimizers | [`src/optimizers`](src/optimizers/README.md) | Derivative-free fit of a scalar `f(x) -> Real`. |

Element formulas are methods on the element type (`hamiltonian`, `element_rhs`, `linear_matrix`, and the map builders). A new element formula goes in `src/physics/hamiltonians.jl`. A matching comparison right-hand side goes in `src/physics/fields.jl`.

## Physics conventions

Taken from `src/core/coordinates.jl` and `src/core/elements.jl`.

The canonical state is `z = (x, p_x, y, p_y, l, δ)`, indices 1 through 6. The source names are `px`, `py`, and `delta`. `p_x = P_x / p_0` and `δ = Δp / p_0`, both divided by the reference momentum `p_0`.

`l = -v Δt`. Positive `l` is an early arrival: `canonical_to_laboratory` maps it to a phase `φ > 0`.

`energy_term(δ, ref)` is the longitudinal term of the s-Hamiltonian in the pair `(l, δ)`. It is not the kinetic energy. On the reference orbit it is zero:

```text
energy_term(δ) = (1/β0) * sqrt((1+δ)^2 + 1/(β0^2 γ0^2)) - 1/β0^2
```

and `∂(energy_term)/∂δ = v(δ) / v_0`.

`RefParticle` is the design particle. Build it with `reference_from_kinetic`. It stores the species, the kinetic energy, `γ0`, `β0`, `p0·c`, and the magnetic rigidity `Bρ`.

Units are metres, radians, MeV, and T·m. Lengths are in m, angles in rad, rest energy, kinetic energy, and `p0·c` in MeV, and `Bρ` in T·m. The normalised gradient `k1 = (∂B_y/∂x) / (Bρ)` is in 1/m². `k1 > 0` focuses in `x`: the quadrupole Hamiltonian is `energy_term(δ) - p_s + (k1/2) (x^2 - y^2)`, so `dp_x/ds = -k1 x`.

## How to check

`symplectic_error` is the symplecticity check `max |Jᵀ S J − S|` for a Jacobian, a Taylor map, a DA map, or one tracking step (`src/core/symplectic.jl`).

A result for an element agrees with `hamiltonian(e, ref, z)`. That call evaluates the same s-Hamiltonian on numbers, on `TaylorN`, and on DA. The fast suite file `test/test_physics_consistency.jl` checks that `element_field` matches `S ∇H` of this Hamiltonian.

Fast tests, from the repository root:

```text
julia --project=. -e "using Pkg; Pkg.test()"
```

`test/` passes or fails. It does not read `benchmarks/`. Optics regression, high-order DA, Aqua, and the optimizer tests run only when `BELL_FULL_TESTS=1`. Comparisons and timings live in `benchmarks/` and are run by hand. A benchmark number comes from a real run and is never fitted.

## Hard rules

- Do not force-push, and do not rewrite the history of `main`.
- Do not touch another person's uncommitted files.
- Do not commit built PDFs.
- Do not bump `version` in `Project.toml` outside a release pull request.
- Record user-visible changes in `CHANGELOG.md` under `[Unreleased]`.
- Use Conventional Commits, as specified in `CONTRIBUTING.md`.
- Write the repository in English: source, documentation, issues, pull requests, and commits.
- Do not mention other accelerator codes anywhere in the repository, issues, pull requests, or commits.
