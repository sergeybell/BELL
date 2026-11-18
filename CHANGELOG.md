# Changelog

All notable changes to Bell are recorded in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and versions follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Detailed historical notes (v0.1–v0.17) stay in
[`docs/HISTORY.md`](docs/HISTORY.md). That file is an archive and is not copied here commit by commit.
From here on, this file is the canonical log.

A short index of the same versions remains in
[`docs/CHANGELOG.md`](docs/CHANGELOG.md).

## [Unreleased]

### Added

- A short README in each source folder under `src/` (the dependency map in `src/README.md` is unchanged), in `ext/`, and in benchmark folders that did not have one. Each note says what is there, what it depends on, what uses it, and where new code of that kind goes.
- `ParticleState` (a canonical 6-vector and a spin slot, empty by default) and `Beam` (a set of such particles). The Thomas–Bargmann–Michel–Telegdi equation is not evaluated: spin is only stored.
- `AbstractMap`, the wrappers `LinearMap` and `DAMap` (`linear_part`, call as a function, `compose` for a pair of wrappers). `LieFactor` is an `AbstractMap` as well.
- Trackers `FieldTracker` (`exact_track_turn`), `LieTracker` (`track` over the factors from `build_factors`), and `MapTracker` (a one-turn function), plus `track!` for a `Beam` and for a mutable 6-vector. The previous tracking functions stay in place.
- Tests are split into files `test/test_*.jl`: structure and coordinate round-trips, an optics regression of `magnetic_2`, and an Aqua check.
- `test/README.md` and `benchmarks/README.md` say what each directory is for. Benchmark numbers come from real runs and are never fitted.
- `test/test_physics_consistency.jl`: `element_field` matches `S ∇H` of the same Hamiltonian on DA of orders 5 and 6, and matches a numerical gradient on `Float64`.
- The top of `docs/HISTORY.md` states that it is the archive through v0.17, and that new entries go in this file.

### Changed

- Descriptions state one model: the element s-Hamiltonian `H`, with `dz/ds = S ∇H`. The field equations are an alternative formulation kept for comparison. Differential algebra and Taylor series are algebras used to compute, not a second physical model. A Taylor map is symplectic only up to its truncation order. A Lie factor applied exactly stays symplectic at any truncation; the current factors carry generators `f3` and `f4` (map order 3). The DA computation still integrates the field equations, and the Lie computation still uses `H`.
- The optics regression no longer reads `benchmarks/`. Its reference numbers are copied into `test/test_optics_regression.jl` from a real run recorded in `benchmarks/twiss/bell/results.csv`. The checks and tolerances are unchanged, and the test still runs only when `BELL_FULL_TESTS=1`.
- The Hénon–Heiles integrator study moved from `test/integrator_study.jl` and `test/output/` to `benchmarks/integrator_study/`. It is run by hand and is not part of `Pkg.test()`.
- CI is a single job: `windows-latest`, Julia 1. The OS and version matrix, coverage collection, the Codecov upload, and the documentation build are removed. The job runs a plain `Pkg.test()`.
- `Pkg.test()` runs the fast suite by default. The optics regression, high-order DA, Aqua, and the optimizer tests run only when `BELL_FULL_TESTS=1`.
- The README is in English and opens with the physics: the element s-Hamiltonian, the equations `dz/ds = S ∇H`, and the variables `z` normalised to the reference particle. Install with `Pkg.develop` on a clone; fast tests are `Pkg.test()`, and full tests set `BELL_FULL_TESTS=1`. Links to `CONTRIBUTING.md` and `AGENTS.md`.
- `LieFactor` is a subtype of `AbstractMap`. Calling `F(z)` does the same as `apply_factor(F, z)`.
- HTML documentation is built from the single tree `docs/src`. The site pages were moved out of `docs/site`. Markdown for the PDF manuals lives in `docs/manual`.
- The architecture page describes the layers `math → core → io → physics → maps/{lie,da} → analysis → tracking interfaces → effects` (there is no `effects` directory yet). Both engines remain: Lie on TaylorSeries, and DA on DiffAlgebra.
- Public names of Bell are listed by an explicit `export` instead of an `@eval export` loop. The set of names is the same; `clear_caches!()` was added.
- The DA map cache (`_DACACHE`) includes the element, `p0`, the order, the integration method itself, `tol`, and `nsteps`. `clear_caches!()` clears the Lie-map, DA-map, and fringe-field caches.
- Plots (`plot_twiss`, `plot_dispersion2`, `plot_chromatic`, `plot_acceleration`, `separatrix!`, `gif_ramp`, `gif_elements`) moved into the `BellPlotsExt` extension. Plots is a weak dependency: the Bell core does not require it, and the methods appear after `using Plots`.
- `hamiltonian(e, ref, z)` evaluates the same s-Hamiltonian on numbers, on `TaylorN`, and on DA. `hamiltonian(e, ref)` still returns the Lie-path polynomial.
- `src/` is laid out in the layers `math`, `core`, `io`, `physics`, `maps`, and `analysis`. Module names and the public API are unchanged.
- `isa` chains in `element_lie_map`, `taylor_map_element`, `element_da_map`, `exact_track_element`, `linear_matrix`, `element_rhs`, `Base.show`, and `_card_line` are replaced by methods. `magnetic_lattice` selects RF elements through `strip_rf`. The behavior is the same; `element_rhs` for `Marker` and `RFCavity` throws `ArgumentError`.
- `track(factors, z; turns)` no longer accepts the ignored arguments `ref` and `path_length`. Symplecticity of a matrix, a `TaylorN` map, a DA map, and a tracking step is one function, `symplectic_error` (the DA case calls `da_symplectic_error`).
- This changelog, the GitHub issue and pull request templates, and the remaining Russian notes in `src/` and `benchmarks/` are in English.

### Fixed

### Removed

- Unused `poisson_bracket`, `lie_series`, `lie_flow`, `vector_field_flow`, and `change_order` were removed from `DiffAlgebra`. `canonical_pairs`, `hamiltonian_vector_field`, and `da_symplectic_error` remain.

## [0.17.0] - 2026-11-02

### Added

- The s-Hamiltonian of `RFCavity`: the thin map (`L = 0`) matches `rf_kick`; the thick map uses the length `L` (`exp(-L:H:)` on the Lie path, drift/kick/drift on the field and DA paths).
- The functions `rf_potential`, `rf_wavenumber`, and `magnetic_lattice`.

### Changed

- Source layout: `src/lie`, `src/da`, `src/common`, `src/longitudinal_models`.
- Element maps (`element_lie_map`, tracking, fringe, DA) no longer contain a separate branch for RF.
