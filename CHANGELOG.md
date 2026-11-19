# Changelog

All notable changes to Bell are recorded in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Foundation examples `reference_and_coordinates`, `lattice_io`, `element_hamiltonians`, and `linear_optics_fodo`. They load `examples/lattices/fodo.bell`; the IO and optics examples also load `magnetic_2.bell`. They show the reference particle, laboratory and canonical coordinates, a structure-file round trip, element s-Hamiltonians, and FODO linear optics. Examples 10–19 are unchanged. ([#57])
- Forward-mode `Dual` numbers in `DiffAlgebra` evaluate the phase-space gradient of a scalar function at the current point, on numbers, on DA, and on `TaylorN`. `symplectic_gradient` applies a pairing matrix `S` (the default is the block pairing `[0 1; -1 0]`). On a DA coordinate jet, `map_vector_field` builds the same vector field with `deriv`, which agrees with the dual partial through order `N-1`. ([#36], [#35])
- A short README in each source folder under `src/` (the dependency map in `src/README.md` is unchanged), in `ext/`, and in benchmark folders that did not have one. Each note says what is there, what it depends on, what uses it, and where new code of that kind goes. ([#26], [#17])
- `test/README.md` and `benchmarks/README.md` say what each directory is for. Benchmark numbers come from real runs and are never fitted. ([#25], [#16])
- `src/README.md` maps the modules from the real `using` and `import` graph, states the split between the Hamiltonian model and the computation, and lists the layering violations. ([#21], [#14])

### Changed

- Examples 10–19 load `examples/lattices/magnetic_2.bell` or `fodo.bell` and the current public names (`path_length`, `laboratory_to_canonical`, `canonical_to_laboratory`, `laboratory_map`, `RFCavity`). Examples 10 and 18 load `Bell.NormalForm`. Comparison tables are no longer part of the examples; example 16 is a thin-RF longitudinal demo. ([#56])
- Root `README.md` and `AGENTS.md` state that magnets and a thick RF cavity are integrated from `H`, and that a thin cavity stays `rf_kick` (no drift/kick/drift as the production path). ([#48])
- Public docstrings in `src/` follow the settled model: magnets and a thick `RFCavity` integrate `dz/ds = S ∇H`, a thin cavity is `rf_kick`, and the field equations stay the comparison formulation. Coordinate notes state `p_x = P_x/p_0`, `δ = Δp/p_0`, and `l = -v Δt`. Version tags and retired paths are gone from those comments. ([#46], [#44])
- `Bell.jl` includes files along math → core → io → physics → maps → analysis → tracking → longitudinal models → optimizers → research, except `Hamiltonians` after `LieAlgebra` and `FringeMaps` after `MapReadings`. `clear_caches!` is a registry in `core/caches.jl`; the facade registers the map caches when it loads. `dispersion_along` carries an RF cavity as a drift by method. The `MapReadings` header states which readings build a map. ([#45], [#43])
- Fringe physics (edge types, profiles, `fringe_hamiltonian`, `FringeModel`) stays in `src/physics/fringe.jl` and does not build maps. Edge maps, the edge-map cache, and `fringe_observables` live in `src/maps/lie/fringe_maps.jl`. Family selection is by method. `fringe_hamiltonian(h0, F, Fp, ref, z)` evaluates the scalar Hamiltonian; the Lie path calls it on `Z`. `element_flow(H, z0; s0, s1)` integrates an s-dependent scalar `H(z, s)` by the dual symplectic gradient at each stage. `TaylorIntegrate` uses that flow. ([#42], [#41])
- A thick `RFCavity` is integrated from its s-Hamiltonian, `dz/ds = S ∇H`, by `element_flow`. `exact_track_element`, `element_da_map`, and `taylor_map_element` use that flow. A thin cavity stays `rf_kick`, and `element_flow` rejects it. The field equations remain for comparison and for tests. ([#40], [#39])
- `element_flow` integrates a magnet from its s-Hamiltonian, `dz/ds = S ∇H`. `exact_track_element`, `element_da_map`, and `taylor_map_element` use that flow for magnets. A marker stays the identity. An RF cavity stays a kick, or drift/kick/drift when it is thick. The field equations remain for comparison and for tests. ([#38], [#37])
- `is_adaptive` selects ODE methods that carry an error estimate. `integrate` uses that trait. Implicit midpoint stops when `errnorm` of the update is within tolerance, including DA and other series. ([#33], [#29])
- `SextupoleFamilies` lives in `src/analysis/`. `place_sextupole_families` checks that optics bounds match `lattice.elements` in length and order. ([#33], [#29])
- `SYMPLECTIC_S` is defined with the coordinates. `LieAlgebra` adds the `TaylorN` method of `symplectic_error`. `Symplectic` no longer imports `LieAlgebra`. ([#33], [#29])
- `insert_cavity` keeps the ring harmonic unless the keyword `harmonic` is set. ([#33], [#29])
- `da_var` and the DA intrinsics keep the coefficient type. A missing monomial is an `ArgumentError`. `abs` of a DA with a zero constant part and higher coefficients throws `DomainError`. ([#33], [#29])
- `canonical_pairs` is no longer part of the public DiffAlgebra API. `closed_orbit` stops at `|Δu| < 1e-12`. ([#33], [#29])
- `RandomStart` rejects `n < 1`. `GridSearch` rejects points of unequal length. ([#33], [#29])
- The energy ramp calls `after_element!`, dispatched on `RFCavity`. ([#33], [#29])
- Notes from 0.1.0 through 0.17.0 are condensed in this file. Each later change cites its pull request. ([#32], [#28])
- `CONTRIBUTING.md` and `AGENTS.md`: a user-visible pull request adds an entry under `[Unreleased]` with its pull request number. `version` in `Project.toml` is bumped, and `[Unreleased]` is renamed, only in a release pull request. ([#32], [#28])
- `AGENTS.md` describes the Hamiltonian model, the `src/` layers, the coordinate conventions, and how to check a result. `CONTRIBUTING.md` is the workflow only: one issue, a branch from `main`, fast tests, and a squash-merge. ([#27], [#19])
- The optics regression does not read `benchmarks/`. Its reference numbers are copied into `test/test_optics_regression.jl` from a real run recorded in `benchmarks/twiss/bell/results.csv`. The checks and tolerances are unchanged, and the test still runs only when `BELL_FULL_TESTS=1`. ([#25], [#16])
- The Hénon–Heiles integrator study moved from `test/integrator_study.jl` and `test/output/` to `benchmarks/integrator_study/`. It is run by hand and is not part of `Pkg.test()`. ([#25], [#16])
- Descriptions state one model: the element s-Hamiltonian `H`, with `dz/ds = S ∇H`. The field equations are an alternative formulation kept for comparison. Differential algebra and Taylor series are algebras used to compute. A Taylor map is symplectic only up to its truncation order. A Lie factor applied exactly stays symplectic at any truncation; the current factors carry generators `f3` and `f4`. The DA path integrates the field equations, and the Lie path uses `H`. ([#24], [#18])
- Repository text outside `docs/` and `examples/` is in English: this changelog, the issue and pull request templates, and the remaining notes in `src/` and `benchmarks/`. ([#23], [#15])
- `Pkg.test()` runs the fast suite. The optics regression, high-order DA, Aqua, and the optimizer tests run only when `BELL_FULL_TESTS=1`. Continuous integration is one job, Julia 1 on `windows-latest`, and runs that fast suite. ([#22], [#12])
- The README is in English and opens with the physics: the element s-Hamiltonian, the equations `dz/ds = S ∇H`, and the variables `z` normalised to the reference particle. Install with `Pkg.develop` on a clone; fast tests are `Pkg.test()`, and full tests set `BELL_FULL_TESTS=1`. ([#20], [#13])

### Fixed

- A soft edge on an `SBend` with `k1 ≠ 0` dropped the body gradient over both fringe half-widths: `Body(−W)` removes the full body Hamiltonian, and the fringe Hamiltonian had no `k1`. The local gradient is now `k1 * F(s)`, so `F = 1` matches the body and the thin-lens strength over the magnet is kept. ([#42], [#41])
- The pole-face wedge of a soft edge is applied at the effective field boundary, local `s = 0`, between the upstream and downstream halves. It had been applied at the field-free end of the fringe (`s = −W` on entrance, `s = +W` on exit). ([#42], [#41])
- Map caches for the linear matrix, the Lie map, the Taylor element map, the DA element map, and the fringe edge map were keyed by `p0` alone. Two reference particles with the same `p0` and a different mass, or the same mass and a different charge, no longer share a cached map: the key is the full `RefParticle`. ([#34], [#30])

### Removed

- `examples/Bell_examples.ipynb`, `examples/_make_notebook.py`, and `examples/common/reference_values.jl`. ([#56])
- `NormalForm.evaluate_map` had no callers. The unused `_conj` helper is gone. `delta_on_torus` remains only as an error that points at `delta_eq_phase_condition`. ([#45], [#43])
- Unused `StructureFile`. ([#33], [#29])
- `hamiltonian_vector_field` had no callers and was removed. `da_symplectic_error` remains. ([#33], [#29])
- Write-only soft-edge RK step counter `_RK_STEPS` (it was stored and never read). ([#34], [#30])
- `docs/HISTORY.md` and the version index `docs/CHANGELOG.md`. ([#32], [#28])

## [0.19.0] - 2026-11-17

### Added

- `ParticleState` (a canonical 6-vector and a spin slot, empty by default) and `Beam` (a set of such particles). The Thomas–Bargmann–Michel–Telegdi equation is not evaluated: spin is only stored. ([#7])
- `AbstractMap`, the wrappers `LinearMap` and `DAMap` (`linear_part`, call as a function, `compose` for a pair of wrappers). `LieFactor` is an `AbstractMap`. Calling `F(z)` does the same as `apply_factor(F, z)`. ([#7])
- Trackers `FieldTracker` (`exact_track_turn`), `LieTracker` (`track` over the factors from `build_factors`), and `MapTracker` (a one-turn function), plus `track!` for a `Beam` and for a mutable 6-vector. The previous tracking functions stay in place. ([#7])
- `clear_caches!()` clears the Lie-map, DA-map, and fringe-field caches. ([#6])

### Changed

- `src/` is laid out in the layers `math`, `core`, `io`, `physics`, `maps`, and `analysis`, with `longitudinal_models`, `optimizers`, and plotting beside that chain. Module names and the public API are unchanged. ([#4])
- `isa` chains in `element_lie_map`, `taylor_map_element`, `element_da_map`, `exact_track_element`, `linear_matrix`, `element_rhs`, `Base.show`, and `_card_line` are replaced by methods. `magnetic_lattice` selects RF elements through `strip_rf`. The behavior is the same; `element_rhs` for `Marker` and `RFCavity` throws `ArgumentError`. ([#5])
- `track(factors, z; turns)` no longer accepts the ignored arguments `ref` and `path_length`. Symplecticity of a matrix, a `TaylorN` map, a DA map, and a tracking step is one function, `symplectic_error` (the DA case calls `da_symplectic_error`). ([#5])
- Public names of Bell are listed by an explicit `export` instead of an `@eval export` loop. The set of names is the same. ([#6])
- The DA map cache (`_DACACHE`) includes the element, `p0`, the order, the integration method itself, `tol`, and `nsteps`. ([#6])
- Plots (`plot_twiss`, `plot_dispersion2`, `plot_chromatic`, `plot_acceleration`, `separatrix!`, `gif_ramp`, `gif_elements`) moved into the `BellPlotsExt` extension. Plots is a weak dependency: the Bell core does not require it, and the methods appear after `using Plots`. ([#6])
- HTML documentation is built from the single tree `docs/src`. The site pages were moved out of `docs/site`. Markdown for the PDF manuals lives in `docs/manual`. The architecture page describes the layers `math → core → io → physics → maps/{lie,da} → analysis → tracking interfaces → effects` (there is no `effects` directory yet). Notes through 0.17 were filed as an archive in `docs/HISTORY.md`. ([#8])

### Removed

- Unused `poisson_bracket`, `lie_series`, `lie_flow`, `vector_field_flow`, and `change_order` were removed from `DiffAlgebra`. `canonical_pairs`, `hamiltonian_vector_field`, and `da_symplectic_error` remain. ([#5])

## [0.18.0] - 2026-11-17

### Added

- `CONTRIBUTING.md`, `AGENTS.md`, and this changelog. ([#1])
- Tests are split into files `test/test_*.jl`: structure and coordinate round-trips, an optics regression of `magnetic_2` against `benchmarks/twiss/bell/results.csv`, and an Aqua check. ([#2])
- `hamiltonian(e, ref, z)` evaluates the same s-Hamiltonian on numbers, on `TaylorN`, and on DA. `hamiltonian(e, ref)` still returns the Lie-path polynomial. ([#3])
- `test/test_physics_consistency.jl`: `element_field` matches `S ∇H` of the same Hamiltonian on DA of orders 5 and 6, and matches a numerical gradient on `Float64`. ([#3])

### Changed

- Test diagnostics are printed only when `BELL_TEST_VERBOSE=1`. ([#2])

## [0.17.0] - 2026-11-02

### Added

- The s-Hamiltonian of `RFCavity`. The thin map (`L = 0`) matches `rf_kick`. The thick map uses the length `L`: `exp(-L:H:)` on the Lie path, and drift/kick/drift on the field and DA paths.
- The functions `rf_potential`, `rf_wavenumber`, and `magnetic_lattice`.

### Changed

- Source layout: `src/lie`, `src/da`, `src/common`, `src/longitudinal_models`. `Smooth` and the acceleration ramp stay as approximate models kept for comparison.
- Element maps (`element_lie_map`, tracking, fringe, and DA) use one path for every element, including `RFCavity`. A storage ring (`sin phi_s = 0`) keeps a constant `p0`.

## [0.16.0] - 2026-11-02

### Added

- Two tracking stacks share elements, coordinates, and integrators. The Lie stack builds a polynomial `hamiltonian` and Dragt–Finn factors. The field stack integrates `dz/ds` from the multipole fields: `FieldTrack` on `Float64` and `DATracking` at any order. `DiffAlgebra` stays pure mathematics.
- Derivative-free optimizers: `NelderMead` (local fit, tolerance on the simplex spread), `RandomStart` (random starts in a box, then Nelder–Mead), and `GridSearch` (a fixed list of points, optional polish). The entry point is `optimize`.

### Changed

- RF is a thin `RFCavity` plus `rf_kick`, on numbers or on DA. `insert_cavity`, `accelerate`, and `exact_track_*` use that kick. On a ramp, `renormalize` runs after each kick.
- Exact tracking lives in `FieldTrack`. `Smooth` remains a separate averaged `(phi, delta)` model.

### Removed

- `LongCavity`, `cavity_field`, `cavity_parameters`, `transit_factor`, and the module `LieFlow` (`lie_flow_map`, `hamiltonian_da`).

## [0.15.0] - 2026-11-02

### Changed

- Public names: `lattice`, `record`, `injection_reference`, `laboratory_particles`, `path_length`, `c_light`, `deuteron`, `proton`, `electron`, and the coordinate indices. Laboratory coordinates are `(x, x', y, y', phi, delta)`.
- `HigherOrder` is `MapReadings` (Newton closed orbit, integer and fractional tune, `chromaticity_series`). The reduced longitudinal model is `Smooth`. `Cavity` is `LongCavity`, and `insert_cavity` requires `in_element`.
- Integrators are one file per method, including implicit midpoint and implicit-midpoint Yoshida. The tracking step selects `method` and keeps `rf_kick`, `renormalize`, and the cavity.
- `symplectic_error` is the check of a matrix, a series map, a DA map, or one tracking step.
- The normal form is `Bell.NormalForm` and is left out of `using Bell`. Sextupole families leave `structure.jl`. Optimizers and plotting have their own folders.
- The fringe model states hard edge and soft edge. The soft edge is `Body(-W) ∘ Soft ∘ Drift(-W)`.

## [0.14.0] - 2026-10-30

### Added

- `chromaticity_series`: the tune series `Q(delta)` from a one-turn DA map (`LieFlow`), to the order of the map.
- `convert_optim(...; close_ring=true)` scales bend angles so the ring closes to `2π` when the total bend already lies within a relative `1e-6` of a full turn.

## [0.13.0] - 2026-10-30

### Added

- `Bell.LieFlow`: an element map is the flow `exp(L :H:)` summed as a Lie series to an arbitrary order. At order `N` the series uses `F = S ∇H` with `H` taken to order `N+1`. For a long or strong element (`L*rho > 3`) the map is obtained by scaling and squaring. The module is removed in 0.16.0.

## [0.12.0] - 2026-10-30

### Added

- `Bell.DiffAlgebra`: truncated power series in any number of variables and any order.
- `Bell.Integrators`: selectable ODE integrators for exact tracking and for checks.
- Neither submodule is re-exported from `Bell`.

## [0.11.0] - 2026-10-30

### Added

- A fringe-field model is chosen when the lattice is built: no edge, a hard edge (linear or nonlinear kick), or a soft edge (Enge, hyperbolic tangent, linear ramp, or a tabulated profile). The edge is integrated by splitting, Magnus, or Taylor.

### Changed

- The repository root is the package. Examples are one folder per example, with `examples/Project.toml`.
- Exported names have docstrings. The documentation site is split into a physics reference, a user's guide, a developer's guide, and an API reference.

### Fixed

- The soft dipole edge Hamiltonian includes the metric term `h*b1`. The cache key of a tabulated profile is the profile content.

## [0.10.0] - 2026-10-28

### Added

- The package is `Bell`. A lattice is a `.bell` structure file with `reference`, `parameters`, `elements`, and `sequence` (`N * cell` repeats a named sequence). `load_structure` builds the `Lattice`. Units in the file are metres, radians, and `1/m²`.
- A converter writes the `.bell` file from an external lattice description. The examples include `lattices/fodo.bell` and `lattices/magnetic_2.bell`.
- Examples `10`–`18` call Bell on a structure file. `examples/Bell_examples.ipynb` runs them.

### Changed

- High-order readings come from the Lie product (`observables`, `ring_lie_factor`). The differential-algebra map stays an independent check.
- The second-order readers are split from 6D tracking, so the reduced `(phi, delta)` model does not own the 6D map.

### Removed

- The earlier `LieRing` package tree.

## [0.9.0] - 2026-10-27

### Added

- The one-turn object is the Lie product `ring_lie_factor`. `concatenate` multiplies Dragt–Finn factors to degree 4: after transport by `R_a`, `f3 = f3a + f3b∘R_a` and `f4 = f4a + f4b∘R_a − ½[f3b∘R_a, f3a]`. Nested brackets are dropped. `one_turn_factor_rf` appends a thin RF kick as one more factor.
- `observables` reads Twiss from the linear factor and the phase slip and chromaticity from the same polynomial (`MapAnalysis`).

### Changed

- The normal form takes that one-turn factor. The homological operator is `F − F ∘ r(μ)` on homogeneous polynomials. `delta_eq` is the average of laboratory `delta` over the torus of the three actions.
- The Beam Physics manual and the programmer manual are HTML. The equilibrium chapter uses a thin RF kick.

## [0.8.0] - 2026-10-27

### Added

- Explicit 6D Dragt–Finn normal form `dragt_finn_nf6d`: fixed point, symplectic eigen-normalisation, and generators through degree 4 (`M = rot ∘ exp(:f3:) ∘ exp(:f4:)`). `delta_eq`, the actions, and the tunes are read from that normal form.

### Changed

- The equilibrium study uses a thin RF kick. The finite-cavity closed-orbit offset recorded in 0.7.0 is left out of that study.

## [0.7.0] - 2026-10-27

### Added

- A 6D Lie normal form of the one-turn map with RF (`normal_form_rf`): linear Courant–Snyder plus synchrotron Twiss, Dragt–Finn to order 3, `delta_eq` in the actions, and `tunes_nf`.
- A finite-length cavity shifts the 6D closed orbit in `delta`. A thin cavity leaves that offset at zero. A separate note records the equilibrium momentum shift and is not part of the Beam Physics manual.

## [0.6.0] - 2026-10-27

### Added

- The lattice code is an installable package, `LieRing`, developed with `Pkg.develop`.
- An equilibrium-momentum study of the magnetic lattice, and a Beam Physics chapter for that shift.

### Fixed

- A sequence whose bends already close the ring, within one degree of a full turn, is not repeated. With the exact deuteron mass the bend sum sits slightly off `2π`.

## [0.5.0] - 2026-10-27

### Changed

- The deuteron rest energy is the CODATA value.
- The model choices of this line are fixed: a chosen injection energy and RF voltage, harmonic, and phase; a longitudinal cavity field; magnets that scale with `p0`; the closed periodic dispersion. Fringe, space charge, and impedances are absent from this model.
- The nonlinear longitudinal term stays in the element Hamiltonian, including the slip coefficients `alpha0`, `alpha1`, and `alpha2`.
- The Beam Physics manual uses canonical variables only. The programmer manual lists the modules and these choices. Both exist as HTML and as PDF.

## [0.4.0] - 2026-10-27

### Added

- A Beam Physics manual (HTML and PDF): the coordinate transform, every element Hamiltonian, the RF model, acceleration, second-order optics, chromaticity, and the optional sextupoles.
- Separatrix and bucket helpers, and a check that the geometric-to-canonical transform is exact on the polynomial and on a large-amplitude round trip.

### Changed

- Element types and lattice-editing helpers live in one elements module. Optimisation is a separate module and is not used by the default lattice. Sextupole strengths default to zero.

## [0.3.0] - 2026-10-27

### Added

- These additions were first committed in the 0.4.0 snapshot.
- Map order 3, with `JET_ORDER = 4` and the Dragt–Finn generator `f4`.
- Second-order optics from the one-turn map: second-order dispersion, momentum compaction and phase slip, and chromaticity.
- A symplectic integrator of the exact Hamiltonian, including a finite cavity, and `renormalize` of the reference after the cavity.
- Optional sextupole families. The default lattice leaves them at zero.
- An acceleration ramp of the reference momentum.

## [0.2.0] - 2026-10-27

### Added

- These additions were first committed in the 0.4.0 snapshot.
- A modular symplectic ring: drift, sector bend without fringe, quadrupole, sextupole, a thin RF cavity, and a marker.
- Hamiltonians in `(x, x', y, y', phi, delta)`, expanded so the maps are order 2, tracked as Dragt–Finn factors `R ∘ exp(:f3:)`.
- An independent Taylor map of the ring, integrated for comparison.
- Linear optics of the example ring: Twiss `beta` and `alpha`, dispersion, tunes, momentum compaction, and circumference.

## [0.1.0] - 2026-10-27

### Added

- Julia examples, kept beside the lattice package: language basics, matrices, a linear FODO cell, dual numbers, Taylor maps, and Poisson brackets.

[Unreleased]: https://github.com/sergeybell/BELL/compare/b1e3056291ecf39565cac308baf85d7d1a553055...HEAD
[0.19.0]: https://github.com/sergeybell/BELL/compare/1ead7e469c024e03f509567f44ec3112ae7b9cb0...b1e3056291ecf39565cac308baf85d7d1a553055
[0.18.0]: https://github.com/sergeybell/BELL/compare/95c5f2f75aa0cc1e642590420cf9e093645296e5...1ead7e469c024e03f509567f44ec3112ae7b9cb0
[0.17.0]: https://github.com/sergeybell/BELL/compare/693edba881116b76e0fd40d52394650502d33bba...95c5f2f75aa0cc1e642590420cf9e093645296e5
[0.16.0]: https://github.com/sergeybell/BELL/compare/dc87414ad315a28f5d2bfe37e779a7b0074b584b...693edba881116b76e0fd40d52394650502d33bba
[0.15.0]: https://github.com/sergeybell/BELL/compare/3b795caba3dfdc6476171323bae745b748f8757d...dc87414ad315a28f5d2bfe37e779a7b0074b584b
[0.14.0]: https://github.com/sergeybell/BELL/compare/9068936f7763848317369d94275c2268dd4c133c...3b795caba3dfdc6476171323bae745b748f8757d
[0.13.0]: https://github.com/sergeybell/BELL/compare/79e23066d9b570d2880ae046fe57f19592c1be8c...9068936f7763848317369d94275c2268dd4c133c
[0.12.0]: https://github.com/sergeybell/BELL/compare/bd712e294501e140252176b4147cbbf6fe9f21bf...79e23066d9b570d2880ae046fe57f19592c1be8c
[0.11.0]: https://github.com/sergeybell/BELL/compare/1959aa8ca46944c734c3b830602e219e70270745...bd712e294501e140252176b4147cbbf6fe9f21bf
[0.10.0]: https://github.com/sergeybell/BELL/compare/23adcb65b955e99d9080b30b0bd3735254990d24...1959aa8ca46944c734c3b830602e219e70270745
[0.9.0]: https://github.com/sergeybell/BELL/compare/37a301d62fe0de7d15d4a6c2bad9fe5b847b9dfd...23adcb65b955e99d9080b30b0bd3735254990d24
[0.8.0]: https://github.com/sergeybell/BELL/compare/f04a411e157f839d52076d53eae0c03e9cabe76a...37a301d62fe0de7d15d4a6c2bad9fe5b847b9dfd
[0.7.0]: https://github.com/sergeybell/BELL/compare/b58606754d73ffdd4a61a59a31ff43d101857485...f04a411e157f839d52076d53eae0c03e9cabe76a
[0.6.0]: https://github.com/sergeybell/BELL/compare/3e65d1a1541e9fe776fd6c8bd6760fe0650c307f...b58606754d73ffdd4a61a59a31ff43d101857485
[0.5.0]: https://github.com/sergeybell/BELL/compare/677763ec0d68c22be8e7aa1fb5927df4a406c929...3e65d1a1541e9fe776fd6c8bd6760fe0650c307f
[0.4.0]: https://github.com/sergeybell/BELL/compare/be991e3e7507497260725351d0059e0f795b7199...677763ec0d68c22be8e7aa1fb5927df4a406c929
[0.3.0]: https://github.com/sergeybell/BELL/compare/be991e3e7507497260725351d0059e0f795b7199...677763ec0d68c22be8e7aa1fb5927df4a406c929
[0.2.0]: https://github.com/sergeybell/BELL/compare/be991e3e7507497260725351d0059e0f795b7199...677763ec0d68c22be8e7aa1fb5927df4a406c929
[0.1.0]: https://github.com/sergeybell/BELL/compare/1abfaae9c6ec2037e98306252ee83721b934be50...be991e3e7507497260725351d0059e0f795b7199
[#1]: https://github.com/sergeybell/BELL/pull/1
[#2]: https://github.com/sergeybell/BELL/pull/2
[#3]: https://github.com/sergeybell/BELL/pull/3
[#4]: https://github.com/sergeybell/BELL/pull/4
[#5]: https://github.com/sergeybell/BELL/pull/5
[#6]: https://github.com/sergeybell/BELL/pull/6
[#7]: https://github.com/sergeybell/BELL/pull/7
[#8]: https://github.com/sergeybell/BELL/pull/8
[#12]: https://github.com/sergeybell/BELL/issues/12
[#13]: https://github.com/sergeybell/BELL/issues/13
[#14]: https://github.com/sergeybell/BELL/issues/14
[#15]: https://github.com/sergeybell/BELL/issues/15
[#16]: https://github.com/sergeybell/BELL/issues/16
[#17]: https://github.com/sergeybell/BELL/issues/17
[#18]: https://github.com/sergeybell/BELL/issues/18
[#19]: https://github.com/sergeybell/BELL/issues/19
[#20]: https://github.com/sergeybell/BELL/pull/20
[#21]: https://github.com/sergeybell/BELL/pull/21
[#22]: https://github.com/sergeybell/BELL/pull/22
[#23]: https://github.com/sergeybell/BELL/pull/23
[#24]: https://github.com/sergeybell/BELL/pull/24
[#25]: https://github.com/sergeybell/BELL/pull/25
[#26]: https://github.com/sergeybell/BELL/pull/26
[#27]: https://github.com/sergeybell/BELL/pull/27
[#28]: https://github.com/sergeybell/BELL/issues/28
[#29]: https://github.com/sergeybell/BELL/issues/29
[#30]: https://github.com/sergeybell/BELL/issues/30
[#32]: https://github.com/sergeybell/BELL/pull/32
[#33]: https://github.com/sergeybell/BELL/pull/33
[#34]: https://github.com/sergeybell/BELL/pull/34
[#35]: https://github.com/sergeybell/BELL/issues/35
[#36]: https://github.com/sergeybell/BELL/pull/36
[#37]: https://github.com/sergeybell/BELL/issues/37
[#38]: https://github.com/sergeybell/BELL/pull/38
[#39]: https://github.com/sergeybell/BELL/issues/39
[#40]: https://github.com/sergeybell/BELL/pull/40
[#41]: https://github.com/sergeybell/BELL/issues/41
[#42]: https://github.com/sergeybell/BELL/pull/42
[#43]: https://github.com/sergeybell/BELL/issues/43
[#44]: https://github.com/sergeybell/BELL/issues/44
[#45]: https://github.com/sergeybell/BELL/pull/45
[#46]: https://github.com/sergeybell/BELL/pull/46
[#47]: https://github.com/sergeybell/BELL/issues/47
[#48]: https://github.com/sergeybell/BELL/pull/48
[#56]: https://github.com/sergeybell/BELL/pull/56
