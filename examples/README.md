# Examples

Teaching order for Bell. Each folder has one `run.jl`. From the repository root:

```text
julia --project=examples examples/01_reference_and_coordinates/run.jl
```

Replace the folder name to run any other example. `[core]` marks the examples that carry the model: the s-Hamiltonian, the linear map, the element maps, the Lie product, the readings of a map that is already built, and the comparison of those constructions.

Figures are written under `<folder>/output/` and are not part of the repository.

`examples/_make_notebook.py` writes `examples/Bell_examples.ipynb` with the captured output of each script. It runs `01`–`12` and `14`–`15` unchanged. For `13_acceleration` the notebook runs `examples/_notebook_13_smoke.jl` (a cavity check, a short fixed-energy tune, and `accelerate` over a few turns). `13_acceleration/run.jl` is the full ramp from 10 MeV to 270 MeV.

| folder | shows | lattices | main API |
| --- | --- | --- | --- |
| `01_reference_and_coordinates` | Reference particle, `energy_term`, laboratory and canonical coordinates, and the polynomial conversion | `fodo.bell` | `reference_from_kinetic`, `describe_reference`, `deuteron`, `proton`, `electron`, `energy_term`, `velocity_ratio`, `laboratory_to_canonical`, `canonical_to_laboratory`, `laboratory_map`, `ring_taylor_map`, `symplectic_error`, `jacobian_of_map` |
| `02_lattice_io` | Load and write a structure file, insert a cavity, drop it, replace the reference | `fodo.bell`, `magnetic_2.bell` | `load_structure`, `write_structure`, `unique_elements`, `kind_string`, `path_length`, `total_bend_angle`, `insert_cavity`, `RFCavity`, `magnetic_lattice`, `with_reference`, `linear_matrix` |
| `03_element_hamiltonians` [core] | s-Hamiltonian of a drift, a quadrupole, a sector bend, a sextupole, a thin cavity, and a thick cavity | `fodo.bell` | `hamiltonian`, `hamiltonian_formula`, `poly_string`, `energy_term`, `curvature`, `Sextupole`, `rf_wavenumber`, `rf_potential` |
| `04_linear_optics` [core] | Element matrices, the one-turn matrix, periodic Twiss, and dispersion. FODO first, then `magnetic_2` | `fodo.bell`, `magnetic_2.bell` | `linear_matrix`, `ring_matrix`, `periodic_optics`, `optics_along`, `Lattice`, `plot_twiss`, `Bell.NormalForm.rotation_residual` |
| `05_element_maps` [core] | Lie series, Taylor map, DA map, and one numeric step of a single element | `fodo.bell` | `lie_map`, `element_lie_map`, `dragt_finn_factor`, `factor_to_map`, `taylor_map_element`, `element_da_map`, `element_flow`, `exact_track_element`, `rf_kick` |
| `06_lie_maps` [core] | Dragt–Finn product on FODO, the Poisson-bracket sign of `concatenate`, then the same product on `magnetic_2` against the Taylor map | `fodo.bell`, `magnetic_2.bell` | `build_factors`, `concatenate`, `product_factors`, `ring_lie_factor`, `apply_factor`, `track`, `pull_generator`, `poisson`, `hom`, `compose_maps`, `ring_taylor_map`, `closed_orbit` |
| `07_map_readings` [core] | Second-order dispersion, momentum compaction, slip factor, and chromaticity | `magnetic_2.bell` | `dispersion_orders`, `dispersion_along`, `path_expansion`, `velocity_expansion`, `observables`, `chromaticity_orders`, `tunes_at_delta`, `measure_tune`, `exact_track_turn`, `plot_chromatic` |
| `08_da_basics` | Differential-algebra arithmetic, then a quadrupole map and a FODO one-turn map | `fodo.bell` | `Bell.DiffAlgebra` (`DASpace`, `da_vars`, `deriv`, `compose`), `map_space`, `element_da_map`, `one_turn_map_da`, `chromaticity_series`, `chromaticity` |
| `09_method_comparison` [core] | Lie path, DA map of `H`, Taylor check, and numeric tracking, on elements, one FODO cell, and the FODO ring | `fodo.bell` | `element_lie_map`, `element_da_map`, `taylor_map_element`, `ring_taylor_map`, `one_turn_map_da`, `ring_lie_factor`, `element_rhs`, `exact_track_element`, `exact_track_turn`, `Bell.Integrators` |
| `10_sextupole_families` | Three sextupole families, the response solve, and a derivative-free fit of the same residual | `magnetic_2.bell` | `place_sextupole_families`, `SEXT_FAMILIES`, `response`, `solve_A`, `solve_B`, `set_family_strengths`, `ring_metrics`, `optimize`, `NelderMead`, `RandomStart`, `GridSearch`, `plot_chromatic`, `plot_dispersion2` |
| `11_trackers_and_beam` | `ParticleState`, `Beam`, and the three trackers on the FODO ring | `fodo.bell` | `ParticleState`, `Beam`, `FieldTracker`, `LieTracker`, `MapTracker`, `LinearMap`, `track!`, `exact_track_turn`, `symplectic_error` |
| `12_rf_longitudinal` | Thin cavity (`rf_kick`), thick cavity integrated from `H`, and the smooth longitudinal model, at fixed energy | `magnetic_2.bell` | `insert_cavity`, `rf_kick`, `magnetic_lattice`, `path_expansion`, `synchrotron_tune`, `LongitudinalParams`, `long_H`, `fixed_points`, `separatrix`, `bucket_height`, `hamiltonian`, `measure_tune` |
| `13_acceleration` | Thin-cavity ramp from 10 MeV to the structure-file energy, with an ensemble | `magnetic_2.bell` | `accelerate`, `renormalize`, `set_cavity_phase`, `revolution_frequency`, `turn_snapshots`, `path_expansion`, `plot_acceleration`, `gif_ramp`, `gif_elements` |
| `14_normal_form` | One-turn Lie factor with a thin cavity, then the homological normal form | `magnetic_2.bell` | `Bell.NormalForm` (`one_turn_factor_rf`, `dragt_finn_nf6d`), `pull_generator`, `rotation_matrix_6`, `ring_taylor_map`, `chromaticity_orders` |
| `15_fringe_fields` | Hard and soft edges, fringe profiles, and the edge integrators on `magnetic_2` | `magnetic_2.bell` | `FringeModel`, `NoEdge`, `HardEdge`, `SoftEdge`, `Enge`, `Tanh`, `LinearRamp`, `Tabulated`, `Splitting`, `Magnus`, `TaylorIntegrate`, `fringe_integral`, `fringe_observables`, `profile_function` |

`lattices/fodo.bell` and `lattices/magnetic_2.bell` live in `examples/lattices/`. `examples/common/print_util.jl` is the shared console layout.
