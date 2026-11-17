# =============================================================================
# Bell.jl - top-level module
# BELL - Beam Evolution in Accelerator Lattice
# =============================================================================
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# Loads every block in dependency order and re-exports public names.
# A script only needs `using Bell`. NormalForm is included as research and
# is not re-exported; load it with `using Bell.NormalForm`.
#
# Layers: math → core → io → physics → maps → analysis →
# longitudinal_models → optimizers. Plotting is ext/BellPlotsExt.jl
# (weak dependency Plots), not a submodule.
# A new `using ..X` follows that direction. Include order is the existing
# load order: Hamiltonians and Symplectic stay after LieAlgebra, Fringe
# after Tracking, SextupoleFamilies after MapReadings.
# =============================================================================
"""
    Bell

BELL - Beam Evolution in Accelerator Lattice.

Two independent physics stacks share elements, coordinates and integrators:

- **Lie path** — polynomial `hamiltonian`, Lie algebra / Dragt–Finn, `Tracking`.
- **Field / DA path** — equations of motion in `Fields`, `FieldTrack` (Float64)
  and `DATracking` (any order N). Thin `RFCavity` uses `rf_kick` in both stacks.

Smooth `(φ, δ)` is a separate averaged model, not 6D RF tracking.
`NormalForm` is research-only (`using Bell.NormalForm`).
Plotting (`plot_twiss` and the other `plot_*` / `gif_*` names) loads with Plots.
"""
module Bell

# pure-mathematics DiffAlgebra / Integrators are not re-exported into Bell
# (generic names): `using Bell.DiffAlgebra`, `using Bell.Integrators`.
include("math/DiffAlgebra.jl")
include("math/integrators/Integrators.jl")

include("core/coordinates.jl")
include("core/elements.jl")
include("core/structure.jl")
include("maps/lie/lie_algebra.jl")
include("physics/hamiltonians.jl")
include("physics/fields.jl")
include("io/optim_converter.jl")
include("maps/lie/linear.jl")
include("core/symplectic.jl")
include("maps/da/field_track.jl")
include("maps/lie/tracking.jl")
include("maps/da/tracking.jl")
include("analysis/map_readings.jl")
include("longitudinal_models/smooth.jl")
include("longitudinal_models/acceleration.jl")
include("physics/fringe.jl")
include("core/sextupole_families.jl")
include("core/caches.jl")
include("core/plotting.jl")
include("optimizers/Optimizers.jl")
include("maps/lie/normal_form.jl")

using .Coordinates, .Elements, .Structure, .LieAlgebra, .Hamiltonians, .Fields,
      .OptimConverter, .Linear, .Symplectic, .FieldTrack, .Tracking, .DATracking,
      .MapReadings, .Smooth, .Acceleration, .Fringe, .SextupoleFamilies, .Optimizers

function __init__()
    # Headless GR for the plotting extension. This has to run before Plots is
    # loaded; BellPlotsExt itself starts only once Plots is already up.
    if !haskey(ENV, "GKSwstype")
        ENV["GKSwstype"] = "100"
    end
    return nothing
end

export AbstractElement, Drift, EdgeMap, EdgeModel, Enge, FringeMethod, FringeModel, FringeProfile,
       GridSearch, HardEdge, JET_ORDER, Lattice, LieFactor, LinearRamp, LongitudinalParams,
       MAP_ORDER, Magnus, Marker, NelderMead, NoEdge, OptiMFile, OptimizationResult, Optimizer,
       Particle, Quadrupole, RFCavity, RandomStart, RefParticle, SBend, SEXT_FAMILIES, SPACE,
       SYMPLECTIC_S, Sextupole, SoftEdge, Splitting, StructureFile, Tabulated, Tanh, TaylorIntegrate,
       VAR_NAMES, Z, accelerate, apply_factor, bucket_height, build_factors, c_light,
       canonical_to_laboratory, chromaticity, chromaticity_orders, chromaticity_series,
       clear_caches!, closed_orbit, coeff, compose_maps, concatenate, convert_optim, curvature,
       describe_reference, deuteron, dispersion_along, dispersion_orders, dragt_finn_factor,
       edge_map, edge_model, electron, elem_length, element_da_map, element_da_maps, element_field,
       element_lie_map, element_maps, element_rhs, energy_term, enge_default, exact_track_element,
       exact_track_turn, exponents_of_order, factor_to_map, find_optics_maxima, fixed_points,
       fringe_hamiltonian, fringe_integral, fringe_observables, fringe_sequence, gif_elements,
       gif_ramp, hamiltonian, hamiltonian_field, hamiltonian_formula, hessian, hom, identity_factor,
       index_delta, index_l, index_px, index_py, index_x, index_y, insert_cavity,
       jacobian_of_function, jacobian_of_map, kind_string, laboratory_map, laboratory_names,
       laboratory_to_canonical, lie_exp, lie_map, linear_matrix, load_structure, long_H,
       longitudinal_hamiltonian_text, magnetic_lattice, magnetic_parameters, map_linear_part,
       map_space, map_string, maxabs, measure_tune, momentum_compaction, observables,
       one_turn_map_da, optics_along, optimize, path_expansion, path_length, periodic_optics,
       place_sextupole_families, plot_acceleration, plot_chromatic, plot_dispersion2, plot_twiss,
       poisson, poly_string, product_factors, profile_function, proton, pull_generator,
       reference_from_kinetic, renormalize, response, revolution_frequency, rf_kick, rf_potential,
       rf_wavenumber, ring_lie_factor, ring_matrix, ring_metrics, ring_taylor_map, rotation_matrix_6,
       separatrix, separatrix!, set_cavity_phase, set_family_strengths, solve_A, solve_B,
       symplectic_error, synchrotron_tune, taylor_map_element, total_bend_angle, track, truncate_to,
       tunes_at_delta, turn_snapshots, unique_elements, velocity_expansion, velocity_ratio,
       with_reference, write_structure

end # module
