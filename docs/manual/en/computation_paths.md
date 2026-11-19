# Computation paths

This chapter names every route that exists now and the file that owns it. It does not change those routes. The target in [Lie factors and differential algebra](lie_and_da.md) is one engine on ``H``; the list below is what the source tree does today.

The layer order is `math`, `core`, `io`, `physics`, `maps`, `analysis`, `tracking`, `longitudinal_models`, `optimizers`, then the normal form as a research module. `math` has no elements and no Hamiltonian. `physics` holds ``H``, the comparison field, and the fringe scalar. `maps` builds maps. `analysis` reads a map that already exists.

## Linear

`src/maps/lie/linear.jl` builds the matrix ``R=\exp(L\,S\,\mathrm{Hess}\,H)`` from the quadratic part of the same ``H``. That matrix is ``\exp(:f_2:)`` in closed form. One-turn optics multiply these matrices, or read ``R`` from the Lie factor. The reading is [Linear optics](linear.md).

## Lie product, Dragt–Finn, BCH

`src/physics/hamiltonians.jl` evaluates ``H`` on the Taylor variables and keeps a polynomial of degree 4 with the constant removed. `src/maps/lie/lie_algebra.jl` is the Poisson bracket, the Lie series, and the product of two Dragt–Finn factors, including the degree-4 Baker–Campbell–Hausdorff term. `src/maps/lie/tracking.jl` builds `element_lie_map` and `ring_lie_factor`, and expands a factor to a polynomial with `factor_to_map`. The factors carry ``f_3`` and ``f_4``. A thin cavity enters as the kick, expanded to the same degree. The physics of the product is [Lie product of the ring](lie.md).

## Flow of ``H``

`src/maps/element_flow.jl` integrates ``\mathrm{d}z/\mathrm{d}s=S\nabla H`` for a magnet, a thick cavity, or an explicit scalar ``H(z,s)``.

* On numbers and on Taylor series, each stage evaluates the symplectic gradient with a dual number (`src/math/dual.jl`).
* On a differential-algebra jet, the vector field is built once by differentiating the polynomial (`deriv` on the identity) and is then composed with the current state. That derivative agrees with the dual partial through order ``N-1``.

A marker, and a magnet of length zero, are the identity. A thin cavity is rejected here. A soft fringe calls the explicit ``H(z,s)`` form, because the fringe scalar depends on ``s``.

## Tracking a point

`src/maps/da/field_track.jl` tracks one particle. The name of the file is historical. A magnet and a thick cavity go through `element_flow`. A thin cavity is `rf_kick`: the exact flow of ``G`` from [Element Hamiltonians](hamiltonians.md), with no series cut. A drift is one step. `exact_track_turn` walks the lattice. `renormalize` installs a new reference after a thin kick on a ramp ([Acceleration](acceleration.md)).

## Differential-algebra maps

`src/maps/da/tracking.jl` integrates the same flow on the identity jet of order ``N`` (`element_da_map`, `one_turn_map_da`). A thin cavity is `rf_kick` on the differential-algebra variables. A thick cavity needs the path length, because the wavenumber does. The cache key is the element, the whole reference particle, the order, the method, the tolerance, and the step count.

## Taylor maps, as a check

`taylor_map_element` and `ring_taylor_map`, in `src/maps/lie/tracking.jl`, integrate the flow in the Taylor algebra and truncate. The resulting polynomial is symplectic only up to the kept degree. It is the check drawn beside the Lie factor, not the object the optics read.

## Field formulas, as a check

`src/physics/fields.jl` is `element_field` and `element_rhs`: the hard-edge ``\mathrm{d}z/\mathrm{d}s`` written from ``h``, ``k_1`` and ``k_2``. The fast tests compare it with ``S\nabla H``. Production tracking does not call it. A cavity has no field-formula method; the cavity is ``H`` or the kick.

## Fringe

`src/physics/fringe.jl` is the scalar: edge type, profile, `fringe_hamiltonian`, and the fringe model. It does not build a map. `src/maps/lie/fringe_maps.jl` builds the edge maps, caches them, and exposes `fringe_observables`. A soft edge integrates ``H(z,s)`` by `element_flow`. The physics is [Fringe fields](fringe.md).

## Trackers

`src/tracking/state.jl` is `ParticleState` and `Beam`: a canonical 6-vector and a spin slot that is only stored. `src/tracking/trackers.jl` selects an engine that already exists.

* `FieldTracker` calls `exact_track_turn`.
* `LieTracker` tracks the factors from the Lie product.
* `MapTracker` applies a one-turn function.

`track!` advances a beam or a mutable 6-vector. No tracker evaluates a spin equation ([Spin](spin.md)).

## Readings and everything downstream

`src/analysis/map_readings.jl` reads coefficients. `observables` builds the one-turn Lie factor of the magnetic lattice and then reads it. `dispersion_along` builds each element Lie map and carries an RF cavity as a drift of the same length. `chromaticity_series` on a lattice builds the one-turn differential-algebra map. Those three call map builders. The coefficient readers do not integrate. The physics is [Readings of the one-turn map](map_readings.md).

`src/analysis/sextupole_families.jl` places and fits families from an optics table that is already built ([Sextupole families](sextupoles.md)).

`src/maps/lie/normal_form.jl` is the normal form. It is loaded last and is not part of the facade: one asks for it by name. The physics is [Normal form](normal_form.md).

`src/longitudinal_models/smooth.jl` is the averaged ``(\varphi,\delta)`` Hamiltonian. It does not track the ring in six dimensions ([Longitudinal motion](longitudinal.md)). `src/longitudinal_models/acceleration.jl` is the ramp: a thin cavity and `renormalize`.

`src/optimizers/` fits a scalar function of several real parameters. It does not build a map.
