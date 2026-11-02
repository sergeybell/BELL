# longitudinal_models

Approximate / expanded longitudinal models, kept for comparison (v0.17).

- `smooth.jl` (`Smooth`): averaged turn-by-turn model in (phi, delta), Q_s formula, buckets.
- `acceleration.jl` (`Acceleration`): energy ramp by thin `rf_kick` + external `renormalize` of the reference after every cavity.

The reference RF model is the `RFCavity` s-Hamiltonian in `src/lie/hamiltonians.jl`; its maps are built like every other element.
