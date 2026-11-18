# longitudinal_models

Approximate longitudinal models, kept for comparison with the element Hamiltonian. `smooth.jl` (`Smooth`) is an averaged turn-by-turn model in `(φ, δ)`. `acceleration.jl` (`Acceleration`) ramps energy with a thin `rf_kick` and `renormalize` after each cavity.

`Smooth` depends on `Coordinates`. `Acceleration` also depends on `Elements`, `Structure`, and `FieldTrack`. `NormalForm` uses `Smooth`. The `Bell` facade re-exports both. Layer edges: [src/README.md](../README.md).

The reference RF model is the `RFCavity` s-Hamiltonian in `src/physics/hamiltonians.jl`. Another comparison model of this kind goes in this folder.
