# Normal forms

Examples: `examples/06_lie_maps/run.jl`, `examples/14_normal_form/run.jl`.
The normal-form names need `using Bell.NormalForm`.

```julia
ref = reference_from_kinetic(deuteron, 270.0)
lat = Lattice(lat0.name, lat0.elements, ref, 5)
cav = RFCavity("RF", 0.2, 5, π)                      # thin RF, stationary bucket
F   = one_turn_factor_rf(lat, cav, ref, path_length(lat))
nf  = dragt_finn_nf6d(F)
nf.mu ./ 2π                                          # (Qx, Qy, Qs)
nf6d_tunes(nf, [Jx, Jy, Js]).Q                       # amplitude-dependent tunes
nf6d_delta_eq(nf, [Jx, Jy, Js])                      # equilibrium momentum offset
nf6d_actions(nf, z)                                  # actions of a tracked point
```

`examples/06_lie_maps` checks the Lie product (`concatenate`) against the Lie series.
Steps of the algorithm: [Normal form (internals)](../dev/normal_form.md).
