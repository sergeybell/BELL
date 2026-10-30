# Tracking and acceleration

Examples: `examples/11_lie_maps`, `examples/12_higher_order`,
`examples/14_rf_acceleration`, `examples/16_blond_longitudinal`.

**Exact 6D tracking** (non-expanded Hamiltonians, the reference):

```julia
z = [1e-3, 0, 0, 0, 0, 1e-4]                 # canonical (x, px, y, py, l, δ)
for n in 1:1000
    z = exact_track_turn(z, lat)
end
```

**Fast symplectic tracking** with element Lie factors:

```julia
factors = build_factors(lat)
z = track(factors, z; turns = 1000)
```

**Acceleration** with a finite cavity (example 14):

```julia
ref_inj = reference_from_kinetic(DEUTERON, 10.0)
lat_rf, _ = insert_cavity(lat, Cavity("RF", 1.0, 0.2, 5, deg2rad(150)))
lat_rf = Lattice(lat_rf.name, lat_rf.elements, ref_inj, 5)
users = [zeros(6), [2e-3, 0, 0, 0, 0, 0]]     # user coords; particle 1 = synchronous
rec, flat = accelerate(lat_rf, ref_inj, users; Ekin_target = 270.0, flat_turns = 1000)
rec.Ekin, rec.u[2]                            # energy per turn, particle 2 per turn
```

The smooth longitudinal model (`LongitudinalParams`, `bucket_height`,
`synchrotron_tune`) is for comparing tracking with the equations; it is not
used by tracking. Algorithms: [Tracking (internals)](../dev/tracking.md),
[RF & longitudinal](../dev/rf.md).
