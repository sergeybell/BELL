# Tracking and acceleration

Examples: `examples/06_lie_maps/run.jl`, `examples/11_trackers_and_beam/run.jl`,
`examples/12_rf_longitudinal/run.jl`, `examples/13_acceleration/run.jl`.

**Point tracking.** A magnet and a thick cavity are [`element_flow`](@ref) of
the s-Hamiltonian. A thin cavity is [`rf_kick`](@ref).

```julia
z = [1e-3, 0.0, 0.0, 0.0, 0.0, 1e-4]    # canonical (x, px, y, py, l, δ)
for _ in 1:1000
    z = exact_track_turn(z, lat)
end
```

**Lie-factor tracking:**

```julia
factors = build_factors(lat)
z = track(factors, z; turns = 1000)
```

**Acceleration** (`examples/13_acceleration/run.jl`). The cavity is a thin
[`RFCavity`](@ref). [`insert_cavity`](@ref) splits a named drift and keeps
[`path_length`](@ref) when the cavity length is zero. Particles are laboratory
coordinates `(x, x', y, y', φ, δ)`. Particle 1 is the synchronous particle.
[`accelerate`](@ref) converts with [`laboratory_to_canonical`](@ref) and
records [`canonical_to_laboratory`](@ref).

```julia
ref_inj = reference_from_kinetic(deuteron, 10.0)
drift = lat.elements[findfirst(e -> e isa Drift, lat.elements)]
lat_rf, _ = insert_cavity(lat, RFCavity("RF", 0.2, 5, deg2rad(150));
                          in_element = drift.name)
lat_rf = Lattice(lat_rf.name, lat_rf.elements, ref_inj, 5)
users = [zeros(6), [2e-3, 0.0, 0.0, 0.0, 0.0, 0.0]]
rec, flat = accelerate(lat_rf, ref_inj, users; Ekin_target = 270.0, flat_turns = 1000)
rec.Ekin, rec.u[2]          # energy per turn, particle 2 per turn
```

The smooth longitudinal model (`LongitudinalParams`, `bucket_height`,
`synchrotron_tune`) is a comparison. [`accelerate`](@ref) does not call it.
Algorithms: [Tracking (internals)](../dev/tracking.md),
[RF Hamiltonian and acceleration](../dev/rf_acceleration.md).
