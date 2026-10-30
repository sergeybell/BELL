# Optics and chromaticity

Examples: `examples/10_linear_optics`, `examples/12_higher_order`,
`examples/13_sextupole_optimization`.

```julia
po  = periodic_optics(ring_matrix(lat), lat.ref, circumference(lat))
po.x.beta, po.x.alpha, po.D4, po.alpha_c          # closed Twiss, dispersion, α_c
opt = optics_along(lat, po)                        # β, α, D, μ along s; total tunes opt.Qx, opt.Qy

obs = observables(lat)                             # one-turn Lie factor and its readings
obs.xi.x, obs.xi.y                                 # (Q0, ξ1, ξ2) per plane
obs.pe.alpha0, obs.pe.alpha1, obs.pe.eta0, obs.pe.eta1
```

Second-order dispersion along the ring: `dispersion_orders(obs.M)` and
`dispersion_along(lat, D1, D2)`. Sextupole correction: example 13.
How these numbers are computed: [Optics (internals)](../dev/optics.md).
