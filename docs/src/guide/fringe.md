# Fringe-field models

Example: `examples/19_fringe_models`.

The lattice is unchanged; the edge model is passed when the ring is built.

```julia
model = FringeModel(dipole = SoftEdge(profile = Enge(gap = 0.05),
                                      method  = Splitting(n = 20, order = 4)),
                    quad = HardEdge(order = 2),
                    override = Dict("BM3" => HardEdge(order = 1)))
F   = ring_lie_factor(lat, model)       # one-turn Lie factor with edges
obs = fringe_observables(lat, model)    # Qx, Qy, ξ1, ξ2, symplecticity, time
```

| choice | options |
|:--|:--|
| model | `NoEdge()` (default), `HardEdge(order = 1, 2, 3; gap, fint)`, `SoftEdge(profile, method)` |
| profile | `Enge(coeffs; gap)`, `Tanh(λ)`, `LinearRamp(width)`, `Tabulated(z, B)` |
| method | `Splitting(n, order = 2/4)`, `Magnus(n, order = 2/4)`, `TaylorIntegrate(tol)` (reference) |

`FringeModel()` reproduces the plain ring exactly. The edge maps themselves are
in the Developer's Guide, [Fringe fields](../dev/fringe.md).
