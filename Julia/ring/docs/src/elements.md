# Element Hamiltonians and Lie maps

*Types: `src/Elements.jl`; Hamiltonians: `src/Hamiltonians.jl`; Lie operators:
`src/LieOperators.jl`; factor-by-factor tracking: `src/Tracking.jl`.*

## Hamiltonians

With ``p_s=\sqrt{(1+\delta)^2-p_x^2-p_y^2}``, curvature ``h=1/\rho`` and
normalised strengths ``k_1=B'/(B\rho)``, ``k_2=B''/(B\rho)``, the
``s``-Hamiltonians are (no fringe fields, hard-edge, field constant in ``s``):

| element | Hamiltonian ``H(z)`` |
|---|---|
| Drift | ``e(\delta)-p_s`` |
| Quadrupole | ``e(\delta)-p_s+\tfrac{k_1}{2}(x^2-y^2)`` |
| Sextupole | ``e(\delta)-p_s+\tfrac{k_2}{6}(x^3-3xy^2)`` |
| Sector bend (with gradient) | ``e(\delta)-(1+hx)p_s+hx+\tfrac{h^2x^2}{2}+\tfrac{k_1}{2}(x^2-y^2)`` |
| Marker | ``0`` |
| Thin RF kick (`RFCavity`) | kick only, see [RF cavity of finite length](rf.md) |
| Finite cavity (`Cavity`) | ``e(\delta)-p_s+A\cos\Phi(l,\delta,s)`` (exact integrator only; as drift in the maps) |

For the bend, ``\partial H/\partial x|_{z=0}=-h+h=0``: the reference arc is a
solution. All Hamiltonians are expanded as `TaylorN` polynomials to
`JET_ORDER = 4` in ``z``, and the constant term is removed.

## Map of an element

For an ``s``-independent ``H`` over length ``L`` the map is the Lie transformation

```math
\mathcal M = \exp(-L:H:),\qquad :f:g = [f,g]=\sum_{i}\left(\frac{\partial f}{\partial q_i}\frac{\partial g}{\partial p_i}-\frac{\partial f}{\partial p_i}\frac{\partial g}{\partial q_i}\right),
```

applied to the coordinate functions as the series
``\exp(:f:)z=\sum_n :f:^n z/n!``, truncated at `MAP_ORDER = 3`.
The map is then Dragt–Finn factorised,

```math
\mathcal M = R\ \circ\ e^{:f_3:}\ \circ\ e^{:f_4:},
```

``R`` the linear (symplectic) matrix and ``f_3,f_4`` homogeneous cubic/quartic
generators. Tracking (`Tracking.jl`) applies ``R`` exactly and each
``e^{:f_k:}`` by integrating the Hamiltonian flow of ``-f_k`` with the implicit
midpoint rule, so the tracking map is symplectic to round-off even though
the generators are truncated.
