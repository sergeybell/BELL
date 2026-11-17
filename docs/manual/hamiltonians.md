# Element Hamiltonians

With ``p_s=\sqrt{(1+\delta)^2-p_x^2-p_y^2}``, curvature ``h=1/\rho`` and
normalised multipoles ``k_1=B'/(B\rho)``, ``k_2=B''/(B\rho)``, the
``s``-Hamiltonians (hard-edge, no fringe fields) are

| element | ``H(z)`` |
|---|---|
| Drift | ``e(\delta)-p_s`` |
| Quadrupole | ``e(\delta)-p_s+\tfrac{k_1}{2}(x^2-y^2)`` |
| Sextupole | ``e(\delta)-p_s+\tfrac{k_2}{6}(x^3-3xy^2)`` |
| Sector bend | ``e(\delta)-(1+hx)p_s+hx+\tfrac12 h^2 x^2+\tfrac{k_1}{2}(x^2-y^2)`` |
| Marker | ``0`` |
| Finite RF cavity | ``e(\delta)-p_s+A\cos\Phi(l,\delta,s)`` (see [Finite-length RF cavity](rf.md)) |

For the bend, ``\partial H/\partial x|_{z=0}=0``: the reference arc is a
solution of Hamilton's equations.

The bend potential is cut at ``hx+\tfrac12 h^2 x^2``. There are no fringe
fields. A gradient in a bend uses the straight-magnet term; curvature
corrections of order ``h k_1 x^3`` are omitted (the lattice used here has
``k_1=0`` in every bend).

The Hamiltonian of each magnet is replaced by its Taylor polynomial of
total degree 4 (square root and energy expanded). That polynomial is the
generator. It is not yet a map.

## Maps from a Hamiltonian

An ``s``-independent Hamiltonian over length ``L`` generates the Lie
transformation ``\exp(-L:H:)``, with

```math
:f:g=[f,g]=\sum_i\Bigl(\partial_{q_i}f\,\partial_{p_i}g-\partial_{p_i}f\,\partial_{q_i}g\Bigr).
```

How the ring is assembled from these transformations is in
[Lie product of the ring](lie.md).