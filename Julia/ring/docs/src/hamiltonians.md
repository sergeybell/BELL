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

## Maps from a Hamiltonian

An ``s``-independent Hamiltonian over length ``L`` generates the Lie map

```math
\mathcal{M}=\exp(-L:H:),\qquad
:f:g=[f,g]=\sum_i\Bigl(\partial_{q_i}f\,\partial_{p_i}g-\partial_{p_i}f\,\partial_{q_i}g\Bigr).
```

Truncation of the Taylor series of ``\mathcal{M}`` at a finite order yields a
polynomial map that is symplectic only up to that order. Factorising

```math
\mathcal{M}=R\circ e^{:f_3:}\circ e^{:f_4:}
```

(Dragt–Finn) and integrating each nonlinear factor by the implicit midpoint
rule restores symplecticity to round-off. An independent construction
integrates Hamilton's equations ``\mathrm{d}z/\mathrm{d}s=S\nabla H`` by a
Runge–Kutta scheme in the truncated polynomial algebra (differential-algebra
map) and is used as a cross-check.