# Lie flow-operator maps (`Bell.LieFlow`)

A second map route to **arbitrary order** ``N``, independent of (and not replacing) the
Dragt–Finn / Lie-factor route of `lie_algebra.jl`, which expands ``H`` to degree 4 and gives
maps of order 3.

```julia
using Bell, Bell.DiffAlgebra
ref = reference_from_kinetic(DEUTERON, 270.0)
M  = lie_flow_map(Quadrupole("QF", 0.3, 1.5), ref, 7)   # Vector{DA}, z = (x, px, y, py, l, δ)
Mr = one_turn_map_flow(lat, 5)                           # ring, order 5
da_symplectic_error(Mr)                                  # all orders, not only linear
```

## Mathematics (DA layer, `src/da/lie_series.jl`)
For an s-independent Hamiltonian ``H`` the map over the length ``L`` is the Lie series
```math
\mathcal M = e^{-L{:}H{:}}\,z = \sum_{n\ge0}\frac{(-L)^n}{n!}{:}H{:}^n z
           = \sum_{n\ge0}\frac{L^n}{n!}(F\cdot\nabla)^n z,
\qquad {:}H{:}g = [H,g],\quad F = S\nabla H ,
```
because ``-{:}H{:}g = [g,H] = F\cdot\nabla g``. `lie_flow(H, L)` sums it term by term
(``u_n = \tfrac{L}{n}F\cdot\nabla u_{n-1}``) until two consecutive terms are below
``\varepsilon\,\|\mathcal M\|`` (or a term vanishes, e.g. the drift).

* **Exactness.** The reference orbit is a solution, so ``H`` has no linear terms and ``F(0)=0``:
  ``F\cdot\nabla`` does not lower the degree. The degree-``N`` coefficients of the map need the
  degree-``(N+1)`` terms of ``H``, therefore ``H`` is built in ``{}_{N+1}D_6`` and the map is kept
  through ``N`` (`change_order`). The result is exact through ``N`` up to roundoff, with the
  full square root ``p_s``, the curvature ``(1+hx)`` and the energy term — nothing is expanded by hand.
* **Scaling and squaring.** The series behaves like ``e^{L\rho}``, ``\rho=\|\partial F/\partial z(0)\|_\infty``.
  For ``L\rho > 1/2`` the map for ``L/2^k`` is computed and squared ``k`` times
  (``\mathcal M_{L} = \mathcal M_{L/2^k}^{\circ 2^k}``, exact for a one-parameter group); if the series
  still does not converge in `nmax` terms, ``k`` is increased.
* `poisson_bracket`, `hamiltonian_vector_field`, `lie_series` (``e^{t{:}H{:}}g`` for any ``g``) and
  `da_symplectic_error` (max coefficient of ``J^TSJ-S`` through order ``N-1``, ``J=\partial\mathcal M/\partial z``)
  are generic (any number of canonical pairs; unpaired variables are parameters).
* `compose` for DA arguments now builds the monomial table in place, skips monomials that do
  not occur and those of degree ``> N`` when ``N(0)=0``, and accumulates in place (same values;
  3–4× faster than the previous generic sum, which is still used for numeric arguments).

## Elements (`src/lie_flow.jl`)
`hamiltonian_da(e, ref, sp)` evaluates the same exact ``H`` as `hamiltonian` (Drift, Cavity
with RF off, Quadrupole, Sextupole, SBend body, Marker) in a DA space of any order.
`lie_flow_map(e, ref, N)` → map in `flow_space(N)` (shared cached space, so maps compose);
`element_flow_maps(lat, N)` caches per (element, p0, N); `one_turn_map_flow(lat, N)` composes them.

**Not included:** SBend edge angles ``e_1,e_2`` and fringe fields (s-dependent: use the fringe
models / an ODE integrator), thin RF cavities (error).

## Tests (`test/runtests.jl`, testset *LieFlow*)
Order 3 equals the Dragt–Finn/Lie route `lie_map(hamiltonian(e), L)` to ``\le 10^{-15}`` for all
four element types; drift map equals ``x + L p_x/p_s``, ``l + L(e'(\delta) - (1+\delta)/p_s)``
exactly; quadrupole linear part equals the cos/cosh matrices to ``10^{-16}``; quadrupole map
equals the DOP853 solution of ``z' = S\nabla H`` (order 5) to ``10^{-12}``; ``J^TSJ - S`` at roundoff
for ``N = 3, 6`` and ``L`` up to 10 m (scaling-and-squaring path); one-turn FODO map at order 3
equals the composed Lie maps to ``10^{-12}``.
Benchmarks against COSY INFINITY: SCT_2 `benchmarks/lie_flow`.