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
  degree-``(N+1)`` terms of ``H``, therefore ``H`` is built in ``{}_{N+1}D_6``, its vector field ``F=S\nabla H`` is kept through ``N``
  (`change_order`) and the series is summed in ``{}_ND_6`` (`vector_field_flow`): the degree-``N`` part of
  ``F\cdot\nabla g`` needs only ``F`` and ``g`` through ``N``. The result is exact through ``N`` up to roundoff, with the
  full square root ``p_s``, the curvature ``(1+hx)`` and the energy term — nothing is expanded by hand.
* **Scaling and squaring.** The series behaves like ``e^{L\rho}``, ``\rho=\|\partial F/\partial z(0)\|_\infty``.
  For ``L\rho > 3`` (at most one digit of cancellation below that) the map for ``L/2^k`` is computed and squared ``k`` times
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

## Higher-order chromaticity: `chromaticity_series`

`chromaticity_series(lat, N)` (or `chromaticity_series(M)` for a ready DA map) reads the series
``Q(\delta) = \sum_n q_n\delta^n`` directly from the one-turn map of order ``N``, with no finite differences:

1. parametric closed orbit ``u(\delta)``: ``\mathcal M_{1..4}(u, 0, \delta) = u`` is solved in the DA by Newton
   iteration with the constant matrix ``(R_0 - I)^{-1}`` (each step gains one order in ``\delta``);
2. the map is shifted to the orbit, ``\tilde{\mathcal M}(z) = \mathcal M(u + z_\perp, l, \delta) - u``, and its linear part
   ``R_{ij}(\delta) = \partial\tilde{\mathcal M}_i/\partial z_j|_{z_\perp=0}`` is a DA in ``\delta``, exact through
   ``\delta^{N-1}`` (the derivative of an order-``N`` map is exact through ``N-1``);
3. ``\cos\mu = \tfrac12\operatorname{tr}R`` is solved for ``\mu(\delta)`` by Newton in the DA, started on the branch
   ``\mu_0\in(0,\pi)`` if ``R_{12}(0)>0`` (else ``2\pi-\arccos``); ``Q = \mu/2\pi``.

It returns ``q_0..q_{N-1}`` and ``\xi_n = n!\,q_n`` (``n = 1..N-1``); planes are treated as uncoupled.
The `Lattice` method drops thin `RFCavity` elements and uses [`one_turn_map_flow`](@ref), i.e. no fringe
fields and no dipole edge angles. The old `chromaticity_orders` (central difference with step `dd` on the
order-3 map) has an error ``\approx q_3^{\rm trunc}\,dd^2`` in ``\xi_1``, where ``q_3^{\rm trunc}`` is the unreliable
cubic coefficient of the truncated map; `chromaticity_series` has no such error.