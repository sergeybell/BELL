# Symplecticity checks and Taylor-map comparison

*`src/Symplectic.jl`, `src/TaylorMap.jl`, `src/LieOperators.jl`, script
`02_lie_maps.jl`.*

A map ``z\mapsto\mathcal M(z)`` is symplectic iff its Jacobian satisfies
``J^TSJ=S`` with ``S=\mathrm{diag}(S_2,S_2,S_2)``, ``S_2=\begin{pmatrix}0&1\\-1&0\end{pmatrix}``.
The code reports ``\max|J^TSJ-S|`` for

* every element Lie map (at ``z=0`` and at a finite amplitude),
* the ring map built from Lie factors, the factor-by-factor tracking map
  (Jacobian computed exactly by running the tracker on `TaylorN`), and the
  exact integrator (including the cavity).

Truncated polynomial maps are symplectic only up to the truncation order (the
error grows like ``|z|^{\rm MAP\_ORDER}``); the Lie-factor tracker and the exact
integrator are symplectic to round-off.

**Taylor-map comparison module.** `TaylorMap.jl` builds the one-turn map
completely independently of the Lie machinery: Hamilton's equations
``dz/ds=S\nabla H`` are integrated with RK4 in the truncated polynomial
algebra (a COSY-style differential-algebra map), element by element, and the
maps are composed. Its coefficients are compared with the Lie-factor product
and with one-turn tracking. Numbers: printed by `02_lie_maps.jl` (see
[Sextupoles and the optional optimiser](sextupoles.md) for the tables reproduced from the
latest run).

Latest run of `02_lie_maps.jl` / `03_second_order.jl`:

| check | value |
|---|---|
| linear ring matrix: Lie product vs product of ``\exp(LSA)`` | ``2.3\times10^{-13}`` |
| linear ring matrix: Taylor/RK4 vs ``\exp(LSA)`` | ``6.9\times10^{-9}`` |
| max 2nd-order coefficient difference Lie vs Taylor (largest coeff 167.5) | ``5.1\times10^{-8}`` |
| symplectic error ring Taylor map at 0 / at finite amplitude | ``6.8\times10^{-13}`` / ``1.4\times10^{-5}`` |
| one turn: Lie tracking vs Taylor map | ``4.7\times10^{-10}`` |
| symplectic error of 1-turn Lie-tracking Jacobian | ``2.4\times10^{-14}`` |
| one turn: exact integrator vs Lie (order 3) / vs Taylor / vs linear | ``2.4\times10^{-8}`` / ``2.4\times10^{-8}`` / ``6.9\times10^{-5}`` |
| symplectic error of exact-integrator one-turn Jacobian | ``1.1\times10^{-13}`` |
| 1000 turns Lie tracking: relative spread of ``J_x``, ``J_y`` | ``2.4\times10^{-3}``, ``2.7\times10^{-3}`` |
