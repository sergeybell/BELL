# Optics: Twiss, dispersion, tunes, chromaticity

```@meta
CurrentModule = Bell
```

Modules: `Linear` (`src/linear.jl`), `HigherOrder` (`src/higher_order.jl`).

## Is the full Hamiltonian built and then truncated for linear optics? — Yes, in this sense:

1. For every element the **full** Hamiltonian (exact ``p_s``, exact ``e(\delta)``,
   bend and multipole terms) is evaluated in `TaylorN` arithmetic and stored as a
   polynomial of degree 4 ([`hamiltonian`](@ref)).
2. **Linear optics uses only its quadratic part.** [`linear_matrix`](@ref) takes
   the Hessian at the origin, ``A_{ij} = \partial^2 H/\partial z_i\partial z_j|_{z=0}``
   ([`hessian`](@ref)), i.e. ``H_2 = \tfrac12 z^\mathsf{T}Az``, and returns the exact
   matrix exponential ``R = \exp(L\,SA)``. The cubic and quartic parts are ignored
   here. There are no textbook element matrices; the closed forms on the
   [Elements](elements.md) page are what this exponential equals.
3. [`ring_matrix`](@ref) multiplies ``R_n\cdots R_1``. Equivalently,
   [`observables`](@ref) builds the full nonlinear one-turn Lie factor
   [`ring_lie_factor`](@ref) and reads the Twiss parameters from its linear part
   `F.R` (identical to `ring_matrix` because concatenation multiplies the ``R``s).

So: the full nonlinear object exists, and the linear optics is its degree-1 (map)
/ degree-2 (Hamiltonian) reading. Chromaticity is **not** a linear-optics
quantity and is read from the nonlinear map (below).

## Linear readings of the one-turn matrix ([`periodic_optics`](@ref))

Per plane, from the 2×2 diagonal block ``m``:

```math
\cos\mu = \tfrac12(m_{11}+m_{22}),\quad \sin\mu = \operatorname{sign}(m_{12})\sqrt{1-\cos^2\mu},\quad
\beta = \frac{m_{12}}{\sin\mu},\quad \alpha = \frac{m_{11}-m_{22}}{2\sin\mu},\quad Q_\text{frac} = \frac{\mu}{2\pi}.
```

Periodic dispersion (fixed point of the 4D map at fixed ``\delta``) and momentum compaction:

```math
D_4 = (D_x, D_x', D_y, D_y') = (I - R_{1:4,1:4})^{-1}R_{1:4,6},\qquad
(l|\delta)_\text{closed} = R_{56} + R_{5,1:4}\cdot D_4,\qquad
\alpha_c = \frac{1}{\gamma_0^2} - \frac{(l|\delta)_\text{closed}}{C}.
```

[`optics_along`](@ref) cuts each element into ``\lceil L/ds\rceil`` slices with
``R = \exp((L/n)SA)``, transports ``(\beta,\alpha)`` with the 2×2 Twiss
transformation, accumulates the phase ``\Delta\mu = \operatorname{atan2}(m_{12},\, m_{11}\beta - m_{12}\alpha)``
(this gives the integer part of the tune), multiplies the dispersion vector
``(D_4, 0, 1)`` by each slice matrix, and reports the closure after one turn.

## Nonlinear readings of the one-turn map (`HigherOrder`)

[`observables`](@ref) does, in this order:

```
F  = ring_lie_factor(lat)          # product of element factors, generators to degree 4
M  = factor_to_map(F)              # polynomial R ∘ e^{:f3:} ∘ e^{:f4:}, cut at degree 3
po = periodic_optics(F.R, ref, C)  # linear part only
pe = path_expansion(M, ref, C)     # α0, α1, η0, η1
xi = chromaticity_orders(M)        # Q0, ξ1, ξ2 per plane
```

### How ξ is extracted ([`chromaticity_orders`](@ref))

Not from a formula like ``-\frac{1}{4\pi}\oint\beta k\,ds``; the tune is
measured on the polynomial map at three momenta and differentiated numerically:

```
for δ in (−dd, 0, +dd)            # dd = 1e-4
    u(δ) = closed_orbit(M, δ)      # Newton: M_{1:4}(u,0,δ) = u, exact polynomial Jacobian, |Δu| < 1e-16
    J    = ∂M/∂z at (u, 0, δ)      # exact derivative of the polynomial
    Q(δ) = tune of 2×2 blocks of J (cos 2πQ = tr/2, sign from J12)
ξ1 = (Q(+dd) − Q(−dd)) / (2 dd)
ξ2 = (Q(+dd) − 2Q(0) + Q(−dd)) / dd²       # Q(δ) ≈ Q0 + ξ1 δ + ξ2 δ²/2
```

Because the Jacobian is taken *on the off-momentum closed orbit*, ``\xi_1``
contains the chromatic terms of the drift/quadrupole/bend Hamiltonians
(e.g. ``-\tfrac12\delta p_x^2``, the ``\delta``-dependence of ``R_x(\delta)``) and the
sextupole feed-down ``k_2D\delta``. All Jacobian entries up to ``\delta^2`` come
from map monomials of degree ≤ 3, so ``\xi_1`` is exact for the degree-3 map and
``\xi_2`` contains every term present in it. The finite-difference error is
``O(dd^2)\sim 10^{-8}`` relative.

### Second-order dispersion and slip ([`dispersion_orders`](@ref), [`path_expansion`](@ref))

Closed orbit ``u(\delta) = D_1\delta + D_2\delta^2``:

```math
D_1 = (I-R_4)^{-1}R_{1:4,6},\qquad D_2 = (I-R_4)^{-1}\,q_{1:4}(w),\quad w = (D_1, 0, 1),
```

with ``q_i`` the quadratic part of ``M_i``. Path length on that orbit:
``l(\delta) = l_1\delta + l_2\delta^2``, ``l_1 = R_{5,1:4}D_1 + R_{56}``,
``l_2 = R_{5,1:4}D_2 + q_5(w)``, and

```math
\alpha_0 = -\frac{l_1}{C} + c_1,\quad \alpha_1 = -\frac{l_2}{C} + c_2,\quad
\eta_0 = -\frac{l_1}{C},\quad \eta_1 = -\frac{l_2}{C} + \frac{l_1}{C}c_1 .
```
