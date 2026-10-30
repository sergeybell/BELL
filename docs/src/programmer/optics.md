# Optics: Twiss, dispersion, chromaticity

```@meta
CurrentModule = Bell
```

Modules: `Linear` (`src/linear.jl`), `HigherOrder` (`src/higher_order.jl`).

## Linear optics

The linear matrix of every element follows from its Hamiltonian:
``R = \exp(L\,S A)`` with ``A = \partial^2 H/\partial z^2`` ([`hessian`](@ref)).
There are no textbook element matrices. [`ring_matrix`](@ref) multiplies the
element matrices, and [`periodic_optics`](@ref) reads the one-turn matrix:

- tunes and Twiss parameters from the 2×2 blocks (``\cos\mu = \tfrac12\mathrm{tr}``,
  ``\beta = R_{12}/\sin\mu``, ``\alpha = (R_{11}-R_{22})/(2\sin\mu)``);
- the periodic dispersion ``D_4 = (I-R_4)^{-1}R_{1:4,6}``;
- the momentum compaction ``\alpha_c = 1/\gamma_0^2 - (l|\delta)_\text{closed}/C``.

[`optics_along`](@ref) transports these functions through slices and
reports the closure after one turn.

```@docs
hessian
linear_matrix
ring_matrix
periodic_optics
optics_along
```

## Readings of the one-turn map

`HigherOrder` reads a polynomial one-turn map; it does not integrate anything.
[`observables`](@ref) is the standard entry point. It builds
[`ring_lie_factor`](@ref), converts it with [`factor_to_map`](@ref), and returns
the periodic optics, path expansion and chromaticities from that one object.

- **Closed orbit at δ.** [`closed_orbit`](@ref) uses Newton iteration with the
  exact polynomial Jacobian.
- **Tunes at δ.** [`tunes_at_delta`](@ref) uses the Jacobian on that orbit.
- **Chromaticity.** [`chromaticity_orders`](@ref) takes central differences of
  Q(δ), so ξ1 and ξ2 include sextupole feed-down on the off-momentum orbit.
- **Dispersion to 2nd order.** [`dispersion_orders`](@ref) solves the fixed-point
  equation order by order in δ.
- **Path length and slip.** [`path_expansion`](@ref) uses the l-row of the map
  on the dispersion orbit together with the velocity expansion.

```@docs
observables
closed_orbit
tunes_at_delta
chromaticity_orders
dispersion_orders
dispersion_along
path_expansion
velocity_expansion
```
