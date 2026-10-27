# Coordinate transformation user <-> canonical

*Implemented in `src/Constants.jl` (`canonical_to_user`, `user_to_canonical`),
used in `src/TaylorMap.jl` (`user_map`) and `src/Acceleration.jl`; checked by
`scripts/06_coordinates.jl`.*

## Exact formulas

Canonical ``\to`` user:

```math
p_s=\sqrt{(1+\delta)^2-p_x^2-p_y^2},\qquad x'=\frac{p_x}{p_s},\quad y'=\frac{p_y}{p_s},\qquad
\phi=\frac{2\pi h}{C}\,\frac{l}{v(\delta)/v_0}.
```

User ``\to`` canonical:

```math
p_x=\frac{x'(1+\delta)}{\sqrt{1+x'^2+y'^2}},\quad p_y=\frac{y'(1+\delta)}{\sqrt{1+x'^2+y'^2}},\qquad
l=\frac{C}{2\pi h}\,\frac{v(\delta)}{v_0}\,\phi .
```

``x``, ``y`` and ``\delta`` are common to both sets. These follow from
``x'=P_x/P_s`` and from ``\phi=-h\omega_0\Delta t``, ``l=-v\Delta t``.

## Applied as full polynomials, not linearised

Both functions are generic Julia code. In `user_map` they are evaluated on the
`TaylorN` variables, which produces the complete Taylor expansion of the
square roots and of ``1/v(\delta)``; the result is truncated only at the map
order (`MAP_ORDER = 3`). The map in user coordinates is

```math
U = T\circ M\circ T^{-1},\qquad T = \text{canonical\_to\_user}.
```

Tracking (`accelerate`, `turn_snapshots`) uses the same functions on
floating-point numbers, i.e. exactly. The code was checked for hidden
linearisations: none were found. (The linear-optics module reports
``D'=dx/d\delta``-slope in ``p_x``; to first order ``p_x=x'``, so this is not an
approximation at the order it is used.)

Examples printed by `06_coordinates.jl` (order 3):

```math
x' = p_x - p_x\delta + \tfrac12 p_x^3 + \tfrac12 p_x p_y^2 + p_x\delta^2+\dots,\qquad
\phi = 0.267624\,l - 0.204533\,l\delta + 0.228642\,l\delta^2+\dots
```

## Symplecticity holds only in canonical variables

``(x,x',y,y',\phi,\delta)`` is **not** a canonical set (``x'`` depends on
``\delta`` and ``p_y``; ``\phi`` is scaled by ``1/v(\delta)``). Therefore the
symplectic condition ``J^TSJ=S`` is checked only for maps in ``z``. The
user-coordinate map is the same physics written in other variables.

## Numeric check (scripts/06_coordinates.jl)

| check | result |
|---|---|
| polynomial round trip ``T(T^{-1}(u))-u``, all orders kept | max coefficient ``\le 1.1\times10^{-16}`` |
| exact functions round trip at ``u=(5\,\mathrm{mm},20\,\mathrm{mrad},5\,\mathrm{mm},20\,\mathrm{mrad},1\,\mathrm{rad},10^{-2})`` | ``0`` |
| order-3 polynomial ``T`` vs exact ``T`` at that point | ``8.7\times10^{-7}`` (truncation) |
| linearised ``T`` vs exact ``T`` at that point | ``7.6\times10^{-3}`` |
| ``\max\vert J^TSJ-S\vert `` canonical one-turn map at ``10^{-3}`` | ``6.8\times10^{-6}`` (order-3 truncation) |
| ``\max\vert J^TSJ-S\vert `` user one-turn map at ``10^{-3}`` | ``0.12`` (non-canonical variables) |

Selected one-turn coefficients (270 MeV, RF off):

| coefficient (canonical / user) | canonical | user |
|---|---|---|
| ``(x\vert p_x)`` / ``(x\vert x')`` | 4.61348258 | 4.61348258 |
| ``(x\vert \delta)`` | 2.42044391e-2 | 2.42044391e-2 |
| ``(l\vert \delta)`` / ``(\phi\vert \delta)`` | 62.6600041 | 16.7693216 |
| ``(l\vert \delta^2)`` / ``(\phi\vert \delta^2)`` | -71.5275409 | -31.9585141 |
| ``(l\vert p_x^2)`` / ``(\phi\vert x'^2)`` | -86.2567313 | -23.0843723 |
| ``(x\vert p_x^2\delta)`` / ``(x\vert x'^2\delta)`` | 343.450443 | 343.541752 |

The last line shows a genuinely nonlinear difference between the systems.
