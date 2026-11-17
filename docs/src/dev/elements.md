# Elements: what map each element produces

```@meta
CurrentModule = Bell
```

This page answers, for every element type, three questions:

1. which **Hamiltonian** the code uses (exact expression, and its expansion as it
   is actually stored),
2. which **objects** are computed from it (6×6 matrix, polynomial map,
   Dragt–Finn generators, exact flow), written out explicitly up to ``\delta^2``
   and beyond where useful,
3. what **adding the element to a lattice** produces as a data structure.

All numbers quoted below were printed by the package itself (deuteron,
``E_\text{kin}=270`` MeV, ``\gamma_0 = 1.14395``, ``\beta_0 = 0.48563``) and agree
with the closed formulas to round-off.

## Variables and the two expansions every element shares

Canonical variables ``z = (x, p_x, y, p_y, l, \delta)``, independent variable
``s``, ``\delta = \Delta p/p_0`` (momentum, not energy), ``l = -v(\delta)\,\Delta t``
(``l>0``: particle ahead of the reference). Every element Hamiltonian contains

```math
p_s = \sqrt{(1+\delta)^2 - p_x^2 - p_y^2},\qquad
e(\delta) = \frac{1}{\beta_0}\sqrt{(1+\delta)^2 + \frac{1}{\beta_0^2\gamma_0^2}} - \frac{1}{\beta_0^2}
```

(`energy_term`). Their Taylor expansions, which the polynomial engine produces
automatically when the exact expression is evaluated in `TaylorN` arithmetic, are

```math
e(\delta) = \delta + \frac{\delta^2}{2\gamma_0^2} - \frac{\beta_0^2\,\delta^3}{2\gamma_0^2}
 + \Big(-\frac{\beta_0^2}{8} + \frac{3\beta_0^4}{4} - \frac{5\beta_0^6}{8}\Big)\delta^4 + O(\delta^5),
```

```math
-p_s = -(1+\delta) + \frac{p_x^2+p_y^2}{2(1+\delta)} + \frac{(p_x^2+p_y^2)^2}{8(1+\delta)^3} + \dots
```

so that the field-free part of every Hamiltonian, after the constant is dropped and
the linear terms cancel (``\delta - \delta``), is

```math
H_\text{drift} = \underbrace{\frac{p_x^2+p_y^2}{2} + \frac{\delta^2}{2\gamma_0^2}}_{H_2}
\;\underbrace{-\,\frac{(p_x^2+p_y^2)\,\delta}{2} - \frac{\beta_0^2\delta^3}{2\gamma_0^2}}_{H_3}
\;+\;\underbrace{\frac{(p_x^2+p_y^2)^2}{8} + \frac{(p_x^2+p_y^2)\,\delta^2}{2} + e_4\,\delta^4}_{H_4}.
```

The derivative ``\partial e/\partial\delta = v(\delta)/v_0`` (`velocity_ratio`), so
``dl/ds`` always contains the exact velocity:
``v/v_0 = 1 + c_1\delta + c_2\delta^2 + c_3\delta^3``,
``c_1 = 1/\gamma_0^2``, ``c_2 = -3\beta_0^2/(2\gamma_0^2)``,
``c_3 = \beta_0^2(5\beta_0^2-1)/(2\gamma_0^2)`` ([`velocity_expansion`](@ref)).

## The pipeline: from element to maps

Nothing is computed when an element is created; an element is a plain immutable
struct. Every map is derived on demand from **one** function,
[`hamiltonian`](@ref), and cached by `(element value, p0)`:

```
element e ──hamiltonian(e, ref)──► H(z)  TaylorN, exact expression expanded to degree JET_ORDER = 4
   │                                 │
   │                                 ├─ linear_matrix(e, ref)   R = exp(L · S · Hess H)      (6×6 Float64)
   │                                 ├─ element_lie_map(e, ref) M = exp(−L :H:) z, cut at degree MAP_ORDER = 3
   │                                 │        └─ dragt_finn_factor(M) → LieFactor(R, f3, f4)
   │                                 └─ taylor_map_element(e)   RK4 of dz/ds = S∇H in TaylorN (cross-check)
   └─ exact_track_element(z, e, ref, C)  integrates the EXACT (non-expanded) vector field
```

* ``R = \exp(L\,SA)`` with ``A_{ij} = \partial^2H/\partial z_i\partial z_j|_0``
  is the exact flow of the quadratic part ``H_2`` (matrix exponential, not a
  truncated series).
* ``M = e^{-L{:}H{:}}z = z - L[H,z] + \tfrac{L^2}{2}[H,[H,z]] - \dots`` with the
  degree-4 ``H``; the series is summed until it converges and the result is cut at
  total degree 3. This is the exact flow of the full ``H`` through degree 3 (all
  terms ``z_iz_jz_k``, e.g. ``x\,\delta^2``, ``p_x\delta^2``, ``\delta^3``).
* The Dragt–Finn factor writes the same map as
  ``\mathcal M = R\circ e^{{:}f_3{:}}\circ e^{{:}f_4{:}}``
  (applied to a point: first ``e^{:f_4:}``, then ``e^{:f_3:}``, then ``R``).

### What "adding an element to a lattice" produces

| step | data structure |
|:--|:--|
| `Drift("D", 0.5)` or a card `D drift L = 0.5` | an immutable `Drift` value (name, L) — no map yet |
| `load_structure` / `Lattice(name, elements, ref, h)` | `Lattice.elements::Vector{AbstractElement}`: the expanded ring in tracking order; repeated cells are the *same* values, so each distinct element is computed once |
| `build_factors(lat)` | `Vector{Any}`: one `LieFactor(name, R, f3, f4, field, field4)` per element, the thin `RFCavity` kept as the struct itself |
| `ring_lie_factor(lat)` | one `LieFactor` of the whole ring (BCH product of the element factors, RF skipped) |
| `ring_matrix(lat)` | ``R_n\cdots R_1`` (6×6) |
| `fringe_sequence(lat, model)` | `Vector{Any}` with `EdgeMap(name_in, parent, :in, model)`, parent, `EdgeMap(name_out, …)` around magnets that have an edge model |

## Drift

**Hamiltonian** (`hamiltonian(::Drift)`): ``H = e(\delta) - p_s``, expansion
``H_2+H_3+H_4`` above.

**Exact flow** (used by `exact_track_element`, one implicit-midpoint Yoshida step;
the field is independent of ``x,y,l`` so the step reproduces the exact solution):

```math
x_f = x + L\frac{p_x}{p_s},\quad y_f = y + L\frac{p_y}{p_s},\quad
l_f = l + L\Big(\frac{v(\delta)}{v_0} - \frac{1+\delta}{p_s}\Big),\quad
p_{x,y},\ \delta\ \text{unchanged}.
```

**Linear matrix** (`linear_matrix`):

```math
R_\text{drift} = \begin{pmatrix}
1 & L & 0 & 0 & 0 & 0\\ 0 & 1 & 0 & 0 & 0 & 0\\
0 & 0 & 1 & L & 0 & 0\\ 0 & 0 & 0 & 1 & 0 & 0\\
0 & 0 & 0 & 0 & 1 & L/\gamma_0^2\\ 0 & 0 & 0 & 0 & 0 & 1 \end{pmatrix}
```

**Polynomial map to degree 3** (`element_lie_map`, exactly what the code stores):

```math
\begin{aligned}
x_f &= x + L p_x - L p_x\delta + L p_x\delta^2 + \tfrac{L}{2}p_x(p_x^2+p_y^2),\\
y_f &= y + L p_y - L p_y\delta + L p_y\delta^2 + \tfrac{L}{2}p_y(p_x^2+p_y^2),\\
l_f &= l + L\big(c_1\delta + c_2\delta^2 + c_3\delta^3\big) - \tfrac{L}{2}(p_x^2+p_y^2) + L(p_x^2+p_y^2)\,\delta,
\end{aligned}
```

i.e. at fixed ``\delta`` the transverse matrix is ``\begin{pmatrix}1 & L/(1+\delta)\\0&1\end{pmatrix}``
expanded to ``\delta^2``.

**Dragt–Finn generators.** All terms of ``H_\text{drift}`` depend only on
``(p_x, p_y, \delta)`` and commute, hence exactly

```math
f_3 = -L\,H_3 = \tfrac{L}{2}(p_x^2+p_y^2)\delta + \tfrac{L\beta_0^2}{2\gamma_0^2}\delta^3,\qquad
f_4 = -L\,H_4 .
```

## Quadrupole

**Hamiltonian**: ``H = e(\delta) - p_s + \tfrac{k_1}{2}(x^2-y^2)``, i.e.
``H_\text{drift} + \tfrac{k_1}{2}(x^2 - y^2)`` (``k_1>0`` focuses in ``x``).

**Linear matrix** (``\omega = \sqrt{k_1}``, ``k_1>0``; for ``k_1<0`` exchange the
planes):

```math
R_x = \begin{pmatrix}\cos\omega L & \frac{\sin\omega L}{\omega}\\ -\omega\sin\omega L & \cos\omega L\end{pmatrix},\quad
R_y = \begin{pmatrix}\cosh\omega L & \frac{\sinh\omega L}{\omega}\\ \omega\sinh\omega L & \cosh\omega L\end{pmatrix},\quad
R_{56} = \frac{L}{\gamma_0^2}.
```

**Chromatic (δ-dependent) part.** The part of ``H`` that is quadratic in
``(x,p_x)`` is, *exactly in* ``\delta``,
``\frac{p_x^2}{2(1+\delta)} + \frac{k_1}{2}x^2``. Hence at fixed ``\delta`` the
transverse matrix is, with ``\kappa = k_1/(1+\delta)``, ``\omega_\delta = \sqrt{\kappa}``,

```math
R_x(\delta) = \begin{pmatrix}\cos\omega_\delta L & \dfrac{\sin\omega_\delta L}{\omega_\delta(1+\delta)}\\[2mm]
-\omega_\delta(1+\delta)\sin\omega_\delta L & \cos\omega_\delta L\end{pmatrix}.
```

The polynomial map stored by the code contains the Taylor coefficients of this
matrix to ``\delta^2`` (terms ``x\delta``, ``p_x\delta``, ``p_x\delta^2``, ``x\delta^2``),
plus the geometric degree-3 terms from ``H_4`` and the ``l`` row as for the drift.
Check printed by the package for ``L=0.21`` m, ``k_1=0.52`` m⁻²:
``(x|x\delta) = 0.0114222``, ``(x|p_x\delta) = -0.208398``,
``(x|p_x\delta^2) = 0.207598``, all equal to the derivatives of ``R_x(\delta)``.

## Sector bend (SBend body)

**Hamiltonian** (curvature ``h = \text{angle}/L``):

```math
H = e(\delta) - (1+hx)\,p_s + hx + \frac{h^2x^2}{2} + \frac{k_1}{2}(x^2-y^2).
```

Expanded (what is stored):

```math
H = \frac{h^2+k_1}{2}x^2 - \frac{k_1}{2}y^2 - h\,x\,\delta + H_{2,\text{drift}}
 + \frac{h\,x}{2}(p_x^2+p_y^2) + H_{3,\text{drift}}
 - \frac{h\,x}{2}(p_x^2+p_y^2)\delta + H_{4,\text{drift}}.
```

``\partial H/\partial x = 0`` at ``z=0`` (checked by the code; a warning is issued
otherwise), so the reference arc is a solution. The pole-face angles `e1, e2`,
`gap`, `fint` are **not** in the body Hamiltonian; edges are added only through a
[`FringeModel`](@ref) (see [Fringe fields](fringe.md)). Without a model the bend is
a pure sector body.

**Linear matrix** (``K = h^2 + k_1 > 0``, ``\omega = \sqrt K``, ``C=\cos\omega L``,
``S=\sin\omega L``; vertical plane: quadrupole with ``-k_1``, a drift if ``k_1=0``):

```math
R = \begin{pmatrix}
C & S/\omega & 0 & 0 & 0 & h(1-C)/K\\
-\omega S & C & 0 & 0 & 0 & hS/\omega\\
0 & 0 & \cdot & \cdot & 0 & 0\\ 0 & 0 & \cdot & \cdot & 0 & 0\\
-hS/\omega & -h(1-C)/K & 0 & 0 & 1 & \dfrac{L}{\gamma_0^2} - \dfrac{h^2(\omega L - S)}{\omega^3}\\
0 & 0 & 0 & 0 & 0 & 1\end{pmatrix}.
```

(All nine non-trivial entries were compared with `linear_matrix` for the
5.33 m, 22.5° bend: identical to 10 digits.) The sign of row 5 follows from
``dl/ds = \partial H/\partial\delta = v/v_0 - (1+hx)(1+\delta)/p_s``: an outward
particle travels a longer path and falls behind.

**Chromatic part at fixed δ**, exact in ``\delta`` to linear order in ``(x,p_x)``:

```math
x' = \frac{p_x}{1+\delta},\qquad p_x' = h\delta - (h^2+k_1)\,x .
```

The stored degree-3 map contains its expansion to ``\delta^2`` together with
the geometric terms ``h x p_x^2``, ``h x p_y^2`` etc.

## Sextupole

**Hamiltonian**: ``H = e(\delta) - p_s + \tfrac{k_2}{6}(x^3 - 3xy^2)``.

**Linear matrix**: identical to the drift (the sextupole term is cubic).

**Map** (degree 3): the leading nonlinear terms are

```math
p_{x,f} = p_x - \frac{k_2L}{2}(x^2-y^2) + O(L^2),\qquad
p_{y,f} = p_y + k_2L\,xy + O(L^2),
```

and the code keeps the full ``O(L^2), O(L^3)`` mixing with the drift (e.g.
``x\,p_x``, ``p_x^2`` terms). Example (``L=0.12``, ``k_2=-3.1``):
``(p_x|x^2) = 0.186 = -k_2L/2``, ``(p_y|xy) = -0.372 = k_2L``.
The generator starts ``f_3 = -\tfrac{k_2L}{6}(x^3-3xy^2) + \dots + f_{3,\text{drift}}``.

## Thin RF cavity (`RFCavity`)

No Hamiltonian; a single kick at the position of the element in the list
([`rf_kick`](@ref), numbers; [`rf_kick_taylor`](@ref), polynomial). For a particle
``(l,\delta)``:

```math
\begin{aligned}
\Delta t &= -\frac{l}{v(\delta)},\qquad \varphi = -h\,\omega_0\,\Delta t,\qquad \omega_0 = \frac{2\pi v_0}{C},\\
E &= \sqrt{p_0^2(1+\delta)^2 + m^2},\qquad E' = E + Z V \sin(\varphi_s + \varphi)\quad[\text{MeV}],\\
\delta' &= \frac{\sqrt{E'^2 - m^2}}{p_0} - 1,\qquad l' = -v(\delta')\,\Delta t,
\end{aligned}
```

``x, p_x, y, p_y`` unchanged. The formulas are exact (no expansion) in tracking.
Linearised: ``(\delta|l) = \dfrac{2\pi h\, ZV\cos\varphi_s}{C\,\beta_0^2E_0}`` (checked numerically).
For ``\sin\varphi_s \neq 0`` the kick has a constant part ``ZV\sin\varphi_s``, so a map
with RF has a non-zero fixed point (found by `fixed_point_6d`).

## Finite cavity (`Cavity`)

* In **maps and linear optics** it is a drift of length ``L`` (RF off):
  `hamiltonian(::Cavity) = hamiltonian(Drift)`.
* In **exact tracking** the s- and t-dependent Hamiltonian is integrated:

```math
H_\text{cav} = e(\delta) - p_s + A\cos\Phi,\qquad
\Phi = \frac{\omega(s - L/2)}{v_0} - \frac{\omega\,l}{v(\delta)} - \varphi_s,\qquad
A = \frac{ZV}{p_0c}\,\frac{c}{\omega L},\quad \omega = h\frac{2\pi v_0}{C}.
```

Details and the kick it produces: [Tracking](tracking.md).

## Marker

Identity map, no factor contribution (``R=I``, ``f_3=f_4=0``).

## Edges

Edges (EdgeIn/EdgeOut) are separate thin items produced by `fringe_sequence`; their
explicit maps are on the [Fringe fields](fringe.md) page.
