# Fringe fields

```@meta
CurrentModule = Bell
```

Module: `Fringe` (`src/fringe.jl`). The lattice stays physical (an [`SBend`](@ref)
carries `e1, e2, gap, fint`); the edge model is a separate object chosen when the
ring is built.

## Template: how an edge is inserted

(How to choose a model: User's Guide, [Fringe-field models](../guide/fringe.md).)


```
fringe_sequence(lat, model):                    # model chosen per element: override[name] > family > NoEdge
    for e in lat.elements:
        m = edge_model(model, e)
        NoEdge  → [e]
        else    → [EdgeMap(e.name*"_in", e, :in, m),  e,  EdgeMap(e.name*"_out", e, :out, m)]

build_factors(lat, model):
    EdgeMap → dragt_finn_factor(edge_map(item, ref))      # a separate thin LieFactor
    element → dragt_finn_factor(element_lie_map(e, ref))  # body unchanged (full length L)
ring_lie_factor(lat, model) = product_factors(all factors)
```

So a bend with an edge model becomes three Lie factors:
**EdgeIn → Body(L) → EdgeOut**. `FringeModel()` (all `NoEdge`) leaves the list
unchanged and reproduces the plain ring exactly. Every edge map is a polynomial
map of degree 3 (`MAP_ORDER`), cached by (parent, side, model, ``p_0``); the model key is `repr(model)` plus the table contents for a `Tabulated` profile.

## Hard edge (`HardEdge(order)`)

Every hard edge is a thin Lie map ``e^{-{:}K{:}}z`` (`lie_map(K, 1)`), hence
symplectic. Notation: ``h`` = curvature of the parent bend, ``\varepsilon`` = `e1`
(entrance) or `e2` (exit), ``g`` = full gap, FINT from the magnet or the override,

```math
\psi = \varepsilon - h\,g\,\mathrm{FINT}\,\frac{1+\sin^2\varepsilon}{\cos\varepsilon}.
```

### Dipole, order 1 (linear)

```math
K = -\frac{h\tan\varepsilon}{2}x^2 + \frac{h\tan\psi}{2}y^2
\quad\Longrightarrow\quad
p_x \to p_x + h\tan\varepsilon\;x,\qquad p_y \to p_y - h\tan\psi\;y .
```

Matrix: ``\begin{pmatrix}1&0\\ h\tan\varepsilon&1\end{pmatrix}`` in ``x``,
``\begin{pmatrix}1&0\\ -h\tan\psi&1\end{pmatrix}`` in ``y``, identity elsewhere.

### Dipole, order 2 and 3 (nonlinear)

The vertical kick uses the effective angle ``\psi \pm x'`` with the exact slope
``t = x' = p_x/p_s``:

```math
K = -\frac{h\tan\varepsilon}{2}x^2 + \frac{h\,y^2}{2}\,T,\qquad
T_\text{in} = \frac{\tan\psi + t}{1 - t\tan\psi} = \tan(\psi + x'),\quad
T_\text{out} = \frac{\tan\psi - t}{1 + t\tan\psi} = \tan(\psi - x').
```

Order 2 keeps ``K`` to degree 3:
``K_3 = -\tfrac{h\tan\varepsilon}{2}x^2 + \tfrac{h\tan\psi}{2}y^2 \pm \tfrac{h}{2\cos^2\psi}\,y^2p_x``
(``+`` entrance, ``-`` exit), which is the term ``h\,y^2p_x/2`` from ``B_s = y\,b'(s)``
of a step field. Order 3 keeps ``K`` to degree 4 (adds ``y^2p_x\delta``, ``y^2p_x^2``
terms of ``T``).

### Quadrupole (order ≥ 2; order 1 = identity)

```math
K = \pm\frac{k_1}{12(1+\delta)}\big(x^3p_x + 3xy^2p_x - 3x^2y\,p_y - y^3p_y\big)\qquad(+\ \text{in},\ -\ \text{out}).
```

``K`` is already of degree 4, so in the degree-4 polynomial space the
``1/(1+\delta)`` factor contributes only its constant 1; the map receives cubic
terms ``\propto k_1 x^3, k_1xy^2, \dots``.

## Soft edge (`SoftEdge(profile, method)`)

### Profiles ``F(z)`` (0 outside, 1 inside, ``z`` along the beam)

| profile | raw formula | half-width ``W_0`` |
|:--|:--|:--|
| `Enge(c; gap=g)` | ``F = \dfrac{1}{1+\exp P(-z/g)},\ P(t) = \sum_{k\ge1} c_k t^{k-1}`` (default `ENGE_DEFAULT`) | first ``z`` (steps ``g/20``, ≤ 40g) where ``F < 10^{-10}`` and ``1-F < 10^{-10}`` |
| `Tanh(λ)` | ``F = \tfrac12\big(1 + \tanh(z/\lambda)\big)`` | ``12\lambda`` |
| `LinearRamp(w)` | ``F = \operatorname{clamp}(z/w + \tfrac12,\ 0,\ 1)`` | ``w/2`` |
| `Tabulated(z_i, B_i)` | ``B/\max\lvert B\rvert``, linear interpolation, mirrored if falling | ``\max\lvert z_i\rvert`` |

[`profile_function`](@ref) then shifts the profile so that ``z=0`` is the effective
field boundary (equal area): it iterates ``c \leftarrow c - G(c)`` with

```math
G(c) = \int_{-L}^{0}F(z+c)\,dz + \int_{0}^{L}\big(F(z+c) - 1\big)\,dz \quad(\text{Simpson, } 4000\text{ intervals}),
```

until ``|G| < 10^{-13}``; the region is ``[-W, W]``, ``W = W_0 + |c|``. At the exit the
profile is reversed, ``F_\text{out}(s) = F(-s)``. [`fringe_integral`](@ref)
``I = \int F(1-F)\,dz`` (``= g\cdot\mathrm{FINT}`` for a dipole).

### Soft-edge Hamiltonian ([`fringe_hamiltonian`](@ref))

In the Frenet frame that follows the field, ``h(s) = h_0F(s)``:

```math
H(s) = e(\delta) - (1+hx)\sqrt{(1+\delta)^2 - (p_x - a_x)^2 - p_y^2} + hx + \frac{h^2x^2}{2}
 + \frac{b_1}{2}(x^2-y^2) + \frac{b_2}{6}(x^3 - 3xy^2) + h\,b_1\Big(\frac{x^3}{3} - \frac{xy^2}{2}\Big),
```

```math
a_x = -\frac{y^2\,h_0F'(s)}{2(1+hx)},\qquad b_1 = h_0F'(s)\,\theta(s),\qquad b_2 = h_0F''(s)\,\theta(s)^2,
```

```math
\theta_\text{in}(s) = h_0\int_{-W}^{s}F,\qquad \theta_\text{out}(s) = -h_0\int_{s}^{W}F
```

(``\theta`` = angle of the orbit to the straight field-free orbit; ``b_1, b_2``
describe a straight EFB seen from the curved frame). ``F'``, ``F''`` are window
averages (central differences over ``\pm w``), so kinked profiles stay integrable.
``F = 1`` gives the body Hamiltonian, ``F = 0`` a drift.

### Methods for the s-dependent ``H`` over ``[-W, W]``

**`Splitting(n, order)`** — ``\Delta s = 2W/n``; for each slice ``[a, a+\Delta s]``:

```
order 2: sub-slices [a, a+Δs]
order 4: sub-slices of lengths W1Δs, W0Δs, W1Δs (W1 = 1/(2−2^{1/3}), W0 = 1−2W1 < 0)
for each sub-slice [a, b]:
    H_k = fringe_hamiltonian with F((a+b)/2), F' = (F(b)−F(a))/(b−a), F'' = (F'(b)−F'(a))/(b−a), θ((a+b)/2)
    M  ← exp(−(b−a) :H_k:) ∘ M        (lie_map, then compose_maps, degree 3)
```

**`Magnus(n, order)`** — on each interval with midpoint ``s_m``:

```math
\text{order 2: } \Omega = \Delta s\,H(s_m);\qquad
\text{order 4: } \Omega = \frac{\Delta s}{2}(H_1 + H_2) - \frac{\sqrt3}{12}\Delta s^2\,[H_1, H_2],\quad
H_{1,2} = H\big(s_m \mp \tfrac{\sqrt3}{6}\Delta s\big),
```

and the interval map is ``e^{-{:}\Omega{:}}`` (`lie_map(Ω, 1)`), composed in order.

**`TaylorIntegrate(tol, nmax)`** — classical RK4 of ``dz/ds = S\nabla H(s)`` in
`TaylorN` (degree 3), starting with 16 steps, doubling the number of steps until the
maximum coefficient change ``< \text{tol}`` or ``n > n_\text{max}``. Not symplectic;
it is the reference for the two methods above.

### Pole-face angle and the wrapper

If ``\varepsilon\ne0`` a linear wedge ``W_\varepsilon = e^{-{:}K_w{:}}``,
``K_w = -h_0\tan\varepsilon\,(x^2-y^2)/2`` (``p_x \mathrel{+}= h_0\tan\varepsilon\,x``,
``p_y \mathrel{-}= h_0\tan\varepsilon\,y``) is applied at the EFB. The soft region
replaces "last ``W`` of drift + first ``W`` of body" of the hard-edge layout, so the
lattice length is unchanged. In order of application to a particle:

```
EdgeIn  :  Drift(−W)  →  [wedge e1]  →  Soft(−W…W)  →  Body(−W)
Body    :  full body map of length L
EdgeOut :  Body(−W)   →  Soft_out(−W…W)  →  [wedge e2]  →  Drift(−W)
```

(``\text{Body}(-W) = e^{+W{:}H_\text{body}{:}}``, ``\text{Drift}(-W)`` likewise; the
magnet must satisfy ``2W < L``.) Implemented for `SBend` only.

## Checks ([`fringe_observables`](@ref))

Tunes, ``\xi_1``, ``\xi_2`` of the ring with the model (via `chromaticity_orders` of
`factor_to_map(ring_lie_factor(lat, model))`), symplecticity error of every edge
map (linear part and Jacobian at a test point) and of the ring ``R``, and the
computing time.


!!! note "Metric term `h\,b_1`"
    In the curved frame the force carries the factor `(1+hx)`, so the gradient term contributes `h b_1 x^2` (horizontal) and `h b_1 xy` (vertical). Since `\int h b_1\,ds = -\tfrac12\int b_2\,ds` (`\theta' = h_0F`), the horizontal sextupole-like term of the straight EFB cancels and the vertical one is `b_2/2`. Without this term (before Oct 2026) SoftEdge differed from COSY FR 3 by `\Delta\xi_x \approx -2.9\cdot10^{-3}`, `\Delta\xi_y \approx +1.0\cdot10^{-3}` at a 10 cm gap; with it the single-dipole second-order map agrees with COSY to `10^{-9}`.
