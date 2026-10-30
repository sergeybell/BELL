# Tracking and integrators

```@meta
CurrentModule = Bell
```

Modules: `Tracking` (`src/tracking.jl`), `ExactIntegrator`
(`src/exact_integrator.jl`), `Acceleration` (`src/acceleration.jl`).

Bell has four ways to move phase-space coordinates through elements. All of them
are generic code: they run on `Float64` points and on `TaylorN` polynomials (the
latter gives exact Jacobians for symplecticity checks).

| route | function | what is integrated | truncation | symplectic | used for |
|:--|:--|:--|:--|:--|:--|
| **exact integrator** | [`exact_track_element`](@ref), [`exact_track_turn`](@ref) | exact (non-expanded) ``H`` incl. finite RF cavity | none (only step error, 4th order) | yes | reference tracking, acceleration, longitudinal studies |
| Lie-factor tracking | [`track`](@ref), [`apply_factor`](@ref) | element factors ``R, f_3, f_4`` | generators to degree 4 (``H`` to degree 4) | yes | fast tracking of the magnetic ring, thin RF |
| Lie series | [`lie_map`](@ref) | ``e^{-L{:}H{:}}z`` | map cut at degree 3 | to degree 3 | element maps |
| RK4 in TaylorN | [`taylor_map_element`](@ref) | ``dz/ds=S\nabla H`` | degree 3 | no | cross-check maps |

## 1. Exact-Hamiltonian integrator (`ExactIntegrator`)

### Equations of motion

For every magnet (drift, quadrupole, sextupole, sector bend are special cases
with ``(h, k_1, k_2)`` taken from the element type, `_params`) the code evaluates
the exact vector field `_field_mag`:

```math
H_\text{mag} = e(\delta) - (1+hx)p_s + hx + \tfrac{h^2x^2}{2} + \tfrac{k_1}{2}(x^2-y^2) + \tfrac{k_2}{6}(x^3-3xy^2),
```

```math
\begin{aligned}
x' &= \frac{(1+hx)\,p_x}{p_s}, & p_x' &= h\,p_s - h - h^2x - k_1x - \tfrac{k_2}{2}(x^2-y^2),\\
y' &= \frac{(1+hx)\,p_y}{p_s}, & p_y' &= k_1 y + k_2\,xy,\\
l' &= \frac{v(\delta)}{v_0} - \frac{(1+hx)(1+\delta)}{p_s}, & \delta' &= 0 .
\end{aligned}
```

No expansion in ``p_x, p_y, \delta`` or ``x`` is made anywhere: ``p_s`` and
``v(\delta)/v_0`` are evaluated as square roots.

For the finite cavity (`_field_cav`), with ``\Phi = \omega(s-L/2)/v_0 - \omega l/v(\delta) - \varphi_s``:

```math
\begin{aligned}
x' &= p_x/p_s, & p_x' &= 0,\qquad y' = p_y/p_s,\qquad p_y' = 0,\\
l' &= \frac{v}{v_0} - \frac{1+\delta}{p_s} - A\sin\Phi\;\frac{\partial\Phi}{\partial\delta}, &
\delta' &= -\frac{\partial H}{\partial l} = -A\,\frac{\omega}{v}\,\sin\Phi
 = -\frac{ZV}{p_0c\,L}\,\frac{c}{v}\,\sin\Phi ,
\end{aligned}
```

with ``\partial\Phi/\partial\delta = \omega\,l\,v_0\,(dv/d\delta)/v^2``. The last
equation is exactly ``dp/ds = (1/v)\,dE/ds`` with the on-axis field
``dE/ds = -(ZV/L)\sin\Phi``. Integrated over the cavity it gives the energy gain
``\Delta E = ZVT\sin(\varphi_s + \varphi)`` with the transit factor
``T = \sin(\theta/2)/(\theta/2)``, ``\theta = \omega L/v_0`` ([`transit_factor`](@ref)).
Transverse RF fields are not modelled.

### Step algorithm

```
exact_track_element(z, e, ref, C; ds):
    Marker                → return z
    Cavity                → (ω, A, v0) = cavity_parameters(e, ref, C)        # recomputed from current ref
                            n = max(4, ceil(L / 0.02 m));  F = _field_cav
    Drift                 → n = 1;               F = _field_mag(h=0,k1=0,k2=0)
    Quadrupole/Sext/SBend → n = ceil(L / ds);    F = _field_mag(h, k1, k2)
    return _integrate(F, z, L, n)

_integrate(F, z, L, n):                       # Yoshida 4th order
    hs = L/n ; s = 0
    W1 = 1/(2 − 2^(1/3)) ; W0 = 1 − 2 W1     # W0 < 0
    repeat n times:
        z = _midpoint(F, z, s, W1·hs) ; s += W1·hs
        z = _midpoint(F, z, s, W0·hs) ; s += W0·hs
        z = _midpoint(F, z, s, W1·hs) ; s += W1·hs
    return z

_midpoint(F, z, s, h):                        # implicit midpoint, order 2, symplectic
    z' = z + h F(z, s + h/2)                  # explicit Euler guess
    repeat ≤ 40:  z' = z + h F((z + z')/2, s + h/2)   until |Δ| ≤ 1e-16
    return z'
```

[`exact_track_turn`](@ref) simply calls `exact_track_element` for
`lat.elements` in order, with a **fixed** reference particle.

## 2. One turn with acceleration ([`accelerate`](@ref))

Energy is gained only with a finite [`Cavity`](@ref) and the exact integrator.
The algorithm, turn by turn:

```
input: lat (contains a Cavity), ref0, particles u_k in user coords (u_1 = synchronous)
C = circumference(lat);  h = lat.harmonic;  ref = ref0
z_k = user_to_canonical(u_k, ref, C, h)
for n = 1 … max_turns:
    for e in lat.elements (in order):
        z_k = exact_track_element(z_k, e, ref, C; ds)          # all particles, current ref
        if e is Cavity:
            (z_k, ref) = renormalize(z_k, ref, 1)              # particle 1 becomes the reference
    record (ref.Ekin, γ0, β0, p0, user coords of all particles)
    if ref.Ekin ≥ Ekin_target (first time):  set φs of every Cavity to phi_flat (default π)
    stop phi_flat turns after that
```

**What is recomputed after each cavity passage** ([`renormalize`](@ref)). Let the
synchronous particle leave the cavity with ``(x_s, p_{xs}, y_s, p_{ys}, l_s, \delta_s)``:

```math
\begin{aligned}
p_0' &= p_0(1+\delta_s),\qquad E_\text{kin}' = \sqrt{p_0'^2 + m^2} - m,\qquad
\texttt{ref}' = \texttt{reference\_from\_kinetic}(E_\text{kin}')\ \Rightarrow\ \gamma_0', \beta_0', B\rho',\\
\delta' &= \frac{1+\delta}{1+\delta_s} - 1,\qquad
p_x' = \frac{p_x - p_{xs}}{1+\delta_s},\qquad p_y' = \frac{p_y - p_{ys}}{1+\delta_s},\qquad
x' = x - x_s,\quad y' = y - y_s,\\
\Delta t &= -\frac{l}{v(\delta)} + \frac{l_s}{v(\delta_s)},\qquad l' = -v(\delta)\,\Delta t .
\end{aligned}
```

The scaling ``p_x/(1+\delta_s)`` is the adiabatic damping. Because the next element
is tracked with ``\texttt{ref}'``, everything that depends on the reference is
rebuilt automatically on the next call: ``\beta_0, \gamma_0`` in ``e(\delta)`` and
``v(\delta)``, the revolution frequency ``\omega_0 = 2\pi\beta_0c/C`` and the RF
frequency ``\omega = h\omega_0`` (the RF follows the beam), and the cavity
amplitude ``A \propto ZV/p_0``. The normalised magnet strengths ``h, k_1, k_2`` are
kept: magnets are assumed to ramp in proportion to ``p_0``.

[`turn_snapshots`](@ref) runs one such turn and records the coordinates after
every element.

## 3. Lie-factor tracking (`Tracking`)

[`build_factors`](@ref) turns each element into one [`LieFactor`](@ref)
``(R, f_3, f_4)`` (see [Elements](elements.md)); a thin `RFCavity` stays in the
list as the struct. [`track`](@ref):

```
track(factors, z; turns, ref, C):
    for turn = 1 … turns, for F in factors (lattice order):
        if F is LieFactor:  z = apply_factor(F, z)
        else (RFCavity):    z = rf_kick(z, F, ref, C)     # exact thin kick, fixed ref
apply_factor(F, z):
    w = z
    if f4 ≠ 0:  w = _flow(field4, w)      # time-1 flow of −S∇f4  = e^{:f4:}
    w = _flow(field, w)                   # time-1 flow of −S∇f3  = e^{:f3:}
    return R · w
_flow(field, w; nsteps = 2, iters = 12):  # implicit midpoint in "time" t ∈ [0,1]
    repeat nsteps: w ← w + ½ field((w + w⁺)/2)  (fixed point iteration, ≤ 12, tol 1e-17)
```

This route is exactly symplectic (each flow is a Hamiltonian flow solved by a
symplectic rule) but the *generators* are truncated: ``f_3, f_4`` represent the
element exactly through degree 3 of the map. The reference energy is fixed;
this route does **not** accelerate (the thin kick changes ``\delta`` of a particle,
but ``p_0`` is never updated).

### Thin RF kick

[`rf_kick`](@ref) applies, at the position of the `RFCavity` in the element list:

```math
\Delta t = -\frac{l}{v(\delta)},\ \ \varphi = \frac{2\pi h v_0}{C}\,\frac{l}{v(\delta)},\ \
E' = \sqrt{p_0^2(1+\delta)^2+m^2} + ZV\sin(\varphi_s+\varphi),\ \
\delta' = \frac{\sqrt{E'^2-m^2}}{p_0}-1,\ \ l' = -v(\delta')\,\Delta t .
```

An energy kick at fixed arrival time is canonical in ``(t,-E)`` and the change of
variables to ``(l,\delta)`` is canonical, so the kick is symplectic.

## 4. Order and accuracy summary

| quantity | exact integrator | Lie-factor tracking |
|:--|:--|:--|
| ``p_s``, ``v(\delta)``, ``e(\delta)`` | exact square roots | expanded to degree 4 in ``H`` |
| longitudinal ``l, \delta`` | exact (no cut in ``\delta``) | same truncation as transverse |
| RF | finite cavity, exact ``\sin\Phi``, t-dependent | thin kick, exact formulas |
| reference energy | updated after every cavity (`renormalize`) | fixed |
| step error | ``O(\Delta s^4)`` | none per factor (flows of 2 midpoint steps) |

## Reference

```@docs
build_factors
ring_lie_factor
element_lie_map
apply_factor
Tracking.apply_nonlinear
Tracking._flow
track
rf_kick
exact_track_element
exact_track_turn
cavity_parameters
transit_factor
renormalize
ExactIntegrator.W1
ExactIntegrator._field_mag
ExactIntegrator._field_cav
ExactIntegrator._midpoint
ExactIntegrator._integrate
```
