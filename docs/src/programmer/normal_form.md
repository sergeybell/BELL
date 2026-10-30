# Normal form

```@meta
CurrentModule = Bell
```

Module: `NormalForm` (`src/normal_form.jl`). Two routes are implemented. Route 1
is the primary one (example 18).

## Route 1: Dragt–Finn 6D normal form ([`dragt_finn_nf6d`](@ref))

### Input: the one-turn factor with RF ([`one_turn_factor_rf`](@ref))

```
Fmag = ring_lie_factor(lat_mag)                 # BCH product of element factors (R, f3, f4)
Mrf  = rf_kick_taylor(cav, ref, C)              # thin kick of rf_kick as a degree-3 polynomial
Frf  = dragt_finn_factor(Mrf)                   # (R_rf, f3_rf, f4_rf)
F    = concatenate(Fmag, Frf)                   # magnetic ring first, then RF
```

`dragt_finn_factor` requires the origin to be a fixed point, so this route is
for ``\varphi_s = 0`` or ``\pi`` (no constant term in the kick). With a polynomial map
as input, `dragt_finn_nf6d(M)` first finds the fixed point ``z^*`` by Newton
([`fixed_point_6d`](@ref)), centres the map ([`center_map_at`](@ref)) and factorises it.

The one-turn map is ``\mathcal M = R\circ e^{{:}f_3{:}}\circ e^{{:}f_4{:}}`` with
``R`` 6×6 symplectic, ``f_3`` cubic, ``f_4`` quartic (map exact to degree 3).

### Steps (`_nf6d_from_factor`)

**Step 1 — linear normalisation** ([`symplectic_eigen_A`](@ref)). Eigen-decompose
``R``; for each eigenvalue ``\lambda`` with ``\operatorname{Im}\lambda>0`` and
eigenvector ``v = a + ib``, scale so that ``a^\mathsf{T}Sb = 1`` (flip ``b`` and
``\mu`` if negative). Plane assignment: ``s`` = pair with the largest
``\delta``-weight, ``x`` = the remaining pair with larger weight in ``(x,p_x)``.
Columns of ``A`` are ``(a, b)`` per plane, ``\mu_k = \arg\lambda_k``. Result:

```math
A^{-1}RA = r(\mu) = \operatorname{diag}\Big(\begin{pmatrix}\cos\mu_k&\sin\mu_k\\-\sin\mu_k&\cos\mu_k\end{pmatrix}\Big),\qquad A^\mathsf{T}SA = S .
```

Coupling and dispersion are included automatically (the full 6×6 ``R`` with RF).

**Step 2 — generators in normal coordinates** ``w``: ``f_{3n}(w) = f_3(Aw)``,
``f_{4n}(w) = f_4(Aw)`` ([`pull_generator`](@ref)).

**Step 3 — order 3.** Solve the homological equation for ``F_3``:

```math
F_3 - F_3\circ r(\mu) = -f_{3n}.
```

The operator ``\mathcal L_k F = F - F\circ r`` is assembled as a matrix on the
monomial basis of degree ``k`` (56 monomials for ``k=3``, 126 for ``k=4``;
[`homological_operator`](@ref)); in the complex basis its eigenvalues are
``1 - e^{-i(\mathbf j-\mathbf n)\cdot\mu}``. It is solved by least squares
(`L3 \ (−f3)`). All cubic terms are removed (no cubic resonance is kept; for
generic tunes none is resonant).

**Step 4 — order 4.** The cubic transformation feeds the quartic order:

```math
f_4' = f_{4n} + \tfrac12[f_{3n}, F_3],\qquad
F_4 - F_4\circ r(\mu) + h_4 = -f_4',\qquad h_4 = \sum_{i\le j} c_{ij}\,J_iJ_j,\quad J_k = \tfrac12(w_{2k-1}^2+w_{2k}^2).
```

``F_4`` and the six coefficients ``c_{ij}`` are solved together, as one linear
system ``[\mathcal L_4\ \ -B](F_4, c) = -f_4'`` with ``B`` the coefficient columns of
``J_iJ_j``, by pseudo-inverse. ``h_4`` is the kept (resonant, phase-independent)
part: the detuning.

**Step 5 — transformation and outputs.**

```math
T = e^{{:}F_3{:}}\circ e^{{:}F_4{:}},\qquad T^{-1} = e^{-{:}F_4{:}}\circ e^{-{:}F_3{:}},\qquad
z = z^* + A\,T(w).
```

* Amplitude-dependent tunes ([`nf6d_tunes`](@ref)):
  ``\mu_k(J) = \mu_k - \partial h_4/\partial J_k = \mu_k - (2c_{kk}J_k + \sum_{j\ne k}c_{kj}J_j)``.
* Actions of a lab point ([`nf6d_actions`](@ref)): ``u = A^{-1}(z-z^*)``,
  ``w = T^{-1}(u)``, ``J_k = (w_{2k-1}^2 + w_{2k}^2)/2``.
* Equilibrium momentum ([`nf6d_delta_eq`](@ref)): with ``U_6 = (A\,T(w))_6``,

```math
\delta_\text{eq}(J) = z^*_6 + \langle U_6\rangle_\theta = z^*_6 + \sum_k d_kJ_k,\qquad
d_k = \tfrac12\Big(\frac{\partial^2U_6}{\partial w_{2k-1}^2} + \frac{\partial^2U_6}{\partial w_{2k}^2}\Big)\Big|_0,
```

  checked by a numerical average over an ``8^3`` torus grid (`d_grid`).
* Quality numbers: ``|A^{-1}RA - r|``, ``|A^\mathsf{T}SA - S|``, ``\max|f_3|`` before
  and after (``f_{3n} + F_3 - F_3\circ r``), ``\max|f_4' + F_4 - F_4\circ r - h_4|``.

## Route 2: map-based RF normal form ([`normal_form_rf`](@ref), older)

```
cav = RFCavity(V·T, h, φs)
M   = one_turn_map_rf(lat, cav)    # RK4 Taylor map of the ring (ring_taylor_map), then rf_kick_taylor
z*  = fixed_point_6d(M); Mc = center_map_at(M, z*)
pe, ξ from the RF-off Taylor map; po from ring_matrix; Qs0 = synchrotron_tune(η0)
A   = linear_A_6d(po, pe, Qs0):  dispersion shear A[1:4,6] = D4, Courant–Snyder blocks
      [√β 0; −α/√β 1/√β], longitudinal block diag((C/2πh)√βs, 1/√βs), βs = Qs0/(h|η0|)
Mn  = A⁻¹ ∘ Mc ∘ A ;  F = dragt_finn_factor(Mn)
δ_eq on 7 tori (Jx, Jy, δm) from ⟨Δφ⟩ = 0 (12³ angle grid, secant iteration);
fit δ_eq = c0 + cx Jx + cy Jy + cs Js
```

Notes (truthful status of this route): the linear ``A`` is approximate (not the
symplectic eigen-normalisation), the cubic generator is not solved by a
homological equation (`resonant_cubic` returns ``h_3 = 0`` and stores
``F_3 = f_3``), and ``h_4`` is taken as ``f_4``. The fitted coefficients are compared
with the analytic counterparts ``-\eta_1/(\eta_0\beta_s)``, ``(2\pi/C)\xi_{x,y}/\eta_0``.
Use route 1 for the normal form itself.

## Reference

```@docs
dragt_finn_nf6d
NF6D
NormalForm._nf6d_from_factor
one_turn_factor_rf
symplectic_eigen_A
homological_operator
nf6d_delta_eq
nf6d_actions
nf6d_tunes
normal_form_rf
NormalFormRF
delta_eq_nf
tunes_nf
nonlinear_normal_form
delta_eq_phase_condition
NormalForm.linear_A_6d
one_turn_map_rf
rf_kick_taylor
fixed_point_6d
center_map_at
evaluate_map
linear_normalizing_matrix
rotation_residual
```
