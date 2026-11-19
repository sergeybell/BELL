# Normal form and the equilibrium momentum

## Normal form of the one-turn factor

The linear part ``R`` of the one-turn factor is put into three rotations by a symplectic matrix ``A``, ``A^{\mathsf T}SA=S``. The generators in the normal coordinates ``w=A^{-1}z`` are the pull-backs ``f(Aw)``.

A rotation acts on a complex monomial
``h_{\mathbf{j}\mathbf{n}}=\prod_k(q_k+ip_k)^{j_k}(q_k-ip_k)^{n_k}``
by the phase ``e^{-i(\mathbf{j}-\mathbf{n})\cdot\boldsymbol{\mu}}``. The homological operator on a generator is the divisor

```math
F-F\circ r(\mu),\qquad
F_{\mathbf{j}\mathbf{n}}
=\frac{-f_{\mathbf{j}\mathbf{n}}}{1-e^{-i(\mathbf{j}-\mathbf{n})\cdot\boldsymbol{\mu}}}.
```

At degree 3, ``\mathbf{j}=\mathbf{n}`` is impossible, so in the absence of a third-order resonance the whole cubic generator is removed. At degree 4 the kernel is the span of products of actions. Those terms stay as the normal-form Hamiltonian ``h_4=\sum_{i\le j}c_{ij}J_i J_j``. The cubic conjugation contributes ``\tfrac12[f_3,F_3]`` to the quartic generator before this step.

The homological operator is applied to the generators. Conjugation of a full polynomial map is not used to build the normal form. The same extraction of generators from a polynomial map that did not start life as a Lie factor is the Dragt–Finn step in [Lie factors and differential algebra](lie_and_da.md).

## The torus

After the normal form there are three amplitudes ``J_x,J_y,J_s`` and three phases. The torus of those actions is the set of all triples of phases at those fixed amplitudes: three circles, not a trajectory and not a list of turns.

One turn adds constant angles to the phases and does not change the amplitudes. Laboratory ``\delta`` is a function of the three phases. The equilibrium momentum is the average of that function over the whole torus at once. Pieces odd in the phases average to zero. What remains is linear in the actions,

```math
\delta_{\mathrm{eq}}=\delta^\ast+d_x J_x+d_y J_y+d_s J_s.
```

For thin RF on the stationary separatrix, ``\delta^\ast`` is a numerical zero. The cubic-RF correction in ``\hat\varphi^2`` is degree ``J^2`` and is not kept.

A long tracking average equals this torus average when the tunes are incommensurate. Tracking does not build the torus.

## Slip formula as a fit of that average

A particle that oscillates in a stationary RF bucket (``\varphi_s=\pi``, below transition) does not have mean relative momentum zero. The RF phase returns to itself after a synchrotron period only if the mean time of flight matches the RF period. With a nonlinear slip factor that condition fixes

```math
\delta_{\mathrm{eq}}
=-\frac{\eta_1\langle\delta^2\rangle+\Lambda}{\eta_0},\qquad
\Lambda=-\frac{2\pi}{C}\bigl(J_x\xi_x+J_y\xi_y\bigr).
```

Here ``\eta_0,\eta_1`` are the time-of-flight slip coefficients of [Readings of the one-turn map](map_readings.md), ``\xi_{x,y}`` the linear chromaticities, and ``J_{x,y}`` the Courant–Snyder actions. For a linear synchrotron oscillation ``\langle\delta^2\rangle=\delta_m^2/2``. A cubic truncation of the RF waveform multiplies the ``\eta_1`` piece by ``1-\hat\varphi^2/32`` and shifts the synchrotron frequency by ``\Omega/\omega_s=1-\hat\varphi^2/16``. That correction is second order in the actions and is not kept in the normal form.

The coefficients ``d_k`` from the torus do not use ``\eta`` or ``\xi``. Substituting ``\delta_m^2=2\sigma_\delta^2 J_s`` into the formula above reproduces them. On the Nuclotron-based lattice at deuteron kinetic energy 270 MeV the two readings agree to relative ``10^{-6}``:

|  | ``d_x`` [1/m] | ``d_y`` [1/m] | ``d_s`` [1/m] |
|---|---|---|---|
| torus average | 0.14665244 | 0.22328941 | ``2.477064\times 10^{-3}`` |
| slip formula | 0.14665257 | 0.22328971 | ``2.477064\times 10^{-3}`` |

Element-by-element integration of the unexpanded Hamiltonians, with the exact thin kick, gives the same ``\delta_{\mathrm{eq}}`` to a few parts in ``10^3`` at ``\delta_m\le 10^{-3}``. The residual grows as ``\delta_m^2`` (the next order, not kept). The formula is a fit of the torus average, not an independent dynamical model.

The mean energy also receives a kinematic piece from ``\langle\delta^2\rangle``:

```math
\Delta E_{\mathrm{eq}}
=\beta_s^2 E_s\Bigl(\delta_{\mathrm{eq}}+\frac{\delta_m^2}{4\gamma_s^2}\Bigr).
```

The RF in this comparison is a thin energy kick at one point of the ring. A finite-length cavity has its own closed-orbit offset in ``\delta``. That offset is the thick-cavity flow, not this average.
