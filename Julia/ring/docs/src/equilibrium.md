# Equilibrium momentum shift

A particle that oscillates in a stationary RF bucket (``\varphi_s=\pi``,
below transition) does not have mean relative momentum zero. The RF phase
returns to itself after a synchrotron period only if the mean time of flight
matches the RF period. With a nonlinear slip factor that condition fixes a
nonzero mean

```math
\delta_{\mathrm{eq}}
=-\frac{\eta_1\langle\delta^2\rangle+\Lambda}{\eta_0},\qquad
\Lambda=-\frac{2\pi}{C}\bigl(J_x\xi_x+J_y\xi_y\bigr).
```

Here ``\eta_0,\eta_1`` are the time-of-flight slip coefficients,
``\xi_{x,y}`` the linear chromaticities, and ``J_{x,y}`` the Courant–Snyder
actions. For a linear synchrotron oscillation ``\langle\delta^2\rangle=\delta_m^2/2``.
A cubic truncation of the RF waveform multiplies the ``\eta_1`` piece by
``1-\hat\varphi^2/32`` and shifts the synchrotron frequency by
``\Omega/\omega_s=1-\hat\varphi^2/16``. That correction is second order in
the actions and is not kept in the normal form below.

The RF is a thin energy kick at one point of the ring. A finite-length
cavity has its own closed-orbit offset in ``\delta``; that is a separate
study.

## From the torus

The one-turn map is the Lie product of the magnetic factors and the thin
kick. After the linear symplectic normalisation and the homological
removal of non-resonant generators, motion is uniform rotation of three
phases at three fixed actions. The torus of those actions is the set of
all triples of phases. Laboratory ``\delta`` is a function on that torus.
``\delta_{\mathrm{eq}}`` is the average of ``\delta`` over the whole torus:

```math
\delta_{\mathrm{eq}}=\delta^\ast+d_x J_x+d_y J_y+d_s J_s.
```

For thin RF at ``\varphi_s=\pi``, ``\delta^\ast=0``. The coefficients ``d_k``
are the quadratic part of the sixth component of the normalising
transformation, averaged over angles. They do not use ``\eta`` or
``\xi``. Substituting ``\delta_m^2=2\sigma_\delta^2 J_s`` into the formula
above reproduces them to relative ``10^{-6}``:

|  | ``d_x`` [1/m] | ``d_y`` [1/m] | ``d_s`` [1/m] |
|---|---|---|---|
| torus average | 0.14665244 | 0.22328941 | ``2.477064\times 10^{-3}`` |
| formula | 0.14665257 | 0.22328971 | ``2.477064\times 10^{-3}`` |

Element-by-element symplectic tracking of the unexpanded Hamiltonians,
with the exact thin kick, gives the same ``\delta_{\mathrm{eq}}`` to a
few parts in ``10^3`` at ``\delta_m\le 10^{-3}``. The residual grows as
``\delta_m^2`` (the next order, not kept). A long time average equals the
torus average when the tunes are incommensurate. The torus itself is not
built from the tracking samples.

The mean energy also receives a kinematic piece from ``\langle\delta^2\rangle``:

```math
\Delta E_{\mathrm{eq}}
=\beta_s^2 E_s\Bigl(\delta_{\mathrm{eq}}+\frac{\delta_m^2}{4\gamma_s^2}\Bigr).
```
