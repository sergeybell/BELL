# Equilibrium momentum shift

In a stationary RF bucket the phase returns to itself after a synchrotron
period, so the turn-averaged relative period lag vanishes,

```math
\langle\tau\rangle=0.
```

With the expansion ``\tau=\eta_0\delta+\eta_1\delta^2+\Lambda`` (``\Lambda`` =
betatron path lengthening, constant on the synchrotron time scale) this fixes
the mean momentum deviation

```math
\delta_{\mathrm{eq}}
=-\frac{1}{\eta_0}\Bigl(\frac{\eta_1\delta_m^2}{2}+\Lambda\Bigr),
\qquad
\delta_m=\frac{\delta_{\max}-\delta_{\min}}{2},
\qquad
\Lambda=-\frac{2\pi}{C}\bigl(J_x\xi_x+J_y\xi_y\bigr).
```

The slip factors ``\eta_{0,1}`` come from the circumference expansion
``C(\delta)/C_0=1+\alpha_0\delta+\alpha_1\delta^2`` and the exact velocity
ratio; ``\xi_{x,y}`` are the chromaticities; ``J_{x,y}`` are Courant–Snyder
actions. A cubic truncation of the RF waveform,
``\sin\varphi\approx\varphi-\varphi^3/6``, multiplies the
``\eta_1\delta_m^2/2`` piece by ``1-\hat\varphi^2/32`` and shifts the
synchrotron frequency by ``\Omega/\omega_s=1-\hat\varphi^2/16``.

The mean energy and Lorentz factor receive both the shift ``\delta_{\mathrm{eq}}``
and a kinematic term from ``\langle\delta^2\rangle``:

```math
\Delta E_{\mathrm{eq}}
=\beta_s^2 E_s\Bigl(\delta_{\mathrm{eq}}+\frac{\delta_m^2}{4\gamma_s^2}\Bigr).
```

These relations are verified by symplectic tracking of the full ring
Hamiltonians (canonical variables) against the formulae above.