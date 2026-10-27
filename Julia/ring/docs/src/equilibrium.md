# Equilibrium momentum shift

In a stationary RF bucket the phase returns to itself after a synchrotron
period only if the mean time of flight matches the RF period. With a
nonlinear slip factor this condition fixes a nonzero mean relative momentum

```math
\delta_{\mathrm{eq}}
= -\frac{\eta_1}{\eta_0}\,\frac{\delta_m^2}{2}
- \frac{\Lambda}{\eta_0},
\qquad
\Lambda = -\frac{2\pi}{C}\bigl(J_x\,\xi_x+J_y\,\xi_y\bigr),
```

where `delta_m` is the synchrotron half-amplitude, `eta_0` and `eta_1` are the
leading slip-factor coefficients, `xi_x,y` the linear chromaticities and
`J_x,y` the Courant-Snyder actions. A cubic truncation of the RF waveform
(`sin phi ~ phi - phi^3/6`) multiplies the `eta_1 delta_m^2/2` piece by
`1 - phi_hat^2/32` and shifts the synchrotron frequency by
`Omega/omega_s = 1 - phi_hat^2/16`.

The mean energy and Lorentz factor receive both the shift `delta_eq` and a
kinematic term from the mean square of delta:

```math
\Delta E_{\mathrm{eq}}
=\beta_s^2 E_s\Bigl(\delta_{\mathrm{eq}}+\frac{\delta_m^2}{4\gamma_s^2}\Bigr).
```

## Numeric check (magnetic ring, deuteron 270 MeV, stationary bucket)

Full-ring symplectic tracking in canonical variables, harmonic least-squares
extraction of `delta_eq` (constant term of a two-harmonic fit of delta over
many synchrotron periods; half-split systematic error at most about 1e-8).
Phase reported in the late-positive convention of the analytic note.

**Longitudinal only** (Lambda = 0), delta_m scan:

| delta_m | delta_eq (track) | -eta1 dm^2/(2 eta0) | rel. |
|---|---|---|---|
| 5.0e-4 | 3.45e-7 | 2.38e-7 | 0.45 |
| 1.0e-3 | 1.06e-6 | 9.50e-7 | 0.11 |
| 2.0e-3 | 3.88e-6 | 3.78e-6 | 0.025 |

**Betatron contribution.** At `delta_m ~ 3e-4` and vertical amplitude 20 mm,
`Lambda_fit / Lambda_form = 1.009`; with sextupoles set to `xi ~ 0` the fitted
Lambda drops by an order of magnitude toward the longitudinal baseline.

**Cubic phase.** For `phi_hat` in [0.30, 0.50], the measured `Omega/omega_s`
agrees with `1 - phi_hat^2/16` to 1e-4, and
`<(delta-deq)^2> / (dm^2/2)` stays within 1% of unity (cubic correction
`1 - phi_hat^2/32` is a few parts in 1e3).

These relations are verified by symplectic tracking of the full ring
Hamiltonians against the formulae above.