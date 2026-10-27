# Assumptions and hypotheses

Every choice that is not a consequence of the lattice file or of the
equations is listed here.

1. **Injection energy 10 MeV** (5 MeV/u): not in `magnetic_2.opt`; chosen as a
   typical linac energy and confirmed by the user. Final 270 MeV (135 MeV/u)
   is from the file header.
2. **RF parameters** (chosen, not from a design): cavity length ``L=1`` m,
   voltage ``V=0.2`` MV, harmonic ``h=5``, ``\phi_s=150^\circ`` during the ramp and
   ``180^\circ`` on the flat top, one cavity in the first `OBM` drift. ``V`` is
   large so that the ramp takes ~2600 turns instead of ``10^5``-``10^6``.
3. **Cavity field** purely longitudinal on axis, sinusoidal in time, uniform
   over ``L``; **no transverse RF fields**, no beam loading.
4. **Magnets ramp with ``p_0``** (normalised strengths fixed); RF frequency
   locked to ``h\omega_0(p_0)`` of the synchronous particle.
5. **No fringe fields**, hard-edge elements, sector bends; **no space charge**,
   no impedances, no radiation, no errors/misalignments.
6. **Test amplitudes**: ``x_0=2`` mm, ``y_0=2`` mm, ``\delta_0=10^{-3}``
   (``x'=y'=0``); ensemble ``\phi\sim N(0,0.6)``, ``\delta\sim N(0,0.25\hat\delta)``,
   60 particles, seed 1.
7. **Optimisation targets A and B** are optional examples; the default
   lattice has sextupoles at zero.
8. **"Number of periods = 2"** in `magnetic_2.opt`: the listed sequence already
   bends by 360 deg, so it is used **once** (as in the COSY deck).
9. **``D(0)=0`` in OptiM** initial conditions vs our closed ``D(0)=0.0111`` m:
   OptiM's value is an input, not the closed solution; all other initial
   Twiss agree to ``10^{-5}``.
10. Mass ``2\times938.272`` MeV (as OptiM), not the deuteron mass 1875.613 MeV.
11. ``\alpha_1`` used for the separatrix is linearly interpolated between the
    maps at 10 and 270 MeV (it turned out equal at both).
12. Map orders: Hamiltonians to order 4, maps to order 3; tracking with the
    exact integrator uses ``ds=0.2`` m sub-steps (4th order).
