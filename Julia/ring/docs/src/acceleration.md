# Acceleration and adiabatic damping

*`src/ExactIntegrator.jl` (`renormalize`), `src/Acceleration.jl`
(`accelerate`, `turn_snapshots`), script `05_rf_acceleration.jl`.*

## Reference momentum update

All particles, including the synchronous particle (#1, initially ``z=0``),
are tracked through the full lattice with the exact integrator. Right after
the cavity the synchronous particle has ``\delta_s\neq0``; it is made the new
reference:

```math
p_0'=p_0(1+\delta_s),\quad 1+\delta'=\frac{1+\delta}{1+\delta_s},\quad
p_x'=\frac{p_x}{1+\delta_s},\ p_y'=\frac{p_y}{1+\delta_s},\quad
\Delta t'=\Delta t-\Delta t_s,\ l'=-v\,\Delta t' .
```

Positions ``x,y`` are unchanged. The magnets are assumed to ramp with ``p_0``
(normalised ``h,k_1,k_2`` fixed) and the RF frequency follows
``h\,\omega_0(p_0)`` (locked to the synchronous particle).

## Adiabatic damping

The physical transverse momentum ``P_x=p_0p_x`` is not changed by a
longitudinal field, so ``p_x\propto1/p_0``. Between kicks the motion is a
betatron oscillation with fixed optics; the Courant–Snyder action in
normalised variables ``J=\tfrac12(\gamma x^2+2\alpha xp_x+\beta p_x^2)`` is an
adiabatic invariant of the *unnormalised* phase space,
``\oint P_x\,dx=p_0\oint p_x\,dx=2\pi p_0J=\text{const}``. Hence
``J\propto1/p_0`` and the amplitude ``A=\sqrt{2\beta J}`` scales as

```math
\frac{A(n)}{A(0)}=\sqrt{\frac{p_{\rm inj}}{p_0(n)}} .
```

Ramp results (``10\to270`` MeV, ``V=0.2`` MV, ``h=5``, ``\phi_s=150^\circ``, then
1000 turns with ``\phi_s=180^\circ``): see the table in
[Longitudinal motion, separatrix and synchrotron tune](longitudinal.md).
