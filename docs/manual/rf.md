# Finite-length RF cavity

The cavity of length ``L`` carries a purely longitudinal on-axis field.
With ``q A_s/p_0=-A\cos(\omega t+\psi)`` and ``t=s/v_0-l/v(\delta)``,

```math
H_{\mathrm{cav}}=e(\delta)-p_s+A\cos\Phi,\qquad
\Phi=\frac{\omega(s-L/2)}{v_0}-\frac{\omega\,l}{v(\delta)}-\phi_s,
```

```math
A=\frac{ZV}{p_0 c}\,\frac{c}{\omega L},\qquad \omega=h\omega_0=\frac{2\pi h v_0}{C}.
```

Hamilton's equation ``\mathrm{d}\delta/\mathrm{d}s=-\partial H/\partial l``
is equivalent to ``\mathrm{d}E/\mathrm{d}s=-(ZV/L)\sin\Phi``. Integrated over
the cavity the energy gain is

```math
\Delta E=ZVT\sin(\phi_s+\phi),\qquad
T=\frac{\sin(\theta/2)}{\theta/2},\quad \theta=\frac{2\pi h L}{C}.
```

``H_{\mathrm{cav}}`` depends on ``s``; it is integrated by a symplectic
``s``-dependent scheme (implicit midpoint, Yoshida 4th order). Transverse RF
fields are omitted. The thin-kick limit ``L\to 0`` recovers an instantaneous
energy change ``ZV\sin(\phi_s+\phi)`` at fixed ``(x,p_x,y,p_y,l)``.