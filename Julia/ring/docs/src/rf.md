# RF cavity of finite length

*`src/ExactIntegrator.jl` (`_field_cav`, `cavity_parameters`,
`transit_factor`), type `Cavity` in `src/Elements.jl`; checks in
`05_rf_acceleration.jl` section 2.*

The cavity is a **finite-length** element of length ``L`` with a purely
longitudinal on-axis field; it is *not* a drift plus a thin kick. With
``q A_s/p_0=-A\cos(\omega t+\psi)``, ``E_s=-\partial A_s/\partial t`` and
``t=s/v_0-l/v(\delta)``, the ``s``-Hamiltonian is

```math
H_{\rm cav}=e(\delta)-p_s+A\cos\Phi,\qquad
\Phi=\frac{\omega(s-L/2)}{v_0}-\frac{\omega\,l}{v(\delta)}-\phi_s,
```

```math
A=\frac{ZV}{p_0c}\,\frac{c}{\omega L},\qquad \omega=h\omega_0=\frac{2\pi h v_0}{C}.
```

Hamilton's equation gives

```math
\frac{d\delta}{ds}=-\frac{\partial H}{\partial l}=-\frac{ZV}{p_0cL}\frac{c}{v}\sin\Phi,
```

i.e. ``dE/ds=-(ZV/L)\sin\Phi``. Integrated over the cavity for a particle
with ``\phi`` the energy gain is

```math
\Delta E=ZVT\sin(\phi_s+\phi),\qquad T=\frac{\sin(\theta/2)}{\theta/2},\ \theta=\frac{\omega L}{v_0}=\frac{2\pi hL}{C}.
```

``H_{\rm cav}`` depends explicitly on ``s``; it is integrated with the implicit
midpoint rule composed to 4th order (Yoshida), which is symplectic for
``s``-dependent Hamiltonians. Transverse RF fields are neglected.

Checks at 270 MeV (``L=1`` m, ``V=0.2`` MV, ``h=5``, ``\phi_s=150^\circ``,
``T=0.99702``):

| check | result |
|---|---|
| synchronous gain vs ``ZVT\sin\phi_s`` | 0.099700856 vs 0.099701840 MeV |
| thin limit ``L=0.1/0.01/0.001`` m vs thin kick (``\phi=0.05``) | ``-6.3\times10^{-9}``, ``-1.5\times10^{-10}``, ``-9.9\times10^{-12}`` in ``\delta`` |
| symplectic error of the cavity Jacobian | ``4.4\times10^{-16}`` |

The thin kick (`RFCavity`, `Tracking.rf_kick`) applies the same energy gain ``ZV\sin(\phi_s+\phi)`` at one point; the finite cavity converges to it as ``L\to0`` (table above).

