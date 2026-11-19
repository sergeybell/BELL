# Coordinates, units, reference particle

## Canonical state

The dynamical state at arc length ``s`` is the six-vector

```math
z=(x,\,p_x,\,y,\,p_y,\,l,\,\delta),
```

with indices 1 through 6. Source names for the momenta and the momentum
deviation are `px`, `py`, and `delta`.

* ``x,y``: transverse displacements from the reference orbit, in metres.
* ``p_x=P_x/p_0``, ``p_y=P_y/p_0``: canonical momenta divided by the reference momentum ``p_0``.
* ``\delta=(p-p_0)/p_0``: relative momentum deviation. This is not an energy offset.
* ``l=-v(\delta)\,\Delta t``: longitudinal coordinate conjugate to ``\delta``. Here ``\Delta t=t-t_0(s)`` is the time lag relative to the reference particle and ``v(\delta)`` is the particle speed. Positive ``l`` is an early arrival.

```math
\frac{v(\delta)}{v_0}
=\frac{1+\delta}{\beta_0\sqrt{(1+\delta)^2+1/(\beta_0\gamma_0)^2}}.
```

The symplectic form on this space is the standard one,

```math
\omega=\mathrm{d}x\wedge\mathrm{d}p_x+\mathrm{d}y\wedge\mathrm{d}p_y+\mathrm{d}l\wedge\mathrm{d}\delta,
```

or ``S=\mathrm{diag}(S_2,S_2,S_2)`` with
``S_2=\begin{pmatrix}0&1\\-1&0\end{pmatrix}``.
Hamilton's equations are ``\mathrm{d}z/\mathrm{d}s=S\nabla H``.
A map ``z\mapsto\mathcal{M}(z)`` is symplectic when its Jacobian satisfies
``J^{\mathsf T}SJ=S``. The scalar ``\max|J^{\mathsf T}SJ-S|`` is the symplecticity check used throughout the package.

## Units

Lengths are in metres, angles in radians, energies in MeV, and magnetic rigidity in T·m.

* Rest energy, kinetic energy, and ``p_0\cdot c`` are in MeV.
* The normalised gradient ``k_1=(\partial B_y/\partial x)/(B\rho)`` is in m``^{-2}``. The sextupole strength ``k_2=(\partial^2 B_y/\partial x^2)/(B\rho)`` is in m``^{-3}``.
* Curvature ``h=1/\rho`` is in m``^{-1}``.
* ``k_1>0`` focuses in ``x``.

## Reference particle

`RefParticle` is the design particle. Build it with `reference_from_kinetic` from a species and a kinetic energy. Do not fill the fields by hand. It stores

* the species (`name`, rest energy, charge number ``Z``, anomalous moment ``G``);
* the kinetic energy that was requested;
* ``\gamma_0``, ``\beta_0``, ``\beta_0\gamma_0``;
* ``p_0\cdot c`` in MeV;
* the magnetic rigidity ``B\rho`` in T·m.

With rest energy ``mc^2`` and kinetic energy ``E_{\mathrm{kin}}``,

```math
\gamma_0=1+\frac{E_{\mathrm{kin}}}{mc^2},\qquad
\beta_0\gamma_0=\sqrt{\gamma_0^2-1},\qquad
p_0 c=mc^2\,\beta_0\gamma_0,
```

and ``B\rho=(p_0 c)/(Z\,e)`` with ``Z`` the signed charge number. Predefined species are the deuteron, the proton, and the electron. The anomalous moment is stored for later spin work and is not used in the orbit.

Every normalised strength (``k_1``, ``k_2``, ``h``, ``p_x``, ``\delta``) refers to this ``p_0``.

## Laboratory coordinates

Plots and laboratory input use

```math
u=(x,\,x',\,y,\,y',\,\varphi,\,\delta),\qquad
x'=\frac{p_x}{p_s},\quad
p_s=\sqrt{(1+\delta)^2-p_x^2-p_y^2}.
```

The RF phase reported with the orbit is

```math
\varphi=\frac{2\pi h}{C}\,\frac{l}{v(\delta)/v_0}.
```

``\varphi>0`` means early arrival: `canonical_to_laboratory` maps positive ``l`` to positive ``\varphi``. The set ``u`` is not canonical. ``x'`` depends on ``\delta`` and on ``p_y``, and ``\varphi`` is scaled by ``1/v(\delta)``. Symplecticity is checked only in ``z``. The conversion ``z\leftrightarrow u`` is exact. A truncated polynomial map applies the same conversion as a Taylor polynomial cut at the map order.
