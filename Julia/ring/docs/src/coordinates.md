# Coordinates and reference particle

*Implemented in `src/Constants.jl`.*

The reference particle has charge ``Z e``, mass ``m``, momentum ``p_0``,
velocity ``v_0=\beta_0 c`` and Lorentz factor ``\gamma_0``. The independent
variable is the arc length ``s`` along the reference orbit.

**User coordinates** (input/output, plots):

```math
u = (x,\ x',\ y,\ y',\ \phi,\ \delta),\qquad x'=\frac{dx}{ds},\quad \delta=\frac{p-p_0}{p_0}.
```

``\delta`` is the relative *momentum* deviation (not energy). ``\phi`` is the RF
phase of arrival relative to the reference particle with the sign convention

```math
\phi = -h\,\omega_0\,\Delta t,\qquad \Delta t = t - t_0(s),
```

so ``\phi>0`` means the particle arrives *early*. ``h`` is the harmonic number and
``\omega_0 = 2\pi v_0/C`` the angular revolution frequency.

**Canonical coordinates** (all maps, Hamiltonians, symplecticity):

```math
z = (x,\ p_x,\ y,\ p_y,\ l,\ \delta),\qquad p_{x,y}=\frac{P_{x,y}}{p_0},\quad l = -v(\delta)\,\Delta t .
```

``(l,\delta)`` is a canonical pair: ``l`` is the path-length-like lag measured
with the particle's own velocity ``v(\delta)``, with

```math
\frac{v(\delta)}{v_0}=\frac{1+\delta}{\sqrt{1+\beta_0^2\left((1+\delta)^2-1\right)}} .
```

The energy term of every Hamiltonian is

```math
e(\delta)=\frac{1}{\beta_0}\sqrt{(1+\delta)^2+\frac{1}{\beta_0^2\gamma_0^2}}-\frac{1}{\beta_0^2},
```

(the constant is irrelevant for the dynamics and makes ``e(0)``-terms vanish after
removing the constant part). Particle data: deuteron, ``Z=1``,
``mc^2 = 2\times 938.272`` MeV (as in OptiM).
