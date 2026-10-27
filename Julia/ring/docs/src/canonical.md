# Canonical phase space

## Variables

The dynamical state at arc length ``s`` is the six-vector

```math
z=(x,\,p_x,\,y,\,p_y,\,l,\,\delta).
```

- ``x,y``: transverse displacements from the reference orbit [m].
- ``p_x=P_x/p_0``, ``p_y=P_y/p_0``: canonical momenta normalised to the
  reference momentum ``p_0``.
- ``\delta=(p-p_0)/p_0``: relative *momentum* deviation.
- ``l=-v(\delta)\,\Delta t``: longitudinal coordinate conjugate to ``\delta``,
  with ``\Delta t=t-t_0(s)`` the time lag relative to the reference particle
  and ``v(\delta)`` the particle speed,

```math
\frac{v(\delta)}{v_0}=\frac{1+\delta}{\beta_0\sqrt{(1+\delta)^2+1/(\beta_0\gamma_0)^2}}.
```

The reference particle has charge ``Ze``, rest energy ``mc^2``, Lorentz factor
``\gamma_0``, velocity ``v_0=\beta_0 c`` and momentum ``p_0``.

The symplectic form on this space is the standard one,

```math
\omega=\mathrm{d}x\wedge\mathrm{d}p_x+\mathrm{d}y\wedge\mathrm{d}p_y+\mathrm{d}l\wedge\mathrm{d}\delta,
```

or in matrix form ``S=\mathrm{diag}(S_2,S_2,S_2)`` with
``S_2=\begin{pmatrix}0&1\\-1&0\end{pmatrix}``. A map ``z\mapsto\mathcal{M}(z)``
is symplectic iff its Jacobian satisfies ``J^\mathsf{T}SJ=S``.

## Output (non-canonical) coordinates

For plots and comparison with laboratory data one uses

```math
u=(x,\,x',\,y,\,y',\,\phi,\,\delta),\qquad
x'=\frac{p_x}{p_s},\quad
p_s=\sqrt{(1+\delta)^2-p_x^2-p_y^2},\quad
\phi=-h\omega_0\Delta t=\frac{2\pi h}{C}\,\frac{l}{v(\delta)/v_0}.
```

``\phi>0`` means early arrival. The set ``u`` is **not** canonical:
``x'`` depends on ``\delta`` and ``p_y``, and ``\phi`` is scaled by
``1/v(\delta)``. Symplecticity is therefore checked only in ``z``. The
conversion ``z\leftrightarrow u`` is exact (nonlinear); truncated polynomial
maps apply it as a full Taylor polynomial to the map order.

## Energy term

Every magnetic Hamiltonian uses

```math
e(\delta)=\frac{1}{\beta_0}\sqrt{(1+\delta)^2+\frac{1}{\beta_0^2\gamma_0^2}}-\frac{1}{\beta_0^2},
```

chosen so that the irrelevant constant vanishes at ``\delta=0`` after removal
of the constant term of ``H``.