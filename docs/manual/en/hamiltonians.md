# Element Hamiltonians

One formula per element, in the canonical ``z`` of [Coordinates](coordinates.md), with the ``\mathrm{energy\_term}`` of [The s-Hamiltonian](hamiltonian_s.md). Write

```math
p_s=\sqrt{(1+\delta)^2-p_x^2-p_y^2}
```

when the transverse vector potential vanishes. Normalised multipoles are ``k_1=B'/(B\rho)`` and ``k_2=B''/(B\rho)``. The body formulas below are hard-edge: the fringe is a separate map, [Fringe fields](fringe.md).

The same expression is the scalar passed to the flow and, on the Lie path, the polynomial of total degree 4 with the constant removed. That polynomial is a generator. It is not yet a map.

## Drift

No field, straight axis, ``a=0``:

```math
H=\mathrm{energy\_term}(\delta)-p_s.
```

A drift of length zero, and a marker, are the identity. The marker has ``H=0``.

## Quadrupole

The normalised gradient ``k_1`` enters as the scalar potential of a straight quadrupole,

```math
H=\mathrm{energy\_term}(\delta)-p_s+\frac{k_1}{2}(x^2-y^2).
```

Then ``\mathrm{d}p_x/\mathrm{d}s=-k_1 x``. A positive ``k_1`` focuses in ``x`` and defocuses in ``y``.

## Sextupole

```math
H=\mathrm{energy\_term}(\delta)-p_s+\frac{k_2}{6}(x^3-3xy^2).
```

The force is quadratic. First-order chromaticity and the second-order dispersion see this term only off axis, where a dispersive orbit gives ``x\sim D\delta``. Family placement is in [Sextupole families](sextupoles.md).

## Sector bend

Curvature ``h=1/\rho`` of the reference arc multiplies the longitudinal momentum and adds the weak-focusing potential that keeps ``x=0`` a solution:

```math
H=\mathrm{energy\_term}(\delta)-(1+hx)\,p_s+hx+\frac12 h^2 x^2+\frac{k_1}{2}(x^2-y^2).
```

On the reference orbit ``\partial H/\partial x=0``, so the design arc satisfies Hamilton's equations. The gradient in a bend uses the straight-magnet term. Curvature corrections of order ``h k_1 x^3`` are omitted. Edge angles are not inside this body Hamiltonian; they belong to the fringe map. A bend map built from the body alone therefore has no pole-face rotation.

## Thick and thin RF cavity

The cavity is the same kind of object as a magnet: one ``H`` when the gap has a length, and the exact flow of the integrated generator when the gap is thin.

### Mode and gauge

A pillbox on the ``\mathrm{TM}_{010}`` mode at ``\omega=h\omega_0``, with ``k=\omega/c`` and ``\omega_0=2\pi v_0/C``, has

```math
E_z=E_{\mathrm{m}} J_0(kr)\sin(\omega t+\varphi_0),\qquad
B_\theta=\frac{E_{\mathrm{m}}}{c} J_1(kr)\cos(\omega t+\varphi_0).
```

The gauge ``\Phi=0`` is a choice. Fields and trajectories do not depend on it; the Hamiltonian changes only by a canonical transformation. With no free charge in the gap, ``\mathbf{E}=-\partial\mathbf{A}/\partial t``, and

```math
A_z=\frac{E_{\mathrm{m}}}{\omega} J_0(kr)\cos(\omega t+\varphi_0)
```

reproduces both ``E_z`` and ``B_\theta``. Near the axis ``J_0(kr)\approx 1-k^2 r^2/4``. The transverse piece of order ``k^2 r^2`` (RF focusing) is dropped. The cavity treated here is longitudinally on axis.

### Transit time and the uniform gap

Substitute ``\omega t=\omega\Delta t+\omega s/v_0`` into the general ``H``. The piece ``ks/\beta_0`` is the particle sliding through a standing wave. Average over a gap of length ``L``. The effective voltage is ``V=E_{\mathrm{m}} L\,T_{\mathrm{tr}}`` with the transit-time factor

```math
T_{\mathrm{tr}}=\frac{\sin(\omega L/(2v_0))}{\omega L/(2v_0)}.
```

For ``L=\lambda/2`` and ``\beta_0=1``, ``T_{\mathrm{tr}}=2/\pi``. The voltage stored on the element is this effective voltage: the transit-time factor is already inside ``V``, not applied again.

Set ``\varphi_s`` so that a particle with ``q=0`` takes the energy ``ZeV\sin\varphi_s`` from a uniform gap, ``\mathrm{d}E/\mathrm{d}s=ZeV\sin(\varphi_s+kq)/L``. In the pair ``(q,P)``,

```math
H_{\mathrm{rf}}=\frac{ZeV}{p_0 c\,k L}\cos(\varphi_s+kq),
\qquad k=\frac{2\pi h\beta_0}{C}.
```

Early arrival is ``\Delta t<0``, hence ``q=-c\Delta t>0``. Passing to ``(l,\delta)`` with ``q=l/\beta(\delta)`` does not change ``H``. Together with the drift part, the thick-cavity Hamiltonian used by the flow is

```math
H=\mathrm{energy\_term}(\delta)-p_s
+\frac{ZV}{p_0 c\,k L}
\cos\!\Big(\varphi_s+\frac{k\,l}{\beta_0\,v(\delta)/v_0}\Big).
```

It does not depend on ``s`` inside the gap. The standing-wave phase was removed by the average that produced ``T_{\mathrm{tr}}``. What remains at finite ``L`` is that the drift and ``H_{\mathrm{rf}}`` act together: they do not commute, and the map is the flow of the sum, not a kick glued to a drift.

``\partial H/\partial\delta`` receives a piece from the phase, because ``l=-v\Delta t`` contains the particle's own speed. That piece is kept. It is written out in [RF motion](rf.md).

### Thin limit

Let ``G=L H_{\mathrm{rf}}``, so

```math
G=\frac{ZV}{p_0 c\,k}\cos(\varphi_s+kq)
```

is independent of ``L``. The map of a gap is ``\exp(-L{:}H{:})`` with ``H=H_{\mathrm{d}}+H_{\mathrm{rf}}`` and ``H_{\mathrm{d}}=\mathrm{energy\_term}-p_s``. As ``L\to 0`` at fixed ``V``, the drift piece is proportional to ``L`` and drops out, while ``L{:}H_{\mathrm{rf}}{:}={:}G{:}`` stays finite. The thin map is the pure kick ``\exp(-{:}G{:})``.

In ``(q,P)`` the generator ``G`` depends only on ``q``, so the flow is closed:

```math
q'=q,\qquad
P'=P+\frac{ZV}{p_0 c}\sin(\varphi_s+kq).
```

Arrival time is fixed and the energy changes by ``ZV\sin(\varphi_s+\varphi)``. Back in ``(l,\delta)``, with no series expansion,

```math
\begin{aligned}
E'&=E+ZV\sin(\varphi_s+\varphi),\\
\delta'&=\sqrt{E'^2-m^2c^4}\,/\,(p_0 c)-1,\\
l'&=l\,\frac{v(\delta')}{v(\delta)}.
\end{aligned}
```

Transverse variables are unchanged. This is the thin-cavity kick. A thin cavity has no continuous Hamiltonian: the flow is asked to reject it, and the kick is applied instead.

### Drift–kick–drift and the BCH remainder

The production map of a thick cavity is the flow of ``H_{\mathrm{d}}+H_{\mathrm{rf}}``. The split "drift ``L/2``, thin kick, drift ``L/2``" is only a comparison. For ``X=-L:H_{\mathrm{d}}:`` and ``Y=-L:H_{\mathrm{rf}}:``,

```math
e^{X/2}e^{Y}e^{X/2}
=\exp\!\Big(X+Y-\tfrac{1}{24}[X,[X,Y]]-\tfrac{1}{12}[Y,[X,Y]]+\cdots\Big).
```

The single commutator ``[X,Y]`` cancels between the two half-drifts. Both ``X`` and ``Y`` are proportional to ``L``, so the first surviving error is cubic:

```math
e^{-(L/2){:}H_{\mathrm{d}}{:}}\,
e^{-L{:}H_{\mathrm{rf}}{:}}\,
e^{-(L/2){:}H_{\mathrm{d}}{:}}
=e^{-L{:}(H_{\mathrm{d}}+H_{\mathrm{rf}}){:}}\bigl(1+O(L^3)\bigr).
```

An asymmetric split (the whole drift, then the kick) would leave an ``O(L^2)`` error. The thin cavity is the ``L\to 0`` limit of the same Hamiltonian system, not a different element.

### Symplecticity of the kick

``\exp(-{:}G{:})`` is a Hamiltonian flow, so it is symplectic in ``(q,P)``. The point transformation to ``(l,\delta)`` is canonical, so the kick is symplectic in ``z`` at any amplitude: ``J^{\mathsf T}SJ=S``. In normalised variables this fails only when the reference momentum itself is replaced after the kick. That step is acceleration, [Acceleration](acceleration.md), not a defect of ``G``.

A cavity with ``\sin\varphi_s\neq 0`` accelerates the reference particle. The Lie product used for a storage-ring factor is written for maps that fix the origin, so the Taylor polynomial of a net-gain cavity has a linear term in ``l`` and is a different problem from the stationary bucket.
