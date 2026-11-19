# RF motion

The cavity Hamiltonian is derived in [Element Hamiltonians](hamiltonians.md). This chapter is the motion that follows from it: the on-axis equations, the one-turn map, and the synchrotron tune. The smooth ring average is [Longitudinal motion](longitudinal.md). A changing reference momentum is [Acceleration](acceleration.md).

## On-axis equations

On axis ``p_s=1+\delta`` and the thick uniform gap of length ``L`` has

```math
H(l,\delta)=\mathrm{energy\_term}(\delta)-(1+\delta)
+\frac{ZV}{p_0 c\,k L}\cos\varphi,
\qquad
\varphi=\varphi_s+\frac{k\,l}{\beta(\delta)},
\qquad
k=\frac{2\pi h\beta_0}{C}.
```

Hamilton's equations are ``\mathrm{d}l/\mathrm{d}s=\partial H/\partial\delta`` and ``\mathrm{d}\delta/\mathrm{d}s=-\partial H/\partial l``.

The energy term contributes ``v/v_0``. The ``-(1+\delta)`` contributes ``-1``. The RF term depends on ``\delta`` only through ``\beta(\delta)`` in the phase. With

```math
\frac{\mathrm{d}\beta}{\mathrm{d}\delta}=\frac{\beta}{(1+\delta)\gamma^2}
```

the pair is

```math
\begin{aligned}
\frac{\mathrm{d}l}{\mathrm{d}s}
&=\frac{\beta(\delta)}{\beta_0}-1
+\frac{ZV}{p_0 c\,L}\,\frac{l\sin\varphi}{\beta\,(1+\delta)\,\gamma^2},\\[2pt]
\frac{\mathrm{d}\delta}{\mathrm{d}s}
&=\frac{ZV}{p_0 c\,L}\,\frac{\sin\varphi}{\beta(\delta)}.
\end{aligned}
```

The second equation is ``\mathrm{d}p/\mathrm{d}s=(1/v)\,\mathrm{d}E/\mathrm{d}s`` with ``\mathrm{d}E/\mathrm{d}s=ZV\sin\varphi/L``. The first term of ``\mathrm{d}l/\mathrm{d}s`` is the drift: a faster particle moves ahead. The last term exists because ``l=-v\Delta t`` contains the particle's own speed. When the field changes ``v``, ``l`` changes even at fixed ``\Delta t``. In the variables below, where the coordinate is the delay itself, that term is absent.

## Phase and energy offset

```math
\varphi=\varphi_s+\frac{2\pi h}{C}\frac{\beta_0}{\beta(\delta)}\,l,
\qquad
\Delta E=E-E_0=p_0 c\Big(\mathcal{E}(\delta)-\frac{1}{\beta_0}\Big).
```

The pair ``(\varphi,\Delta E)`` is the pair ``(q,P)`` up to constant factors, so it is canonical up to the constant multiplier ``k\,p_0 c`` on the Hamiltonian. On axis

```math
\frac{\mathrm{d}\varphi}{\mathrm{d}s}
=\frac{2\pi h}{C}\Big(1-\frac{\beta_0}{\beta(\delta)}\Big),
\qquad
\frac{\mathrm{d}(\Delta E)}{\mathrm{d}s}
=\frac{ZV}{L}\sin\varphi.
```

The energy kick needs no further conversion. The nonlinearity of ``\delta(\Delta E)`` sits in the phase equation. In ``(l,\delta)`` the kick itself is the nonlinear map.

## One turn

Follow the particle once per turn, at a fixed azimuth. The one-turn map is the composition of every element. With a single thin cavity and no acceleration (``\sin\varphi_s=0``, or ``ZV\sin\varphi_s`` exactly cancelling a constant loss), the same map is applied every turn. The ring is autonomous: a closed orbit, tunes, a normal form and a bucket are properties of that one map.

Without a series expansion the thin kick is the map in [Element Hamiltonians](hamiltonians.md), and the rest of the ring changes the arrival through the time of flight ``T(\delta)``:

```math
\begin{aligned}
\mathcal{E}_{n+1}&=\mathcal{E}(\delta_n)+\frac{ZV}{p_0 c}S_n,\\
\delta_{n+1}&=\sqrt{\mathcal{E}_{n+1}^2-1/(\beta_0^2\gamma_0^2)}-1,\\
l_{n+1}&=\frac{\beta_{n+1}}{\beta_n}\,l_n
-C_0\frac{\beta_{n+1}}{\beta_0}\Big(\frac{T_{n+1}}{T_0}-1\Big),
\end{aligned}
```

with ``S_n=\sin\varphi_n-\sin\varphi_s``. The first piece of ``l_{n+1}`` recomputes ``l=\beta q`` at fixed ``\Delta t``. The second piece is the slip around the ring.

To first order in ``(\delta,l)``, with the linear slip ``\eta_0`` of [Readings of the one-turn map](map_readings.md),

```math
\delta_{n+1}=\delta_n+\mathcal{V}\,S_n,
\qquad
l_{n+1}=l_n-C_0\eta_0\,\delta_{n+1},
\qquad
\mathcal{V}=\frac{ZV}{\beta_0^2 E_0}.
```

The quadratic slip ``\eta_1`` enters at the next order, in ``T/T_0-1=\eta_0\delta+\eta_1\delta^2+\cdots``. The coefficient in front of ``\delta^2`` depends on whether the coordinate is ``l``, ``\delta`` or ``\Delta E``. The root ``T(\delta)=T_0`` does not.

## Synchrotron tune

Linearise the one-turn map: a thin cavity and a ring that advances ``\Delta l=-\eta_0 C_0\delta`` per turn. The kick and the trace give

```math
\cos 2\pi Q_s
=1-\frac{2\pi h\,\eta_0\,ZV\cos\varphi_s}{2\beta_0^2 E_0}.
```

For small ``Q_s``,

```math
Q_s=\sqrt{\frac{h\,\eta_0\,ZV\cos\varphi_s}{2\pi\beta_0^2 E_0}},
\qquad
\eta_0=\alpha_0-\frac{1}{\gamma_0^2}.
```

Stability of the origin requires ``\eta_0\cos\varphi_s>0``. Above transition ``\varphi_s=0`` in this sine convention; below transition ``\varphi_s=\pi``. Comparing ``Q_s`` from the trace of the longitudinal block with this square root checks the sign and the normalisation of ``H_{\mathrm{rf}}``. They differ by a term of order ``Q_s^2``.

A thick cavity is not replaced by this kick inside the tracker. The flow integrates ``H`` over ``L``. The tune formula is the thin, one-kick-per-turn reading of the same ``H``.
