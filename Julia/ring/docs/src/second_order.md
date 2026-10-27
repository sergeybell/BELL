# Second-order dispersion and momentum compaction

## Closed orbit in ``\delta``

```math
z_\perp(\delta)=D_1\delta+D_2\delta^2+O(\delta^3),
```

from ``M_\perp(z_\perp(\delta),\delta)=z_\perp(\delta)`` order by order.

## Path length and ``\alpha_i``

On the closed orbit the geometric circumference admits the expansion

```math
\frac{C(\delta)}{C_0}=1+\alpha_0\delta+\alpha_1\delta^2+\alpha_2\delta^3+\cdots.
```

With the canonical coordinate ``l`` related to time of flight by
``l=-v\Delta t``, one turn gives
``\Delta l=-C(\delta)+C_0\,v(\delta)/v_0``, hence

```math
\alpha_0=-\frac{l_1}{C}+\frac{1}{\gamma_0^2},\qquad
\alpha_1=-\frac{l_2}{C}+c_2,\qquad
\alpha_2=-\frac{l_3}{C}+c_3,
```

where ``v/v_0-1=c_1\delta+c_2\delta^2+c_3\delta^3`` and
``c_1=1/\gamma_0^2``. Equivalently, ``\alpha_i`` may be fitted from
``C(\delta)/C_0`` evaluated on closed orbits at several ``\delta``.

## Slip factor

The relative change of revolution period is
``\Delta T/T_0=\eta_0\delta+\eta_1\delta^2+\eta_2\delta^3+\cdots``.
Directly from time of flight,
``\eta_0=-l_1/C=\alpha_0-1/\gamma_0^2``. The BLonD convention (Equations of
Motion, eq. 15) expresses higher ``\eta_i`` through ``\alpha_i``:

```math
\begin{aligned}
\eta_0&=\alpha_0-\frac1{\gamma_s^2},\\
\eta_1&=\frac{3\beta_s^2}{2\gamma_s^2}+\alpha_1-\alpha_0\eta_0,\\
\eta_2&=-\frac{\beta_s^2(5\beta_s^2-1)}{2\gamma_s^2}+\alpha_2-2\alpha_0\alpha_1
+\frac{\alpha_1}{\gamma_s^2}+\alpha_0^2\eta_0-\frac{3\beta_s^2\alpha_0}{2\gamma_s^2}.
\end{aligned}
```

The two expansions of ``\eta_1`` (TOF identity vs BLonD) differ by terms of
order ``\alpha_0^2``; both are consistent when used with the matching drift
formula (see below).

## BLonD longitudinal drift (comparison)

BLonD updates arrival time after an energy kick by (eq. 20, solver `exact`)

```math
\Delta t\leftarrow\Delta t+T_0\Biggl[\Bigl(1+\alpha_0\delta+\alpha_1\delta^2+\alpha_2\delta^3\Bigr)
\frac{1+\Delta E/E_s}{1+\delta}-1\Biggr],
```

with ``\delta`` computed exactly from ``\Delta E``, and by the linearisation
(eq. 19, solver `simple`)

```math
\Delta t\leftarrow\Delta t+T_0\,\eta_0\,\frac{\Delta E}{\beta_s^2 E_s}.
```

Element-by-element integration of the ring Hamiltonians reproduces the
`exact` formula at the level ``\sim 10^{-6}`` relative in ``\Delta\phi`` per
turn when ``\alpha_{0,1,2}`` are taken from the lattice; the `simple` formula
deviates at the level of ``\alpha_1\delta`` (a few ``10^{-3}`` at
``|\delta|=10^{-3}``). Thus ``\alpha_1`` *is* the leading nonlinear
longitudinal term in the BLonD sense; the ring model includes it
automatically through the exact element Hamiltonians (to all orders in the
exact integrator, to map order in the truncated maps).

![Second-order dispersion along the ring.](assets/sextupoles_dispersion2.png)