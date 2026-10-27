# Longitudinal motion, separatrix and synchrotron tune

*`src/SecondOrder.jl` (`LongitudinalParams`, `long_H`, `fixed_points`,
`separatrix`, `bucket_height`, `synchrotron_tune`), plots in
`src/Plotting.jl`, script `05_rf_acceleration.jl`.*

## Turn-by-turn map and smooth Hamiltonian

With one cavity per turn, the motion of ``(\phi,\delta)`` per turn ``n`` is

```math
\phi_{n+1}=\phi_n-2\pi h\left(\eta_0\delta_{n+1}+\eta_1\delta_{n+1}^2\right),\qquad
\delta_{n+1}=\delta_n+A\left(\sin(\phi_s+\phi_n)-\sin\phi_s\right),\quad A=\frac{ZVT}{\beta_0^2E_0}.
```

(The tracking itself does not use this map; it integrates the exact
Hamiltonians element by element. The map and the Hamiltonian below are used
for formulas and for drawing the separatrix.) Treating ``n`` as continuous,

```math
H(\phi,\delta)=-2\pi h\left(\frac{\eta_0\delta^2}{2}+\frac{\eta_1\delta^3}{3}\right)+A\left(\cos(\phi_s+\phi)+\phi\sin\phi_s\right).
```

Fixed points on ``\delta=0``: stable ``\phi=0``, unstable
``\phi_u=\pi-2\phi_s``; ``\eta_1`` adds a second pair at
``\delta\approx-\eta_0/\eta_1`` (far outside the bucket here, 0.69 at 10 MeV).
The **separatrix** is the level line ``H(\phi,\delta)=H(\phi_u,0)``; the bucket
is the region around ``(0,0)`` inside it. Stability requires
``\eta_0\cos\phi_s>0``: below transition (``\eta_0<0``) and with our sign of
``\phi`` this means ``90^\circ<\phi_s\le180^\circ``; acceleration uses
``\phi_s=150^\circ`` (gain ``ZVT\sin\phi_s>0``), the flat top ``\phi_s=180^\circ``.

## Synchrotron tune

Small amplitude:

```math
Q_s^{\rm smooth}=\sqrt{\frac{h\,\eta_0\,ZVT\cos\phi_s}{2\pi\beta_0^2E_0}},\qquad
\cos(2\pi Q_s^{\rm kick})=1-\frac{(2\pi Q_s^{\rm smooth})^2}{2}\ \ \text{(one kick per turn)}.
```

| ``E`` [MeV], ``\phi_s=180^\circ`` | tracking | smooth | kick |
|---|---|---|---|
| 10 | 0.0784653 | 0.0777018 | 0.0784950 |
| 270 | 0.0129407 | 0.0129376 | 0.0129411 |

At 10 MeV ``Q_s`` is large and the discrete-kick formula is needed (1% vs
0.04%).

## Ramp 10 -> 270 MeV

Parameters: ``L=1`` m, ``V=0.2`` MV, ``h=5`` (``f_{RF}`` 1.313 -> 6.200 MHz),
``\phi_s=150^\circ``; gain 0.099702 MeV/turn; 270 MeV reached at turn 2608,
then 1000 turns at ``\phi_s=180^\circ``. Bucket half-height at injection
``+0.02425/-0.02369`` (asymmetric because of ``\eta_1``), stationary bucket at
270 MeV ``0.00976``.

| result | value |
|---|---|
| tracked | 64 particles x 3608 turns (330 s) |
| mean gain per turn | 0.099698 MeV (expected 0.099702) |
| x amplitude 2 mm -> | 0.8701 mm (ratio 0.43506) |
| y amplitude 2 mm -> | 0.8658 mm (ratio 0.43291) |
| prediction ``\sqrt{p_{inj}/p_{fin}}`` | 0.43142 |
| ensemble particles lost from the bucket | 0 of 60 |
| longitudinal invariant (turn 300 -> 2493) | 1.0000 -> 0.9508 |

The first ~200 turns are not adiabatic (gain = 1% of ``E_{kin}`` per turn at
10 MeV), which explains most of the 5% drift of the longitudinal invariant.

Ensemble: 60 particles with ``\phi\sim N(0,0.6)`` rad,
``\delta\sim N(0,0.25\,\hat\delta_{\rm bucket})``, kept if inside the injection
separatrix, on their dispersion orbit, ``x'=y'=0``. This is deliberately
mismatched to the bucket so that filamentation is visible. ``\eta_1`` and
``\alpha_1`` enter the separatrix; ``\alpha_0=0.23047028``,
``\alpha_1=0.33900771`` from the ring map (energy independent, as expected
for scaled magnets).

![Acceleration overview (05_rf_acceleration.jl). Top left: reference kinetic energy vs turn; dashed line: start of the flat top. Top right: longitudinal phase space at the end of the ramp (turn 2608); black: analytic separatrix H = H(phi_u,0) of the smooth Hamiltonian with eta0, eta1 at that energy and phi_s = 150 deg; red cross: unstable fixed point; green: stable fixed point; blue dots: the 60 ensemble particles; orange: test particle #4 (delta0 = 1e-3) over the last 300 turns (small, near the centre). Bottom left: betatron amplitudes of the x and y test particles vs turn and the prediction A0 sqrt(p_inj/p0(n)). Bottom right: synchrotron tune; lines: smooth and one-kick formulas evaluated at every turn (step at turn 2608 where phi_s changes to 180 deg); markers: tracking, tune fitted to delta of particle #4 in 200-turn windows moved by 50 turns, plotted at the window centre (windows containing the ramp end are skipped).](assets/acceleration.png)

**Animations** (HTML only; the PDF shows the static figure):
`assets/longitudinal_ramp.gif` - ensemble every 20 turns during the ramp
and flat top, with the separatrix of the current energy and ``\phi_s``; right
panel with ``\delta`` divided by the current bucket half-height (removes the
adiabatic shrinking and shows the filamentation).
`assets/longitudinal_elements.gif` - one turn at injection, frame = exit of
each element; right panel: change since the turn start (only the cavity
changes ``\delta``; ``\phi`` changes in every element through ``\eta``-like
path/velocity terms).

![Animation: ensemble during the ramp (GIF)](assets/longitudinal_ramp.gif)

![Animation: one turn element by element (GIF)](assets/longitudinal_elements.gif)
