# LieRing - symplectic ring modelling on Lie algebra (phase 1)

Own modular Julia program (folder `ring`, own `Project.toml`).  All physics is
derived here from element Hamiltonians in OUR coordinates; COSY implementation
was not consulted.  COSY/OptiM numbers are used only as external checks
(`scripts/reference_values.jl`).

Run (from this folder):

    julia --project=. scripts/01_linear_optics.jl   # phase 1: linear optics + plot
    julia --project=. scripts/02_lie_maps.jl        # 2nd-order Lie maps, tracking, Taylor-map check

## 1. Module list (`src/`)

| file | module | role |
|---|---|---|
| `Constants.jl` | Constants | c, Particle (deuteron), RefParticle (gamma0, beta0, p0, Brho), coordinate definitions, exact canonical <-> user transform |
| `Jets.jl` | Jets | TaylorSeries JetSpace: 6 variables, order 3 (H) / map order 2; truncation, printing |
| `Elements.jl` | Elements | Drift, SBend, Quadrupole, Sextupole, RFCavity (thin), Marker, Lattice |
| `Hamiltonians.jl` | Hamiltonians | exact H of each element, expanded by TaylorN arithmetic to order 3 |
| `LieOperators.jl` | LieOperators | Poisson bracket, exp(:f:) Lie series, element Lie map exp(-L:H:), Dragt-Finn factor R o exp(:f3:), concatenation |
| `OptiMParser.jl` | OptiMParser | reader of OptiM `.opt` (variables, header, Initial, sequence, element list) |
| `LinearOptics.jl` | LinearOptics | R = exp(L S A), one-turn matrix, closed Twiss, D, D', tunes, alpha_c, propagation + closure check |
| `Symplectic.jl` | Symplectic | `J^T S J - S`; exact Jacobians of maps and of tracking code (TaylorN AD) |
| `Tracking.jl` | Tracking | factor-by-factor symplectic tracking (implicit midpoint per factor), thin RF kick |
| `TaylorMap.jl` | TaylorMap | INDEPENDENT truncated Taylor map of the ring (RK4 in TaylorN arithmetic), user-coordinate map, chromaticity, alpha_0/alpha_1 |
| `NormalForm.jl` | NormalForm | linear normalisation (done); nonlinear normal form = stub (phase 4) |
| `Plotting.jl` | Plotting | beta / D / layout PNG |
| `LieRing.jl` | LieRing | includes everything, re-exports |

## 2. Coordinates and conventions

Independent variable: arc length `s` along the reference orbit.

User coordinates `u = (x, x', y, y', phi, delta)`:

* `x, y` [m]; `x` radially outward in bends.
* `x' = px/ps`, `y' = py/ps`, `ps = sqrt((1+delta)^2 - px^2 - py^2)` (outside bends = dx/ds).
* `phi = -h * omega0 * dt` [rad], `dt = t - t_ref(s)`, `omega0 = 2 pi beta0 c / C`.
  **Sign: late particle (dt > 0) has phi < 0.**  h = 1 when no RF is present.
* `delta = (p - p0)/p0` (momentum, not energy).

Canonical coordinates used inside all maps: `z = (x, px, y, py, l, delta)`,
`px = Px/p0`, `l = -v(delta) dt` [m].

Derivation: in the s-Hamiltonian with longitudinal pair `(t, -E)`, apply
`F2(t, delta; s) = -(t - s/v0) E(delta)/p0`; then `l = dF2/ddelta = -(t - s/v0) v(delta)`
and H gains `dF2/ds = E(delta)/(p0 v0)`, giving the energy term

    e(delta) = sqrt((1+delta)^2 + 1/(beta0 gamma0)^2)/beta0 - 1/beta0^2
             = delta + delta^2/(2 gamma0^2) - delta^3/(2 beta0^2 gamma0^2) + ...

Exact relations: `phi = 2 pi h l / (C v/v0)`, `x' = px/ps`.
Hamilton: `dz/ds = [z, H] = S grad H`, `[f,g] = sum_k f_q g_p - f_p g_q`,
element map `z(L) = exp(-L :H:) z`.

## 3. Element Hamiltonians (exact; expanded to order 3 by the code)

    Drift        H = e(delta) - ps
    Quadrupole   H = e(delta) - ps + k1/2 (x^2 - y^2),            k1 = G/(B rho)
    Sextupole    H = e(delta) - ps + k2/6 (x^3 - 3 x y^2),        k2 = B''/(B rho)
    Sector bend  H = e(delta) - (1+hx) ps + h x + h^2 x^2/2 + k1/2 (x^2 - y^2),  h = angle/L
    RF (thin)    E -> E + Z V sin(phi_s + phi) at fixed dt;  delta, l recomputed (canonical)

Bend: `-(1+hx) a_s = h x + h^2 x^2/2` reproduces the uniform design field in
the curved frame; `dH/dx = 0` at z = 0 (reference arc is a solution).
No fringe fields (hard edge, sector).  The phase-1 maps are order 2 in z,
so phi/l contains delta^2, x^2, x'^2, x*delta, ... terms (printed by script 02).

## 4. Roadmap

1. **Linear optics** (DONE): parser, exp(L S A) element matrices, closed Twiss,
   D, D', tunes, alpha_c, closure check, plot, comparison with OptiM/COSY.
2. **2nd-order Lie maps + tracking** (DONE, first version): H3 -> Lie maps,
   Dragt-Finn f3 per element and for the ring, symplecticity per element/ring,
   factor-by-factor symplectic tracking.  Next: sextupole families from file,
   chromaticity correction, 3rd order (JET_ORDER = 4).
3. **Taylor-map comparison** (DONE, first version): independent RK4/DA map,
   agreement with Lie product, user-coordinate map, chromaticity, alpha_0/1.
   Next: map-coefficient comparison with COSY output maps after an explicit,
   documented coordinate conversion.
4. **Normal form** (DONE): one-turn object is the Lie product of Dragt-Finn
   factors (BCH to f4). Homological equation on the generators
   (divisor 1 - exp(-i(j-n)·μ)). δ_eq is the torus average of laboratory δ.
5. **RF / synchrotron averaging**: thin RF in the one-turn product at φ_s=π;
   finite cavity and ramp remain a separate path (origin-shifting maps).

## 5. Phase-1 results (magnetic_2.opt) - see script output for full precision

* Lattice: 84 elements, 16 SBend (5.329 m, 22.5 deg), 40 quads, 10 sextupoles (all S = 0),
  C = 117.388296 m. Note: the .opt header comment says "8 octants x 12 dipoles, 1.44 m",
  but the actual lattice of magnetic_2.opt is 2 arcs x 8 bends of 5.33 m.
* "Number of periods=2" in the file, but the listed sequence already bends 360 deg ->
  used once (as in the COSY deck).
* Closed optics at s = 0: beta_x = 19.367522 m, beta_y = 4.054470 m (OptiM Initial 19.3675 / 4.05447),
  alpha and D' agree with OptiM to ~1e-10; D_x = 0.01109 m (OptiM Initial lists 0.0).
  CLOSED check: all closure errors < 1e-12.
* Tunes Qx = 2.46172014, Qy = 2.46133733 (COSY wp_natural: frac 0.44628 / 0.47668 - differ by ~0.015, open).
* alpha_c = 0.23047028 (also from (1/C) int h D ds), COSY alpha0 = 0.23036260.
* Natural chromaticity (order-2 map): xi_x = -1.4615, xi_y = -2.2255
  (COSY chrom_check_I0_FR0: -1.8408 / -1.8477; textbook -(1/4pi) int beta K ds gives -2.81 / -2.93,
  which omits the bend/dispersion terms that the exact H contains). Open item.
* alpha_1 = 0.339 (no sextupoles); Lie product vs RK4 Taylor map agree to 5e-8 in 2nd-order coefficients;
  1-turn tracking Jacobian symplectic error 1.4e-14 (truncated Taylor map at the same point: 2e-4).

## 6. Phase 2 (done) - second order, sextupoles, finite RF, acceleration

New modules: `SecondOrder.jl` (D2, alpha0/1, eta0/1, xi1/xi2, longitudinal H), `ExactIntegrator.jl`
(exact-H symplectic integrator, implicit midpoint + Yoshida-4, finite cavity, renormalisation),
`Sextupoles.jl` (families at maxima, response, solutions A/B), `Acceleration.jl` (ramp, tune fit),
`Plotting2.jl`. Jets now JET_ORDER = 4 / MAP_ORDER = 3; LieFactor has f3 and f4.
Scripts: `03_second_order.jl`, `04_sextupoles.jl`, `05_rf_acceleration.jl`.

Conventions: dC/C = alpha0 d + alpha1 d^2, dT/T0 = eta0 d + eta1 d^2 (eta0 = alpha0 - 1/g^2,
eta1 = alpha1 - alpha0/g^2 + 1.5 b^2/g^2 + 1/g^4); Q(d) = Q0 + xi1 d + xi2 d^2/2 (xi2 = d2Q/dd2).

Finite cavity Hamiltonian (own derivation, on-axis longitudinal field E_s = -(V/L) sin Phi):

    H = e(delta) - ps + A cos(Phi),  Phi = omega (s - L/2)/v0 - omega l/v(delta) - phi_s,
    A = (Z V/(p0 c)) c/(omega L),  omega = h omega0
    -> energy gain Z V T sin(phi_s + phi), T = sin(th/2)/(th/2), th = 2 pi h L / C.

Acceleration: after the cavity the synchronous particle becomes the new reference
(p0' = p0(1+delta_s), delta' = (1+delta)/(1+delta_s) - 1, px' = px/(1+delta_s), dt' = dt - dt_s);
magnets ramp with p0, RF frequency = h f0(p0).  Injection 10 MeV (5 MeV/u) is an ASSUMPTION
(not in magnetic_2.opt).

Longitudinal motion (turn n):  dphi/dn = -2 pi h (eta0 d + eta1 d^2),
dd/dn = (Z V T/(b^2 E))(sin(phi_s+phi) - sin phi_s);  Q_s = sqrt(h eta0 Z V T cos phi_s/(2 pi b^2 E)),
stable iff eta0 cos(phi_s) > 0 (our phi sign) -> below transition phi_s in (90,180) deg.

Next steps: coefficient-level comparison with COSY after a coordinate dictionary,
transverse RF / fringe fields, spin.


## Status (v0.9)
The ring is a product of Lie factors through the normal form. See HISTORY.md.
Docs: docs/build/index.html and docs/pdf/LieRing_Beam_Physics_Manual.pdf.
