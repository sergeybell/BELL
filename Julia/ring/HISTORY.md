# HISTORY.md - version notes for Julia/ring (LieRing)

Local tags and commits live in `C:\GIT_REPOS\BELL` (never pushed by this
work). Phases 1 and 2 were not snapshotted separately; their content is
summarised here so a future reader / AI can reconstruct the thought process
without digging through chat. The code as committed under tag `v0.4` contains
everything of phases 1-3.

## v0.1 - Julia learning examples (sibling folder)

See commit / tag `v0.1` of `Julia/` (examples 01-06, not this `ring/` tree):
hello, matrices, FODO linear optics, dual numbers, Taylor maps, Poisson
brackets. Goal: learn Julia and DA / Lie tools before writing the ring code.

## v0.2 - Phase 1: modular Lie ring + linear optics (historical)

**Goal.** Own modular Julia program for symplectic ring modelling built on Lie
algebra, in a new folder `Julia/ring` with its own `Project.toml`. All physics
derived from Hamiltonians in coordinates `(x,x',y,y',phi,delta)`; COSY
*implementation* forbidden, COSY / OptiM *numbers* used only as external
checks. Particle: deuteron. Lattice: `magnetic_2.opt` (Nuclotron-based).

**What was built.**

* Modules: `Constants`, `Jets`, `Elements`, `Hamiltonians`, `LieOperators`,
  `OptiMParser`, `LinearOptics`, `Symplectic`, `Tracking`, `TaylorMap`,
  `NormalForm` (stub), `Plotting`, main `LieRing`.
* Element types: Drift, sector SBend (no fringe), Quadrupole, Sextupole,
  thin RFCavity (placeholder), Marker.
* Hamiltonians expanded to order 3 (maps order 2); Dragt-Finn factors
  `R o exp(:f3:)`; factor-by-factor symplectic tracking; independent DA
  (RK4) Taylor map of the ring for comparison.
* Linear optics: closed Twiss beta/alpha, dispersion D/D', tunes, alpha_c,
  circumference, layout plot `output/magnetic_2_twiss.png`, "closed" check.
* Scripts: `01_linear_optics.jl`, `02_lie_maps.jl`, `reference_values.jl`.

**Key results vs OptiM / COSY.**

| quantity | LieRing | reference | note |
|---|---|---|---|
| C | 117.388296 m | | |
| beta_x(0), beta_y(0) | 19.3675, 4.05447 m | OptiM same | agree to 5 digits |
| D'(0) | -5.0616e-4 | OptiM | agree |
| D(0) | 0.011094 m | OptiM 0 | OptiM input, not closed |
| Qx, Qy | 2.46172, 2.46134 | COSY frac 0.446 / 0.477 | differ by ~0.015 |
| alpha_c | 0.230470 | COSY 0.230363 | close |

**Open issues left at end of phase 1.** Tunes differ from COSY by ~0.015
(suspect: fringe / element model in the COSY deck). Sextupoles in the file
are all zero. Finite RF and acceleration not yet done. Coordinate transform
user <-> canonical was exact but not yet documented or checked as a
polynomial.

## v0.3 - Phase 2: second order, chromaticity, sextupoles, RF + ramp (historical)

**Goal.** Second-order optics, chromaticities, optional sextupole families,
finite-length RF cavity, acceleration 10 -> 270 MeV with adiabatic damping.

**What was added.**

* Orders raised: `JET_ORDER=4`, `MAP_ORDER=3`; `f4` in Dragt-Finn.
* `SecondOrder.jl`: D2, alpha0/1, eta0/1, xi1/2 from the map; longitudinal H text.
* `ExactIntegrator.jl`: Yoshida 4th-order implicit midpoint of the exact H,
  including finite cavity with transit-time factor; `renormalize` for the
  reference after the cavity (= adiabatic damping).
* `Sextupoles.jl` (later split): place families at D / beta_x / beta_y maxima,
  response, example solutions A/B.
* `Acceleration.jl`, `Plotting2.jl`; scripts 03, 04, 05.
* Injection energy ASSUMED 10 MeV (not in the .opt); RF: L=1 m, V=0.2 MV,
  h=5, phi_s=150 deg / 180 deg flat.

**Key numbers (natural lattice, sextupoles 0).**

| | value |
|---|---|
| D2(0), D2'(0) | -10.055 m, 0.459 |
| alpha0, alpha1 | 0.230470, 0.339008 |
| eta0, eta1 (270 MeV) | -0.533784, 1.017209 |
| xi1_x, xi1_y | -1.4615, -2.2255 |
| xi2_x, xi2_y | 8.555, 8.104 |
| cavity gain vs ZVT sin(phi_s) | 0.099701 vs 0.099702 MeV |
| thin-limit delta (L=0.001 m) | diff 1e-11 |
| Qs track / formula (270 MeV, phi_s=180) | 0.012941 / 0.012938 |
| amp damping 2 mm -> | 0.87 mm; pred sqrt(p_inj/p)=0.431 |

Optional A (xi1=0, alpha1=0) and B (xi1=0, min rms D2) computed; default
lattice left at k2=0. First ~200 turns of the ramp not adiabatic (1% gain
per turn at 10 MeV).

**Assumptions made explicit later (see Beam Physics Manual).** Magnets ramp
with p0; RF frequency locked; no transverse RF fields; no fringe / space
charge; periods=2 interpretation; OptiM D(0)=0.

## v0.4 - Phase 3: restructure + Beam Physics Manual (this commit)

**Code.** All element types + lattice-editing helpers in `Elements.jl`.
Optimisation moved to `Optimizers.jl` (NOT used by default; only
`04_optional_sextupole_optimization.jl`). `Plotting2` merged into
`Plotting`. Sextupoles default to zero. Separatrix / bucket helpers in
`SecondOrder`; GIFs of longitudinal phase space (ensemble during ramp;
element-by-element one turn). New script `06_coordinates.jl` proves the
user <-> canonical transform is exact nonlinear (polynomial residual 0,
large-amplitude round trip 0, linearised residual 7.6e-3).

**Docs.** Beam Physics Manual (physics only; programmer stub): HTML via
Documenter.jl in `docs/build/`, PDF via xelatex in
`docs/pdf/LieRing_Beam_Physics_Manual.pdf`. Includes a dedicated section on
the coordinate transform, every element Hamiltonian, RF model, acceleration,
second-order optics, chromaticity, sextupoles (optional), Assumptions, Code
map, and captioned figures / GIFs.

**How to open.** HTML: open `docs/build/index.html`. PDF: open
`docs/pdf/LieRing_Beam_Physics_Manual.pdf`. Rebuild: `julia --project=docs
docs/make.jl` and `julia docs/make_pdf.jl`.

## v0.5 - Canonical manuals, BLonD drift, fixed GIF scales

**Accepted design choices** (removed from "hypotheses"): injection 10 MeV;
RF L/V/h/phi_s as chosen; longitudinal cavity field only; magnets ramp with
p0; no fringe/space charge/impedances; test amplitudes; optimiser optional;
periods=2 is a file error (use 360 deg once); D(0)=0.011 m is the exact
closed solution; mass = CODATA deuteron 1875.61294257 MeV.

**GIF right panel (longitudinal_elements).** Was auto-zoomed to
1.2 * max particle change over the turn (arbitrary). Now fixed:
  Delta phi axis = +/- 1.5 |2 pi h eta0| d_ens   (full-turn slip scale)
  Delta delta axis = +/- 1.5 A                   (RF kick amplitude)
Left panel: absolute (phi,delta) with analytic separatrix, bucket-height axes.

**alpha_1 / BLonD.** Nonlinear longitudinal term IS included (exact element
H, or alpha0/1/2 in BLonD exact drift). BLonD formulas (eqs 19-20, eta 15)
implemented as blond_exact_drift / blond_simple_drift / blond_alphas.
Magnetic one-turn: exact vs BLonD exact rel err ~2e-6; simple ~1e-3 at
|delta|=1e-3. With RF, 32 turns: max |dphi|=1.7e-4, |ddelta|=1.5e-5.

**Manuals rewritten.** Beam Physics: impersonal, canonical variables only,
no Julia file references. Programmer Manual: modules, scripts, design
choices. Both HTML (Documenter) and PDF (xelatex).

## v0.6 - Reusable package + equilibrium-energy app

**Package.** `Julia/ring` is now Julia package `LieRing` v0.6.0 (UUID in
Project.toml). Install with `Pkg.develop(path=...)`. README shows a few-line
user lattice in code and OptiM load. Scripts under `scripts/` remain demos.

**App.** `Julia/apps/equilibrium_energy` depends on LieRing via develop.
Studies delta_eq from note ch.1-3 (magnetic): tests A (dm scan), B (Lambda
vs xi J), C (cubic phase). Note phase convention in app output (phi>0 late).
Stationary bucket phi_s=pi, 270 MeV, magnetic_2, exact integrator.

**Parser fix.** With CODATA mass, bend is 360.08 deg (not machine-eps to 2pi);
periods=2 no longer doubles the sequence (use once if within 1 deg of 360).

**Normal form.** Plan documented in NormalForm.jl for 6D NF with RF ->
delta_eq(Jx,Jy,Js); not implemented yet.

**Docs.** Beam Physics chapter "Equilibrium momentum shift".

### v0.6b - equilibrium_energy app results

App `Julia/apps/equilibrium_energy`: tests A/B/C vs note ch.1-3.
delta_eq from harmonic LS fit (constant + 1st + 2nd harmonic); note phase
convention (phi>0 late). At dm=2e-3, rel err vs -eta1 dm^2/(2 eta0) ~ 2.5%;
Y-plane Lambda at 20 mm within 1%; Omega/omega_s matches 1-phi^2/16 to 1e-4.
Beam Physics chapter `equilibrium.md` updated with the table.
