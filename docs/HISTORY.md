# HISTORY.md - version notes for Bell (was Julia/ring, LieRing)

Local tags and commits live in `C:\GIT_REPOS\BELL` (never pushed by this
work). Phases 1 and 2 were not snapshotted separately; their content is
summarised here so a future reader / AI can reconstruct the thought process
without digging through chat. The code as committed under tag `v0.4` contains
everything of phases 1-3. From v0.10 the program is the Julia package `Bell`
in `Julia/`, not the folder `Julia/ring`. This is not a product release.

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

## v0.7 - 6D Lie normal form with RF + equilibrium_energy research note (2026-10-27)

**Goal.** Make the Lie / map machinery visible (not only tracking): Dragt-Finn
normal form of the one-turn map *with RF*, derive `delta_eq(J)`, explain the
A-test residual, and write a *separate* Russian research note (not in the Beam
Manual).

**LieRing changes.**

* `NormalForm.jl` rewritten: `one_turn_map_rf`, `rf_kick_taylor`,
  `fixed_point_6d`, `normal_form_rf` (linear CS + synchrotron Twiss,
  Dragt-Finn to order 3), `delta_eq_phase_condition` / `delta_eq_nf`,
  `tunes_nf`, `evaluate_map`.
* NF coefficients of `delta_eq = c_x J_x + c_y J_y + c_s J_s` match the note
  formulas to relative `< 1e-5`:
  NF `(0.146652, 0.223289, 392.952)` vs note `(0.146653, 0.223290, 392.950)`.
* Thin-RF fixed point `delta* ~ 0`; finite cavity `delta*(V=0.2 MV, L=1 m)
  = +1.073e-7`, linear in V, grows with L (TOF / cavity energy-gain asymmetry).

**App `equilibrium_energy`.**

* `scripts/run_v07.jl`: four-way `delta_eq` (a note, b NF, c exact tracking,
  d truncated Taylor) + residual V/L scans + long-term mean Taylor vs Lie.
* Separate note `apps/equilibrium_energy/note/main.pdf` (Russian, 5 pp.) with
  captioned figures; not part of the Beam Manual.

**Key residual origin.** The A-test floor `~1e-7` independent of `dm` is the
6D closed-orbit `delta*` of the *finite* RF cavity, not a tracking error.
Tracking `delta_eq ≈ note + delta*`; thin-RF NF / Taylor match the bare note.

**Local only.** Commits + tag `v0.7`; no push.

## v0.8 - thin RF everywhere, explicit 6D Dragt-Finn normal form, models M1-M4 (2026-10-27)

**Goal.** The v0.7 research note was too thin. Make RF a thin kick in every
model of `apps/equilibrium_energy`, remove the finite-cavity `delta*`
contribution from numbers and text (one sentence only: finite cavity is a
separate future study), derive `delta_eq` from an explicit normal form, define
the models, explain every number, and develop the time-evolution / averaging
figures.

**LieRing (`NormalForm.jl`).** New explicit order-by-order normal form
`dragt_finn_nf6d(M)` (+ `nf6d_delta_eq`, `nf6d_actions`, `nf6d_tunes`,
`symplectic_eigen_A`, type `NF6D`):

* fixed point z*, symplectic eigen-normalisation A (A^T S A = S,
  A^-1 R A = rotations; residual 1e-13);
* Dragt-Finn `M_n = rot o e^{:f3:} e^{:f4:}`; order 3: homological equation
  `f3' = f3 + F3 - F3 o rot = 0` solved as a 56x56 linear system built by
  exact map composition (f3 0.26 -> 3e-16); order 4: F4 removes all
  non-resonant quartic terms, `h4 = sum c_ij J_i J_j` (residual 1.5e-13);
* `delta_eq = z*_6 + d.J`, `d_k = (1/2)(d_qq + d_pp) U6 |_0`,
  `U6 = (A T3 T4)_6`, checked by 8^3 torus grid average.
* Result: `d = (0.14665244, 0.22328941, 2.4770639e-3)` vs formula
  `(2pi xi_x/(C eta0), 2pi xi_y/(C eta0), -eta1 sigma_delta^2/eta0)`
  `= (0.14665257, 0.22328971, 2.4770642e-3)`: rel. diff <= 1.4e-6.
  Old `normal_form_rf` / phase-condition fit kept (unused by v0.8 app).

**Constants check.** `DEUTERON` mass 1875.61294257 MeV = CODATA 2018
(m_d c^2 = 1875.612 942 57(57) MeV); c exact; G_d = -0.142987272 (CODATA 2018
-0.1429872722) - nothing changed. The OptiM file mass 1876.544 (= 2*938.272)
is parsed but not used for the reference particle (documented in the note).

**App `equilibrium_energy`.** `scripts/run_v08.jl` (thin RF, kick at s=0):
M1 analytic formula, M2 normal form, M3 exact symplectic element-by-element
tracking + exact thin kick, M4 iterated order-3 Taylor map. Test A (5 dm,
12000 turns), Qs(J) check, Test B (x/y, 5-20 mm, 6000 turns), fast betatron
averaging of path length (within turn, turn-by-turn), slow synchrotron
averaging (running mean, period mean, LS fit, convergence vs N, spectra),
M3 vs M4. Numbers/tables written to `note/gen/*.tex` and `\input` by the
note (no hand-copied numbers). Note rewritten (16 pp., Russian).

Key numbers (thin RF): dm=1e-3: M1 0.9520e-6, M2 0.9546e-6, M3 0.9499e-6,
M4 0.9519e-6; dm=2e-3: M1 3.789e-6, M2 3.831e-6, M3 3.737e-6, M4 3.770e-6.
Test B: M1/M2/M3 agree to 0.05%. Qs(J): NF vs tracking < 1e-6.
No constant offset any more (v0.7 `~1e-7` was the finite-cavity fixed point).

**Local only.** Commit + tag `v0.8`; no push.

## v0.9 - one-turn object is the Lie product (2026-10-27)

**Lie algebra through the normal form.** `concatenate` multiplies Dragt-Finn
factors to degree 4: cubics add after transport by R, quartics receive
`-½[f3b ∘ Ra, f3a]` (sign fixed on pairs of elements against the Lie series,
residual `~1e-16`). Nested brackets are dropped. `ring_lie_factor` is that
product over the lattice; `one_turn_factor_rf` appends the thin kick as one
more factor.

**Normal form on generators.** `dragt_finn_nf6d` takes the one-turn factor.
The homological operator is `F − F ∘ r(μ)` on homogeneous polynomials
(eigenvalues `1 − exp(−i(j−n)·μ)`). Map conjugation is not used to solve it.
`δ_eq` is the average of laboratory `δ` over the torus of the three actions.
On `magnetic_2` with thin RF: `f3` after the cubic step `4e-16`; `d` matches
the time-of-flight formula to relative `1e-6`, same digits as v0.8.

**Checks.** Polynomial substitution (`ring_taylor_map`) and the unexpanded
Hamiltonian integrator remain independent. They do not define the one-turn
object.

**Docs.** Beam Physics chapter on the Lie product and the torus; equilibrium
chapter uses thin RF; research note models M1–M4 rewritten accordingly.
BLonD slip formulae are not part of this path.

Scripts: `08_lie_product.jl`, `09_lie_normal_form.jl`. App `run_v08.jl` builds
M2 from `one_turn_factor_rf`.

## v0.9 - Lie product as the ring object; manuals in HTML

**Object.** The one-turn map is the Dragt–Finn product `ring_lie_factor`.
`observables(lat)` reads Twiss from `F.R`, `η` and `ξ` from the same
polynomial. `MapAnalysis.jl` holds those readers; `TaylorMap.jl` is RK4 only
(a check). Linear optics is the closed form of `exp(:f2:)`, not a second
theory of the lattice.

**App.** `equilibrium_energy` moved to `C:\GIT_REPOS\SCT_2\Julia\equilibrium_energy`.
It still depends on this package via `Pkg.develop`. `Julia/apps` is removed.

**Docs.** Open `index.html` or `docs/build/index.html` in a browser (HTML site,
like BLonD). PDFs: `docs/pdf/LieRing_Beam_Physics_Manual.pdf` and
`LieRing_Programmer_Manual.pdf`. Architecture figures sit on the home page.

## v0.10 - Bell package; lattice in its own file; examples notebook

**Layout.** The ring experiment `Julia/ring` (LieRing) is gone. The program is
the Julia package `Bell` in `Julia/src`: one file per block. A concrete ring
is a structure file under `Julia/lattices/`, not an OptiM object and not a
MAD-X `LINE=(quad, L=..., k1=...)`. OptiM is a converter: it writes the
`.bell` file; `load_structure` builds the `Lattice`.

**Structure file.** Three parts plus a sequence of names:

* `reference` — particle, kinetic energy, harmonic
* `parameters` — named numbers (the optimizer changes these)
* `elements` — one card per named magnet / drift / cavity
* `sequence` — names only; `N * cell` repeats a named sequence

Units are SI (`m`, `rad`, `1/m^2`). FODO is `lattices/fodo.bell`. The
Nuclotron-based ring is `lattices/magnetic_2.bell`, converted from
`magnetic_2.opt`.

**Blocks.** Coordinates, element types, structure parser, Hamiltonians, Lie
algebra, linear optics, higher-order readers, normal form, Lie tracking,
exact-H integrator, reduced longitudinal motion, acceleration, optimizers,
plotting, symplectic checks. Differential algebra (`ring_taylor_map`) stays
an independent check: high-order numbers come from the Lie product
(`observables` / `ring_lie_factor`), not from RK4 in Taylor arithmetic.
`SecondOrder` was split so 6D tracking does not live in the reduced
`(φ, δ)` module.

**Examples.** Numbered scripts are `01`–`06` (Julia language, no Bell) and
`10`–`18` (Bell + structure file). There are no `07`–`09`: those names were
the old `Julia/ring` scripts that became `16`–`18`. Helpers, not examples:

* `print_util.jl` — `banner` / `section` for `10`–`18`
* `reference_values.jl` — OptiM path, COSY numbers for comparison tables
* `linear_optics.jl` — leftover 2×2 transfer-matrix helpers from the
  tutorial FODO; `03_fodo_linear.jl` is self-contained and does not include it

BLonD `(dt, dE)` formulae stay in `16_blond_longitudinal.jl`, not in `src/`.
Notebook `examples/Bell_examples.ipynb` runs every numbered example in order
(executed: 18 code cells, 0 errors).

**Block coverage.** Every `src/` block is called from `10`–`18`. The notebook
table lists which example uses which file. FODO through
`load_structure("lattices/fodo.bell")` is a notebook cell, not a separate `07`.

**Docs.** `docs/architecture.md` (English) says how the blocks connect, where
orders come from, how Lie tracking differs from the exact integrator, and
what the structure file is. Physics chapters remain under `docs/src/`.
Programmer notes: `docs/src/programmer/`.

**Downstream.** `C:\GIT_REPOS\SCT_2\Julia\equilibrium_energy` should
`Pkg.develop` `Julia/`, not `Julia/ring`.

**Not a product.** Package version is `0.10.0`. There is no 1.0 release.

## v0.11 - Fringe-field models (2026-10-30)

**Commits.** `0bf92c1`, `fdb71ef`, `dd9d21a`, `3d067b4`, `0d6454c`, `bd712e2`.

**Code.** New block `src/fringe.jl`: the fringe-field model of each magnet is
chosen when the lattice is built. Models: `NoEdge`, hard edge (linear edge
kick or nonlinear), soft edge (Enge, Tanh, linear ramp, tabulated profile),
integrated by splitting, Magnus or Taylor. Example `19_fringe_models` runs
every model on `magnetic_2`: gap -> 0 limit, convergence, symplecticity.

**Physics fix (`bd712e2`).** The soft dipole edge Hamiltonian now includes
the metric term `h*b1`; the cache key of `Tabulated` profiles is content based.

**Docs.** Beam Physics chapter on fringe fields; Programmer notes for `fringe.jl`.

**Minor (no own version), 2026-10-28..30.** Repository root is the package
(`0531680`, `5109955`); one folder per example and `examples/Project.toml`
(`20c180e`, `0349f24`, `8653af4`, `7360c8b`); new README / .gitignore
(`814259d`); corrected `magnetic_2` lattice from SCT_2 (`55fd17f`);
docstrings for every exported name and LaTeX formulas (`5687cb2`, `a06600e`);
Programmer Manual as HTML site with explicit algorithms (`36b642a`,
`4af0096`); docs site restructured into Physics / Guide / Developer / API
(`62f273f`); PDF rebuild (`d944ac2`).

## v0.12 - Differential algebra and selectable integrators (2026-10-30)

**Commit.** `79e2306`.

**Code.** Pure-mathematics submodules, not re-exported into `Bell`:

* `Bell.DiffAlgebra` (`src/da/`) - general truncated power series
  (differential algebra) in any number of variables and any order.
* `Bell.Integrators` (`src/integrators/`) - selectable ODE integrators
  used by the exact-H tracking and checks.

Use with `using Bell.DiffAlgebra`, `using Bell.Integrators`.

## v0.13 - Lie flow operator to arbitrary order (2026-10-30)

**Commits.** `457b5b5`, `9068936`.

**Code.** `Bell.LieFlow`: element maps as the flow operator `exp(L :H:)`
summed as a Lie series to arbitrary order N, with a Lie-series layer in
`DiffAlgebra`. The series at order N is built from `F = S grad H` with `H`
taken to order N+1; for long / strong elements (`L*rho > 3`) the map is
obtained by scaling and squaring.

**Docs.** Developer notes `docs/site/dev/lie_flow.md`.

## v0.14 - Chromaticity series from the one-turn DA map (2026-10-30)

**Commit.** `846891a`.

**Code.** `chromaticity_series`: exact series `Q(delta)` of the tunes from a
one-turn DA map (LieFlow), to the order of the map.

**Minor.** `convert_optim(...; close_ring=true)` scales bend angles so that
the ring closes exactly to `2 pi` (`26e6cea`).

**Package version** `0.14.0` (not a product release).

## v0.15 - Core review: names, one scheme, research out of the export (2026-11-02)

**Code.** Public names: `lattice`, `record`, `injection_reference`,
`laboratory_particles`, `path_length`, `c_light`, `deuteron` / `proton` /
`electron`, `index_x`…, `position_indices` / `momentum_indices`. Laboratory
coordinates `(x, x', y, y', φ, δ)`. `Cavity` is `LongCavity`. `insert_cavity`
requires `in_element`.

`differential_algebra.jl` is split: integrator, series, Hamiltonian live in
their places. Taylor maps of `H` sit in tracking. `HigherOrder` is
`MapReadings` (Newton closed orbit, integer vs fractional tune,
`chromaticity_series`). `longitudinal.jl` is `Smooth`. Integrators are one
file per method, including implicit midpoint and implicit-midpoint Yoshida.
Exact tracking keeps `rf_kick`, `renormalize` and the cavity; the step is
`method`. HardEdge vs SoftEdge (`Body(-W)∘Soft∘Drift(-W)`) is written in
fringe; `enge_default` is an unverified set. `symplectic_error` is the single
check of a matrix, a series map, a DA map, or a tracking step.

Sextupole families leave `structure.jl`. `src/optimizers/` holds `GridSearch`
and `RandomStart` stubs. Plotting moves to `src/plotting/`. Normal form moves
to `src/research/` and is not re-exported by `using Bell`.

**Tests / CI.** Mathematical pendulum on every integrator
(`benchmarks/pendulum`). GitHub Actions runs `Pkg.test`. Package version
`0.15.0` (not a product release). Examples are not edited in this branch.

## v0.16 - Lie path vs field/DA path; thin RF only (2026-11-02)

**Code.** Two stacks that share elements, coordinates and integrators but do
not cross:

- Lie: polynomial `hamiltonian` → `lie_algebra` / Dragt–Finn → `Tracking`.
- Field/DA: `src/da/fields.jl` (`dz/ds` from multipole fields) → `FieldTrack`
  (Float64) and `src/da/tracking.jl` (`DATracking`, any order N). No
  `hamiltonian_da`. DiffAlgebra stays pure maths; fields and DA tracking sit
  beside it under `src/da/`.

Removed: `LongCavity`, `cavity_field` / `cavity_parameters` / `transit_factor`,
module `LieFlow` (`lie_flow_map`, `hamiltonian_da`). RF is only thin
`RFCavity` + `rf_kick` (numbers or DA). `insert_cavity` / `accelerate` /
`exact_track_*` use the thin kick; `renormalize` after each kick on a ramp.
`exact_integrator.jl` is `field_track.jl` (`FieldTrack`). Smooth stays a
separate averaged `(φ, δ)` model.

**Tests.** DATracking replaces the LieFlow testset (field map vs Lie at order
3, analytic drift/quad, symplecticity, chromaticity_series, thin RF turn).
Package version `0.16.0`. Examples still not edited.

**Optimizers (same branch).** `NelderMead` — local derivative-free fit with
`tol` on the simplex spread in `f`. `RandomStart` — MultiStart: random
initial points in a box, local Nelder–Mead from each. `GridSearch` — fixed
list of points, optional polish. Entry point `optimize`.


**Note (stale references).** Removed with LongCavity/LieFlow: `docs/site/dev/lie_flow.md` (replaced by `dev/da_tracking.md`); `ExactIntegrator` references in the docs now point to `FieldTrack`. Example 14 and `benchmarks/twiss/bell/linear_routes.jl` were moved to the thin `RFCavity` / `one_turn_map_da` path.

## v0.17 - RF cavity s-Hamiltonian; lie/ da/ common/ longitudinal_models/ (2026-11-02)

**Code.** `RFCavity` has an s-Hamiltonian like every element:
`H = energy_term(delta) - ps + G/L`, `G = rf_potential(cav, ref, C) = ZV/(p0 c k) cos(phi_s + k l/(beta0 v/v0))`,
`k = rf_wavenumber = 2 pi h beta0 / C`. New field `L` (`RFCavity(name, V, h, phi_s, L)`; 4-argument form
is thin, `L = 0`). Thin map `exp(-:G:)` equals `rf_kick` (Taylor coefficients to 2e-15); thick map
`exp(-L:H:)` (Lie) or drift/kick/drift splitting (field/DA). `Tracking.build_factors`, `track`,
`ring_lie_factor`, `taylor_map_element`, `element_maps`, `fringe` factors and `DATracking.element_da_map(s)`
have no special RF branch; `element_lie_map(e, ref, path_length)`. New `magnetic_lattice` (RF removed /
thick -> drift), used by MapReadings. Storage ring (`sin phi_s = 0`): constant p0, turn-independent
symplectic one-turn map. Acceleration in H (p0(s), `Delta H = (p0'/p0)[x px + y py + l(1+delta)]`) is a
documented TODO in `docs/site/dev/rf_acceleration.md` (with a summary of the COSY RF procedure); the
`accelerate` path (`rf_kick` + `renormalize`) is unchanged.

**Layout.** `src/common/` (coordinates, elements, structure, optim_converter, map_readings, symplectic,
sextupole_families), `src/lie/` (lie_algebra, hamiltonians, linear, tracking, fringe, normal_form),
`src/da/` (DiffAlgebra, lie_series, fields, field_track, tracking), `src/longitudinal_models/`
(smooth, acceleration: approximate models kept for comparison). Exports unchanged except the additions
`rf_potential`, `rf_wavenumber`, `magnetic_lattice`.

**Tests.** "RF Hamiltonian": thin limit vs `rf_kick`, thick -> thin O(L^2), symplecticity, Q_s vs Smooth
(1.5e-5), 10-turn Lie vs exact field tracking. Package version `0.17.0`.
