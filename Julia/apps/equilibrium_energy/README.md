# EquilibriumEnergy — study of δ_eq (equilibrium momentum shift)

Depends on LieRing (`Pkg.develop`). Physics: note *equilibrium-energy-v3*
chapters 1–3 (magnetic ring: base, betatron, cubic phase).

## Phase convention in THIS app (note convention)

```
phi_note > 0  <=>  particle is LATE
phi_note = - phi_LieRing
(LieRing uses phi = -h ω0 Δt with phi>0 early.)
```

## Run

```
julia --project=. setup.jl          # once: develop LieRing + instantiate
julia --project=. scripts/run_all.jl
```

Outputs under `output/`: `A_delta_eq.png`, `B_lambda.png`, `C_cubic_phase.png`, `run.log`.

Lattice: `magnetic_2.opt`, deuteron, 270 MeV, stationary bucket φ_s = π
(LieRing; late sync in the note), exact integrator in canonical variables.
δ_eq is the constant term of a two-harmonic least-squares fit of δ(n).
## v0.8 (thin RF everywhere)

```
julia --project=. scripts/run_v08.jl     # ~7 min; log in output/run_v08.log
cd note; latexmk -xelatex main.tex       # note/main.pdf
```

Models: M1 analytic formula, M2 6D Dragt-Finn normal form (`dragt_finn_nf6d`),
M3 exact symplectic element-by-element tracking + exact thin RF kick at s=0,
M4 iterated order-3 Taylor map. All numbers in the note come from
`note/gen/*.tex`, written by `run_v08.jl`. `run_v07.jl` (finite cavity) is
kept only for history.