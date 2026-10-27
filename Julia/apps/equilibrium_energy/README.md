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