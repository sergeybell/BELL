# Code map

| file | physics implemented |
|---|---|
| `src/Constants.jl` | constants, deuteron, reference particle, ``v(\delta)``, ``e(\delta)``, exact user <-> canonical transform |
| `src/Jets.jl` | truncated polynomial space (TaylorSeries `TaylorN`), orders |
| `src/Elements.jl` | **all** element types (Drift, SBend, Quadrupole, Sextupole, thin RFCavity, finite Cavity, Marker), Lattice, lattice-editing helpers (sextupole placement, cavity insertion) |
| `src/Hamiltonians.jl` | element Hamiltonians as polynomials |
| `src/LieOperators.jl` | Poisson bracket, ``\exp(:f:)``, Dragt–Finn factors ``R,f_3,f_4`` |
| `src/OptiMParser.jl` | reader of OptiM `.opt` files |
| `src/LinearOptics.jl` | linear maps, closed Twiss, dispersion, tunes, ``\alpha_c`` |
| `src/Symplectic.jl` | ``J^TSJ-S`` checks |
| `src/Tracking.jl` | factor-by-factor Lie tracking, thin RF kick |
| `src/TaylorMap.jl` | independent DA (RK4) one-turn Taylor map, closed orbit, tunes vs ``\delta`` |
| `src/NormalForm.jl` | linear normalisation (nonlinear normal form: planned) |
| `src/SecondOrder.jl` | ``D_2``, ``\alpha_{0,1}``, ``\eta_{0,1}``, ``\xi_{1,2}``, longitudinal Hamiltonian, separatrix, ``Q_s`` |
| `src/ExactIntegrator.jl` | exact-H symplectic integrator, finite RF cavity, reference renormalisation |
| `src/Acceleration.jl` | ramp loop, tune fit, element-by-element snapshots |
| `src/Plotting.jl` | all figures and GIFs |
| `src/Optimizers.jl` | OPTIONAL sextupole response and example solutions |
| `scripts/01_linear_optics.jl` | linear optics + OptiM/COSY comparison |
| `scripts/02_lie_maps.jl` | Lie maps, symplecticity, Lie vs Taylor, tracking |
| `scripts/03_second_order.jl` | ``D_2``, ``\alpha/\eta``, chromaticities |
| `scripts/04_optional_sextupole_optimization.jl` | optional optimiser |
| `scripts/05_rf_acceleration.jl` | RF cavity, ramp, separatrix, GIFs |
| `scripts/06_coordinates.jl` | coordinate transformation checks |
