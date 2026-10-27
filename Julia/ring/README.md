# LieRing

Modular Julia package for **symplectic ring modelling** with Lie maps and
exact Hamiltonians. Coordinates are canonical
`(x, p_x, y, p_y, l, δ)`; user angles `(x', φ)` are for input/output only.

Install (from a clone):

```julia
using Pkg
Pkg.develop(path = raw"C:\GIT_REPOS\BELL\Julia\ring")   # or your path
using LieRing
```

## Minimal example: define a lattice in code

```julia
using LieRing

# Reference: deuteron at 270 MeV kinetic
ref = reference_from_kinetic(DEUTERON, 270.0)
brho = ref.brho_Tm

# FODO cell (lengths in m, normalised strengths)
L_d = 2.0; L_q = 0.4; L_b = 1.5
k1  = 0.3          # 1/m^2
h   = 0.05         # 1/m  (sector bend)
elems = [
    Drift("D1", L_d),
    Quadrupole("QF", L_q,  k1),
    Drift("D2", L_d),
    SBend("B",  L_b, h * L_b, 0.0),   # angle = h*L, k1=0
    Drift("D3", L_d),
    Quadrupole("QD", L_q, -k1),
    Drift("D4", L_d),
    SBend("B2", L_b, h * L_b, 0.0),
]
lat = Lattice("fodo", elems, ref, 1)   # harmonic placeholder

C  = circumference(lat)
R  = ring_matrix(lat)
po = periodic_optics(R, ref, C)
println("Qx = ", po.x.Q, "  Qy = ", po.y.Q, "  alpha_c = ", po.alpha_c)
```

## From an OptiM file

```julia
of  = parse_optim(raw"path\to\lattice.opt")   # uses CODATA particle mass
lat = of.lattice
po  = periodic_optics(ring_matrix(lat), lat.ref, circumference(lat))
```

## Typical workflow

| task | entry points |
|------|----------------|
| Optics | `ring_matrix`, `periodic_optics`, `optics_along` |
| Maps | `element_lie_map`, `build_factors`, `ring_taylor_map` |
| 2nd order | `path_expansion`, `blond_alphas`, `chromaticity_orders`, `dispersion_orders` |
| Tracking | `exact_track_element`, `exact_track_turn`, `track` (Lie factors) |
| RF / ramp | `Cavity`, `insert_cavity`, `accelerate`, `renormalize` |
| Optional sextupole fit | `Optimizers` (`solve_A`, `solve_B`) — not used by default |

Physics: see `docs/build/index.html` (Beam Physics Manual).  
Code map: Programmer Manual in the same docs tree.

## Applications

Applications live **outside** this package and depend on it via `Pkg.develop`, e.g.

`C:\GIT_REPOS\BELL\Julia\apps\equilibrium_energy\`

## Scripts in this repo

Under `scripts/` — regression / demos for the Nuclotron-based `magnetic_2`
lattice (not required to use the package).
