# LieRing

Modular Julia package for **symplectic ring modelling** with Lie maps and
exact Hamiltonians. Coordinates are canonical
`(x, p_x, y, p_y, l, δ)`; user angles `(x', φ)` are for input/output only.

**Open the manuals in a browser** (same idea as BLonD HTML docs):
[`index.html`](index.html) or [`docs/build/index.html`](docs/build/index.html).

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
| Lie product | `element_lie_map`, `build_factors`, `ring_lie_factor`, `one_turn_factor_rf` |
| One-turn readings | `observables` (Twiss, eta, xi from the Lie factor) |
| Normal form | `dragt_finn_nf6d`, `nf6d_delta_eq`, `nf6d_actions` |
| Check (DA map) | `ring_taylor_map` |
| Slip / chromaticity | `path_expansion`, `chromaticity_orders` |
| Tracking | `exact_track_element`, `exact_track_turn`, `track` (factor flows) |
| RF / ramp | `RFCavity`, `Cavity`, `insert_cavity`, `accelerate`, `renormalize` |

Physics: open [`docs/build/index.html`](docs/build/index.html) in a browser.  
Code map: [Programmer Manual](docs/build/programmer/index.html) in the same HTML tree. PDFs are in `docs/pdf/`.

## Applications

Applications live **outside** this package and depend on it via `Pkg.develop`, e.g.

`C:\GIT_REPOS\SCT_2\Julia\equilibrium_energy\`

## Scripts in this repo

Under `scripts/` — regression / demos for the Nuclotron-based `magnetic_2`
lattice (not required to use the package).
