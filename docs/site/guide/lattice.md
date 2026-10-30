# Loading a lattice

Examples: `examples/10_linear_optics`, `examples/15_coordinates`.

```julia
using Bell
lat = load_structure("examples/lattices/magnetic_2.bell")   # Lattice: elements, ref, harmonic
lat.ref.Ekin_MeV, circumference(lat), total_bend_angle(lat)
```

A `.bell` file has four sections (`reference`, `parameters`, `elements`,
`sequence`); cards look like `QF quadrupole L = 0.21 k1 = 0.52` or
`BM sbend L = 5.33 angle = 22.5 deg e1 = 0 gap = 0.05 fint = 0.5`.
An OptiM file is converted once with `convert_optim(opt_file, bell_file)`.

**Another energy:** `with_reference(lat, reference_from_kinetic(DEUTERON, 10.0))`
(normalised strengths kept, maps rebuilt automatically).

**Edits** return a new lattice: `insert_cavity(lat, Cavity("RF", 1.0, 0.2, 5, phi_s))`
returns `(lattice, index)`; `place_sextupole_families` / `set_family_strengths`
add sextupoles (example 13).

**Coordinates:** tracking uses canonical ``z = (x, p_x, y, p_y, l, \delta)``;
convert with `user_to_canonical(u, ref, C, h)` / `canonical_to_user` from user
``(x, x', y, y', \varphi, \delta)`` (example 15).
