# Loading a lattice

Examples: `examples/01_reference_and_coordinates/run.jl`,
`examples/02_lattice_io/run.jl`.

```julia
using Bell
lat = load_structure("examples/lattices/magnetic_2.bell")   # elements, ref, harmonic
lat.ref.Ekin_MeV, path_length(lat), total_bend_angle(lat)
```

A `.bell` file has four sections (`reference`, `parameters`, `elements`,
`sequence`). Cards look like `QF quadrupole L = 0.21 k1 = 0.52` or
`BM sbend L = 5.33 angle = 22.5 deg e1 = 0 gap = 0.05 fint = 0.5`.
A `.opt` file is converted once with `convert_optim(opt_file, bell_file)`.
[`load_structure`](@ref) is what builds the [`Lattice`](@ref).

**Another energy:** [`with_reference`](@ref)`(lat, reference_from_kinetic(deuteron, 10.0))`
keeps the elements and the normalised strengths. Maps are built with the new
reference.

**Edits** return a new lattice. [`insert_cavity`](@ref) splits the named drift
into drift, a thin [`RFCavity`](@ref), and drift. The ring harmonic stays
`lat.harmonic` unless the keyword `harmonic` is set.

```julia
drift = lat.elements[findfirst(e -> e isa Drift, lat.elements)]
lat_rf, index = insert_cavity(lat, RFCavity("RF", 0.2, 5, phi_s);
                              in_element = drift.name)
```

A thick cavity is `RFCavity(name, V, h, phi_s, L)` with `L > 0`. It is
integrated from its s-Hamiltonian; see
[RF Hamiltonian and acceleration](../dev/rf_acceleration.md).
`place_sextupole_families` and `set_family_strengths` add
sextupoles (`examples/10_sextupole_families/run.jl`).

**Coordinates.** Tracking uses the canonical state
``z = (x, p_x, y, p_y, l, \delta)``, with ``l = -v\Delta t``. Laboratory
coordinates are ``(x, x', y, y', \varphi, \delta)``. Positive ``l`` is an
early arrival: [`canonical_to_laboratory`](@ref) maps it to ``\varphi > 0``.

```julia
C = path_length(lat)
h = lat.harmonic
z = laboratory_to_canonical(u, lat.ref, C, h)
u = canonical_to_laboratory(z, lat.ref, C, h)
```

`examples/01_reference_and_coordinates/run.jl` prints both maps, including the
polynomial form on the coordinate jet.
