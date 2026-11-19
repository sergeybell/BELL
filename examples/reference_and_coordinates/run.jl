# reference_and_coordinates.jl
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# Reference particle, energy_term, and the exact change of coordinates
# between the laboratory 6-vector and the canonical 6-vector, on
# lattices/fodo.bell. The one-turn Taylor map is rewritten with laboratory_map.
#
#   julia --project=examples examples/reference_and_coordinates/run.jl

using Bell, Printf

include(joinpath(@__DIR__, "..", "common", "print_util.jl"))

banner("Reference particle and laboratory / canonical coordinates (FODO)")

lat = load_structure(joinpath(@__DIR__, "..", "lattices", "fodo.bell"))
ref = lat.ref
C = path_length(lat)
h = lat.harmonic


section("1. Reference particle stored on the lattice",
        "load_structure builds RefParticle from the reference block.",
        "reference_from_kinetic(deuteron, Ekin) rebuilds the same particle.")

println(describe_reference(ref))
println()
built = reference_from_kinetic(deuteron, ref.Ekin_MeV)
@printf("  rebuilt deuteron:  beta = %.10f   p0*c = %.6f MeV   Brho = %.8f T m\n",
        built.beta, built.p0_MeV, built.brho_Tm)
@printf("  difference vs file: |dbeta| = %.3e   |dp0| = %.3e   |dBrho| = %.3e\n",
        abs(built.beta - ref.beta), abs(built.p0_MeV - ref.p0_MeV), abs(built.brho_Tm - ref.brho_Tm))
@printf("  path length C = %.6f m    harmonic h = %d    elements = %d\n",
        C, h, length(lat.elements))


section("2. Same kinetic energy, three species",
        "gamma, beta, p0 and Brho follow the species. Electron charge is negative, so Brho is negative.")

@printf("  %-10s %10s %12s %12s %14s %14s\n",
        "species", "Ekin MeV", "gamma", "beta", "p0*c MeV", "Brho T m")
for (name, species) in (("deuteron", deuteron), ("proton", proton), ("electron", electron))
    r = reference_from_kinetic(species, ref.Ekin_MeV)
    @printf("  %-10s %10.1f %12.6f %12.8f %14.4f %14.6f\n",
            name, r.Ekin_MeV, r.gamma, r.beta, r.p0_MeV, r.brho_Tm)
end


section("3. energy_term and velocity_ratio",
        "energy_term(delta) is the longitudinal piece of the s-Hamiltonian.",
        "On the reference orbit it is 0.  d(energy_term)/d(delta) = v(delta)/v0.")

@printf("  energy_term(0)     = %.3e\n", energy_term(0.0, ref))
@printf("  velocity_ratio(0)  = %.16f\n", velocity_ratio(0.0, ref))
println()
@printf("  %12s %16s %16s %16s %12s\n",
        "delta", "energy_term", "v/v0", "dE/ddelta", "|diff|")
for delta in (0.0, 1e-3, 1e-2, 0.1)
    eps = 1e-7
    slope = (energy_term(delta + eps, ref) - energy_term(delta - eps, ref)) / (2eps)
    ratio = velocity_ratio(delta, ref)
    @printf("  %12.4e %16.8e %16.8f %16.8f %12.2e\n",
            delta, energy_term(delta, ref), ratio, slope, abs(slope - ratio))
end


section("4. Coordinate names and indices",
        "Canonical z = (x, px, y, py, l, delta).  Laboratory u = (x, x', y, y', phi, delta).",
        "l = -v dt.  Positive l is an early arrival, and canonical_to_laboratory maps it to phi > 0.")

@printf("  %5s  %-12s  %-12s  %s\n", "index", "canonical", "laboratory", "constant")
rows = (
    (index_x, "index_x"),
    (index_px, "index_px"),
    (index_y, "index_y"),
    (index_py, "index_py"),
    (index_l, "index_l"),
    (index_delta, "index_delta"),
)
for (i, constant) in rows
    @printf("  %5d  %-12s  %-12s  %s\n", i, VAR_NAMES[i], laboratory_names[i], constant)
end


section("5. Point conversion, round trip",
        "u = (5 mm, 20 mrad, -3 mm, -10 mrad, 0.4 rad, delta = 1e-2).",
        "px = x' (1+delta) / sqrt(1 + x'^2 + y'^2),   l = phi C (v/v0) / (2 pi h).")

u = [5e-3, 0.02, -3e-3, -0.01, 0.4, 1e-2]
z = laboratory_to_canonical(u, ref, C, h)
back = canonical_to_laboratory(z, ref, C, h)

@printf("  %-18s %16s %16s %16s\n", "component", "laboratory", "canonical", "round trip")
for i in 1:6
    @printf("  %-18s %16.8e %16.8e %16.8e\n",
            laboratory_names[i] * "  " * VAR_NAMES[i], u[i], z[i], back[i])
end
println()
@printf("  max |u - canonical_to_laboratory(laboratory_to_canonical(u))| = %.3e\n",
        maximum(abs, back - u))


section("6. One-turn map in both coordinate systems",
        "M is the canonical Taylor map of the 16-cell line.",
        "laboratory_map rewrites it as U = T o M o Tinv.",
        "The laboratory Jacobian is not symplectic: (x, x', y, y', phi, delta) is not a canonical set.")

M = ring_taylor_map(lat)
U = laboratory_map(M, ref, C, h)

pairs = (
    ("x | px   vs  x | x'", 1, [0, 1, 0, 0, 0, 0]),
    ("x | delta", 1, [0, 0, 0, 0, 0, 1]),
    ("px | delta  vs  x' | delta", 2, [0, 0, 0, 0, 0, 1]),
    ("l | delta  vs  phi | delta", 5, [0, 0, 0, 0, 0, 1]),
    ("l | delta^2  vs  phi | delta^2", 5, [0, 0, 0, 0, 0, 2]),
    ("px | px delta  vs  x' | x' delta", 2, [0, 1, 0, 0, 0, 1]),
)
@printf("  %-36s %16s %16s\n", "coefficient", "canonical", "laboratory")
for (label, i, exponents) in pairs
    @printf("  %-36s %16.8e %16.8e\n", label, coeff(M[i], exponents), coeff(U[i], exponents))
end

println()
z0 = [1e-3, 1e-3, 1e-3, 1e-3, 0.0, 1e-3]
@printf("  symplectic error of canonical Jacobian at z0   : %.3e\n",
        symplectic_error(jacobian_of_map(M, z0)))
@printf("  symplectic error of laboratory Jacobian at z0  : %.3e\n",
        symplectic_error(jacobian_of_map(U, z0)))
note("The canonical Taylor map is symplectic up to truncation order $MAP_ORDER. The laboratory set is not canonical.")

println()
println("Done.")
