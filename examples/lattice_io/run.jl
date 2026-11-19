# lattice_io.jl
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# Load a structure file, write it back, and edit the Lattice:
# path length, bend angle, unique elements, a thin cavity, the magnetic
# lattice, and a new reference particle. Lattices: fodo.bell, magnetic_2.bell.
#
#   julia --project=examples examples/lattice_io/run.jl

using Bell, Printf

include(joinpath(@__DIR__, "..", "common", "print_util.jl"))

banner("Lattice IO: load, write, cavity, reference (FODO and magnetic_2)")

fodo_file = joinpath(@__DIR__, "..", "lattices", "fodo.bell")
mag_file = joinpath(@__DIR__, "..", "lattices", "magnetic_2.bell")
fodo = load_structure(fodo_file)
mag = load_structure(mag_file)


section("1. Two lattices from examples/lattices")

function summarize(lat, file)
    println("  file : ", abspath(file))
    @printf("  name : %s    harmonic : %d    elements : %d    distinct : %d\n",
            lat.name, lat.harmonic, length(lat.elements), length(unique_elements(lat)))
    @printf("  path length         : %.7f m\n", path_length(lat))
    @printf("  total bending angle : %.9f rad   (%.6f deg)\n",
            total_bend_angle(lat), rad2deg(total_bend_angle(lat)))
    @printf("  reference           : %s, Ekin = %.4f MeV, beta = %.6f\n",
            lat.ref.particle.name, lat.ref.Ekin_MeV, lat.ref.beta)
    println()
end

summarize(fodo, fodo_file)
summarize(mag, mag_file)
@printf("  FODO bending angle / (2 pi) = %.6f\n", total_bend_angle(fodo) / (2pi))
@printf("  magnetic_2 bending angle / (2 pi) = %.6f\n", total_bend_angle(mag) / (2pi))


section("2. Distinct elements and kind_string")

for (title, lat) in (("FODO", fodo), ("magnetic_2", mag))
    println("  ", title)
    for e in unique_elements(lat)
        n = count(isequal(e), lat.elements)
        @printf("    %-8s x%-4d %-12s  %s\n", e.name, n, kind_string(e), string(e))
    end
    println()
end


section("3. write_structure and load the file back",
        "load_structure stores the evaluated numbers. The written file has no parameters block.",
        "Loading that file builds a Lattice whose name is the new file name.")

outdir = joinpath(@__DIR__, "output")
mkpath(outdir)
written = joinpath(outdir, "fodo_roundtrip.bell")
write_structure(written, fodo; extra_header = "# numeric round trip of lattices/fodo.bell")
back = load_structure(written)
text = read(written, String)

@printf("  written : %s\n", abspath(written))
@printf("  parameters block in the written file : %s\n",
        occursin(r"(?m)^parameters\s*$", text) ? "yes" : "no")
@printf("  name : %s -> %s\n", fodo.name, back.name)
@printf("  elements : %d -> %d\n", length(fodo.elements), length(back.elements))
@printf("  path length difference      : %.3e m\n", path_length(back) - path_length(fodo))
@printf("  bend angle difference       : %.3e rad\n", total_bend_angle(back) - total_bend_angle(fodo))
@printf("  |delta p0|                  : %.3e MeV\n", abs(back.ref.p0_MeV - fodo.ref.p0_MeV))
n_mismatch = count(p -> !isequal(p[1], p[2]), zip(fodo.elements, back.elements))
@printf("  elements that differ        : %d\n", n_mismatch)
n_mismatch == 0 || error("write_structure / load_structure changed an element")
back.harmonic == fodo.harmonic || error("harmonic was not preserved")


section("4. insert_cavity, then magnetic_lattice",
        "The cavity harmonic is the cavity's own h. The ring harmonic stays unless harmonic= is set.",
        "magnetic_lattice drops a thin cavity. The two drift halves stay, so the path length is unchanged.")

drift = fodo.elements[findfirst(e -> e isa Drift, fodo.elements)]
cav = RFCavity("RF", 0.05, 5, 0.0)
lat_rf, icav = insert_cavity(fodo, cav; in_element = drift.name)
placed = lat_rf.elements[icav]

@printf("  split drift %s (L = %.4f m) at index %d\n", drift.name, drift.L, icav)
println("  ", lat_rf.elements[icav - 1])
println("  ", placed)
println("  ", lat_rf.elements[icav + 1])
println()
@printf("  cavity harmonic = %d    ring harmonic = %d (was %d)\n",
        placed.harmonic, lat_rf.harmonic, fodo.harmonic)
@printf("  elements : %d -> %d\n", length(fodo.elements), length(lat_rf.elements))
@printf("  path length difference : %.3e m\n", path_length(lat_rf) - path_length(fodo))
@printf("  kind_string(cavity) = %s    L = %.3f m\n", kind_string(placed), placed.L)

lat_mag = magnetic_lattice(lat_rf)
@printf("  after magnetic_lattice: elements = %d, RF cavities = %d, path length difference = %.3e m\n",
        length(lat_mag.elements),
        count(e -> e isa RFCavity, lat_mag.elements),
        path_length(lat_mag) - path_length(lat_rf))
@printf("  bend angle unchanged: %.3e rad\n", total_bend_angle(lat_mag) - total_bend_angle(fodo))


section("5. with_reference",
        "Normalised strengths stay on the elements. Brho changes with the species,",
        "so the physical gradient k1*Brho changes. The drift R56 changes with gamma.")

ref_p = reference_from_kinetic(proton, fodo.ref.Ekin_MeV)
lat_p = with_reference(fodo, ref_p)
iq = findfirst(e -> e isa Quadrupole, fodo.elements)
qf = fodo.elements[iq]
id = findfirst(e -> e isa Drift, fodo.elements)
drift0 = fodo.elements[id]

@printf("  %-12s %12s %14s %14s %14s\n",
        "reference", "beta", "Brho T m", "k1 1/m^2", "k1*Brho T/m")
for (label, lat) in (("deuteron", fodo), ("proton", lat_p))
    @printf("  %-12s %12.6f %14.6f %14.6f %14.6f\n",
            label, lat.ref.beta, lat.ref.brho_Tm, qf.k1, qf.k1 * lat.ref.brho_Tm)
end
println()
Rd = linear_matrix(drift0, fodo.ref)
Rp = linear_matrix(drift0, lat_p.ref)
@printf("  drift %s R56 (deuteron) = %.8f m\n", drift0.name, Rd[5, 6])
@printf("  drift %s R56 (proton)   = %.8f m\n", drift0.name, Rp[5, 6])
@printf("  same element object     : %s\n", lat_p.elements[iq] === qf ? "yes" : "no")
@printf("  quadrupole k1 unchanged : %.6f 1/m^2\n", lat_p.elements[iq].k1)

println()
println("Done.")
