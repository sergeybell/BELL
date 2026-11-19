# 16_blond_longitudinal/run.jl
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# Thin RF cavity on magnetic_2: path expansion (alpha, eta) and element-by-element
# tracking in a stationary bucket. The folder name is unchanged.
#
#   julia --project=examples examples/16_blond_longitudinal/run.jl

using Bell, Printf, LinearAlgebra, TaylorSeries

include(joinpath(@__DIR__, "..", "common", "print_util.jl"))

banner("Thin RF cavity and longitudinal motion")

lat = load_structure(joinpath(@__DIR__, "..", "lattices", "magnetic_2.bell"))
ref = lat.ref
C = path_length(lat)
h = 5
V = 0.2
phis = pi

@printf("  lattice %s, C = %.6f m, E_kin = %.1f MeV\n", lat.name, C, ref.Ekin_MeV)
@printf("  thin cavity: V = %.3f MV, h = %d, phi_s = %.1f deg\n", V, h, rad2deg(phis))


section("1. Momentum compaction and slip from the one-turn map",
        "Magnets only. dC/C = alpha0 d + alpha1 d^2; dT/T0 = eta0 d + eta1 d^2.")

t0 = time()
M = ring_taylor_map(lat)
@printf("  Taylor map in %.1f s\n", time() - t0)

pe = path_expansion(M, ref, C)
qs = synchrotron_tune(ref, pe.eta0, V, h, phis)

println()
@printf("  alpha0 = %.8f   alpha1 = %.8f   gamma_t = %.6f\n", pe.alpha0, pe.alpha1, pe.gamma_t)
@printf("  eta0   = %.8f   eta1   = %.8f   (TOF identity: %.8f, diff %.1e)\n",
        pe.eta0, pe.eta1, pe.eta1_formula, pe.eta1 - pe.eta1_formula)
@printf("  smooth synchrotron tune, stationary bucket: Q_s = %.6f\n", qs)


section("2. Stationary bucket: thin kick plus magnetic tracking",
        "Particle on the linear dispersion orbit at delta = 1e-3, phi = 0.",
        "Fixed energy: the synchronous phase is 180 deg, and the reference is not rebuilt.")

idrift = something(findfirst(e -> e isa Drift && startswith(e.name, "OBM"), lat.elements),
                   findfirst(e -> e isa Drift, lat.elements))
cav = RFCavity("RF", V, h, phis)
lat_rf, icav = insert_cavity(lat, cav; in_element = lat.elements[idrift].name, harmonic = h)

@printf("  cavity in the middle of drift %s (element %d)\n", lat.elements[idrift].name, icav)
@printf("  path length unchanged: %.6f m\n", path_length(lat_rf))

po = periodic_optics(ring_matrix(lat), ref, C)
d0 = 1e-3
u0 = [po.D4[1] * d0, po.D4[2] * d0, 0.0, 0.0, 0.0, d0]

function track_stationary(lat_rf, ref, C, h, u0, nturns)
    z = laboratory_to_canonical(u0, ref, C, h)
    phis_tr = Float64[]
    dels_tr = Float64[]
    for _ in 1:nturns
        for e in lat_rf.elements
            z = exact_track_element(z, e, ref, C; ds = 0.2)
        end
        uu = canonical_to_laboratory(z, ref, C, h)
        push!(phis_tr, uu[5])
        push!(dels_tr, uu[6])
    end
    return phis_tr, dels_tr
end

const N_TURNS = 128
t0 = time()
phis_tr, dels_tr = track_stationary(lat_rf, ref, C, h, u0, N_TURNS)

println()
@printf("  %d turns in %.1f s\n", N_TURNS, time() - t0)
@printf("  %6s %14s %14s\n", "turn", "phi [rad]", "delta")
for n in (1, 8, 16, 32, 64, 128)
    @printf("  %6d %14.6e %14.6e\n", n, phis_tr[n], dels_tr[n])
end

qm = measure_tune(dels_tr)
println()
@printf("  tune from delta over %d turns = %.6f\n", N_TURNS, qm)
@printf("  smooth Q_s                     = %.6f   (rel. diff %.3e)\n", qs, qm / qs - 1)
println("  The smooth tune is the small-amplitude formula; the tracked particle has a finite delta.")

println()
println("Done.")
