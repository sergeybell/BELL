# rf_longitudinal/run.jl
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# Thin and thick RF on lattices/magnetic_2.bell, at fixed energy.
# The thin cavity is rf_kick, placed by insert_cavity. The thick cavity is
# integrated from its s-Hamiltonian. The smooth (phi, delta) model is the
# averaged comparison: synchrotron_tune, LongitudinalParams, long_H,
# fixed_points, separatrix, bucket_height.
#
#   julia --project=examples examples/rf_longitudinal/run.jl

using Bell, Printf

include(joinpath(@__DIR__, "..", "common", "print_util.jl"))

const V = 0.2
const HARM = 5
const PHIS = pi
const DS = 0.5
const N_THIN = 48
const N_THICK = 8
const STARTED = time()

function track_turns(lattice, z0, nturns; ds)
    z = copy(z0)
    phis_tr = Float64[]
    dels_tr = Float64[]
    ref = lattice.ref
    C = path_length(lattice)
    h = lattice.harmonic
    for _ in 1:nturns
        for e in lattice.elements
            z = exact_track_element(z, e, ref, C; ds = ds)
        end
        u = canonical_to_laboratory(z, ref, C, h)
        push!(phis_tr, u[5])
        push!(dels_tr, u[6])
    end
    return z, phis_tr, dels_tr
end

function with_thick_cavity(lattice, drift_name, cav, harmonic)
    i = findfirst(e -> e.name == drift_name && e isa Drift, lattice.elements)
    i === nothing && error("no drift named $drift_name")
    host = lattice.elements[i]
    host.L > cav.L || error("cavity length $(cav.L) m does not fit in $(host.name)")
    gap = host.L - cav.L
    elems = AbstractElement[]
    for (j, e) in enumerate(lattice.elements)
        if j == i
            push!(elems, Drift(host.name * "_u", gap / 2))
            push!(elems, cav)
            push!(elems, Drift(host.name * "_d", gap / 2))
        else
            push!(elems, e)
        end
    end
    return Lattice(lattice.name * "+thickRF", elems, lattice.ref, harmonic), i + 1
end

function print_turn_table(phis_tr, dels_tr, turns)
    @printf("  %6s %14s %14s\n", "turn", "phi [rad]", "delta")
    for n in turns
        @printf("  %6d %14.6e %14.6e\n", n, phis_tr[n], dels_tr[n])
    end
    return nothing
end


banner("Thin and thick RF, and the smooth longitudinal model")

lat = load_structure(joinpath(@__DIR__, "..", "lattices", "magnetic_2.bell"))
ref = lat.ref
C = path_length(lat)


section("1. Magnetic ring",
        "Compaction is read from the magnetic lattice, before a cavity is added.",
        "alpha and eta come from the one-turn Lie map. The linear matrix gives",
        "the same eta0 as eta_slip of periodic_optics.")

println(describe_reference(ref))
println()
@printf("  %s    elements = %d    C = %.6f m    file harmonic = %d\n",
        lat.name, length(lat.elements), C, lat.harmonic)

t_map = @elapsed begin
    magnets = magnetic_lattice(lat)
    M = factor_to_map(ring_lie_factor(magnets))
end
pe = path_expansion(M, ref, C)
po = periodic_optics(ring_matrix(lat), ref, C)
qs = synchrotron_tune(ref, pe.eta0, V, HARM, PHIS)

@printf("  one-turn Lie map in %.2f s\n", t_map)
println()
@printf("  alpha0 = %.8f    alpha1 = %.8f    gamma_t = %.6f\n", pe.alpha0, pe.alpha1, pe.gamma_t)
@printf("  eta0   = %.8f    eta1   = %.8f    eta1 from alphas = %.8f    diff %.1e\n",
        pe.eta0, pe.eta1, pe.eta1_formula, pe.eta1 - pe.eta1_formula)
@printf("  linear eta_slip = %.8f    alpha_c = %.8f    |eta0 - eta_slip| = %.1e\n",
        po.eta_slip, po.alpha_c, abs(pe.eta0 - po.eta_slip))
@printf("  V = %.3f MV    h = %d    phi_s = %.1f deg\n", V, HARM, rad2deg(PHIS))
@printf("  smooth synchrotron tune Q_s = %.6f\n", qs)
@printf("  stable when eta0 * cos(phi_s) > 0:  product = %.6e\n", pe.eta0 * cos(PHIS))


section("2. Thin cavity: insert_cavity and rf_kick",
        "insert_cavity splits a drift and leaves the path length unchanged.",
        "A thin cavity has L = 0. Its map is rf_kick, not a continuous Hamiltonian.",
        "laboratory_to_canonical takes (x, x', y, y', phi, delta).")

idrift = something(findfirst(e -> e isa Drift && startswith(e.name, "OBM"), lat.elements),
                   findfirst(e -> e isa Drift, lat.elements))
drift_name = lat.elements[idrift].name
thin = RFCavity("RF", V, HARM, PHIS)
lat_thin, icav = insert_cavity(lat, thin; in_element = drift_name, harmonic = HARM)

@printf("  cavity %s in the middle of %s (element %d), L = %.3f m\n",
        thin.name, drift_name, icav, thin.L)
@printf("  path length %.9f m -> %.9f m\n", C, path_length(lat_thin))

u_lab = [1.0e-3, 2.0e-4, -5.0e-4, -1.0e-4, 0.25, 8.0e-4]
z_can = laboratory_to_canonical(u_lab, ref, C, HARM)
u_back = canonical_to_laboratory(z_can, ref, C, HARM)
println()
@printf("  laboratory u     = %s\n", join([@sprintf("% .4e", v) for v in u_lab], " "))
@printf("  canonical z      = %s\n", join([@sprintf("% .4e", v) for v in z_can], " "))
@printf("  |round trip - u| = %.3e\n", maximum(abs, u_back - u_lab))

z_kick = rf_kick(z_can, thin, ref, C)
z_elem = exact_track_element(z_can, thin, ref, C)
t_sym = @elapsed err_kick = symplectic_error(w -> rf_kick(w, thin, ref, C), z_can)
println()
@printf("  |rf_kick - exact_track_element| = %.3e\n", maximum(abs, z_kick - z_elem))
@printf("  delta: %.6e -> %.6e\n", z_can[6], z_kick[6])
@printf("  symplectic error of rf_kick = %.3e    (%.3f s)\n", err_kick, t_sym)


section("3. Smooth model in (phi, delta)",
        "LongitudinalParams builds eta0, eta1, and the voltage coefficient from",
        "the magnetic alpha0, alpha1. long_H is that averaged Hamiltonian.",
        "It is not an element s-Hamiltonian and it does not track the ring.")

for line in split(chomp(longitudinal_hamiltonian_text()), '\n')
    println("  ", line)
end
p = LongitudinalParams(ref, pe.alpha0, pe.alpha1, V, HARM, PHIS)
fp = fixed_points(p)
Hs, inside = separatrix(p)
height = bucket_height(p)
@printf("  h = %d    eta0 = %.8f    eta1 = %.8f    A = %.6e    phi_s = %.4f rad\n",
        p.h, p.eta0, p.eta1, p.A, p.phi_s)
@printf("  fixed points: stable phi = %.4f rad, unstable phi = %.4f rad, delta2 = %.6f\n",
        fp.stable, fp.unstable, fp.delta2)
@printf("  long_H(0, 0) = %.6e    long_H(unstable, 0) = %.6e    separatrix H = %.6e\n",
        long_H(p, 0.0, 0.0), long_H(p, fp.unstable, 0.0), Hs)
@printf("  bucket half-height  +%.6e / %.6e\n", height[1], height[2])
@printf("  inside(0, 0) = %s    inside(0, 1.5 * half-height) = %s\n",
        inside(0.0, 0.0), inside(0.0, 1.5 * height[1]))

phis_acc = deg2rad(150.0)
p_acc = LongitudinalParams(ref, pe.alpha0, pe.alpha1, V, HARM, phis_acc)
qs_acc = synchrotron_tune(ref, pe.eta0, V, HARM, phis_acc)
bh_acc = bucket_height(p_acc)
println()
@printf("  accelerating phase %.1f deg: Q_s = %.6f    half-height +%.6e / %.6e\n",
        rad2deg(phis_acc), qs_acc, bh_acc[1], bh_acc[2])
note("No energy ramp is tracked. The numbers above are the smooth model at fixed energy.")


section("4. Stationary bucket, thin cavity, $(N_THIN) turns",
        "One particle on the linear dispersion orbit at delta = 1e-3, phi = 0.",
        "exact_track_element through the ring. The reference energy is not rebuilt.")

d0 = 1.0e-3
z0 = [po.D4[1] * d0, po.D4[2] * d0, 0.0, 0.0, 0.0, d0]
@printf("  canonical start = %s\n", join([@sprintf("% .4e", v) for v in z0], " "))
t_thin = @elapsed z_thin, phi_thin, del_thin = track_turns(lat_thin, z0, N_THIN; ds = DS)
println()
@printf("  %d turns, ds = %.2f m, in %.2f s\n", N_THIN, DS, t_thin)
print_turn_table(phi_thin, del_thin, (1, 8, 16, 32, N_THIN))
qm = measure_tune(del_thin)
println()
@printf("  tune from delta over %d turns = %.6f\n", N_THIN, qm)
@printf("  smooth Q_s                     = %.6f    (tracked / smooth - 1 = %.3e)\n",
        qs, qm / qs - 1)
note("The smooth tune is the small-amplitude formula. The tracked delta is finite.")


section("5. Thick cavity",
        "Same voltage and phase, L > 0. The drift is shortened so the path length",
        "stays the circumference. exact_track_element integrates dz/ds = S grad H.",
        "A few turns show the same stationary bucket; this is not an acceleration ramp.")

thick = RFCavity("RFthick", V, HARM, PHIS, 0.40)
lat_thick, ithick = with_thick_cavity(lat, drift_name, thick, HARM)
C_thick = path_length(lat_thick)
z_thick_1 = exact_track_element(z_can, thick, ref, C; ds = 0.05)
H_at = hamiltonian(thick, ref, C, z_can)

@printf("  %s    index %d\n", lat_thick.elements[ithick], ithick)
@printf("  path length %.9f m    (magnetic %.9f m)\n", C_thick, C)
@printf("  H(probe) = %.8f\n", H_at)
@printf("  one pass, probe particle:  thin delta %.6e -> %.6e    thick %.6e -> %.6e\n",
        z_can[6], z_kick[6], z_can[6], z_thick_1[6])
println("  thick pass minus thin kick (the thick map also drifts its length):")
for i in 1:6
    @printf("    %-5s % .3e\n", VAR_NAMES[i], z_thick_1[i] - z_kick[i])
end
t_sym_th = @elapsed err_thick = symplectic_error(
    w -> exact_track_element(w, thick, ref, C; ds = 0.1), z_can)
@printf("  symplectic error of one thick pass (ds = 0.1 m) = %.3e    (%.3f s)\n",
        err_thick, t_sym_th)

t_thick = @elapsed z_end, phi_th, del_th = track_turns(lat_thick, z0, N_THICK; ds = DS)
println()
@printf("  %d turns, ds = %.2f m, in %.2f s\n", N_THICK, DS, t_thick)
print_turn_table(phi_th, del_th, (1, 4, N_THICK))
@printf("\n  canonical end = %s\n", join([@sprintf("% .4e", v) for v in z_end], " "))

@printf("\n  wall time %.2f s\n", time() - STARTED)
println("Done.")
