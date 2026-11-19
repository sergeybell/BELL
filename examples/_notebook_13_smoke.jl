# Notebook smoke for 13_acceleration.
#
# 13_acceleration/run.jl is the full ramp: 10 MeV to 270 MeV, an ensemble,
# a flat top, and the figures. This script keeps the thin-cavity checks,
# a short fixed-energy tune, and accelerate over a few turns with two particles.
#
# The notebook generator sets BELL_EXAMPLES_DIR. From a shell whose working
# directory is examples/:
#
#   julia --project=. _notebook_13_smoke.jl

using Bell, Printf

examples_dir = get(ENV, "BELL_EXAMPLES_DIR", pwd())
include(joinpath(examples_dir, "common", "print_util.jl"))

banner("Notebook smoke of the acceleration example")
note("13_acceleration/run.jl ramps from 10 MeV to 270 MeV with an ensemble and writes the figures.")
note("This script checks the thin cavity, tracks a short fixed-energy sample, and calls accelerate for a few turns.")

lat270 = load_structure(joinpath(examples_dir, "lattices", "magnetic_2.bell"))
ref270 = lat270.ref
mass = ref270.mass_MeV
E_inj = 10.0
V, h, phis = 0.2, 5, deg2rad(150.0)
ref_inj = reference_from_kinetic(deuteron, E_inj)
cav = RFCavity("RF", V, h, phis)
idrift = something(findfirst(e -> e isa Drift && startswith(e.name, "OBM"), lat270.elements),
                   findfirst(e -> e isa Drift, lat270.elements))
lat_rf270, icav = insert_cavity(lat270, cav; in_element = lat270.elements[idrift].name)
C = path_length(lat_rf270)
lat_rf = Lattice(lat_rf270.name, lat_rf270.elements, ref_inj, h)


section("1. Injection reference and one thin kick")

println(describe_reference(ref_inj))
println()
@printf("  cavity in %s at index %d\n", lat270.elements[idrift].name, icav)
@printf("  V = %.3f MV    h = %d    phi_s = %.1f deg    C = %.6f m\n",
        V, h, rad2deg(phis), C)
f0 = revolution_frequency(ref_inj, C)
@printf("  injection revolution frequency = %.4f kHz    RF frequency = %.4f MHz\n",
        f0 / 1e3, h * f0 / 1e6)

z = rf_kick(zeros(6), cav, ref270, C)
dE = V * sin(phis)
dE_sync = sqrt((ref270.p0_MeV * (1 + z[6]))^2 + mass^2) - sqrt(ref270.p0_MeV^2 + mass^2)
println()
@printf("  synchronous energy gain at 270 MeV = %.6f MeV\n", dE_sync)
@printf("  Z V sin(phi_s)                    = %.6f MeV\n", dE)
err = symplectic_error(w -> rf_kick(w, cav, ref270, C), zeros(6))
@printf("  symplectic error of rf_kick       = %.3e\n", err)


section("2. Fixed energy at injection, 32 turns",
        "Phase is set to 180 deg. The full script uses 1024 turns at injection and at 270 MeV.")

function fixed_energy_sample(lattice, z0, ref0, path, nturns)
    z = copy(z0)
    ref = ref0
    dd = Float64[]
    for _ in 1:nturns
        for e in lattice.elements
            z = exact_track_element(z, e, ref, path; ds = 0.5)
            if e isa RFCavity
                wrapped, ref = renormalize([z], ref, 1)
                z = wrapped[1]
            end
        end
        push!(dd, z[6])
    end
    return dd
end

n_tune = 32
Lfix = Lattice("fix", set_cavity_phase(lat_rf, pi).elements, ref_inj, h)
po = periodic_optics(ring_matrix(Lfix), ref_inj, C)
d0 = 1e-4
z_tune = [po.D4[1] * d0, po.D4[2] * d0, 0.0, 0.0, 0.0, d0]
t0 = time()
dd = fixed_energy_sample(Lfix, z_tune, ref_inj, C, n_tune)
t_tune = time() - t0
eta0 = po.alpha_c - 1 / ref_inj.gamma^2
qs = synchrotron_tune(ref_inj, eta0, V, h, pi)
qm = measure_tune(dd)
@printf("  %d turns, ds = 0.5 m, in %.1f s\n", n_tune, t_tune)
@printf("  alpha_c = %.8f    eta0 = %.8f\n", po.alpha_c, eta0)
@printf("  Q_s from the sample = %.6f    smooth Q_s = %.6f\n", qm, qs)
note("A few dozen turns show the call. The smooth tune is the small-amplitude formula.")


section("3. Smooth bucket and a short accelerate",
        "alpha1 is left at 0 here. The full script reads alpha0 and alpha1 from the one-turn map.",
        "The target is 0.3 MeV above injection, with two particles and a flat top of 2 turns.")

p = LongitudinalParams(ref_inj, po.alpha_c, 0.0, V, h, phis)
bh = bucket_height(p)
fp = fixed_points(p)
Hs, inside = separatrix(p)
@printf("  half-height +%.4e / %.4e    unstable phi = %.4f rad\n", bh[1], bh[2], fp.unstable)
@printf("  long_H(0, 0) = %.6e    separatrix H = %.6e    inside(0, 0) = %s\n",
        long_H(p, 0.0, 0.0), Hs, inside(0.0, 0.0))

E_stop = E_inj + 0.3
users = [
    zeros(6),
    [po.D4[1] * 1e-3, po.D4[2] * 1e-3, 0.0, 0.0, 0.0, 1e-3],
]
t0 = time()
rec, flat = accelerate(lat_rf, ref_inj, users; Ekin_target = E_stop, flat_turns = 2, ds = 0.5)
t_acc = time() - t0
@printf("  accelerate towards %.1f MeV: %d recorded turns in %.1f s\n",
        E_stop, rec.turn[end], t_acc)
@printf("  flat top starts at turn %d    final E_kin = %.4f MeV\n", flat, rec.Ekin[end])
@printf("  particle 2 final delta = %.6e\n", rec.u[2][end][6])
note("plot_acceleration, gif_ramp, and gif_elements are written by 13_acceleration/run.jl.")

println()
println("Done.")
