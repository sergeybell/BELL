# =============================================================================
# 07_blond_longitudinal.jl - optional BLonD-like exact drift vs element tracking
# =============================================================================
#
# WHAT
#   1. Extract alpha0, alpha1, alpha2 (and BLonD eta0/1/2) from our ring map
#      via blond_alphas (fit of C(d)/C0 on closed orbits).
#   2. BLonD Equations of Motion (https://blond-admin.github.io/BLonD/modules/
#      equations_of_motion.html) exact drift (eq. 20):
#        dt' = dt + T0 [(1 + a0 d + a1 d^2 + a2 d^3) * (1+dE/Es)/(1+d) - 1]
#      with exact delta from energy; simple drift (eq. 19 / solver=simple):
#        dt' = dt + T0 eta0 * dE/(beta^2 Es)
#   3. Compare one-turn Delta phi, Delta delta for magnetic lattice (no RF):
#      element-by-element exact integrator  vs  BLonD exact / simple drift
#      (particles on their dispersion orbit, several delta).
#   4. With RF: thin kick + BLonD exact drift vs full element tracking over
#      a few turns (synchrotron oscillation of a test particle).
#
# WHY
#   Shows that alpha_1 (and alpha_2) ARE the nonlinear longitudinal terms in
#   the BLonD sense, and that our element Hamiltonians reproduce the same
#   time-of-flight when alpha_i are taken from the lattice.
# =============================================================================
include(joinpath(@__DIR__, "..", "src", "LieRing.jl"))
using .LieRing, Printf, LinearAlgebra
include(joinpath(@__DIR__, "reference_values.jl"))
section(t, lines...) = (println("\n--- ", t, " ", "-"^max(0, 66 - length(t))); foreach(l -> println("  > ", l), lines))
println("="^72); println(" LieRing - BLonD-like longitudinal drift vs element tracking"); println("="^72)

of = parse_optim(OPTIM_FILE); lat = of.lattice; ref = lat.ref; C = circumference(lat); h = 5
@printf("  mass used = %.8f MeV (CODATA deuteron); OptiM header Mass = %.3f MeV (ignored)\n", ref.mass_MeV, of.mass_MeV)

section("1. Momentum compaction from the ring map (BLonD alpha_0,1,2)")
t0 = time(); M = ring_taylor_map(lat); @printf("  Taylor map in %.1f s\n", time() - t0)
pe = path_expansion(M, ref, C)
mc = blond_alphas(M, ref, C)
@printf("  path_expansion (order 2):  alpha0 = %.8f  alpha1 = %.8f\n", pe.alpha0, pe.alpha1)
@printf("  blond_alphas fit:   alpha0 = %.8f  alpha1 = %.8f  alpha2 = %.8f\n", mc.alpha0, mc.alpha1, mc.alpha2)
@printf("  BLonD eta (eq.15):         eta0 = %.8f  eta1 = %.8f  eta2 = %.8f\n", mc.eta0, mc.eta1, mc.eta2)
@printf("  path_expansion eta:        eta0 = %.8f  eta1 = %.8f  (TOF identity)\n", pe.eta0, pe.eta1)
@printf("  |alpha0 fit - path| = %.2e  |alpha1 fit - path| = %.2e\n", mc.alpha0 - pe.alpha0, mc.alpha1 - pe.alpha1)

section("2. One-turn magnetic drift: element tracking vs BLonD",
        "Particles on the closed orbit at fixed delta (x = D delta, ...). No RF.",
        "Report Delta phi = phi_out - phi_in and the BLonD prediction from dt.")
po = periodic_optics(ring_matrix(lat), ref, C)
@printf("  %10s %14s %14s %14s %12s %12s\n", "delta", "Dphi exact", "Dphi BLonD", "Dphi simple", "rel BLonD", "rel simple")
for d in (-2e-3, -1e-3, -5e-4, 5e-4, 1e-3, 2e-3)
    u = [po.D4[1]*d, po.D4[2]*d, 0.0, 0.0, 0.0, d]   # user; phi=0
    z = user_to_canonical(u, ref, C, h)
    z2 = exact_track_turn(z, lat, ref, C; ds = 0.2)
    u2 = canonical_to_user(z2, ref, C, h)
    dphi_ex = u2[5] - u[5]
    dt0, dE = phi_delta_to_dt_dE(0.0, d, ref, C, h)
    dt_b = blond_exact_drift(dt0, dE, ref, C, mc.alpha0, mc.alpha1, mc.alpha2)
    dt_s = blond_simple_drift(dt0, dE, ref, C, mc.eta0)
    dphi_b = dt_dE_to_phi_delta(dt_b, dE, ref, C, h)[1]
    dphi_s = dt_dE_to_phi_delta(dt_s, dE, ref, C, h)[1]
    @printf("  %10.4e %14.6e %14.6e %14.6e %12.2e %12.2e\n",
            d, dphi_ex, dphi_b, dphi_s, dphi_b/dphi_ex - 1, dphi_s/dphi_ex - 1)
end

section("3. With RF: thin kick + BLonD exact drift vs full tracking (32 turns)",
        "Stationary bucket phi_s = 180 deg, V = 0.2 MV, L_cav ignored in BLonD",
        "(instantaneous kick); our tracking uses the finite cavity.")
V, phis = 0.2, pi
cav = Cavity("RF", 1.0, V, h, phis)
lat_rf, _ = insert_cavity(lat, cav)
lat_rf = Lattice(lat_rf.name, lat_rf.elements, ref, h)
function compare_rf(lat_rf, ref, C, h, mc, d0, V, phis, po; turns = 32)
    u = [po.D4[1]*d0, po.D4[2]*d0, 0.0, 0.0, 0.0, d0]
    z = user_to_canonical(u, ref, C, h)
    dt, dE = phi_delta_to_dt_dE(0.0, d0, ref, C, h)
    phis_tr = Float64[]; dels_tr = Float64[]; phis_b = Float64[]; dels_b = Float64[]
    # Fixed energy, stationary bucket: do NOT renormalise (test particle != sync).
    for n in 1:turns
        for e in lat_rf.elements
            z = exact_track_element(z, e, ref, C; ds = 0.2)
        end
        uu = canonical_to_user(z, ref, C, h)
        push!(phis_tr, uu[5]); push!(dels_tr, uu[6])
        w0 = 2pi * revolution_frequency(ref, C)
        phi = -h * w0 * dt
        dE = dE + ref.particle.charge * V * sin(phis + phi) - ref.particle.charge * V * sin(phis)
        dt = blond_exact_drift(dt, dE, ref, C, mc.alpha0, mc.alpha1, mc.alpha2)
        pb, db = dt_dE_to_phi_delta(dt, dE, ref, C, h)
        push!(phis_b, pb); push!(dels_b, db)
    end
    return phis_tr, phis_b, dels_tr, dels_b
end
phis_tr, phis_b, dels_tr, dels_b = compare_rf(lat_rf, ref, C, h, mc, 1e-3, V, phis, po)
@printf("  %6s %12s %12s %12s %12s\n", "turn", "phi track", "phi BLonD", "d track", "d BLonD")
for n in (1, 4, 8, 16, 32)
    @printf("  %6d %12.5f %12.5f %12.5e %12.5e\n", n, phis_tr[n], phis_b[n], dels_tr[n], dels_b[n])
end
@printf("  max |phi_track - phi_BLonD| over 32 turns = %.3e rad\n", maximum(abs, phis_tr .- phis_b))
@printf("  max |delta_track - delta_BLonD|           = %.3e\n", maximum(abs, dels_tr .- dels_b))
println("  (small residual expected: finite cavity length + transverse coupling vs thin kick + 1D drift)")
println("\nDone.")
