# 16_blond_longitudinal.jl
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# BLonD-like (dt, dE) drift compared with element-by-element tracking.
# The BLonD formulae live here, not in the 6D ring path.
#
#   julia --project=examples examples/16_blond_longitudinal/run.jl

using Bell, Printf, LinearAlgebra, TaylorSeries

include(joinpath(@__DIR__, "..", "common", "print_util.jl"))
include(joinpath(@__DIR__, "..", "common", "reference_values.jl"))

function blond_alphas(M, ref::RefParticle, C::Real;
                      deltas = [-2e-3, -1e-3, -5e-4, 5e-4, 1e-3, 2e-3])
    ys = Float64[]
    xs0 = Float64[]
    xs1 = Float64[]
    xs2 = Float64[]
    for d in deltas
        u = closed_orbit(M, d)
        z = [u; 0.0; d]
        zout = [evaluate(M[i], z) for i in 1:6]
        Dl = zout[5] - z[5]
        vr = velocity_ratio(d, ref)
        Cd_over_C = vr - Dl / C
        push!(ys, Cd_over_C - 1)
        push!(xs0, d)
        push!(xs1, d^2)
        push!(xs2, d^3)
    end
    A = hcat(xs0, xs1, xs2)
    a = A \ ys
    alpha0, alpha1, alpha2 = a[1], a[2], a[3]
    g2 = ref.gamma^2
    b2 = ref.beta^2
    eta0 = alpha0 - 1 / g2
    eta1 = 1.5 * b2 / g2 + alpha1 - alpha0 * eta0
    eta2 = -b2 * (5b2 - 1) / (2g2) + alpha2 - 2 * alpha0 * alpha1 + alpha1 / g2 +
           alpha0^2 * eta0 - 1.5 * b2 * alpha0 / g2
    return (alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2,
            eta0 = eta0, eta1 = eta1, eta2 = eta2)
end

function blond_exact_drift(dt::Real, dE::Real, ref::RefParticle, C::Real,
                           alpha0, alpha1, alpha2; length_ratio = 1.0)
    Es = ref.gamma * ref.mass_MeV
    E = Es + dE
    p0 = ref.p0_MeV
    p = sqrt(max(E^2 - ref.mass_MeV^2, 0.0))
    delta = p / p0 - 1
    T0 = C / (ref.beta * C_LIGHT)
    path = 1 + alpha0 * delta + alpha1 * delta^2 + alpha2 * delta^3
    vel_ratio_inv = (1 + dE / Es) / (1 + delta)
    return dt + length_ratio * T0 * (path * vel_ratio_inv - 1)
end

function blond_simple_drift(dt::Real, dE::Real, ref::RefParticle, C::Real, eta0;
                            length_ratio = 1.0)
    Es = ref.gamma * ref.mass_MeV
    T0 = C / (ref.beta * C_LIGHT)
    delta = dE / (ref.beta^2 * Es)
    return dt + length_ratio * T0 * eta0 * delta
end

function phi_delta_to_dt_dE(phi, delta, ref::RefParticle, C::Real, h::Real)
    w0 = 2pi * revolution_frequency(ref, C)
    dt = -phi / (h * w0)
    Es = ref.gamma * ref.mass_MeV
    p = ref.p0_MeV * (1 + delta)
    E = sqrt(p^2 + ref.mass_MeV^2)
    return dt, E - Es
end

function dt_dE_to_phi_delta(dt, dE, ref::RefParticle, C::Real, h::Real)
    w0 = 2pi * revolution_frequency(ref, C)
    phi = -h * w0 * dt
    Es = ref.gamma * ref.mass_MeV
    E = Es + dE
    p = sqrt(max(E^2 - ref.mass_MeV^2, 0.0))
    return phi, p / ref.p0_MeV - 1
end

banner("BLonD-like longitudinal drift vs element tracking")

of = convert_optim(OPTIM_FILE, STRUCTURE_FILE)
lat = load_structure(STRUCTURE_FILE)
ref = lat.ref
C = circumference(lat)
h = 5

@printf("  mass used = %.8f MeV (CODATA deuteron); OptiM header Mass = %.3f MeV (ignored)\n",
        ref.mass_MeV, of.mass_MeV)


section("1. Momentum compaction from the ring map (BLonD alpha_0,1,2)")

t0 = time()
M = ring_taylor_map(lat)
@printf("  Taylor map in %.1f s\n", time() - t0)

pe = path_expansion(M, ref, C)
mc = blond_alphas(M, ref, C)

println()
@printf("  path_expansion (order 2):  alpha0 = %.8f  alpha1 = %.8f\n", pe.alpha0, pe.alpha1)
@printf("  blond_alphas fit:   alpha0 = %.8f  alpha1 = %.8f  alpha2 = %.8f\n",
        mc.alpha0, mc.alpha1, mc.alpha2)
@printf("  BLonD eta (eq.15):         eta0 = %.8f  eta1 = %.8f  eta2 = %.8f\n",
        mc.eta0, mc.eta1, mc.eta2)
@printf("  path_expansion eta:        eta0 = %.8f  eta1 = %.8f  (TOF identity)\n",
        pe.eta0, pe.eta1)
@printf("  |alpha0 fit - path| = %.2e  |alpha1 fit - path| = %.2e\n",
        mc.alpha0 - pe.alpha0, mc.alpha1 - pe.alpha1)


section("2. One-turn magnetic drift: element tracking vs BLonD",
        "Particles on the closed orbit at fixed delta (x = D delta, ...). No RF.")

po = periodic_optics(ring_matrix(lat), ref, C)

@printf("  %10s %14s %14s %14s %12s %12s\n",
        "delta", "Dphi exact", "Dphi BLonD", "Dphi simple", "rel BLonD", "rel simple")
for d in (-2e-3, -1e-3, -5e-4, 5e-4, 1e-3, 2e-3)
    u = [po.D4[1]*d, po.D4[2]*d, 0.0, 0.0, 0.0, d]
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
        "Stationary bucket phi_s = 180 deg, V = 0.2 MV.",
        "BLonD uses an instantaneous kick; tracking uses the finite cavity.")

V, phis = 0.2, pi
cav = Cavity("RF", 1.0, V, h, phis)
lat_rf, _ = insert_cavity(lat, cav)
lat_rf = Lattice(lat_rf.name, lat_rf.elements, ref, h)

function compare_rf(lat_rf, ref, C, h, mc, d0, V, phis, po; turns = 32)
    u = [po.D4[1]*d0, po.D4[2]*d0, 0.0, 0.0, 0.0, d0]
    z = user_to_canonical(u, ref, C, h)
    dt, dE = phi_delta_to_dt_dE(0.0, d0, ref, C, h)
    phis_tr = Float64[]
    dels_tr = Float64[]
    phis_b = Float64[]
    dels_b = Float64[]
    # Fixed energy, stationary bucket: do NOT renormalise.
    for n in 1:turns
        for e in lat_rf.elements
            z = exact_track_element(z, e, ref, C; ds = 0.2)
        end
        uu = canonical_to_user(z, ref, C, h)
        push!(phis_tr, uu[5])
        push!(dels_tr, uu[6])
        w0 = 2pi * revolution_frequency(ref, C)
        phi = -h * w0 * dt
        dE = dE + ref.particle.charge * V * sin(phis + phi) - ref.particle.charge * V * sin(phis)
        dt = blond_exact_drift(dt, dE, ref, C, mc.alpha0, mc.alpha1, mc.alpha2)
        pb, db = dt_dE_to_phi_delta(dt, dE, ref, C, h)
        push!(phis_b, pb)
        push!(dels_b, db)
    end
    return phis_tr, phis_b, dels_tr, dels_b
end

phis_tr, phis_b, dels_tr, dels_b = compare_rf(lat_rf, ref, C, h, mc, 1e-3, V, phis, po)

@printf("  %6s %12s %12s %12s %12s\n", "turn", "phi track", "phi BLonD", "d track", "d BLonD")
for n in (1, 4, 8, 16, 32)
    @printf("  %6d %12.5f %12.5f %12.5e %12.5e\n",
            n, phis_tr[n], phis_b[n], dels_tr[n], dels_b[n])
end

println()
@printf("  max |phi_track - phi_BLonD| over 32 turns = %.3e rad\n",
        maximum(abs, phis_tr .- phis_b))
@printf("  max |delta_track - delta_BLonD|           = %.3e\n",
        maximum(abs, dels_tr .- dels_b))
println("  (small residual expected: finite cavity length + transverse coupling vs thin kick + 1D drift)")

println()
println("Done.")
