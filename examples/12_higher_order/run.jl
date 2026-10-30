# 12_higher_order.jl
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# D2, alpha0/alpha1, eta0/eta1, xi1/xi2 on the natural lattice (sextupoles 0).
# Cross-check: exact-H tracking vs the one-turn map.
#
#   julia --project=examples examples/12_higher_order/run.jl

using Bell, Printf, LinearAlgebra, TaylorSeries

include(joinpath(@__DIR__, "..", "common", "print_util.jl"))
include(joinpath(@__DIR__, "..", "common", "reference_values.jl"))

banner("Second-order dispersion, compaction, slip, chromaticity")

of = convert_optim(OPTIM_FILE, STRUCTURE_FILE)
lat = load_structure(STRUCTURE_FILE)
ref = lat.ref
C = circumference(lat)

@printf("  lattice %s, C = %.6f m, E_kin = %.1f MeV, gamma0 = %.8f, beta0 = %.8f\n",
        lat.name, C, ref.Ekin_MeV, ref.gamma, ref.beta)


section("1. One-turn Taylor map (order 3)",
        "Element maps by RK4 in truncated polynomials, composed.")

t0 = time()
M = ring_taylor_map(lat)
@printf("  built in %.1f s\n", time() - t0)


section("2. Second-order dispersion",
        "(I-R4) D2 = Q4(D1,0,1): delta^2 part of the closed orbit.")

d = dispersion_orders(M)
@printf("  s=0:  D = % .8f m   D' = % .8e   D2 = % .8f m   D2' = % .8e\n",
        d.D1[1], d.D1[2], d.D2[1], d.D2[2])

al = dispersion_along(lat, d.D1, d.D2)
@printf("  along ring (element exits): max D = %.5f m, max|D2| = %.5f m, D2 range [%.5f, %.5f] m\n",
        maximum(al.D), maximum(abs, al.D2), minimum(al.D2), maximum(al.D2))
@printf("  closure: D2(C) - D2(0) = %.2e, D2'(C) - D2'(0) = %.2e\n",
        al.D2[end] - al.D2[1], al.D2p[end] - al.D2p[1])

# Same orbit pushed by RK4 element maps (not the definition of D2(s)).
function _d2_rk4(lat, D1, D2)
    dvar = Z[6]
    z = [[D1[i] * dvar + D2[i] * dvar^2 for i in 1:4]; zero(dvar); dvar]
    D2s = [D2[1]]
    for Me in element_maps(lat)
        z = [truncate_to(evaluate(m, z), 2) for m in Me]
        push!(D2s, coeff(z[1], [0, 0, 0, 0, 0, 2]))
    end
    return D2s
end

println()
@printf("  vs RK4 element maps: max |D2 Lie minus RK4| = %.2e\n",
        maximum(abs, al.D2 .- _d2_rk4(lat, d.D1, d.D2)))


section("3. Momentum compaction and slip factor to 2nd order",
        "dC/C = alpha0 d + alpha1 d^2 ;  dT/T0 = eta0 d + eta1 d^2.")

p = path_expansion(M, ref, C)
c1, c2, c3 = velocity_expansion(ref)

@printf("  l1 = (l|d)_closed = %.8f m,  l2 = (l|d^2)_closed = %.8f m\n", p.l1, p.l2)
@printf("  v/v0 - 1 = %.8f d + %.8f d^2  (check vs exact at d=1e-3: %.2e)\n",
        c1, c2, velocity_ratio(1e-3, ref) - 1 - (c1 * 1e-3 + c2 * 1e-6))

println()
@printf("  alpha0 = %.8f   alpha1 = %.8f   gamma_t = %.6f\n", p.alpha0, p.alpha1, p.gamma_t)
@printf("  eta0   = %.8f   eta1   = %.8f   (identity formula: %.8f, diff %.1e)\n",
        p.eta0, p.eta1, p.eta1_formula, p.eta1 - p.eta1_formula)
@printf("  delta_t (2nd fixed point of longitudinal motion, alpha-bucket) ~ -eta0/eta1 = %.4f\n",
        -p.eta0 / p.eta1)

obs = observables(lat)
println()
@printf("  observables (Lie factor): eta0 = %.8f  eta1 = %.8f  (Taylor map eta0 = %.8f, eta1 = %.8f)\n",
        obs.pe.eta0, obs.pe.eta1, p.eta0, p.eta1)
@printf("  observables Twiss from F.R: Qx = %.8f  Qy = %.8f  alpha_c = %.8f\n",
        obs.po.x.mu / 2pi, obs.po.y.mu / 2pi, obs.po.alpha_c)


section("4. Cross-check by exact-H tracking",
        "Closed orbit at fixed delta (Newton), dC/C from l, quadratic fit.")

Jmap(dl) = [evaluate(differentiate(M[i], j), [closed_orbit(M, dl); 0.0; dl]) for i in 1:4, j in 1:4]

function exact_orbit(dl)
    u = closed_orbit(M, dl)
    J = Jmap(dl) - I
    for _ in 1:8
        z = exact_track_turn([u; 0.0; dl], lat)
        du = -J \ (z[1:4] - u)
        u += du
        maximum(abs, du) < 1e-15 && break
    end
    return u
end

dls = [-1.5e-3, -1e-3, -0.5e-3, 0.0, 0.5e-3, 1e-3, 1.5e-3]
xs = Float64[]
dCs = Float64[]
for dl in dls
    u = exact_orbit(dl)
    z = exact_track_turn([u; 0.0; dl], lat)
    push!(xs, u[1])
    push!(dCs, -z[5] / C + velocity_ratio(dl, ref) - 1)
end

V = hcat(dls, dls .^ 2, dls .^ 3)
cx = V \ xs
cc = V \ dCs

@printf("  fit x_co = a d + b d^2 + c d^3:    a = %.8f (map D %.8f),  b = %.6f (map D2 %.6f)\n",
        cx[1], d.D1[1], cx[2], d.D2[1])
@printf("  fit dC/C = a d + b d^2 + c d^3:    a = %.8f (map alpha0 %.8f), b = %.6f (map alpha1 %.6f)\n",
        cc[1], p.alpha0, cc[2], p.alpha1)


section("5. Chromaticity, first and second order",
        "Map: Q(d) from the Jacobian on the off-momentum fixed point.",
        "Tracking: exact H, 256 turns of 1e-5 m oscillation, linear-tune fit.")

ch = chromaticity_orders(M)
@printf("  map:  Qx0 = %.8f  xi1_x = %.6f  xi2_x = d2Qx/dd2 = %.4f\n", ch.x...)
@printf("        Qy0 = %.8f  xi1_y = %.6f  xi2_y = d2Qy/dd2 = %.4f\n", ch.y...)

qxt = Float64[]
qyt = Float64[]
for dl in dls
    u = exact_orbit(dl)
    z = [u; 0.0; dl] .+ [1e-5, 0, 1e-5, 0, 0, 0]
    xx = Float64[]
    yy = Float64[]
    for n in 1:256
        z = exact_track_turn(z, lat)
        push!(xx, z[1] - u[1])
        push!(yy, z[3] - u[3])
    end
    push!(qxt, measure_tune(xx))
    push!(qyt, measure_tune(yy))
end

V2 = hcat(ones(length(dls)), dls, dls .^ 2, dls .^ 3)
fx = V2 \ qxt
fy = V2 \ qyt

println()
@printf("  tracking fit: Qx0 = %.8f  xi1_x = %.6f  xi2_x = %.4f\n", fx[1], fx[2], 2fx[3])
@printf("                Qy0 = %.8f  xi1_y = %.6f  xi2_y = %.4f\n", fy[1], fy[2], 2fy[3])

println()
@printf("  %-8s %13s %13s %13s %13s %13s\n",
        "delta", "Qx track", "Qx map", "Qx map-poly", "Qy track", "Qy map")
for (k, dl) in enumerate(dls)
    tm = tunes_at_delta(M, dl)
    @printf("  %-8.1e %13.9f %13.9f %13.9f %13.9f %13.9f\n",
            dl, qxt[k], tm[1], ch.x[1] + ch.x[2] * dl + ch.x[3] * dl^2 / 2, qyt[k], tm[2])
end

outdir = joinpath(@__DIR__, "output")
mkpath(outdir)
png = plot_chromatic(
    [(title = "natural lattice (sextupoles 0)",
      d = dls, qx_track = qxt, qy_track = qyt,
      qx_poly = x -> ch.x[1] + ch.x[2] * x + ch.x[3] * x^2 / 2,
      qy_poly = x -> ch.y[1] + ch.y[2] * x + ch.y[3] * x^2 / 2)],
    joinpath(outdir, "chromaticity_natural.png"))
println()
println("  written: ", abspath(png))


section("6. Longitudinal motion from eta0, eta1")

println(longitudinal_hamiltonian_text())
@printf("  example at 270 MeV, h=5, V=0.2 MV, T=1, phi_s=180 deg: Q_s = %.6f\n",
        synchrotron_tune(ref, p.eta0, 0.2, 5, pi))


section("7. Exact integrator vs Lie tracking vs linear map (one turn, no RF)")

z0 = [1e-3, 1e-4, 1e-3, -1e-4, 0.0, 1e-3]
ze = exact_track_turn(z0, lat)
zl = track(build_factors(lat), z0)
zm = [evaluate(m, z0) for m in M]
zr = ring_matrix(lat) * z0

@printf("  |exact - Lie(order 3)| = %.2e,  |exact - Taylor map| = %.2e,  |exact - linear| = %.2e\n",
        maximum(abs, ze - zl), maximum(abs, ze - zm), maximum(abs, ze - zr))

Je = jacobian_of_function(z -> exact_track_turn(z, lat), z0)
@printf("  symplectic error of exact-integrator one-turn Jacobian at z0: %.2e\n",
        symplectic_error(Je))


section("8. Notes",
        "COSY comparison items (from phase 1) remain open:",
        "tunes differ from COSY by ~0.015; natural chromaticity differs (COSY -1.84/-1.85);",
        "sextupole strengths in magnetic_2.opt are all zero (natural lattice above).")

@printf("  COSY natural chromaticity (reference only): xi_x = %.6f, xi_y = %.6f\n",
        COSY_REF.chrom_x[1], COSY_REF.chrom_y[1])

println()
println("Done.")
