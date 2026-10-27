# =============================================================================
# run_v07.jl - four-way comparison of delta_eq + residual study + map vs Lie
# =============================================================================
using EquilibriumEnergy
using LieRing, Printf, Plots, Statistics, LinearAlgebra
gr()

const OPTIM = raw"C:\GIT_REPOS\SCT\OptiM\magnetic\magnetic_2.opt"
const OUT = joinpath(@__DIR__, "..", "output")
const NOTE = joinpath(@__DIR__, "..", "note")
mkpath(OUT); mkpath(joinpath(NOTE, "figs"))

section(t) = println("\n", "="^72, "\n ", t, "\n", "="^72)

section("0. Setup + normal form (thin RF effective V*T)")
lat, ref, C, h, V, Ttf, pe, xi, po, info = setup_ring(OPTIM)
of = parse_optim(OPTIM)
lat_mag = Lattice(of.lattice.name, of.lattice.elements, ref, h)
nf = normal_form_rf(lat_mag, ref, C, h, V, Float64(pi); T = Ttf)
@printf("  C=%.6f  eta0=%.6f eta1=%.6f  xi=(%.4f,%.4f)  Qs=%.6f\n",
        C, pe.eta0, pe.eta1, xi.x[2], xi.y[2], info.Qs)
@printf("  NF: delta*(thin)=%.3e  beta_s=%.4e\n", nf.zstar[6], nf.beta_s)
@printf("  NF coeffs cx,cy,cs = %.6e, %.6e, %.6e\n", nf.c_x, nf.c_y, nf.c_s)
@printf("  note             = %.6e, %.6e, %.6e\n", nf.c_x_note, nf.c_y_note, nf.c_s_note)

# ---- Residual: finite-cavity fixed point ----
section("1. Residual origin: 6D fixed point of finite cavity")
function exact_fp(lat, ref, C; z0 = zeros(6), niter = 30)
    z = Float64.(z0); Id = Matrix{Float64}(I, 6, 6)
    for _ in 1:niter
        f0 = exact_track_turn(z, lat, ref, C) .- z
        J = zeros(6, 6); eps = 1e-8
        for j in 1:6
            zp = copy(z); zp[j] += eps; zm = copy(z); zm[j] -= eps
            J[:, j] = (exact_track_turn(zp, lat, ref, C) .- exact_track_turn(zm, lat, ref, C)) / (2eps) .- Id[:, j]
        end
        dz = -J \ f0; z .+= dz
        maximum(abs, dz) < 1e-14 && break
    end
    return z
end
zstar_f = exact_fp(lat, ref, C)
@printf("  finite cavity L=1 V=0.2: delta*=%+.6e  |z_perp|=%.3e\n",
        zstar_f[6], hypot(zstar_f[1], zstar_f[2], zstar_f[3], zstar_f[4]))

rows_res = []
for Vv in (0.05, 0.1, 0.2, 0.4)
    cav = Cavity("RF", 1.0, Vv, h, Float64(pi))
    latv, _ = insert_cavity(lat_mag, cav)
    latv = Lattice(latv.name, latv.elements, ref, h)
    zs = exact_fp(latv, ref, C)
    push!(rows_res, (V = Vv, Lc = 1.0, dstar = zs[6]))
    @printf("  V=%.2f  delta*=%+.6e\n", Vv, zs[6])
end
for Lc in (0.1, 0.5, 1.0, 2.0)
    cav = Cavity("RF", Lc, 0.2, h, Float64(pi))
    try
        latv, _ = insert_cavity(lat_mag, cav)
        latv = Lattice(latv.name, latv.elements, ref, h)
        zs = exact_fp(latv, ref, circumference(latv))
        push!(rows_res, (V = 0.2, Lc = Lc, dstar = zs[6]))
        @printf("  Lc=%.2f delta*=%+.6e\n", Lc, zs[6])
    catch e
        @printf("  Lc=%.2f FAIL %s\n", Lc, e)
    end
end

# ---- A: four methods ----
section("2. Test A: delta_eq three+map ways (zero betatron)")
dms = [0.5e-3, 0.75e-3, 1e-3, 1.5e-3, 2e-3]
turns_A = max(800, round(Int, 40 / info.Qs))
rows_A = []
# Taylor map with thin RF for (d)
Mthin = one_turn_map_rf(lat_mag, RFCavity("t", V * Ttf, h, Float64(pi)), ref, C)
for dm0 in dms
    u0 = [po.D4[1] * dm0, po.D4[2] * dm0, 0.0, 0.0, 0.0, dm0]
    # (a) note
    # measure dm from tracking for fair compare
    dser, pser = track_longitudinal(lat, ref, C, h, u0; turns = turns_A)
    m = measure_synchrotron(dser, pser; Qs = info.Qs)
    dmf = m.delta_m_fit
    a_note = formula_A(pe.eta0, pe.eta1, dmf)
    a_note_fp = a_note + zstar_f[6]   # include finite-cavity closed orbit
    # (b) NF (thin) + finite-cavity FP offset
    b_nf = delta_eq_nf(nf; delta_m = dmf) + zstar_f[6]
    # (c) tracking
    c_tr = m.delta_eq
    # (d) iterate truncated Taylor map (thin RF, start on dispersion orbit)
    z = user_to_canonical(u0, ref, C, h)
    N = turns_A
    ds_t = zeros(N)
    for n in 1:N
        z = evaluate_map(Mthin, z)
        ds_t[n] = z[6]
    end
    md = measure_synchrotron(ds_t, fill(0.0, N); Qs = info.Qs)  # phis unused for deq fit
    d_tay = md.delta_eq
    push!(rows_A, (dm = dmf, a = a_note, a_fp = a_note_fp, b = b_nf, c = c_tr, d = d_tay,
                   err = m.err_split, phat = m.phi_hat))
    @printf("  dm=%.4e  note=%+.4e note+FP=%+.4e  NF+FP=%+.4e  track=%+.4e  Taylor=%+.4e  err=%.2e\n",
            dmf, a_note, a_note_fp, b_nf, c_tr, d_tay, m.err_split)
end

# Plot A: delta_eq vs dm^2
pA = plot(xlabel = "δ_m²", ylabel = "δ_eq", legend = :bottomright,
          title = "δ_eq vs δ_m² (longitudinal only)", size = (700, 480))
scatter!(pA, [r.dm^2 for r in rows_A], [r.c for r in rows_A], label = "(c) exact tracking", ms = 6)
scatter!(pA, [r.dm^2 for r in rows_A], [r.a_fp for r in rows_A], label = "(a) note + δ*", ms = 5, markershape = :diamond)
plot!(pA, [r.dm^2 for r in rows_A], [r.b for r in rows_A], label = "(b) NF map ⟨Δφ⟩=0 + δ*", lw = 2)
plot!(pA, [r.dm^2 for r in rows_A], [r.d for r in rows_A], label = "(d) truncated Taylor iter.", lw = 1, ls = :dash)
slope = -pe.eta1 / (2 * pe.eta0)
plot!(pA, [0.0, maximum(r.dm^2 for r in rows_A)],
      [zstar_f[6], zstar_f[6] + slope * maximum(r.dm^2 for r in rows_A)],
      label = "slope −η₁/(2η₀) + δ*", ls = :dot, lw = 2)
savefig(pA, joinpath(OUT, "v07_A_delta_eq.png"))
savefig(pA, joinpath(NOTE, "figs", "fig_A_delta_eq.png"))
println("  wrote fig A")

# ---- B: Lambda ----
section("3. Test B: Lambda vs J (NF + tracking + note)")
dmB = 3e-4
amps = [5e-3, 10e-3, 15e-3, 20e-3]
turns_B = max(1000, round(Int, 50 / info.Qs))
rows_B = []
for (plane, getJ, makeu) in (
        (:X, A -> action_from_amp(A, po.x.beta), (A, dm) -> [po.D4[1]*dm + A, po.D4[2]*dm, 0.0, 0.0, 0.0, dm]),
        (:Y, A -> action_from_amp(A, po.y.beta), (A, dm) -> [po.D4[1]*dm, po.D4[2]*dm, A, 0.0, 0.0, dm]))
    for A in amps
        J = getJ(A)
        u = makeu(A, dmB)
        dser, pser = track_longitudinal(lat, ref, C, h, u; turns = turns_B)
        m = measure_synchrotron(dser, pser; Qs = info.Qs)
        xi_x = plane == :X ? xi.x[2] : 0.0
        xi_y = plane == :Y ? xi.y[2] : 0.0
        Jx = plane == :X ? J : 0.0
        Jy = plane == :Y ? J : 0.0
        Lam_note = formula_B_Lambda(C, Jx, xi.x[2], Jy, xi.y[2])
        deq_note = -(pe.eta1 * m.delta_m_fit^2 / 2 + Lam_note) / pe.eta0 + zstar_f[6]
        deq_nf = delta_eq_nf(nf; Jx = Jx, Jy = Jy, delta_m = m.delta_m_fit) + zstar_f[6]
        Lam_fit = -pe.eta0 * (m.delta_eq - zstar_f[6]) - pe.eta1 * m.delta_m_fit^2 / 2
        push!(rows_B, (plane = plane, A = A, J = J, Lam_note = Lam_note, Lam_fit = Lam_fit,
                       deq_tr = m.delta_eq, deq_nf = deq_nf, deq_note = deq_note))
        @printf("  %s A=%.0fmm J=%.3e  Lam_note=%+.4e Lam_fit=%+.4e  deq_tr=%+.4e NF=%+.4e\n",
                plane, A*1e3, J, Lam_note, Lam_fit, m.delta_eq, deq_nf)
    end
end
pB = plot(xlabel = "J [m]", ylabel = "Λ", legend = :topleft,
          title = "Λ vs action (note / tracking)", size = (700, 480))
for pl in (:X, :Y)
    rs = filter(r -> r.plane == pl, rows_B)
    scatter!(pB, [r.J for r in rs], [r.Lam_fit for r in rs], label = "Λ_fit track ($pl)", ms = 5)
    plot!(pB, [r.J for r in rs], [r.Lam_note for r in rs], label = "−(2π/C) J ξ ($pl)", lw = 2)
end
savefig(pB, joinpath(OUT, "v07_B_lambda.png"))
savefig(pB, joinpath(NOTE, "figs", "fig_B_lambda.png"))

# ---- Residual plot ----
section("4. Residual plots")
pR1 = plot(xlabel = "V [MV]", ylabel = "δ*", legend = :topleft,
           title = "Closed-orbit δ* vs RF voltage (L=1 m)", size = (600, 420))
rv = filter(r -> r.Lc == 1.0, rows_res)
scatter!(pR1, [r.V for r in rv], [r.dstar for r in rv], label = "finite cavity FP", ms = 6)
plot!(pR1, [r.V for r in rv], [r.dstar for r in rv], label = nothing, lw = 1)
savefig(pR1, joinpath(OUT, "v07_residual_V.png"))
savefig(pR1, joinpath(NOTE, "figs", "fig_residual_V.png"))
pR2 = plot(xlabel = "L_cav [m]", ylabel = "δ*", legend = :topleft,
           title = "Closed-orbit δ* vs cavity length (V=0.2 MV)", size = (600, 420))
rl = filter(r -> r.V == 0.2 && r.Lc != 1.0 || (r.V == 0.2 && r.Lc == 1.0), rows_res)
# unique by Lc
seen = Dict{Float64,Float64}()
for r in rows_res
    r.V == 0.2 && (seen[r.Lc] = r.dstar)
end
Lcs = sort(collect(keys(seen)))
scatter!(pR2, Lcs, [seen[L] for L in Lcs], label = "finite cavity FP", ms = 6)
plot!(pR2, Lcs, [seen[L] for L in Lcs], label = nothing)
savefig(pR2, joinpath(OUT, "v07_residual_L.png"))
savefig(pR2, joinpath(NOTE, "figs", "fig_residual_L.png"))

# ---- Symplecticity: long-term drift of mean delta ----
section("5. Truncated Taylor vs factorised Lie: long-term δ mean")
dm0 = 1e-3
u0 = [po.D4[1]*dm0, po.D4[2]*dm0, 0.0, 0.0, 0.0, dm0]
z0 = user_to_canonical(u0, ref, C, h)
Nlong = 5000
factors = build_factors(lat_mag)
cav_t = RFCavity("t", V * Ttf, h, Float64(pi))
mean_win = Int[]
mean_t = Float64[]; mean_l = Float64[]
let zt = copy(z0), zl = copy(z0), acc_t = 0.0, acc_l = 0.0
    for n in 1:Nlong
        zt = evaluate_map(Mthin, zt)
        zl = track(factors, zl; turns = 1, ref = ref, C = C)
        zl = rf_kick(zl, cav_t, ref, C)
        acc_t += zt[6]; acc_l += zl[6]
        if n % 500 == 0
            push!(mean_win, n)
            push!(mean_t, acc_t / n)
            push!(mean_l, acc_l / n)
        end
    end
end
@printf("  after %d turns: <δ>_Taylor=%.6e  <δ>_Lie+kick=%.6e  (note+FP=%.6e)\n",
        Nlong, mean_t[end], mean_l[end], formula_A(pe.eta0, pe.eta1, dm0) + zstar_f[6])
pS = plot(xlabel = "turns", ylabel = "running ⟨δ⟩", legend = :right,
          title = "Long-term mean δ: truncated Taylor vs Lie factors", size = (700, 420))
plot!(pS, mean_win, mean_t, label = "truncated Taylor map", lw = 2)
plot!(pS, mean_win, mean_l, label = "Lie factors + thin RF", lw = 2)
hline!(pS, [formula_A(pe.eta0, pe.eta1, dm0) + zstar_f[6]], label = "note+δ*", ls = :dash)
savefig(pS, joinpath(OUT, "v07_symplectic_drift.png"))
savefig(pS, joinpath(NOTE, "figs", "fig_symplectic_drift.png"))

# Save summary table
open(joinpath(OUT, "v07_summary.txt"), "w") do io
    println(io, "Phase v0.7 summary")
    println(io, "Finite-cavity delta* (V=0.2,L=1) = ", zstar_f[6])
    println(io, "NF coeffs match note: cx,cy,cs")
    println(io, "A table:")
    for r in rows_A
        @printf(io, "dm=%.4e note+FP=%+.6e NF+FP=%+.6e track=%+.6e Taylor=%+.6e\n",
                r.dm, r.a_fp, r.b, r.c, r.d)
    end
end
println("Done. Outputs in ", OUT, " and ", NOTE)
