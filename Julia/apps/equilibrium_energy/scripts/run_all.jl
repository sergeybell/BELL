# =============================================================================
# run_all.jl - Tests A, B, C for delta_eq (note ch.1-3) on magnetic_2
# =============================================================================
# Phase in tables/plots: NOTE convention (phi>0 = LATE) = -LieRing phi.
# delta_eq from harmonic LS fit (not naive mean).
# julia --project=. scripts/run_all.jl
# =============================================================================
using EquilibriumEnergy
using LieRing, Printf, Plots, Statistics, LinearAlgebra
gr()

const OPTIM = raw"C:\GIT_REPOS\SCT\OptiM\magnetic\magnetic_2.opt"
const OUT = joinpath(@__DIR__, "..", "output")
mkpath(OUT)

section(t) = println("\n", "="^72, "\n ", t, "\n", "="^72)

section("0. Setup (270 MeV, stationary bucket, phi_s = pi LieRing = late sync in note)")
lat, ref, C, h, V, Ttf, pe, xi, po, info = setup_ring(OPTIM)
@printf("  mass = %.8f MeV (CODATA); OptiM header Mass = %.3f (ignored)\n", ref.mass_MeV, info.mass_file)
@printf("  C = %.6f m  beta=%.6f gamma=%.6f\n", C, ref.beta, ref.gamma)
@printf("  eta0 = %.8f  eta1 = %.8f  alpha0 = %.8f  alpha1 = %.8f\n", pe.eta0, pe.eta1, pe.alpha0, pe.alpha1)
@printf("  xi1_x = %.6f  xi1_y = %.6f\n", xi.x[2], xi.y[2])
@printf("  Qs_smooth = %.6f  (phi_s = 180 deg);  eta0*cos(phi_s) = %.6f  (need >0 for stability in energy form)\n",
        info.Qs, pe.eta0 * cos(Float64(pi)))
@printf("  Twiss s=0: bx=%.3f ax=%.4e  by=%.3f ay=%.4e  D=%.5f\n",
        po.x.beta, po.x.alpha, po.y.beta, po.y.alpha, po.D4[1])
# Verify below transition + stationary: LieRing synchrotron_tune finite
@assert pe.eta0 < 0 "expected below transition"
@assert isfinite(info.Qs) && info.Qs > 0 "Qs should be real and positive at phi_s=pi"
println("  OK: below transition, stationary bucket at phi_s=pi is stable.")

# -----------------------------------------------------------------------------
section("A. Longitudinal only: delta_m scan, zero betatron (on dispersion orbit)")
dms = [0.5e-3, 0.75e-3, 1.0e-3, 1.5e-3, 2.0e-3]
# turns: many synchrotron periods; Qs~0.013 => T~77 turns; use ~40 periods
turns_A = max(800, round(Int, 40 / info.Qs))
@printf("  tracking %d turns (~%.0f synchrotron periods)\n", turns_A, turns_A * info.Qs)
rows_A = []
for dm0 in dms
    # start at phi=0 (both conventions), delta=dm0 on dispersion orbit, x'=y'=0
    u0 = [po.D4[1]*dm0, po.D4[2]*dm0, 0.0, 0.0, 0.0, dm0]
    dser, pser = track_longitudinal(lat, ref, C, h, u0; turns = turns_A)
    m = measure_synchrotron(dser, pser; Qs = info.Qs)
    dmf = m.delta_m_fit
    pred = formula_A(pe.eta0, pe.eta1, dmf)
    dE_m = dE_eq(ref, m.delta_eq, dmf)
    dE_p = dE_eq(ref, pred, dmf)
    push!(rows_A, (dm0 = dm0, m = m, pred = pred, dE_m = dE_m, dE_p = dE_p))
    @printf("  dm0=%.2e  dm_fit=%.4e dm_hr=%.4e  <d2>=%.4e  deq=%+.4e form=%+.4e rel=%+.2e err=%.2e phat=%.4f\n",
            dm0, dmf, m.delta_m, m.delta2, m.delta_eq, pred, m.delta_eq/pred - 1, m.err_split, m.phi_hat)
end

pA1 = plot(xlabel = "delta_m", ylabel = "delta_eq", legend = :bottomright,
           title = "A: equilibrium shift (no betatron)")
scatter!(pA1, [r.m.delta_m_fit for r in rows_A], [r.m.delta_eq for r in rows_A],
         label = "tracking <delta>", ms = 5)
plot!(pA1, [r.m.delta_m_fit for r in rows_A], [r.pred for r in rows_A],
      label = "-eta1 dm^2 / (2 eta0)", lw = 2)
pA2 = plot(xlabel = "delta_m", ylabel = "<(delta-deq)^2> / (dm^2/2)", legend = :bottomright,
           title = "A: <delta^2> check (expect ~1)")
scatter!(pA2, [r.m.delta_m_fit for r in rows_A],
         [r.m.delta2 / (r.m.delta_m_fit^2 / 2) for r in rows_A], label = "tracking", ms = 5)
hline!(pA2, [1.0], label = "1", ls = :dash)
savefig(plot(pA1, pA2, layout = (1, 2), size = (1200, 450)), joinpath(OUT, "A_delta_eq.png"))
println("  wrote ", joinpath(OUT, "A_delta_eq.png"))

# -----------------------------------------------------------------------------
section("B. Betatron contribution Lambda")
# Choose dm small so Lambda stands out; A_x, A_y large
dmB = 3e-4
amps = [5e-3, 10e-3, 15e-3, 20e-3]   # 5..20 mm
turns_B = max(1000, round(Int, 50 / info.Qs))
@printf("  dm0 = %.2e  turns = %d;  Lambda = -2pi/C (Jx xi_x + Jy xi_y)\n", dmB, turns_B)
# baseline longitudinal-only prediction at this dm (after measuring dm)
u00 = [po.D4[1]*dmB, po.D4[2]*dmB, 0.0, 0.0, 0.0, dmB]
d0, p0 = track_longitudinal(lat, ref, C, h, u00; turns = turns_B)
m0 = measure_synchrotron(d0, p0; Qs = info.Qs)
deq_long = formula_A(pe.eta0, pe.eta1, m0.delta_m_fit)
@printf("  baseline (J=0): dm_fit=%.4e deq_meas=%+.4e deq_form=%+.4e err=%.2e\n", m0.delta_m_fit, m0.delta_eq, deq_long, m0.err_split)

rows_Bx = []; rows_By = []
for A in amps
    Jx = action_from_amp(A, po.x.beta)
    # x oscillation about dispersion orbit: x = D*dm + A, x'= D'*dm (alpha~0)
    u = [po.D4[1]*dmB + A, po.D4[2]*dmB, 0.0, 0.0, 0.0, dmB]
    dser, pser = track_longitudinal(lat, ref, C, h, u; turns = turns_B)
    m = measure_synchrotron(dser, pser; Qs = info.Qs)
    Lam = formula_B_Lambda(C, Jx, xi.x[2], 0.0, xi.y[2])
    deq_pred = -(pe.eta1 * m.delta_m_fit^2 / 2 + Lam) / pe.eta0
    Lam_fit = -pe.eta0 * m.delta_eq - pe.eta1 * m.delta_m_fit^2 / 2
    push!(rows_Bx, (A = A, J = Jx, m = m, Lam = Lam, Lam_fit = Lam_fit, deq_pred = deq_pred))
    @printf("  X A=%.0f mm  J=%.3e  Lam_form=%+.4e  Lam_fit=%+.4e  rel=%+.2e  deq_meas=%+.4e deq_pred=%+.4e\n",
            A*1e3, Jx, Lam, Lam_fit, Lam_fit/Lam - 1, m.delta_eq, deq_pred)
end
for A in amps
    Jy = action_from_amp(A, po.y.beta)
    u = [po.D4[1]*dmB, po.D4[2]*dmB, A, 0.0, 0.0, dmB]
    dser, pser = track_longitudinal(lat, ref, C, h, u; turns = turns_B)
    m = measure_synchrotron(dser, pser; Qs = info.Qs)
    Lam = formula_B_Lambda(C, 0.0, xi.x[2], Jy, xi.y[2])
    deq_pred = -(pe.eta1 * m.delta_m_fit^2 / 2 + Lam) / pe.eta0
    Lam_fit = -pe.eta0 * m.delta_eq - pe.eta1 * m.delta_m_fit^2 / 2
    push!(rows_By, (A = A, J = Jy, m = m, Lam = Lam, Lam_fit = Lam_fit, deq_pred = deq_pred))
    @printf("  Y A=%.0f mm  J=%.3e  Lam_form=%+.4e  Lam_fit=%+.4e  rel=%+.2e  deq_meas=%+.4e deq_pred=%+.4e\n",
            A*1e3, Jy, Lam, Lam_fit, Lam_fit/Lam - 1, m.delta_eq, deq_pred)
end

# Control: sextupoles to xi~0 (optional Optimizers solve_A)
section("B-control. Sextupoles to xi~0 (Optimizers.solve_A) -> Lambda ~ 0")
lat_mag0 = Lattice(lat.name, [e isa Cavity ? Drift(e.name, e.L) : e for e in lat.elements], ref, h)
# remove cavity for placement (use original magnetic lattice from file)
of0 = parse_optim(OPTIM)
lat_m = Lattice(of0.lattice.name, of0.lattice.elements, ref, h)
po_m = periodic_optics(ring_matrix(lat_m), ref, C)
opt_m = optics_along(lat_m, po_m; ds = 0.05)
lat_s, fam_idx, _, _ = place_sextupole_families(lat_m, opt_m)
base, cols = response(lat_s, fam_idx)
kA = solve_A(base, cols)
lat_corr = set_family_strengths(lat_s, fam_idx, kA)
# check xi
Mc = ring_taylor_map(lat_corr)
xic = chromaticity_orders(Mc)
@printf("  after solve_A: xi1_x = %.4e  xi1_y = %.4e  (target 0)\n", xic.x[2], xic.y[2])
# rebuild RF lattice with corrected magnets
cav = Cavity("RF", info.Lc, V, h, Float64(pi))
lat_c, _ = insert_cavity(lat_corr, cav)
lat_c = Lattice(lat_c.name, lat_c.elements, ref, h)
A_test = 15e-3
Jx_t = action_from_amp(A_test, po.x.beta)
u = [po.D4[1]*dmB + A_test, po.D4[2]*dmB, 0.0, 0.0, 0.0, dmB]
# dispersion of corrected lattice may differ slightly — use po_m / closed from lat_corr
po_c = periodic_optics(ring_matrix(lat_corr), ref, circumference(lat_corr))
u = [po_c.D4[1]*dmB + A_test, po_c.D4[2]*dmB, 0.0, 0.0, 0.0, dmB]
dser, pser = track_longitudinal(lat_c, ref, circumference(lat_c), h, u; turns = turns_B)
m = measure_synchrotron(dser, pser; Qs = info.Qs)
pec = path_expansion(Mc, ref, circumference(lat_corr))
Lam_nat = formula_B_Lambda(C, Jx_t, xi.x[2], 0.0, xi.y[2])
Lam_cor = formula_B_Lambda(circumference(lat_corr), Jx_t, xic.x[2], 0.0, xic.y[2])
Lam_fit = -pec.eta0 * m.delta_eq - pec.eta1 * m.delta_m_fit^2 / 2
@printf("  A_x=15 mm: Lam_form(natural xi)=%+.4e  Lam_form(corr xi)=%+.4e  Lam_fit=%+.4e\n",
        Lam_nat, Lam_cor, Lam_fit)
@printf("  deq_meas=%+.4e  (expect ~ longitudinal only when xi~0)\n", m.delta_eq)

pB = plot(xlabel = "A [mm]", ylabel = "Lambda", legend = :topleft, title = "B: Lambda vs amplitude")
scatter!(pB, [r.A*1e3 for r in rows_Bx], [r.Lam_fit for r in rows_Bx], label = "Lam_fit (X)", ms = 5)
plot!(pB, [r.A*1e3 for r in rows_Bx], [r.Lam for r in rows_Bx], label = "-(2pi/C) Jx xi_x", lw = 2)
scatter!(pB, [r.A*1e3 for r in rows_By], [r.Lam_fit for r in rows_By], label = "Lam_fit (Y)", ms = 5, shape = :diamond)
plot!(pB, [r.A*1e3 for r in rows_By], [r.Lam for r in rows_By], label = "-(2pi/C) Jy xi_y", lw = 2, ls = :dash)
savefig(pB, joinpath(OUT, "B_lambda.png"))
println("  wrote ", joinpath(OUT, "B_lambda.png"))

# -----------------------------------------------------------------------------
section("C. Cubic phase: larger phi_hat ~ 0.3..0.5")
# Lower V to increase phi_hat = w_rf |eta0| dm / w_s ; w_s^2 ~ V, so phi_hat ~ 1/sqrt(V)
# Or larger dm. Target phat 0.3..0.5 with still |eta1 dm/eta0|<<1
# At V=0.2, phat ~ (from A) was small. Use V_c = 0.05 MV and dm up to 3e-3
Vc = 0.05
lat_c2, ref2, C2, h2, _, T2, pe2, _, po2, info2 = setup_ring(OPTIM; V = Vc)
@printf("  V = %.3f MV  Qs = %.5f\n", Vc, info2.Qs)
turns_C = max(1200, round(Int, 50 / info2.Qs))
# choose dm so phat ~ w_rf|eta0|dm/ws ; phat ≈ 2pi h |eta0| dm / (2pi Qs) = h|eta0|dm/Qs
# phat ≈ h*|eta0|*dm/Qs
dms_C = Float64[]
for target in (0.30, 0.35, 0.40, 0.45, 0.50)
    dm = target * info2.Qs / (h * abs(pe2.eta0))
    push!(dms_C, dm)
end
function run_C(lat_c2, ref2, C2, h2, po2, pe2, info2, dms_C, turns_C)
    rows = []
    for dm0 in dms_C
        u0 = [po2.D4[1]*dm0, po2.D4[2]*dm0, 0.0, 0.0, 0.0, dm0]
        dser, pser = track_longitudinal(lat_c2, ref2, C2, h2, u0; turns = turns_C)
        m = measure_synchrotron(dser, pser; Qs = info2.Qs)
        dmf = m.delta_m_fit
        fac = formula_C_factor(m.phi_hat)
        pred = formula_C(pe2.eta0, pe2.eta1, dmf, m.phi_hat)
        pred0 = formula_A(pe2.eta0, pe2.eta1, dmf)
        ratio = m.delta2 / (dmf^2 / 2)
        Om_ratio = m.Qs_meas / info2.Qs
        push!(rows, (dm0 = dm0, m = m, fac = fac, pred = pred, pred0 = pred0, ratio = ratio, Om_ratio = Om_ratio))
        @printf("  dm0=%.3e dm_fit=%.3e phat=%.3f  <d2>/(dm2/2)=%.4f (1-p2/32=%.4f)  deq=%+.3e form=%+.3e form0=%+.3e  Qs_m/Qs0=%.4f (1-p2/16=%.4f) err=%.2e\n",
                dm0, dmf, m.phi_hat, ratio, fac, m.delta_eq, pred, pred0, Om_ratio, formula_C_Omega(m.phi_hat), m.err_split)
    end
    return rows
end
rows_C = run_C(lat_c2, ref2, C2, h2, po2, pe2, info2, dms_C, turns_C)

pC1 = plot(xlabel = "phi_hat (note, late+)", ylabel = "<(d-deq)^2> / (dm^2/2)",
           legend = :bottomleft, title = "C: waveform factor")
scatter!(pC1, [r.m.phi_hat for r in rows_C], [r.ratio for r in rows_C], label = "tracking", ms = 5)
plot!(pC1, [r.m.phi_hat for r in rows_C], [r.fac for r in rows_C], label = "1 - phi_hat^2/32", lw = 2)
pC2 = plot(xlabel = "phi_hat", ylabel = "delta_eq", legend = :bottomright, title = "C: delta_eq")
scatter!(pC2, [r.m.phi_hat for r in rows_C], [r.m.delta_eq for r in rows_C], label = "tracking", ms = 5)
plot!(pC2, [r.m.phi_hat for r in rows_C], [r.pred for r in rows_C], label = "with (1-p^2/32)", lw = 2)
plot!(pC2, [r.m.phi_hat for r in rows_C], [r.pred0 for r in rows_C], label = "without cubic", lw = 1, ls = :dash)
pC3 = plot(xlabel = "phi_hat", ylabel = "Qs_meas / Qs_smooth", legend = :bottomleft, title = "C: frequency")
scatter!(pC3, [r.m.phi_hat for r in rows_C], [r.Om_ratio for r in rows_C], label = "tracking", ms = 5)
plot!(pC3, [r.m.phi_hat for r in rows_C], [formula_C_Omega(r.m.phi_hat) for r in rows_C],
      label = "1 - phi_hat^2/16", lw = 2)
savefig(plot(pC1, pC2, pC3, layout = (1, 3), size = (1500, 420)), joinpath(OUT, "C_cubic_phase.png"))
println("  wrote ", joinpath(OUT, "C_cubic_phase.png"))

# -----------------------------------------------------------------------------
section("Summary tables")
println("A:")
@printf("  %10s %12s %12s %12s %10s %10s\n", "dm_meas", "deq_meas", "deq_form", "rel.err", "err_split", "phat")
for r in rows_A
    @printf("  %10.4e %12.4e %12.4e %10.2e %10.2e %10.4f\n",
            r.m.delta_m_fit, r.m.delta_eq, r.pred, r.m.delta_eq/r.pred-1, r.m.err_split, r.m.phi_hat)
end
println("B (X):")
@printf("  %8s %12s %12s %10s\n", "A_mm", "Lam_form", "Lam_fit", "rel")
for r in rows_Bx
    @printf("  %8.1f %12.4e %12.4e %10.2e\n", r.A*1e3, r.Lam, r.Lam_fit, r.Lam_fit/r.Lam-1)
end
println("Done. Outputs in ", OUT)
