# =============================================================================
# 05_rf_acceleration.jl - finite-length RF cavity and acceleration 10 -> 270 MeV
# =============================================================================
#
# WHAT
#   1. Parameters: injection energy (NOT in magnetic_2.opt -> ASSUMED 5 MeV/u,
#      i.e. 10 MeV for the deuteron, the usual Nuclotron linac energy - flagged),
#      final 270 MeV (135 MeV/u, header of the .opt), cavity L = 1 m, V = 0.2 MV,
#      h = 5, phi_s = 150 deg (below transition: eta0 < 0 requires cos(phi_s) < 0
#      in our convention phi = -h w0 dt).  V is chosen large so that the ramp
#      takes ~2600 turns (a real ramp has ~1e5-1e6 turns; physics is the same
#      as long as the ramp is adiabatic: T_s << ramp time).
#   2. Cavity checks: energy gain of the synchronous particle vs Z V T sin(phi_s)
#      (transit-time factor T), thin limit vs Tracking.rf_kick, symplecticity
#      of the cavity map (exact Jacobian via TaylorN).
#   3. Synchrotron tune at FIXED energy (phi_s = 180 deg, stationary bucket) at
#      10 MeV and 270 MeV: tracking vs Q_s formula with eta0 (and eta1 noted).
#   4. Ramp: particles  #1 reference/synchronous, #2 x0 = 2 mm, #3 y0 = 2 mm,
#      #4 delta0 = 1e-3 (started on its dispersion orbit); then 1000 turns flat
#      top (phi_s = 180 deg).  Energy vs turn, phi-delta, adiabatic damping of
#      x, y amplitudes vs A0 sqrt(p_inj/p0), synchrotron tune vs formula, and
#      the longitudinal adiabatic invariant  (area in (dt, dE) conserved).
#
# ASSUMPTIONS: magnets ramp with p0; RF frequency = h f0(p0) locked to the
#   synchronous particle; cavity field purely longitudinal on axis; no space
#   charge, no fringe fields; lattice = natural magnetic_2 (sextupoles 0).
#
# HOW TO RUN:  julia --project=. scripts/05_rf_acceleration.jl   (a few minutes)
# WHAT YOU SEE: numbers + output/acceleration.png
# =============================================================================
include(joinpath(@__DIR__, "..", "src", "LieRing.jl"))
using .LieRing, Printf, LinearAlgebra, TaylorSeries
include(joinpath(@__DIR__, "reference_values.jl"))
section(t, lines...) = (println("\n--- ", t, " ", "-"^max(0, 66 - length(t))); foreach(l -> println("  > ", l), lines))
println("="^72); println(" LieRing phase 2 - finite RF cavity, synchrotron motion, acceleration"); println("="^72)

of = parse_optim(OPTIM_FILE); lat270 = of.lattice; ref270 = lat270.ref
mass = ref270.mass_MeV
E_inj = 10.0          # ASSUMPTION (not in the .opt file): 5 MeV/u
E_fin = of.energy_MeV # 270 MeV from the .opt header
Lc, V, h, phis = 1.0, 0.2, 5, deg2rad(150.0)

section("1. Parameters")
ref_inj = reference_from_kinetic(DEUTERON, E_inj; mass_MeV = mass)
cav = Cavity("RF", Lc, V, h, phis)
lat_rf270, icav = insert_cavity(lat270, cav)
C = circumference(lat_rf270)
lat_rf = Lattice(lat_rf270.name, lat_rf270.elements, ref_inj, h)
T = transit_factor(cav, ref270, C)
@printf("  injection E_kin = %.1f MeV (ASSUMED, 5 MeV/u; not given in magnetic_2.opt), final %.1f MeV (135 MeV/u)\n", E_inj, E_fin)
@printf("  cavity: L = %.2f m in the first OBM drift (element %d), V = %.3f MV, h = %d, phi_s = %.1f deg\n", Lc, icav, V, h, rad2deg(phis))
@printf("  C = %.6f m (unchanged), transit factor T = sin(th/2)/(th/2), th = 2 pi h L/C = %.5f rad -> T = %.6f\n", C, 2pi * h * Lc / C, T)
for (nm, rf) in (("injection", ref_inj), ("final", ref270))
    f0 = revolution_frequency(rf, C)
    @printf("  %-9s: beta = %.6f gamma = %.6f p0 = %.4f MeV/c Brho = %.5f T m  f0 = %.4f kHz  f_RF = %.4f MHz\n",
            nm, rf.beta, rf.gamma, rf.p0_MeV, rf.brho_Tm, f0 / 1e3, h * f0 / 1e6)
end
dE = V * T * sin(phis)
@printf("  energy gain per turn Z V T sin(phi_s) = %.6f MeV -> about %d turns for the ramp\n", dE, ceil(Int, (E_fin - E_inj) / dE))

section("2. Cavity checks", "Single passes through the cavity at 270 MeV (no renormalisation).")
z = exact_track_element(zeros(6), cav, ref270, C)
@printf("  synchronous particle: delta_out = %.10e  -> dE = %.9f MeV   vs Z V T sin(phi_s) = %.9f MeV\n",
        z[6], (sqrt((ref270.p0_MeV * (1 + z[6]))^2 + mass^2) - sqrt(ref270.p0_MeV^2 + mass^2)), dE)
lt = 0.05 * C * velocity_ratio(0.0, ref270) / (2pi * h)        # phi = +0.05 rad
zin = [0.0, 0.0, 0.0, 0.0, lt, 2e-4]
thin = rf_kick(zin, RFCavity("thin", V, h, phis), ref270, C)
for Lx in (0.1, 0.01, 0.001)
    zc = exact_track_element(zin, Cavity("c", Lx, V, h, phis), ref270, C)
    @printf("  thin limit L = %-6g m: delta_out = %.12e  (thin kick %.12e, diff %.2e)\n", Lx, zc[6], thin[6], zc[6] - thin[6])
end
Jc = jacobian_of_function(w -> exact_track_element(w, cav, ref270, C), [1e-3, 1e-4, 1e-3, 1e-4, lt, 2e-4])
@printf("  symplectic error of cavity map Jacobian: %.2e\n", symplectic_error(Jc))

section("3. Synchrotron tune at fixed energy (phi_s = 180 deg)",
        "Q_s = sqrt(h eta0 Z V T cos(phi_s) / (2 pi beta^2 E)), eta0 = alpha0 - 1/gamma^2.")
function fixed_energy_qs(rf; turns = 1024, d0 = 1e-4)
    L = Lattice("fix", set_cavity_phase(lat_rf, pi).elements, rf, h)
    po = periodic_optics(ring_matrix(L), rf, C)
    zs = [zeros(6), [po.D4[1] * d0, po.D4[2] * d0, 0, 0, 0, d0]]
    ref = rf; dd = Float64[]
    for n in 1:turns
        for e in L.elements
            zs = [exact_track_element(w, e, ref, C; ds = 0.2) for w in zs]
            e isa Cavity && ((zs, ref) = renormalize(zs, ref, 1))
        end
        push!(dd, zs[2][6])
    end
    eta0 = po.alpha_c - 1 / rf.gamma^2
    return measure_tune(dd), synchrotron_tune(rf, eta0, V, h, pi; T = T), eta0, po.alpha_c
end
qs_rows = []
for rf in (ref_inj, ref270)
    qm, qf, eta0, ac = fixed_energy_qs(rf)
    push!(qs_rows, (rf.Ekin_MeV, qm, qf, eta0))
    qd = acos(1 - 2pi^2 * qf^2) / 2pi   # one-kick-per-turn map: cos(2 pi Qs) = 1 - (2 pi Qs_smooth)^2 / 2
    @printf("  E = %6.1f MeV: alpha_c = %.8f eta0 = %.8f  Q_s tracking = %.7f  smooth formula = %.7f (diff %.2e)  discrete-kick formula = %.7f (diff %.2e)\n",
            rf.Ekin_MeV, ac, eta0, qm, qf, qm / qf - 1, qd, qm / qd - 1)
end

section("4. Longitudinal Hamiltonian parameters (eta0, eta1) vs energy",
        "alpha0, alpha1 from the order-3 ring Taylor map at 10 and 270 MeV (RF off).",
        "eta0 = alpha0 - c1, eta1 = alpha1 - alpha0 c1 - c2 + c1^2, (c1, c2) = velocity expansion",
        "at the energy of interest; alpha1 is interpolated linearly in E_kin between the two maps.")
pe = Dict{Float64,Any}()
for rf in (ref_inj, ref270)
    Lm = Lattice("m", lat_rf.elements, rf, h)
    pe[rf.Ekin_MeV] = path_expansion(ring_taylor_map(Lm), rf, C)
    @printf("  E = %6.1f MeV: alpha0 = %.8f alpha1 = %.8f eta0 = %.8f eta1 = %.8f\n", rf.Ekin_MeV,
            pe[rf.Ekin_MeV].alpha0, pe[rf.Ekin_MeV].alpha1, pe[rf.Ekin_MeV].eta0, pe[rf.Ekin_MeV].eta1)
end
a0 = pe[E_inj].alpha0
alpha1_at(E) = pe[E_inj].alpha1 + (pe[E_fin].alpha1 - pe[E_inj].alpha1) * (E - E_inj) / (E_fin - E_inj)
function long_params(E, ph)
    rf = reference_from_kinetic(DEUTERON, E; mass_MeV = mass)
    c1, c2 = velocity_expansion(rf)
    a1 = alpha1_at(E)
    LongitudinalParams(h, a0 - c1, a1 - a0 * c1 - c2 + c1^2, V * T / (rf.beta^2 * rf.gamma * mass), ph)
end
lp_inj = long_params(E_inj, phis)
bh = bucket_height(lp_inj)
fp = fixed_points(lp_inj)
@printf("  injection bucket (phi_s = %.0f deg): phi_u = %.4f rad, half-height +%.4e / %.4e, 2nd fixed-point pair at delta = %.3f\n",
        rad2deg(phis), fp.unstable, bh[1], bh[2], fp.delta2)
lp_top = long_params(E_fin, pi)
@printf("  stationary bucket at 270 MeV (phi_s = 180 deg): half-height %.4e\n", bucket_height(lp_top)[1])

section("5. Acceleration ramp with an ensemble",
        "Particles: #1 synchronous, #2 x0 = 2 mm, #3 y0 = 2 mm, #4 delta0 = 1e-3 (on its dispersion orbit),",
        "#5.. ensemble (x = y = 0 on the dispersion orbit): phi ~ N(0, 0.6 rad), delta ~ N(0, 0.25 bucket height),",
        "kept only if inside the injection separatrix -> a MISMATCHED, filling distribution that filaments.")
using Random
Random.seed!(1)
po_inj = periodic_optics(ring_matrix(lat_rf), ref_inj, C)
_, inside = separatrix(lp_inj)
ensemble = Vector{Vector{Float64}}()
while length(ensemble) < 60
    f = 0.6 * randn(); d = 0.25 * bh[1] * randn()
    inside(f, d) || continue
    push!(ensemble, [po_inj.D4[1] * d, po_inj.D4[2] * d, 0, 0, f, d])
end
users = [zeros(6), [2e-3, 0, 0, 0, 0, 0], [0, 0, 2e-3, 0, 0, 0], [po_inj.D4[1] * 1e-3, po_inj.D4[2] * 1e-3, 0, 0, 0, 1e-3]]
append!(users, ensemble)
# (b) one turn, element by element, at injection (before the ramp)
zs0 = [user_to_canonical(u, ref_inj, C, h) for u in users]
snaps, enames, svals, _, _ = turn_snapshots(lat_rf, ref_inj, zs0)
t0 = time()
rec, flat = accelerate(lat_rf, ref_inj, users; Ekin_target = E_fin, flat_turns = 1000)
@printf("  tracked %d particles x %d turns in %.1f s; 270 MeV reached at turn %d (E_kin = %.6f MeV)\n",
        length(users), rec.turn[end], time() - t0, flat, rec.Ekin[flat+1])
@printf("  mean energy gain per turn during ramp: %.6f MeV (expected %.6f)\n", (rec.Ekin[flat+1] - rec.Ekin[1]) / flat, dE)
nlost = count(j -> abs(rec.u[j][end][6]) > 2 * bucket_height(lp_top)[1], 5:length(users))
println("  ensemble particles outside 2x the final bucket height at the end: $nlost of $(length(users) - 4)")
po = periodic_optics(ring_matrix(lat270), ref270, C)
bx, ax, by, ay = po.x.beta, po.x.alpha, po.y.beta, po.y.alpha
amp(x, xp, b, a) = sqrt(b * ((1 + a^2) / b * x^2 + 2a * x * xp + b * xp^2))
Ax = [amp(u[1], u[2], bx, ax) for u in rec.u[2]]
Ay = [amp(u[3], u[4], by, ay) for u in rec.u[3]]
pred = sqrt.(rec.p0[1] ./ rec.p0)
@printf("  x amplitude: start %.4f mm, end %.4f mm, ratio %.5f;  y: %.4f -> %.4f mm, ratio %.5f\n",
        Ax[1] * 1e3, Ax[end] * 1e3, Ax[end] / Ax[1], Ay[1] * 1e3, Ay[end] * 1e3, Ay[end] / Ay[1])
@printf("  adiabatic-damping prediction sqrt(p_inj/p_final) = sqrt(%.4f/%.4f) = %.5f\n", rec.p0[1], rec.p0[end], pred[end])
# Q_s: formula every turn (smooth curve), tracking in sliding windows (markers)
phis_turn(n) = n > flat ? pi : phis
qs_smooth = Float64[]; qs_kick = Float64[]
for (n, E) in zip(rec.turn, rec.Ekin)
    lp = long_params(E, phis_turn(n))
    a = 2pi * h * lp.eta0 * lp.A * cos(lp.phi_s)   # (2 pi Qs)^2 for the smooth Hamiltonian
    q = a > 0 ? sqrt(a) / 2pi : NaN
    push!(qs_smooth, q); push!(qs_kick, a > 0 && a < 4 ? acos(1 - a / 2) / 2pi : NaN)
end
dser = [u[6] for u in rec.u[4]]
W = 200; step = 50
qs_t = Float64[]; qs_m = Float64[]
for a in 1:step:length(dser)-W
    b = a + W - 1
    (a <= flat + 1 <= b) && continue           # do not mix ramp and flat top in one window
    push!(qs_t, rec.turn[div(a + b, 2)]); push!(qs_m, measure_tune(dser[a:b]))
end
println("  synchrotron tune: formula every turn; tracking = fit of x_{n+1}+x_{n-1} = 2cos(2piQ)x_n + b")
println("  over sliding windows of $W turns (step $step) of delta of particle #4; marker at the window centre.")
@printf("  %8s %10s %12s %12s %12s\n", "turn", "E_kin", "Q_s track", "Q_s smooth", "Q_s kick")
for k in 1:6:length(qs_t)
    n = Int(qs_t[k])
    @printf("  %8d %10.3f %12.6f %12.6f %12.6f\n", n, rec.Ekin[n+1], qs_m[k], qs_smooth[n+1], qs_kick[n+1])
end# longitudinal adiabatic invariant: area in (dt, dE) ~ phi_hat * delta_hat * beta^2 E /(h w0) conserved
function amps(a, b)
    ph = maximum(abs(rec.u[4][n][5]) for n in a:b); dl = maximum(abs(rec.u[4][n][6]) for n in a:b)
    n = div(a + b, 2)
    rf = reference_from_kinetic(DEUTERON, rec.Ekin[n]; mass_MeV = mass)
    w0 = 2pi * revolution_frequency(rf, C)
    return ph, dl, ph * dl * rf.beta^2 * rf.gamma * mass / (h * w0)
end
println("  longitudinal adiabatic invariant: area of the (dt, dE) ellipse ~ phi_hat * dE_hat/(h w0),")
println("  dE_hat = beta^2 E delta_hat; amplitudes = max over ~1.5 synchrotron periods around the turn:")
@printf("  %8s %10s %10s %12s %14s\n", "turn", "E_kin", "phi_hat", "delta_hat", "invariant")
inv0 = 0.0
for (n0, w) in ((1, 20), (300, 120), (1000, 100), (2000, 110), (flat - 115, 115))
    ai = amps(n0, n0 + w)
    n0 == 300 && (global inv0 = ai[3])
    rel = inv0 > 0 ? @sprintf("(/turn-300 value: %.4f)", ai[3] / inv0) : ""
    @printf("  %8d %10.3f %10.5f %12.4e %14.5e %s\n", n0, rec.Ekin[n0], ai[1], ai[2], ai[3], rel)
end
println("  (the first ~200 turns are NOT adiabatic: at 10 MeV the gain of 0.1 MeV/turn is 1% of E_kin per turn)")
outdir = joinpath(@__DIR__, "..", "output"); mkpath(outdir)
lp_end = long_params(rec.Ekin[flat+1], phis)
ens = 5:length(users)
bend = bucket_height(lp_end)[1]
info = (flat_start = flat, n_end = flat, lp_end = lp_end, dmax_end = 1.6 * bend, philim = (-2.6, 1.6),
        ens_phi = [rec.u[j][flat+1][5] for j in ens], ens_delta = [rec.u[j][flat+1][6] for j in ens],
        tr_phi = [u[5] for u in rec.u[4][flat-299:flat+1]], tr_delta = [u[6] for u in rec.u[4][flat-299:flat+1]],
        Ax = Ax, Ay = Ay, Ax_pred = Ax[1] .* pred, qs_smooth = qs_smooth, qs_kick = qs_kick,
        qs_turn = qs_t, qs_meas = qs_m, win = W)
png = plot_acceleration(rec, info, joinpath(outdir, "acceleration.png"))
println("  written: ", abspath(png))
section("6. Animated GIFs")
frames = collect(0:20:rec.turn[end])
lps = [long_params(rec.Ekin[n+1], phis_turn(n)) for n in frames]
g1 = gif_ramp(rec, frames, lps, joinpath(outdir, "longitudinal_ramp.gif"); ens = ens, philim = (-pi - 0.5, pi + 0.5), dlim = 1.4 * bh[1])
println("  (a) ", abspath(g1), "  (every 20th turn; left: absolute, right: delta / current bucket height)")
g2 = gif_elements(snaps, enames, svals, lp_inj, joinpath(outdir, "longitudinal_elements.gif"); ens = ens,
                  philim = (-pi - 0.5, pi + 0.5), dlim = 1.4 * bh[1])
println("  (b) ", abspath(g2), "  (one turn at injection, frame = exit of each element)")
println("\nDone.")