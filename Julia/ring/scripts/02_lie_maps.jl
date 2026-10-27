# =============================================================================
# 02_lie_maps.jl - phase 2 preview: 2nd-order Lie generators, tracking, checks
# =============================================================================
#
# WHAT
#   (phase 2: JET_ORDER = 4 -> H = H2+H3+H4, maps exact to order 3, f3 AND f4 factors)
#   1. For every distinct element: its Hamiltonian H = H2 + H3 (printed),
#      the order-2 Lie map exp(-L:H:) z, the consistency of its linear part
#      with exp(L S A) (script 01), the Dragt-Finn generator f3, the
#      reconstruction error, and symplecticity:
#        (a) truncated Taylor map Jacobian at the origin and at amplitude,
#        (b) factor-by-factor (implicit-midpoint) Jacobian at amplitude.
#   2. Ring: factorised product R_ring o exp(:f3_ring:), compared with the
#      INDEPENDENT Taylor map (RK4 in TaylorN arithmetic, TaylorMap.jl).
#   3. Chromaticity, alpha_0, alpha_1 from the order-2 ring map;
#      second-order path-length / phase coefficients in USER coordinates
#      (phi|delta^2), (phi|x^2), (phi|x'^2), ... .
#   4. Tracking: one turn factor-by-factor vs ring Taylor map; 1000 turns
#      with Courant-Snyder action monitor; symplecticity of the 1-turn
#      tracking Jacobian at amplitude.
#   5. Comparison with COSY output numbers (chromaticity, alpha_0).
#
# WHY
#   Shows that the Lie machinery is correct (two independent ways agree)
#   and that only the factorised tracking is symplectic at amplitude.
#
# HOW TO RUN
#   julia --project=. scripts/02_lie_maps.jl        (about a minute)
#
# WHAT YOU SEE
#   Per-element block, ring block, tables; all errors printed as numbers.
#   Taylor-map symplectic errors at amplitude are ~|z|^2 (NOT zero, expected),
#   factorised tracking errors are ~1e-15.
# =============================================================================
include(joinpath(@__DIR__, "..", "src", "LieRing.jl"))
using .LieRing, Printf, LinearAlgebra, TaylorSeries
include(joinpath(@__DIR__, "reference_values.jl"))

section(t, lines...) = (println("\n--- ", t, " ", "-"^max(0, 66 - length(t))); foreach(l -> println("  > ", l), lines))

println("="^72)
println(" LieRing - 2nd-order Lie maps, symplectic tracking, Taylor-map check")
println("="^72)
of = parse_optim(OPTIM_FILE); lat = of.lattice; ref = lat.ref; C = circumference(lat); hRF = lat.harmonic
z_amp = [3e-3, 2e-4, 2e-3, 2e-4, 1e-3, 2e-3]      # test point at amplitude (canonical)
println("  test amplitude z_amp (x,px,y,py,l,delta) = ", z_amp)

section("1. Elements: H, Lie map, f3, symplecticity",
        "Lie map = exp(-L:H:) z to order 2; f3 = Dragt-Finn cubic generator.")
@printf("  %-5s %10s %10s %10s %10s %10s %10s\n", "elem", "|lin-expLSA|", "recon.err", "TM sym@0", "TM sym@amp", "LF sym@amp", "nterms f3")
for e in unique_elements(lat)
    H = hamiltonian(e, ref)
    M = element_lie_map(e, ref)
    F = dragt_finn_factor(M; name = e.name)
    lin = maximum(abs, map_linear_part(M) - linear_matrix(e, ref))
    rec = maximum(maxabs(a - b) for (a, b) in zip(factor_to_map(F), M))
    s0 = symplectic_error(jacobian_of_map(M, zeros(6)))
    sa = symplectic_error(jacobian_of_map(M, z_amp))
    sf = symplectic_error(jacobian_of_function(z -> apply_factor(F, z), z_amp))
    nt = count(!iszero, vcat([F.f3[3].coeffs]...))
    @printf("  %-5s %10.1e %10.1e %10.1e %10.1e %10.1e %10d\n", e.name, lin, rec, s0, sa, sf, nt)
end
for e in unique_elements(lat)
    e isa SBend || e isa Quadrupole || e isa Sextupole || continue
    e.name in ("BM", "QF", "SD") || continue
    println("\n  ", e.name, ": ", hamiltonian_formula(e))
    println("    H (order<=4) = ", poly_string(hamiltonian(e, ref); digits = 5))
    F = dragt_finn_factor(element_lie_map(e, ref))
    println("    f3           = ", poly_string(F.f3; digits = 5))
end
println("\n  Drift H shows the exact kinematic terms: delta^2/(2 gamma0^2) and -delta^3/(2 beta0^2 gamma0^2), px^2 delta/2 ...")
println("    H_drift = ", poly_string(hamiltonian(Drift("D", 1.0), ref); digits = 5))

section("2. Ring map: factorised Lie product vs independent Taylor map",
        "Lie: F_ring = concatenation of all element factors (f3 transported by R).",
        "Taylor: RK4 integration of dz/ds = S grad H in TaylorN arithmetic, composed.")
factors = build_factors(lat)
Fring = foldl((a, b) -> concatenate(a, b), factors[2:end]; init = factors[1])
Mlie = factor_to_map(Fring)
t0 = time()
Mtm = ring_taylor_map(lat)
@printf("  Taylor map built in %.1f s\n", time() - t0)
R = ring_matrix(lat)
@printf("  |R(Lie) - R(exp LSA product)|      = %.2e\n", maximum(abs, Fring.R - R))
@printf("  |R(Taylor/RK4) - R(exp LSA product)| = %.2e\n", maximum(abs, map_linear_part(Mtm) - R))
d2 = maximum(maxabs(TaylorN(a[2], JET_ORDER) - TaylorN(b[2], JET_ORDER)) for (a, b) in zip(Mlie, Mtm))
m2 = maximum(maxabs(TaylorN(a[2], JET_ORDER)) for a in Mtm)
@printf("  max |2nd-order coeff difference| = %.2e  (largest 2nd-order coeff %.3e)\n", d2, m2)
@printf("  symplectic error of ring Taylor map Jacobian: at 0 %.1e, at z_amp %.1e\n",
        symplectic_error(jacobian_of_map(Mtm, zeros(6))), symplectic_error(jacobian_of_map(Mtm, z_amp)))

section("3. Second-order quantities from the ring map",
        "Chromaticity from tunes on the off-momentum closed orbit; dC/C = a0 d + a1 d^2.")
xi = chromaticity(Mtm)
xil = chromaticity(Mlie)
a0, a1 = momentum_compaction(Mtm, ref, C)
po = periodic_optics(R, ref, C)
@printf("  natural chromaticity  xi_x = %.6f   xi_y = %.6f   (Taylor map)\n", xi[1], xi[2])
@printf("                        xi_x = %.6f   xi_y = %.6f   (Lie product map)\n", xil[1], xil[2])
@printf("  alpha_0 = %.8f (linear optics: %.8f)   alpha_1 = %.6f\n", a0, po.alpha_c, a1)
U = user_map(Mtm, ref, C, hRF)
un = USER_NAMES
println("  selected one-turn coefficients in USER coordinates (x,x',y,y',phi,delta), h = $hRF:")
for (i, ex, lab) in [(5, [0,0,0,0,0,1], "(phi|delta)"), (5, [0,0,0,0,0,2], "(phi|delta^2)"),
                     (5, [2,0,0,0,0,0], "(phi|x^2)"), (5, [0,2,0,0,0,0], "(phi|x'^2)"),
                     (5, [0,0,0,2,0,0], "(phi|y'^2)"), (5, [1,0,0,0,0,1], "(phi|x delta)"),
                     (1, [0,0,0,0,0,1], "(x|delta)"), (1, [0,0,0,0,0,2], "(x|delta^2)"),
                     (1, [2,0,0,0,0,0], "(x|x^2)"), (2, [1,0,0,0,0,1], "(x'|x delta)")]
    @printf("    %-15s = % .6e\n", lab, coeff(U[i], ex))
end

section("4. Tracking", "Factor-by-factor symplectic tracking (implicit midpoint per factor).")
u0 = [2e-3, 0.0, 1e-3, 0.0, 0.0, 1e-3]
z0 = user_to_canonical(u0, ref, C, hRF)
println("  start (user) x=2 mm, y=1 mm, delta=1e-3  -> canonical ", round.(z0, sigdigits = 6))
z1 = track(factors, z0)
z1tm = [evaluate(m, z0) for m in Mtm]
@printf("  one turn: max |Lie tracking - Taylor map| = %.2e (expected O(|z|^3))\n", maximum(abs, z1 - z1tm))
J1 = jacobian_of_function(z -> track(factors, z), z0)
@printf("  symplectic error of 1-turn tracking Jacobian at z0: %.2e (Taylor map: %.2e)\n",
        symplectic_error(J1), symplectic_error(jacobian_of_map(Mtm, z0)))
# Courant-Snyder action about the off-momentum closed orbit
co = closed_orbit(Mtm, z0[6])
Jx(w) = (x = w[1] - co[1]; p = w[2] - co[2]; (po.x.beta * p^2 + 2po.x.alpha * x * p + (1 + po.x.alpha^2) / po.x.beta * x^2) / 2)
Jy(w) = (x = w[3] - co[3]; p = w[4] - co[4]; (po.y.beta * p^2 + 2po.y.alpha * x * p + (1 + po.y.alpha^2) / po.y.beta * x^2) / 2)
w = copy(z0); jx = Float64[]; jy = Float64[]
t0 = time()
for n in 1:1000
    global w = track(factors, w)
    push!(jx, Jx(w)); push!(jy, Jy(w))
end
@printf("  1000 turns in %.1f s; action Jx: mean %.4e, rel. spread %.2e; Jy: mean %.4e, rel. spread %.2e\n",
        time() - t0, sum(jx) / 1000, (maximum(jx) - minimum(jx)) / (sum(jx) / 1000),
        sum(jy) / 1000, (maximum(jy) - minimum(jy)) / (sum(jy) / 1000))
println("  (spread is due to the nonlinear (2nd-order) terms and d-dependence of Twiss, not to loss of symplecticity)")
@printf("  final delta = %.15e (must equal start %.15e: magnetic lattice)\n", w[6], z0[6])

section("5. Comparison with COSY output numbers",
        "COSY chromaticity convention (dQ/d of which energy variable) not verified - see notes.")
@printf("  %-12s %14s %14s %12s  %s\n", "quantity", "LieRing", "COSY", "diff", "source")
for (n, a, (b, src)) in [("xi_x", xi[1], COSY_REF.chrom_x), ("xi_y", xi[2], COSY_REF.chrom_y),
                         ("alpha_0", a0, COSY_REF.alpha0)]
    @printf("  %-12s %14.8f %14.8f %12.3e  %s\n", n, a, b, a - b, src)
end
println("\nDone.")
