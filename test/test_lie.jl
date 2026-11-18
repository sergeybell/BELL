# RF s-Hamiltonian (thin limit = rf_kick). One testset, so it stays in this file
# rather than a separate test_rf.jl.
using LinearAlgebra
using TaylorSeries

@testset "RF Hamiltonian" begin
    ref = reference_from_kinetic(deuteron, 270.0)
    cell = AbstractElement[Quadrupole("QF", 0.2, 1.2), Drift("D1", 0.5), SBend("B", 1.5, pi / 8), Drift("D2", 0.5),
                           Quadrupole("QD", 0.4, -1.2), Drift("D3", 0.5), SBend("B", 1.5, pi / 8), Drift("D4", 0.5)]
    hrf = 4
    mag = Lattice("ring8", AbstractElement[e for _ in 1:8 for e in cell], ref, hrf)
    C = path_length(mag)
    eta0 = -ring_lie_factor(mag).R[5, 6] / C
    phis = eta0 < 0 ? Float64(pi) : 0.0
    V = 0.2
    cav = RFCavity("RF", V, hrf, phis)
    cn(p) = norm(p, Inf)
    Mlie = element_lie_map(cav, ref, C)
    Mk = [truncate_to(c, MAP_ORDER) for c in rf_kick(collect(Z), cav, ref, C)]
    dthin = maximum(cn(Mlie[i] - Mk[i]) for i in 1:6)
    vprintf("  thin exp(-:G:) vs rf_kick Taylor (order %d): %.3e\n", MAP_ORDER, dthin)
    @test dthin < 1e-13
    z = [1e-3, 2e-4, -1e-3, 1e-4, 0.05, 1e-3]
    zk = rf_kick(z, cav, ref, C); zl = [evaluate(m, z) for m in Mlie]
    vprintf("  thin point: |lie - rf_kick| = %.3e (order-3 truncation)\n", maximum(abs.(zl - zk)))
    for L in (0.1, 0.01, 0.001)
        thick = RFCavity("RFL", V, hrf, phis, L)
        Mt = element_lie_map(thick, ref, C)
        D2 = element_lie_map(Drift("d", L / 2), ref)
        Ms = compose_maps(D2, compose_maps(Mlie, D2))
        d = maximum(cn(Mt[i] - Ms[i]) for i in 1:6)
        zf = exact_track_element(z, thick, ref, C; ds = L / 50)
        zt = [evaluate(m, z) for m in Mt]
        vprintf("  L = %.0e: |thick Lie - D(L/2) kick D(L/2)| = %.3e ; |thick field - thin kick| = %.3e ; symp err thick = %.2e\n",
                L, d, maximum(abs.(zf - zk)), symplectic_error(Mt))
        L == 0.001 && @test d < 1e-6
    end
    vprintf("  symplectic_error thin RF map = %.2e\n", symplectic_error(Mlie))
    @test symplectic_error(Mlie) < 1e-12
    lat, _ = insert_cavity(mag, cav; in_element = "D1")
    F = ring_lie_factor(lat)
    R = F.R
    Qs_map = acos(clamp((R[5, 5] + R[6, 6]) / 2, -1, 1)) / 2pi
    Qs_sm = synchrotron_tune(ref, eta0, V, hrf, phis)
    vprintf("  eta0 = %.5f  phi_s = %.3f  Q_s map = %.6e  Q_s smooth = %.6e  rel = %.2e\n", eta0, phis, Qs_map, Qs_sm, Qs_map / Qs_sm - 1)
    @test isapprox(Qs_map, Qs_sm; rtol = 2e-2)
    fs = build_factors(lat)
    z0 = [1e-4, 0.0, 1e-4, 0.0, 0.1, 1e-4]
    z10 = track(fs, z0; turns = 10)
    w = z0; for _ in 1:10; w = apply_factor(F, w); end
    zex = z0; for _ in 1:10; zex = exact_track_turn(zex, lat; ds = 0.1); end
    vprintf("  10 turns: |element factors - one-turn factor| = %.2e, |Lie - exact field| = %.2e\n", maximum(abs.(z10 - w)), maximum(abs.(z10 - zex)))
    @test maximum(abs.(z10 - zex)) < 1e-6
end

@testset "closed_orbit stops at 1e-12" begin
    # Contracting transverse map: unique fixed point (0.2, 0, 0, 0), independent of δ.
    M = [0.1 + 0.5 * Z[1], 0.3 * Z[2], 0.4 * Z[3], 0.2 * Z[4], Z[5], Z[6]]
    u = closed_orbit(M, 1e-3)
    @test isapprox(u[1], 0.2; atol = 1e-12)
    @test norm(u[2:4]) < 1e-12
end

@testset "implicit midpoint on TaylorN" begin
    h = 0.05
    z1 = Bell.Integrators.implicit_midpoint((z, s) -> [-z[1]], [0.2 + Z[1]], 0.0, h)
    exact = (1 - h / 2) / (1 + h / 2) * 0.2
    @test abs(constant_term(z1[1]) - exact) < 1e-12
end
