using LinearAlgebra
using Bell.DiffAlgebra, Bell.Integrators

@testset "DATracking (field maps on DA)" begin
    tcoef(p, ex) = Bell.LieAlgebra.TaylorSeries.getcoeff(p, ex)
    maxdiff(M, ML) = maximum(abs(M[i].c[m] - tcoef(ML[i], monomial_exponents(M[1].sp, m))) for i in 1:6, m in 1:nmonomials(M[1].sp))
    ref = reference_from_kinetic(deuteron, 270.0)
    els = [Drift("D", 1.3), Quadrupole("Q", 0.4, 1.7), Sextupole("S", 0.3, 15.0), SBend("B", 1.2, 0.3, 0.2)]
    @testset "agrees with Lie path at order 3" begin
        for e in els
            Mda = element_da_map(e, ref, 3; tol = 1e-14)
            ML = lie_map(hamiltonian(e, ref), elem_length(e))
            @test maxdiff(Mda, ML) < 1e-10
        end
    end
    @testset "analytic drift and quadrupole" begin
        N = 6; sp = map_space(N); X = da_vars(sp); L = 2.0
        M = element_da_map(Drift("D", L), ref, N; tol = 1e-14)
        ps = sqrt((1 + X[6])^2 - X[2]^2 - X[4]^2)
        de = deriv(energy_term(X[6], ref), 6)
        @test coeffnorm(M[1] - (X[1] + L * X[2] / ps)) < 1e-12
        @test coeffnorm(M[3] - (X[3] + L * X[4] / ps)) < 1e-12
        @test coeffnorm(truncate_da(M[5] - (X[5] + L * (de - (1 + X[6]) / ps)), N - 1)) < 1e-11
        @test all(coeffnorm(M[i] - X[i]) < 1e-14 for i in (2, 4, 6))
        k = 1.7; Lq = 0.6; w = sqrt(k) * Lq
        R = linear_part(element_da_map(Quadrupole("Q", Lq, k), ref, N; tol = 1e-14))
        @test norm(R[1:2, 1:2] - [cos(w) sin(w)/sqrt(k); -sqrt(k)*sin(w) cos(w)]) < 1e-12
        @test norm(R[3:4, 3:4] - [cosh(w) sinh(w)/sqrt(k); sqrt(k)*sinh(w) cosh(w)]) < 1e-12
    end
    @testset "Float64 field track matches DA constant part" begin
        q = Quadrupole("Q", 0.5, 2.0)
        z0 = [1e-3, 0.0, -2e-3, 0.0, 0.0, 1e-4]
        zf = exact_track_element(z0, q, ref, 100.0; ds = 0.01, method = DormandPrince853())
        # re-integrate with tight tol for fair comparison
        zf = Bell.Integrators.integrate(element_rhs(q, ref), collect(z0), 0.0, 0.5;
                                         method = DormandPrince853(), tol = 1e-14).z
        M = element_da_map(q, ref, 5; tol = 1e-14)
        ze = [m(z0) for m in M]
        @test maximum(abs(ze[i] - zf[i]) for i in 1:6) < 1e-10
    end
    @testset "symplecticity of DA maps" begin
        for N in (3, 6), L in (0.5, 3.0)
            for e in (Quadrupole("Q", L, 1.2), SBend("B", L, 0.2L, 0.1), Sextupole("S", L, 5.0))
                M = element_da_map(e, ref, N; tol = 1e-14)
                sc = maximum(coeffnorm, M)
                @test da_symplectic_error(M) < 1e-10 * max(1.0, sc^2)
            end
        end
    end
    @testset "compose and one-turn; chromaticity_series" begin
        A = element_da_map(Quadrupole("Q", 0.4, 1.7), ref, 5; tol = 1e-14)
        B = element_da_map(SBend("B", 1.2, 0.3), ref, 5; tol = 1e-14)
        C = compose(A, B)
        @test maximum(coeffnorm(C[i] - A[i](B)) for i in 1:6) < 1e-13
        lat = Lattice("fodo", AbstractElement[Quadrupole("QF", 0.3, 1.5), Drift("D1", 1.0), SBend("B", 1.2, 0.3),
                      Drift("D2", 1.0), Quadrupole("QD", 0.3, -1.5), Sextupole("S", 0.2, 3.0)], ref, 1)
        M = one_turn_map_da(lat, 3; cache = false, tol = 1e-14)
        ML = collect(Z); for e in lat.elements; ML = compose_maps(lie_map(hamiltonian(e, ref), elem_length(e)), ML); end
        @test maxdiff(M, ML) < 1e-9
        @test da_symplectic_error(one_turn_map_da(lat, 5; tol = 1e-14)) < 1e-9
        cell = AbstractElement[Quadrupole("QF", 0.2, 1.2), Drift("D1", 0.5), SBend("B", 1.5, pi / 8), Drift("D2", 0.5),
                               Quadrupole("QD", 0.4, -1.2), Drift("D3", 0.5), SBend("B", 1.5, pi / 8), Drift("D4", 0.5),
                               Sextupole("S", 0.1, 4.0), Quadrupole("QF", 0.2, 1.2)]
        ring = Lattice("ring8", AbstractElement[e for _ in 1:8 for e in cell], ref, 1)
        c5 = chromaticity_series(ring, 5); c7 = chromaticity_series(ring, 7)
        @test length(c5.xi_x) == 4 && length(c7.xi_y) == 6
        @test isapprox(c5.xi_x, c7.xi_x[1:4]; rtol = 1e-8) && isapprox(c5.xi_y, c7.xi_y[1:4]; rtol = 1e-8)
        vprintln("  chromaticity_series ring8: Q = ", (c7.x[1], c7.y[1]), " xi_x = ", c7.xi_x[1:3], " xi_y = ", c7.xi_y[1:3])
    end
    @testset "thin RF kick on Float64 and DA" begin
        Lpath = 100.0
        cav = RFCavity("RF", 0.01, 1, 0.0)
        z0 = [0.0, 0.0, 0.0, 0.0, 0.01, 0.0]
        z1 = rf_kick(z0, cav, ref, Lpath)
        @test z1[1:4] == z0[1:4]
        @test z1[6] != z0[6] || z1[5] != z0[5]
        lat = Lattice("kick", AbstractElement[Drift("D", 1.0), cav, Drift("D2", 1.0)], ref, 1)
        zt = exact_track_turn(z0, lat; ds = 0.5)
        @test length(zt) == 6
    end
end
