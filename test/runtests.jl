# test/runtests.jl - tests of the pure-mathematics modules DiffAlgebra and Integrators
#   julia --project=. -e "using Pkg; Pkg.test()"     (or: julia --project=. test/runtests.jl)
using Test, LinearAlgebra
using Bell
using Bell.DiffAlgebra, Bell.Integrators

order_estimate(e1, e2, r = 2) = log(e1 / e2) / log(r)

@testset "DiffAlgebra" begin
    sp = DASpace(3, 10); x, y, z = da_vars(sp)
    @test nmonomials(sp) == binomial(13, 10)
    @testset "series" begin
        sp1 = DASpace(1, 15); t = da_var(sp1, 1)
        @test maximum(abs(coefficient(exp(t), [k]) - 1 / factorial(big(k))) for k in 0:15) < 1e-16
        @test maximum(abs(coefficient(sin(t), [k]) - (isodd(k) ? (-1)^((k - 1) ÷ 2) / factorial(big(k)) : 0)) for k in 0:15) < 1e-16
        @test maximum(abs(coefficient(log(1 + t), [k]) - (-1)^(k + 1) / k) for k in 1:15) < 1e-15
        @test maximum(abs(coefficient(atan(t), [k]) - (isodd(k) ? (-1)^((k - 1) ÷ 2) / k : 0)) for k in 0:15) < 1e-15
        @test maximum(abs(coefficient(sqrt(1 + t), [k]) - binomial(big(1) / 2, k)) for k in 0:15) < 1e-15
    end
    @testset "identities" begin
        a = 0.3 + x + 0.5y * z - 0.2x^2; b = 1.7 + y - x * z
        @test coeffnorm(sin(a)^2 + cos(a)^2 - 1) < 1e-14
        @test coeffnorm(exp(log(b)) - b) < 1e-13
        @test coeffnorm(sqrt(b)^2 - b) < 1e-14
        @test coeffnorm(exp(a + b) - exp(a) * exp(b)) < 1e-12
        @test coeffnorm(tan(atan(a)) - a) < 1e-13
        @test coeffnorm(cosh(a)^2 - sinh(a)^2 - 1) < 1e-13
        @test coeffnorm(b * inv(b) - 1) < 1e-14
        @test coeffnorm(b^2.5 - b^2 * sqrt(b)) < 1e-13
        @test coeffnorm(a / b * b - a) < 1e-14
    end
    @testset "calculus and evaluation" begin
        p = 3x^3 * y + 2y * z^2 - x
        @test coeffnorm(deriv(p, 1) - (9x^2 * y - 1)) < 1e-15
        @test coeffnorm(deriv(antideriv(p, 3), 3) - p) < 1e-15
        @test abs(p([0.2, -0.7, 1.1]) - (3 * 0.2^3 * -0.7 + 2 * -0.7 * 1.1^2 - 0.2)) < 1e-15
        sp1 = DASpace(1, 12); t = da_var(sp1, 1)
        @test abs(exp(sin(t))([0.01]) - exp(sin(0.01))) < 1e-15
        @test linear_part([x + y^2, 2y - x * z, z]) == [1.0 0 0; 0 2 0; 0 0 1]
    end
    @testset "composition and inverse" begin
        M = [x + 0.3y^2 - 0.1x * z + 0.05sin(y)^3, y - 0.2x^2 + 0.1z^3, 0.5z + x * y]
        Mi = inverse_map(M)
        I1 = compose(M, Mi); I2 = compose(Mi, M)
        @test maximum(coeffnorm(I1[i] - da_vars(sp)[i]) for i in 1:3) < 1e-14
        @test maximum(coeffnorm(I2[i] - da_vars(sp)[i]) for i in 1:3) < 1e-14
        N = [sin(x) + y, z^2 + x, cos(y) - 1]
        f = exp(x) * y + z
        @test coeffnorm(compose([f], N)[1] - (exp(N[1]) * N[2] + N[3])) < 1e-14
    end
end

@testset "Integrators" begin
    ho(z, s) = [z[2], -z[1]]
    err(m, n) = norm(integrate(ho, [1.0, 0.0], 0.0, 2pi; method = m, nsteps = n).z - [1.0, 0.0])
    @testset "fixed-step orders" begin
        for (m, n, o) in [(RK4(), 40, 4), (DormandPrince5(), 20, 5), (DormandPrince853(), 6, 8),
                          (GraggBulirschStoer(4), 6, 8), (GraggBulirschStoer(5), 4, 10)]
            p = order_estimate(err(m, n), err(m, 2n))
            println("  ", m, ": order ", round(p, digits = 2))
            @test abs(p - o) < 0.6
        end
        flows = SplitFlows((z, h) -> [z[1] + h * z[2], z[2]], (z, h) -> [z[1], z[2] - h * z[1]])
        for o in (2, 4, 6)
            e(n) = norm(integrate(flows, [1.0, 0.0], 0.0, 2pi; method = Yoshida(o), nsteps = n).z - [1.0, 0.0])
            p = order_estimate(e(40), e(80)); println("  Yoshida($o): order ", round(p, digits = 2))
            @test abs(p - o) < 0.6
        end
    end
    @testset "adaptive: Kepler e = 0.6, one period" begin
        ecc = 0.6; z0 = [1 - ecc, 0.0, 0.0, sqrt((1 + ecc) / (1 - ecc))]
        kep(z, s) = (r3 = (z[1]^2 + z[2]^2)^1.5; [z[3], z[4], -z[1] / r3, -z[2] / r3])
        for (m, tol, lim) in [(DormandPrince5(), 1e-12, 1e-9), (DormandPrince853(), 1e-13, 1e-10), (GraggBulirschStoer(8), 1e-14, 1e-11)]
            r = integrate(kep, z0, 0.0, 2pi; method = m, tol = tol)
            println("  Kepler ", m, ": err ", norm(r.z - z0), " steps ", r.nsteps, " nfev ", r.nfev)
            @test norm(r.z - z0) < lim
        end
    end
    @testset "DA initial condition: x' = x^2" begin
        # x(t) = x0/(1 - t x0), x0 = a + e:  coefficient of e^k = t^(k-1)/(1 - t a)^(k+1), k >= 1
        N = 10; sp = DASpace(1, N); a = 1.0; t1 = 0.5
        x0 = da_var(sp, 1; at = a)
        for (m, tol) in [(DormandPrince853(), 1e-15), (GraggBulirschStoer(8), 1e-15)]
            r = integrate((z, s) -> [z[1]^2], [x0], 0.0, t1; method = m, tol = tol)
            ex = [k == 0 ? a / (1 - t1 * a) : t1^(k - 1) / (1 - t1 * a)^(k + 1) for k in 0:N]
            rel = maximum(abs(coefficient(r.z[1], [k]) - ex[k+1]) / abs(ex[k+1]) for k in 0:N)
            println("  DA ODE ", m, ": max rel coeff err ", rel)
            @test rel < 1e-12
        end
    end
end
@testset "LieFlow (flow-operator element maps)" begin
    tcoef(p, ex) = Bell.LieAlgebra.TaylorSeries.getcoeff(p, ex)
    maxdiff(M, ML) = maximum(abs(M[i].c[m] - tcoef(ML[i], monomial_exponents(M[1].sp, m))) for i in 1:6, m in 1:nmonomials(M[1].sp))
    ref = reference_from_kinetic(DEUTERON, 270.0)
    els = [Drift("D", 1.3), Quadrupole("Q", 0.4, 1.7), Sextupole("S", 0.3, 15.0), SBend("B", 1.2, 0.3, 0.2)]
    @testset "= Dragt-Finn/Lie route at order 3" begin
        for e in els
            @test maxdiff(lie_flow_map(e, ref, 3), lie_map(hamiltonian(e, ref), elem_length(e))) < 1e-13
        end
    end
    @testset "analytic drift and quadrupole" begin
        N = 6; sp = flow_space(N); X = da_vars(sp); L = 2.0
        M = lie_flow_map(Drift("D", L), ref, N)
        ps = sqrt((1 + X[6])^2 - X[2]^2 - X[4]^2)
        de = deriv(energy_term(X[6], ref), 6)
        @test coeffnorm(M[1] - (X[1] + L * X[2] / ps)) < 1e-14
        @test coeffnorm(M[3] - (X[3] + L * X[4] / ps)) < 1e-14
        @test coeffnorm(truncate_da(M[5] - (X[5] + L * (de - (1 + X[6]) / ps)), N - 1)) < 1e-13
        @test all(coeffnorm(M[i] - X[i]) == 0 for i in (2, 4, 6))
        k = 1.7; Lq = 0.6; w = sqrt(k) * Lq
        R = linear_part(lie_flow_map(Quadrupole("Q", Lq, k), ref, N))
        @test norm(R[1:2, 1:2] - [cos(w) sin(w)/sqrt(k); -sqrt(k)*sin(w) cos(w)]) < 1e-14
        @test norm(R[3:4, 3:4] - [cosh(w) sinh(w)/sqrt(k); sqrt(k)*sinh(w) cosh(w)]) < 1e-14
    end
    @testset "force form z' = S grad H (DOP853)" begin
        N = 5; q = Quadrupole("Q", 0.5, 2.0); M = lie_flow_map(q, ref, N)
        F = hamiltonian_vector_field(hamiltonian_da(q, ref, flow_space(N + 1)))
        z0 = da_vars(flow_space(N + 1))
        r = Bell.Integrators.integrate((z, s) -> [f(z) for f in F], z0, 0.0, 0.5; method = DormandPrince853(), tol = 1e-14).z
        @test maximum(coeffnorm(M[i] - change_order(r[i], flow_space(N))) for i in 1:6) < 1e-11
    end
    @testset "symplecticity vs order and length (scaling and squaring)" begin
        for N in (3, 6), L in (0.5, 3.0, 10.0)
            for e in (Quadrupole("Q", L, 1.2), SBend("B", L, 0.2L, 0.1), Sextupole("S", L, 5.0))
                M = lie_flow_map(e, ref, N)
                sc = maximum(coeffnorm, M)
                @test da_symplectic_error(M) < 1e-13 * sc^2
            end
        end
    end
    @testset "compose fast path = generic evaluation; ring map" begin
        sp = flow_space(5); X = da_vars(sp)
        A = lie_flow_map(Quadrupole("Q", 0.4, 1.7), ref, 5); B = lie_flow_map(SBend("B", 1.2, 0.3), ref, 5)
        C = compose(A, B)
        @test maximum(coeffnorm(C[i] - A[i](B)) for i in 1:6) < 1e-14
        lat = Lattice("fodo", AbstractElement[Quadrupole("QF", 0.3, 1.5), Drift("D1", 1.0), SBend("B", 1.2, 0.3),
                      Drift("D2", 1.0), Quadrupole("QD", 0.3, -1.5), Sextupole("S", 0.2, 3.0)], ref, 1)
        M = one_turn_map_flow(lat, 3; cache = false)
        ML = collect(Z); for e in lat.elements; ML = compose_maps(lie_map(hamiltonian(e, ref), elem_length(e)), ML); end
        @test maxdiff(M, ML) < 1e-12
        @test da_symplectic_error(one_turn_map_flow(lat, 6)) < 1e-11
    end
    @testset "chromaticity_series (exact, no finite differences)" begin
        cell = AbstractElement[Quadrupole("QF", 0.2, 1.2), Drift("D1", 0.5), SBend("B", 1.5, pi / 8), Drift("D2", 0.5),
                               Quadrupole("QD", 0.4, -1.2), Drift("D3", 0.5), SBend("B", 1.5, pi / 8), Drift("D4", 0.5),
                               Sextupole("S", 0.1, 4.0), Quadrupole("QF", 0.2, 1.2)]
        ring = Lattice("ring8", AbstractElement[e for _ in 1:8 for e in cell], ref, 1)
        c5 = chromaticity_series(ring, 5); c7 = chromaticity_series(ring, 7)
        @test length(c5.xi_x) == 4 && length(c7.xi_y) == 6
        @test isapprox(c5.xi_x, c7.xi_x[1:4]; rtol = 1e-10) && isapprox(c5.xi_y, c7.xi_y[1:4]; rtol = 1e-10)
        MDF = factor_to_map(ring_lie_factor(ring)); sp = flow_space(3)
        MD = [DA(sp, [tcoef(MDF[i], monomial_exponents(sp, m)) for m in 1:nmonomials(sp)]) for i in 1:6]
        cd = chromaticity_series(MD)
        @test isapprox(cd.xi_x, c7.xi_x[1:2]; rtol = 1e-10) && isapprox(cd.xi_y, c7.xi_y[1:2]; rtol = 1e-10)
        q0 = tunes_at_delta(MDF, 0.0)
        @test abs(q0[1] - c7.x[1]) < 1e-12 && abs(q0[2] - c7.y[1]) < 1e-12
        h = chromaticity_orders(MDF; dd = 1e-5)
        @test abs(h.x[2] - c7.xi_x[1]) < 1e-5 * abs(c7.xi_x[1]) && abs(h.y[2] - c7.xi_y[1]) < 1e-5 * abs(c7.xi_y[1])
        println("  chromaticity_series ring8: Q = ", (c7.x[1], c7.y[1]), " xi_x = ", c7.xi_x[1:3], " xi_y = ", c7.xi_y[1:3])
    end
end

@testset "convert_optim close_ring" begin
    opt = tempname() * ".opt"
    write(opt, "OptiM\nEnergy[MeV]=270   Mass[MeV]=1875.612945\nbegin lattice. Number of periods=1\n" *
               join(fill("O1 BM", 16), " ") * "\nend lattice\nbegin list\nO1 L[cm]=10\n" *
               "BM L[cm]=532.94148 B[kG]=2.56106599 G[kG/cm]=0\nend list\n")
    for (cr, tol) in ((true, 1e-15), (false, Inf))
        bell = tempname() * ".bell"
        m = convert_optim(opt, bell; close_ring = cr)
        lat = load_structure(bell)
        tot = sum(e -> e isa SBend ? e.angle : 0.0, lat.elements)
        cr ? (@test abs(tot - 2pi) < 1e-13 && abs(lat.elements[2].angle - pi / 8) < 1e-15) :
             (@test 1e-9 < abs(tot - 2pi) < 1e-6)
        @test any(w -> occursin("close_ring", w), m.warnings) == cr
    end
end
