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