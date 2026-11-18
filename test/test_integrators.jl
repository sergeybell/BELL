using LinearAlgebra
using Bell.DiffAlgebra, Bell.Integrators

order_estimate(e1, e2, r = 2) = log(e1 / e2) / log(r)

@testset "Integrators" begin
    ho(z, s) = [z[2], -z[1]]
    err(m, n) = norm(integrate(ho, [1.0, 0.0], 0.0, 2pi; method = m, nsteps = n).z - [1.0, 0.0])
    @testset "fixed-step orders" begin
        for (m, n, o) in [(RK4(), 40, 4), (DormandPrince5(), 20, 5), (DormandPrince853(), 6, 8),
                          (GraggBulirschStoer(4), 6, 8), (GraggBulirschStoer(5), 4, 10)]
            p = order_estimate(err(m, n), err(m, 2n))
            vprintln("  ", m, ": order ", round(p, digits = 2))
            @test abs(p - o) < 0.6
        end
        flows = SplitFlows((z, h) -> [z[1] + h * z[2], z[2]], (z, h) -> [z[1], z[2] - h * z[1]])
        for o in (2, 4, 6)
            e(n) = norm(integrate(flows, [1.0, 0.0], 0.0, 2pi; method = Yoshida(o), nsteps = n).z - [1.0, 0.0])
            p = order_estimate(e(40), e(80)); vprintln("  Yoshida($o): order ", round(p, digits = 2))
            @test abs(p - o) < 0.6
        end
    end
    @testset "adaptive: Kepler e = 0.6, one period" begin
        ecc = 0.6; z0 = [1 - ecc, 0.0, 0.0, sqrt((1 + ecc) / (1 - ecc))]
        kep(z, s) = (r3 = (z[1]^2 + z[2]^2)^1.5; [z[3], z[4], -z[1] / r3, -z[2] / r3])
        for (m, tol, lim) in [(DormandPrince5(), 1e-12, 1e-9), (DormandPrince853(), 1e-13, 1e-10), (GraggBulirschStoer(8), 1e-14, 1e-11)]
            r = integrate(kep, z0, 0.0, 2pi; method = m, tol = tol)
            vprintln("  Kepler ", m, ": err ", norm(r.z - z0), " steps ", r.nsteps, " nfev ", r.nfev)
            @test norm(r.z - z0) < lim
        end
    end
    # Order-10 DA ODE: full suite only.
    if full_tests()
        @testset "DA initial condition: x' = x^2" begin
            N = 10; sp = DASpace(1, N); a = 1.0; t1 = 0.5
            x0 = da_var(sp, 1; at = a)
            for (m, tol) in [(DormandPrince853(), 1e-15), (GraggBulirschStoer(8), 1e-15)]
                r = integrate((z, s) -> [z[1]^2], [x0], 0.0, t1; method = m, tol = tol)
                ex = [k == 0 ? a / (1 - t1 * a) : t1^(k - 1) / (1 - t1 * a)^(k + 1) for k in 0:N]
                rel = maximum(abs(coefficient(r.z[1], [k]) - ex[k+1]) / abs(ex[k+1]) for k in 0:N)
                vprintln("  DA ODE ", m, ": max rel coeff err ", rel)
                @test rel < 1e-12
            end
        end
    end
end

@testset "pendulum on every integrator" begin
    pend(z, s) = [z[2], -sin(z[1])]
    energy(z) = z[2]^2 / 2 - cos(z[1])
    theta = 0.3
    z0 = [theta, 0.0]
    E0 = energy(z0)
    period = 2pi * (1 + theta^2 / 16 + 11 * theta^4 / 3072)
    flows = SplitFlows((z, h) -> [z[1] + h * z[2], z[2]],
                       (z, h) -> [z[1], z[2] - h * sin(z[1])])
    function check(name, z, nfev)
        dE = abs(energy(z) - E0)
        dist = hypot(z[1] - z0[1], z[2] - z0[2])
        vprintln("  pendulum ", name, ": dE ", dE, " dist ", dist, " nfev ", nfev)
        @test dE < 1e-6
        @test dist < 2e-4
    end
    for (name, method, n) in (("RK4", RK4(), 200),
                              ("ImplicitMidpoint", ImplicitMidpoint(), 2000),
                              ("ImplicitYoshida(4)", ImplicitYoshida(4), 200))
        r = integrate(pend, z0, 0.0, period; method, nsteps = n)
        check(name, r.z, r.nfev)
    end
    r = integrate(flows, z0, 0.0, period; method = Yoshida(4), nsteps = 200)
    check("Yoshida(4) split", r.z, r.nfev)
    for (name, method, tol) in (("DormandPrince5", DormandPrince5(), 1e-12),
                                ("DormandPrince853", DormandPrince853(), 1e-13),
                                ("GraggBulirschStoer", GraggBulirschStoer(6), 1e-12))
        r = integrate(pend, z0, 0.0, period; method, tol)
        check(name, r.z, r.nfev)
    end
end
