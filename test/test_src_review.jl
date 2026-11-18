using Bell.DiffAlgebra, Bell.Integrators

@testset "is_adaptive" begin
    @test is_adaptive(DormandPrince5())
    @test is_adaptive(DormandPrince853())
    @test is_adaptive(GraggBulirschStoer())
    @test !is_adaptive(RK4())
    @test !is_adaptive(ImplicitMidpoint())
    @test !is_adaptive(ImplicitYoshida(4))
    @test !is_adaptive(Leapfrog())
    @test !is_adaptive(Yoshida(4))
    @test_throws ErrorException integrate((z, s) -> [-z[1]], [1.0], 0.0, 1.0; method = RK4())
end

@testset "implicit midpoint stops on errnorm" begin
    sp = DASpace(1, 4)
    z0 = [da_var(sp, 1; at = 0.1)]
    nfev = Ref(0)
    f(z, s) = (nfev[] += 1; [-z[1]])
    h = 0.2
    z = implicit_midpoint(f, z0, 0.0, h)
    exact = (1 - h / 2) / (1 + h / 2) * 0.1
    @test abs(constant_part(z[1]) - exact) < 1e-12
    @test nfev[] < 40
end

@testset "DiffAlgebra coefficient type and missing monomial" begin
    sp = DASpace(1, 3)
    x32 = da_var(sp, 1; at = 0.5f0)
    @test eltype(x32.c) == Float32
    @test eltype(exp(x32).c) == Float32
    @test eltype(sin(x32).c) == Float32
    @test_throws DomainError abs(da_var(sp, 1))
    @test constant_part(abs(da_const(sp, 0.0))) == 0.0
    @test constant_part(abs(da_const(sp, -2.0))) == 2.0

    a = da_const(DASpace(2, 2), 1.0)
    @test_throws ArgumentError coefficient(a, [-1, 2])
    wide = zeros(Int, 16)
    wide[1] = -1
    wide[2] = 2
    @test_throws ArgumentError coefficient(da_const(DASpace(16, 2), 1.0), wide)
end

@testset "optimizer guards" begin
    @test_throws ErrorException RandomStart(0; lower = [0.0], upper = [1.0])
    @test_throws ErrorException GridSearch([[0.0, 0.0], [1.0]])
end
