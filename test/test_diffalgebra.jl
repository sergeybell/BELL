using Bell.DiffAlgebra

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
