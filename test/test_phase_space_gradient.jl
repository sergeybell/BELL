# Dual phase-space gradient against a central difference, TaylorN, and DA.
using Bell.DiffAlgebra: DASpace, Dual, coeffnorm, constant_part, da_const, da_vars, deriv,
                        identity_map, map_vector_field, phase_space_gradient, symplectic_gradient,
                        truncate_da, vector_field_on_identity
using TaylorSeries: constant_term, differentiate

_grad_ref = reference_from_kinetic(deuteron, 270.0)

function quadratic(z)
    return 0.5 * z[1]^2 - 2 * z[1] * z[3] + z[2]^2 + 0.25 * z[4]^2 - z[5] + 3 * z[6]
end

function quadratic_gradient(z)
    return [z[1] - 2 * z[3], 2 * z[2], -2 * z[1], 0.5 * z[4], -1.0, 3.0]
end

function quadrupole_like(z)
    x, px, y, py, _, d = z
    en = energy_term(d, _grad_ref)
    ps = sqrt((1 + d)^2 - px^2 - py^2)
    return en - ps + (1.7 / 2) * (x^2 - y^2)
end

function central_gradient(f, z; step = 1.0e-6)
    g = similar(z, Float64)
    for j in eachindex(z)
        zp = collect(Float64, z)
        zm = collect(Float64, z)
        zp[j] += step
        zm[j] -= step
        g[j] = (f(zp) - f(zm)) / (2 * step)
    end
    return g
end

@testset "Phase-space gradient" begin
    z0 = [1.0e-3, -2.0e-4, 3.0e-3, 1.5e-4, -5.0e-3, 2.0e-3]

    @testset "Dual arithmetic" begin
        x = Dual(2.0, 1.0)
        @test x + 3.0 == Dual(5.0, 1.0)
        @test 3.0 + x == Dual(5.0, 1.0)
        @test x - 3.0 == Dual(-1.0, 1.0)
        @test 3.0 - x == Dual(1.0, -1.0)
        @test -x == Dual(-2.0, -1.0)
        @test 4 * x == Dual(8.0, 4.0)
        @test x * Dual(3.0, 1.0) == Dual(6.0, 2.0 + 3.0)
        @test x / 4 == Dual(0.5, 0.25)
        @test Dual(3.0, 1.0) / Dual(2.0, 1.0) == Dual(1.5, -0.25)
        @test inv(Dual(4.0, 1.0)) == Dual(0.25, -1 / 16)
        @test sqrt(Dual(4.0, 1.0)) == Dual(2.0, 0.25)
        @test x^2 == Dual(4.0, 4.0)
        @test x^3 == Dual(8.0, 12.0)
        @test x^0 == Dual(1.0, 0.0)
        @test x^-1 == inv(x)
        @test exp(Dual(0.5, 1.0)) ≈ Dual(exp(0.5), exp(0.5))
        @test log(Dual(2.0, 1.0)) ≈ Dual(log(2.0), 0.5)
        @test sin(Dual(0.4, 1.0)) ≈ Dual(sin(0.4), cos(0.4))
        @test cos(Dual(0.4, 1.0)) ≈ Dual(cos(0.4), -sin(0.4))
        @test phase_space_gradient(z -> 1.0, [0.2, -0.3]) == [0.0, 0.0]
    end

    @testset "Float64 quadratic" begin
        g = phase_space_gradient(quadratic, z0)
        @test g ≈ quadratic_gradient(z0) atol = 1e-12
        @test g ≈ central_gradient(quadratic, z0) atol = 1e-8
        gS = symplectic_gradient(quadratic, z0)
        @test gS ≈ symplectic_gradient(quadratic, z0; S = SYMPLECTIC_S) atol = 1e-12
        v = [sum(SYMPLECTIC_S[i, j] * g[j] for j in 1:6) for i in 1:6]
        @test symplectic_gradient(quadratic, z0; S = SYMPLECTIC_S) ≈ v atol = 1e-12
    end

    @testset "Float64 quadrupole-like" begin
        g = phase_space_gradient(quadrupole_like, z0)
        gfd = central_gradient(quadrupole_like, z0)
        @test maximum(abs, g) > 1e-6
        @test g[5] == 0.0
        @test g ≈ gfd atol = 1e-8
        v = symplectic_gradient(quadrupole_like, z0; S = SYMPLECTIC_S)
        vfd = [sum(SYMPLECTIC_S[i, j] * gfd[j] for j in 1:6) for i in 1:6]
        @test v ≈ vfd atol = 1e-8
    end

    @testset "default pairing and a passed S" begin
        for k in 1:6
            col = let k = k
                symplectic_gradient(z -> z[k], zeros(6))
            end
            @test col == SYMPLECTIC_S[:, k]
        end
        S = zeros(2, 2)
        S[1, 1] = 3.0
        S[2, 2] = 4.0
        @test symplectic_gradient(z -> z[1] + 2 * z[2], [0.2, 0.3]; S = S) == [3.0, 8.0]
        @test_throws ArgumentError symplectic_gradient(z -> z[1], [0.0])
        @test_throws ArgumentError symplectic_gradient(quadratic, z0; S = zeros(5, 5))
    end

    @testset "TaylorN" begin
        g = phase_space_gradient(quadratic, Z)
        H = quadratic(Z)
        for j in 1:6
            @test maxabs(g[j] - differentiate(H, j)) < 1e-12
        end
        gq = phase_space_gradient(quadrupole_like, Z)
        Hq = quadrupole_like(Z)
        g0 = phase_space_gradient(quadrupole_like, zeros(6))
        for j in 1:6
            @test maxabs(truncate_to(gq[j] - differentiate(Hq, j), JET_ORDER - 1)) < 1e-10
            @test abs(constant_term(gq[j]) - g0[j]) < 1e-12
        end
    end

    @testset "DA on identity" begin
        sp = DASpace(6, 4)
        z = identity_map(sp)
        g = phase_space_gradient(quadratic, z)
        H = quadratic(z)
        for j in 1:6
            @test coeffnorm(g[j] - deriv(H, j)) < 1e-12
        end
        gq = phase_space_gradient(quadrupole_like, z)
        Hq = quadrupole_like(z)
        for j in 1:6
            @test coeffnorm(truncate_da(gq[j] - deriv(Hq, j), sp.no - 1)) < 1e-10
        end
        V = map_vector_field(quadrupole_like, z; S = SYMPLECTIC_S)
        Vid = vector_field_on_identity(quadrupole_like, sp; S = SYMPLECTIC_S)
        Vdual = [sum(SYMPLECTIC_S[i, j] * gq[j] for j in 1:6) for i in 1:6]
        for i in 1:6
            @test V[i] == Vid[i]
            @test coeffnorm(truncate_da(V[i] - Vdual[i], sp.no - 1)) < 1e-10
        end
    end

    @testset "DA at a shifted jet" begin
        sp = DASpace(6, 4)
        z = da_vars(sp; at = z0)
        g = phase_space_gradient(quadrupole_like, z)
        g0 = phase_space_gradient(quadrupole_like, z0)
        H = quadrupole_like(z)
        for j in 1:6
            @test abs(constant_part(g[j]) - g0[j]) < 1e-12
            @test coeffnorm(truncate_da(g[j] - deriv(H, j), sp.no - 1)) < 1e-10
        end
        zbad = identity_map(DASpace(2, 3))
        zbad[1] = zbad[1]^2
        @test_throws ArgumentError map_vector_field(z -> z[1], zbad)
    end

    @testset "coordinate on the identity" begin
        sp = DASpace(2, 2)
        V = vector_field_on_identity(z -> z[1], sp)
        @test coeffnorm(V[1]) == 0.0
        @test V[2] == da_const(sp, -1.0)
        @test all(v -> coeffnorm(v) == 0.0, vector_field_on_identity(z -> 1.0, sp))
    end
end
