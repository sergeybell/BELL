# One magnetic formulation: element_field == S ∇H of hamiltonian(e, ref, z).
using Random
using TaylorSeries
using Bell.DiffAlgebra

tcoef(p, ex) = getcoeff(p, ex)

function field_from_da_hamiltonian(H)
    g = [deriv(H, j) for j in 1:6]
    return [sum(SYMPLECTIC_S[i, j] * g[j] for j in 1:6) for i in 1:6]
end

@testset "Physics consistency" begin
    ref = reference_from_kinetic(deuteron, 270.0)
    els = [
        Drift("D", 1.3),
        Quadrupole("Q", 0.4, 1.7),
        Sextupole("S", 0.3, 15.0),
        SBend("B", 1.2, 0.3, 0.2),
    ]

    @testset "DA S∇H matches element_field through order $N" for N in (5, 6)
        # deriv drops degree N+1, so H is built one order higher than the jet we keep.
        sp = DASpace(6, N + 1)
        z = da_vars(sp)
        for e in els
            H = hamiltonian(e, ref, z)
            Fh = field_from_da_hamiltonian(H)
            h, k1, k2 = magnetic_parameters(e)
            Ff = element_field(z, ref, h, k1, k2)
            for i in 1:6
                @test coeffnorm(truncate_da(Ff[i] - Fh[i], N)) < 1e-12
            end
        end
    end

    @testset "TaylorN jet matches element_field through order 3" begin
        sp = DASpace(6, 3)
        z = da_vars(sp)
        for e in els
            Ft = hamiltonian_field(hamiltonian(e, ref))
            h, k1, k2 = magnetic_parameters(e)
            Ff = element_field(z, ref, h, k1, k2)
            for i in 1:6
                err = maximum(abs(Ff[i].c[m] - tcoef(Ft[i], monomial_exponents(sp, m)))
                              for m in 1:nmonomials(sp))
                @test err < 1e-12
            end
        end
    end

    @testset "Float64 central difference" begin
        z0 = 1e-3 .* randn(MersenneTwister(270), 6)
        eps = 1e-6
        for e in els
            g = zeros(6)
            for j in 1:6
                zp = copy(z0)
                zm = copy(z0)
                zp[j] += eps
                zm[j] -= eps
                g[j] = (hamiltonian(e, ref, zp) - hamiltonian(e, ref, zm)) / (2eps)
            end
            Fnum = [sum(SYMPLECTIC_S[i, j] * g[j] for j in 1:6) for i in 1:6]
            h, k1, k2 = magnetic_parameters(e)
            Ff = element_field(z0, ref, h, k1, k2)
            @test maximum(abs, Ff) > 1e-8
            @test maximum(abs, Fnum .- Ff) < 1e-7
        end
    end

    @testset "scalar type and Lie polynomial" begin
        z0 = zeros(6)
        for e in els
            @test hamiltonian(e, ref, z0) isa Float64
            @test hamiltonian(e, ref, z0) ≈ -1
            H = @test_nowarn hamiltonian(e, ref)
            @test H isa TaylorN
            @test constant_term(H) == 0
        end
        @test hamiltonian(Marker("M"), ref, z0) == 0
        @test iszero(hamiltonian(Marker("M"), ref))
        sp = DASpace(6, 2)
        @test coeffnorm(hamiltonian(Marker("M"), ref, da_vars(sp))) == 0
        thin = RFCavity("thin", 0.01, 1, 0.0)
        @test_throws ErrorException hamiltonian(thin, ref)
        @test_throws ErrorException hamiltonian(thin, ref, 20.0)
        @test_throws ErrorException hamiltonian(thin, ref, z0)
        thick = RFCavity("RF", 0.2, 1, 0.0, 0.5)
        @test (@test_nowarn hamiltonian(thick, ref, 40.0)) isa TaylorN
        @test hamiltonian(thick, ref, 40.0, z0) isa Float64
        @test hamiltonian(thick, ref, 40.0, da_vars(DASpace(6, 3))) isa DA
    end
end
