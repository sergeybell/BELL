# Magnet flow dz/ds = S ∇H against the field equations kept for comparison.
using Bell.DiffAlgebra, Bell.Integrators

_coeff_gap(A, B) = maximum(coeffnorm(A[i] - B[i]) for i in eachindex(A))

@testset "element_flow" begin
    ref = reference_from_kinetic(deuteron, 270.0)
    magnets = (
        Drift("D", 1.3),
        Quadrupole("Q", 0.4, 1.7),
        SBend("B", 1.2, 0.3, 0.2),
    )

    @testset "marker and zero length are a copy" begin
        z0 = [0.1, -0.2, 0.3, -0.4, 0.5, -0.6]
        for e in (Marker("M"), Drift("Z", 0.0), Quadrupole("Q0", 0.0, 2.0), SBend("B0", 0.0, 0.0))
            zf = element_flow(e, ref, z0; method = RK4(), nsteps = 4)
            @test zf == z0
            @test !(zf === z0)
            zf[1] = 0.0
            @test z0[1] == 0.1
        end
        @test_throws ArgumentError element_flow(RFCavity("RF", 0.01, 1, 0.0), ref, z0)
    end

    @testset "DA orders 5 and 6 match the field equations" begin
        method = RK4()
        nsteps = 4
        for N in (5, 6), e in magnets
            sp = map_space(N)
            z0 = collect(da_vars(sp))
            Mh = element_flow(e, ref, z0; method = method, nsteps = nsteps)
            Mf = integrate_fixed(element_rhs(e, ref), collect(da_vars(sp)), 0.0, elem_length(e),
                                 nsteps; method = method).z
            @test _coeff_gap(Mh, Mf) < 1e-8
        end
    end

    @testset "Float64 track matches the field equations" begin
        z0 = [1e-3, -2e-4, 1.5e-3, 1e-4, -1e-3, 2e-3]
        method = ImplicitYoshida(4)
        ds = 0.05
        for e in magnets
            zf = exact_track_element(z0, e, ref, 20.0; ds = ds, method = method)
            n = e isa Drift ? 1 : max(1, ceil(Int, elem_length(e) / ds))
            zg = integrate_fixed(element_rhs(e, ref), collect(z0), 0.0, elem_length(e), n;
                                 method = method).z
            @test maximum(abs, zf .- zg) < 1e-8
            zh = element_flow(e, ref, z0; method = method, ds = ds)
            nh = max(1, ceil(Int, elem_length(e) / ds))
            zhg = integrate_fixed(element_rhs(e, ref), collect(z0), 0.0, elem_length(e), nh;
                                  method = method).z
            @test maximum(abs, zh .- zhg) < 1e-8
        end
    end

    @testset "s-dependent Hamiltonian" begin
        z0 = [0.2, 0.0, 0.0, 0.0, 0.0, 0.0]
        H(z, s) = s * z[1]^2 / 2
        zf = element_flow(H, z0; s0 = 0.0, s1 = 1.0, method = RK4(), nsteps = 4)
        @test zf[1] ≈ 0.2 atol = 1e-12
        @test zf[2] ≈ -0.1 atol = 1e-12
        @test zf[3:6] ≈ zeros(4) atol = 1e-12

        e = Quadrupole("Q", 0.3, 0.5)
        z1 = [1e-3, -2e-4, 1.2e-3, 1e-4, 0.0, 1e-3]
        Ha(z, s) = hamiltonian(e, ref, z)
        @test element_flow(Ha, z1; s0 = 0.0, s1 = e.L, method = RK4(), nsteps = 6) ≈
              element_flow(e, ref, z1; method = RK4(), nsteps = 6) atol = 1e-10
    end

    @testset "DA map symplecticity at orders 4 and 5" begin
        for N in (4, 5), e in magnets
            M = element_da_map(e, ref, N; method = DormandPrince853(), tol = 1e-12)
            sc = maximum(coeffnorm, M)
            @test symplectic_error(M) < 1e-8 * max(1.0, sc^2)
        end
    end
end
