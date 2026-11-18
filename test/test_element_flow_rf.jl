# Thick RF is the flow of H. Thin RF is the kick, checked against the flow of G.
using Bell.DiffAlgebra, Bell.Integrators

_rf_coeff_gap(A, B) = maximum(coeffnorm(A[i] - B[i]) for i in eachindex(A))

# The previous production step for a thick cavity: drift / kick / drift.
function _strang_thick(z, e, ref, path_length; ds, method)
    n = max(1, ceil(Int, e.L / ds))
    hs = e.L / n
    kick = RFCavity(e.name, e.voltage_MV / n, e.harmonic, e.phi_s, 0.0)
    fd = element_rhs(Drift(e.name, hs), ref)
    w = collect(z)
    for _ in 1:n
        w = integrate_fixed(fd, w, 0.0, hs / 2, 1; method).z
        w = rf_kick(w, kick, ref, path_length)
        w = integrate_fixed(fd, w, 0.0, hs / 2, 1; method).z
    end
    return w
end

@testset "element_flow RF" begin
    ref = reference_from_kinetic(deuteron, 270.0)
    C = 20.0
    points = (
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1e-3, -2e-4, 1.5e-3, 1e-4, 0.03, 2e-3],
        [0.0, 0.0, 0.0, 0.0, -0.08, -4e-3],
        [2e-3, 1e-4, -1e-3, 2e-4, 0.12, 1e-2],
    )

    @testset "thin kick matches the flow of the generator" begin
        thin = RFCavity("RF", 0.3, 2, 0.0)
        ex = try
            element_flow(thin, ref, points[1])
            nothing
        catch e
            e
        end
        @test ex isa ArgumentError
        @test occursin("rf_kick", ex.msg)
        @test occursin("rf_potential", ex.msg)
        for (phi, V) in ((0.0, 0.3), (Float64(pi), 0.3), (pi / 2, 0.5))
            cav = RFCavity("RF", V, 2, phi)
            for z0 in points
                zk = rf_kick(z0, cav, ref, C)
                @test exact_track_element(z0, cav, ref, C) == zk
                rhs = (w, _) -> symplectic_gradient(u -> rf_potential(cav, ref, C; z = u), w; S = SYMPLECTIC_S)
                zg = integrate(rhs, z0, 0.0, 1.0; method = DormandPrince853(), tol = 1e-14).z
                @test maximum(abs, zk .- zg) < 1e-12
            end
        end
        sp = map_space(3)
        M = element_da_map(thin, ref, 3; path_length = C)
        @test _rf_coeff_gap(M, rf_kick(collect(da_vars(sp)), thin, ref, C)) == 0
    end

    @testset "thick flow approaches drift-kick-drift" begin
        thick = RFCavity("RFL", 1.5, 2, 0.0, 0.8)
        z0 = [1e-3, -2e-4, 1.5e-3, 1e-4, 0.05, 8e-3]
        missing = try
            element_flow(thick, ref, z0; method = RK4(), nsteps = 2)
            nothing
        catch e
            e
        end
        @test missing isa ArgumentError
        @test occursin("path_length", missing.msg)

        method = ImplicitYoshida(4)
        gaps = Float64[]
        for ds in (0.4, 0.1, 0.025)
            znew = exact_track_element(z0, thick, ref, C; ds = ds, method = method)
            n = max(1, ceil(Int, thick.L / ds))
            zflow = element_flow(thick, ref, z0; method = method, nsteps = n, path_length = C)
            @test znew == zflow
            zold = _strang_thick(z0, thick, ref, C; ds = ds, method = method)
            push!(gaps, maximum(abs, znew .- zold))
        end
        @test gaps[1] > gaps[2]
        @test gaps[2] > gaps[3]
        @test gaps[3] < 1e-9

        N = 3
        damethod = RK4()
        dagaps = Float64[]
        for n in (2, 4, 8)
            M = element_da_map(thick, ref, N; method = damethod, nsteps = n, path_length = C)
            if n == 2
                Mf = element_flow(thick, ref, collect(da_vars(map_space(N)));
                                  method = damethod, nsteps = n, path_length = C)
                @test _rf_coeff_gap(M, Mf) == 0
            end
            ds = thick.L / n
            Mold = _strang_thick(collect(da_vars(map_space(N))), thick, ref, C; ds = ds, method = damethod)
            push!(dagaps, _rf_coeff_gap(M, Mold))
        end
        @test dagaps[1] > dagaps[2]
        @test dagaps[2] > dagaps[3]
        @test dagaps[3] < 1e-5
    end
end
