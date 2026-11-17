using LinearAlgebra
using Bell.DiffAlgebra

@testset "tracking interfaces" begin
    ref = reference_from_kinetic(deuteron, 270.0)
    lat = Lattice("fodo", AbstractElement[
        Quadrupole("QF", 0.3, 1.5), Drift("D1", 1.0), SBend("B", 1.2, 0.3),
        Drift("D2", 1.0), Quadrupole("QD", 0.3, -1.5), Sextupole("S", 0.2, 3.0),
    ], ref, 1)
    z0 = [1e-3, -2e-4, 1.5e-3, 1e-4, 2e-3, 1e-3]
    spin = (0.0, 0.0, 1.0)

    @testset "ParticleState and Beam" begin
        owned = copy(z0)
        p = ParticleState(owned, spin)
        owned[1] = 5.0
        @test p.z == z0
        @test p.spin === (0.0, 0.0, 1.0)
        @test ParticleState(z0).spin === nothing
        @test_throws ArgumentError ParticleState([1.0, 2.0])
        beam = Beam([z0, 2 .* z0])
        @test length(beam.particles) == 2
        @test beam.particles[1].z == z0
        @test beam.particles[2].z ≈ 2 .* z0
        @test beam.particles[1].z !== z0
    end

    @testset "map wrappers" begin
        R = Matrix{Float64}(I, 6, 6)
        R[1, 2] = 0.4
        S = Matrix{Float64}(I, 6, 6)
        S[3, 4] = -0.2
        A = LinearMap(R)
        B = LinearMap(S)
        @test linear_part(A) == R
        @test linear_part(A) !== R
        @test A(z0) ≈ R * z0
        @test compose(A, B)(z0) ≈ A(B(z0))
        @test_throws ArgumentError LinearMap(zeros(5, 5))
        @test LieFactor <: AbstractMap
        sp = map_space(2)
        X = da_vars(sp)
        comps = copy(X)
        comps[1] = X[1] + 0.2 * X[2]
        D = DAMap(comps)
        @test D isa AbstractMap
        @test linear_part(D)[1, 2] ≈ 0.2
        @test D(z0)[1] ≈ z0[1] + 0.2 * z0[2]
        @test D(z0)[2] ≈ z0[2]
        @test compose(D, DAMap(X))(z0)[1] ≈ D(z0)[1]
        @test_throws ArgumentError DAMap(X[1:5])
    end

    @testset "FieldTracker matches exact_track_turn" begin
        tracker = FieldTracker()
        @test tracker.ds == 0.1
        @test tracker.method isa Bell.Integrators.ImplicitYoshida
        @test tracker.method.order == 4
        beam = Beam([ParticleState(z0, spin), ParticleState(0.5 .* z0)])
        z1 = exact_track_turn(z0, lat; ds = tracker.ds, method = tracker.method)
        z2 = exact_track_turn(0.5 .* z0, lat; ds = tracker.ds, method = tracker.method)
        @test track!(beam, lat, tracker; turns = 1) === beam
        @test beam.particles[1].z == z1
        @test beam.particles[2].z == z2
        @test beam.particles[1].spin == spin
        @test z0[1] == 1e-3
        z = copy(z0)
        @test track!(z, lat, tracker; turns = 2) === z
        w = exact_track_turn(exact_track_turn(z0, lat; ds = tracker.ds, method = tracker.method),
                              lat; ds = tracker.ds, method = tracker.method)
        @test z == w
    end

    @testset "LieTracker matches track" begin
        factors = build_factors(lat)
        @test factors[1] isa AbstractMap
        @test factors[1](z0) == apply_factor(factors[1], z0)
        tracker = LieTracker(factors)
        z_lie = track(factors, z0; turns = 1)
        @test track!(copy(z0), lat, LieTracker(lat)) == z_lie
        beam = Beam(ParticleState(z0, spin))
        track!(beam, lat, tracker; turns = 1)
        @test beam.particles[1].z == z_lie
        @test beam.particles[1].spin == spin
        z = copy(z0)
        track!(z, lat, tracker)
        @test z == z_lie
    end

    @testset "MapTracker applies a one-turn map" begin
        R = Matrix{Float64}(I, 6, 6)
        R[1, 2] = 1.0
        tracker = MapTracker(LinearMap(R))
        beam = Beam(ParticleState(z0, spin))
        track!(beam, lat, tracker; turns = 2)
        @test beam.particles[1].z ≈ R * (R * z0)
        @test beam.particles[1].spin == spin
        fn = MapTracker(z -> z .+ 1)
        w = zeros(6)
        track!(w, lat, fn; turns = 3)
        @test w == fill(3.0, 6)
    end
end
