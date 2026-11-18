# Fast end-to-end: a two-particle beam through a drift and a thin RF cavity.
# FieldTracker, LieTracker, and MapTracker each advance one turn.
# The spin slot is stored only and stays as it was.

@testset "trackers end to end" begin
    ref = reference_from_kinetic(proton, 100.0)
    lat = Lattice("tiny", AbstractElement[
        Drift("D", 0.4),
        RFCavity("RF", 0.05, 1, 0.0),
    ], ref, 1)
    z0 = [1e-4, -2e-5, 1.5e-4, 1e-5, 2e-4, 1e-4]
    spin = (0.0, 1.0, 0.0)

    function fresh_beam()
        return Beam([ParticleState(z0, spin), ParticleState(0.5 .* z0)])
    end

    function advanced(beam)
        z = beam.particles[1].z
        @test length(z) == 6
        @test all(isfinite, z)
        @test all(isfinite, beam.particles[2].z)
        @test z != z0
        @test beam.particles[1].spin == spin
        @test beam.particles[2].spin === nothing
        @test z0 == [1e-4, -2e-5, 1.5e-4, 1e-5, 2e-4, 1e-4]
        return nothing
    end

    @testset "FieldTracker" begin
        beam = fresh_beam()
        @test track!(beam, lat, FieldTracker(ds = 0.4); turns = 1) === beam
        advanced(beam)
    end

    @testset "LieTracker" begin
        beam = fresh_beam()
        @test track!(beam, lat, LieTracker(lat); turns = 1) === beam
        advanced(beam)
    end

    @testset "MapTracker" begin
        beam = fresh_beam()
        one_turn = LinearMap(ring_matrix(lat))
        @test track!(beam, lat, MapTracker(one_turn); turns = 1) === beam
        advanced(beam)
    end
end
