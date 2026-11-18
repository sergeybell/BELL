# Soft-edge gradient, pole-face placement, and model dispatch.
# Profiles are linear ramps so the suite stays fast.

_lin_gap(A, B) = maximum(abs, map_linear_part(A) - map_linear_part(B))

function _placed_edge(e, side, model, ref, place)
    F0, W, _ = profile_function(model.profile)
    h0 = curvature(e)
    F = side === :in ? F0 : (s -> F0(-s))
    E = Bell.Fringe._Edge(h0, F, W, side)
    k1 = e.k1
    method = model.method
    Sfull = Bell.FringeMaps._soft(E, method, ref, k1, -W, W)
    Sleft = Bell.FringeMaps._soft(E, method, ref, k1, -W, 0.0)
    Sright = Bell.FringeMaps._soft(E, method, ref, k1, 0.0, W)
    ang = side === :in ? e.e1 : e.e2
    Wg = lie_map(-h0 * tan(ang) * (Z[1]^2 - Z[3]^2) / 2, 1.0)
    S = if place === :start
        compose_maps(Sfull, Wg)
    elseif place === :end
        compose_maps(Wg, Sfull)
    else
        compose_maps(Sright, compose_maps(Wg, Sleft))
    end
    Bm = lie_map(hamiltonian(e, ref), -W)
    Dm = lie_map(hamiltonian(Drift("_", 1.0), ref), -W)
    return side === :in ? compose_maps(Bm, compose_maps(S, Dm)) : compose_maps(Dm, compose_maps(S, Bm))
end

@testset "fringe" begin
    ref = reference_from_kinetic(proton, 200.0)

    @testset "family and override dispatch" begin
        ramp = SoftEdge(profile = LinearRamp(0.2), method = Splitting(n = 2, order = 2))
        model = FringeModel(dipole = HardEdge(), quad = ramp, sext = NoEdge(),
                            override = Dict("Q2" => HardEdge(order = 1)))
        @test edge_model(model, SBend("B", 1.0, 0.1)) == HardEdge()
        @test edge_model(model, Quadrupole("Q", 0.2, 1.0)) == ramp
        @test edge_model(model, Sextupole("S", 0.2, 1.0)) == NoEdge()
        @test edge_model(model, Drift("D", 0.2)) == NoEdge()
        @test edge_model(model, Quadrupole("Q2", 0.2, 1.0)) == HardEdge(order = 1)
        seq = fringe_sequence(Lattice("L", AbstractElement[Drift("D", 0.2), SBend("B", 1.0, 0.1)], ref, 1),
                              FringeModel())
        @test length(seq) == 2
        edged = fringe_sequence(Lattice("L", AbstractElement[SBend("B", 1.0, 0.1)], ref, 1),
                                FringeModel(dipole = HardEdge(order = 1)))
        @test length(edged) == 3
        @test edged[1] isa EdgeMap && edged[3] isa EdgeMap
    end

    @testset "fringe Hamiltonian matches the body at F = 1" begin
        e = SBend("B", 1.2, 0.25, 0.4)
        z = [0.01, -0.002, 0.015, 0.001, -0.02, 0.003]
        h0 = curvature(e)
        Hd = hamiltonian(Drift("D", 1.0), ref, z)
        @test fringe_hamiltonian(h0, 1.0, 0.0, ref, z; k1 = e.k1) ≈ hamiltonian(e, ref, z)
        @test fringe_hamiltonian(h0, 0.0, 0.0, ref, z) ≈ Hd
        dH = fringe_hamiltonian(h0, 1.0, 0.0, ref, z; k1 = 0.4) -
             fringe_hamiltonian(h0, 1.0, 0.0, ref, z; k1 = 0.0)
        @test dH ≈ 0.4 * (z[1]^2 - z[3]^2) / 2
        # The edge field multiplies the body gradient by F, so F = 0 drops k1.
        E0 = Bell.Fringe._Edge(0.0, _ -> 0.0, 0.2, :in)
        @test fringe_hamiltonian(E0, ref, 0.0, z; k1 = e.k1) ≈ Hd
    end

    @testset "soft edge keeps the body gradient" begin
        # Body(-W) removes the full body Hamiltonian. The soft Hamiltonian has to
        # put k1*F back; otherwise the map is a quadrupole of length L-2W between drifts.
        L, k1, width = 1.0, 0.01, 0.2
        e = SBend("B", L, 0.0, k1)
        model = SoftEdge(profile = LinearRamp(width), method = Splitting(n = 8, order = 2))
        _, W, _ = profile_function(model.profile)
        Min = edge_map(EdgeMap("B_in", e, :in, model), ref)
        Mb = element_lie_map(e, ref)
        Mout = edge_map(EdgeMap("B_out", e, :out, model), ref)
        R = map_linear_part(compose_maps(Mout, compose_maps(Mb, Min)))
        Rb = map_linear_part(Mb)
        Qshort = lie_map(hamiltonian(e, ref), L - 2W)
        D = lie_map(hamiltonian(Drift("d", 1.0), ref), W)
        Rbug = map_linear_part(compose_maps(D, compose_maps(Qshort, D)))
        @test abs(R[2, 1] - Rb[2, 1]) < abs(R[2, 1] - Rbug[2, 1])
        @test abs(R[4, 3] - Rb[4, 3]) < abs(R[4, 3] - Rbug[4, 3])
        @test abs(R[2, 1] - Rb[2, 1]) < 5e-4
    end

    @testset "pole face sits on the effective field boundary" begin
        # The kick belongs at local s = 0. Before the fix the entrance wedge ran
        # at s = -W and the exit wedge at s = +W, both on the field-free end.
        model = SoftEdge(profile = LinearRamp(0.4), method = Splitting(n = 4, order = 2))
        cases = ((SBend("Bin", 2.0, 0.4, 0.0, 0.18, 0.0, 0.0, 0.0), :in),
                 (SBend("Bout", 2.0, 0.4, 0.0, 0.0, -0.12, 0.0, 0.0), :out))
        for (e, side) in cases
            M = edge_map(EdgeMap(e.name * "_edge", e, side, model), ref)
            @test _lin_gap(M, _placed_edge(e, side, model, ref, :efb)) < 1e-12
            @test _lin_gap(M, _placed_edge(e, side, model, ref, :start)) > 1e-4
            @test _lin_gap(M, _placed_edge(e, side, model, ref, :end)) > 1e-4
            @test symplectic_error(map_linear_part(M)) < 1e-8
        end
    end

    @testset "TaylorIntegrate follows the s-dependent flow" begin
        e = SBend("B", 1.5, 0.2, 0.0)
        profile = LinearRamp(0.3)
        ms = SoftEdge(profile = profile, method = Splitting(n = 12, order = 4))
        mt = SoftEdge(profile = profile, method = TaylorIntegrate(tol = 1e-4, nmax = 16))
        Rs = map_linear_part(edge_map(EdgeMap("B_in", e, :in, ms), ref))
        Rt = map_linear_part(edge_map(EdgeMap("B_in", e, :in, mt), ref))
        @test Rs ≈ Rt atol = 1e-3
    end

    @testset "tabulated profiles do not share a cache key" begin
        z = collect(range(-0.1, 0.1, length = 5))
        p1 = Tabulated(z, [0.0, 0.2, 0.5, 0.8, 1.0])
        p2 = Tabulated(z, [0.0, 0.0, 0.0, 0.2, 1.0])
        @test Bell.FringeMaps._model_key(SoftEdge(profile = p1)) !=
              Bell.FringeMaps._model_key(SoftEdge(profile = p2))
    end

    @testset "unsupported elements" begin
        @test_throws ErrorException edge_map(EdgeMap("S_in", Sextupole("S", 0.2, 1.0), :in, HardEdge()), ref)
        quad = Quadrupole("Q", 0.3, 1.0)
        @test_throws ErrorException edge_map(
            EdgeMap("Q_in", quad, :in, SoftEdge(profile = LinearRamp(0.1), method = Splitting(n = 2, order = 2))),
            ref)
        short = SBend("B", 0.05, 0.0, 0.0)
        @test_throws ErrorException edge_map(
            EdgeMap("B_in", short, :in, SoftEdge(profile = LinearRamp(0.2), method = Splitting(n = 2, order = 2))),
            ref)
    end
end
