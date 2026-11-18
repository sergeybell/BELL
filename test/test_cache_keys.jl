# Physics caches key the full RefParticle. Same p0 with another mass or charge must not collide.
using LinearAlgebra
using Bell.DiffAlgebra, Bell.Integrators

# reference_from_kinetic recomputes p0·c and does not land two species on one Float64.
# Build the second reference from the shared p0·c; β0 and γ0 still follow that momentum and the mass.
function ref_at_p0(species::Particle, p0::Real)
    m = species.mass_MeV
    p0 = Float64(p0)
    E = sqrt(p0 * p0 + m * m)
    gamma = E / m
    bg = p0 / m
    beta = bg / gamma
    return RefParticle(species, m, E - m, gamma, beta, bg, p0, p0 * 1e6 / c_light / species.charge)
end

poly_gap(A, B) = maximum(maxabs(A[i] - B[i]) for i in eachindex(A))
da_gap(A, B) = maximum(coeffnorm(A[i] - B[i]) for i in eachindex(A))

@testset "Cache keys use the full reference particle" begin
    clear_caches!()
    pref = reference_from_kinetic(proton, 200.0)
    dref = ref_at_p0(deuteron, pref.p0_MeV)
    @test pref.p0_MeV == dref.p0_MeV
    @test pref.mass_MeV != dref.mass_MeV
    @test pref.particle.charge == dref.particle.charge
    @test pref.brho_Tm == dref.brho_Tm
    @test pref.beta != dref.beta

    drift = Drift("D", 1.0)
    quad = Quadrupole("Q", 0.4, 1.5)

    @testset "mass at fixed p0" begin
        for e in (drift, quad)
            Rp = linear_matrix(e, pref)
            Rd = linear_matrix(e, dref)
            @test maximum(abs, Rp - Rd) > 1e-6
            Sp = Bell.Linear._SA(e, pref)
            @test Sp === Bell.Linear._SA(e, pref)
            @test Sp !== Bell.Linear._SA(e, dref)
            @test Rp ≈ exp(elem_length(e) * Sp)

            Lp = element_lie_map(e, pref)
            Ld = element_lie_map(e, dref)
            @test poly_gap(Lp, Ld) > 1e-8
            @test Lp === element_lie_map(e, pref)
            @test Lp !== Ld
        end

        Mp = element_da_map(drift, pref, 2; method = RK4(), nsteps = 4)
        Md = element_da_map(drift, dref, 2; method = RK4(), nsteps = 4)
        @test da_gap(Mp, Md) > 1e-8

        latp = Lattice("p", AbstractElement[drift], pref, 1)
        latd = Lattice("d", AbstractElement[drift], dref, 1)
        Cp = element_da_maps(latp, 2; cache = true, method = RK4(), nsteps = 4)
        Cd = element_da_maps(latd, 2; cache = true, method = RK4(), nsteps = 4)
        @test da_gap(Cp[1], Cd[1]) > 1e-8
        @test Cp[1] === element_da_maps(latp, 2; cache = true, method = RK4(), nsteps = 4)[1]
        @test Cp[1] !== Cd[1]

        Tp = taylor_map_element(drift, pref; ds = 0.5)
        Td = taylor_map_element(drift, dref; ds = 0.5)
        @test poly_gap(Tp, Td) > 1e-8
        Ep = element_maps(latp; ds = 0.5)
        Ed = element_maps(latd; ds = 0.5)
        @test poly_gap(Ep[1], Ed[1]) > 1e-8
        @test Ep[1] === element_maps(latp; ds = 0.5)[1]
        @test Ep[1] !== Ed[1]

        # This hard edge does not use the reference, so the polynomials agree.
        # The cache must still keep a separate entry for each reference particle.
        edge = EdgeMap("Q_in", quad, :in, HardEdge(order = 2))
        Fp = edge_map(edge, pref)
        Fd = edge_map(edge, dref)
        @test poly_gap(Fp, Fd) == 0
        @test Fp === edge_map(edge, pref)
        @test Fp !== Fd
    end

    @testset "charge at fixed mass" begin
        clear_caches!()
        m = proton.mass_MeV
        r1 = reference_from_kinetic(Particle("ion", m, 1.0, 0.0), 200.0)
        r2 = reference_from_kinetic(Particle("ion", m, 2.0, 0.0), 200.0)
        @test r1.p0_MeV == r2.p0_MeV
        @test r1.mass_MeV == r2.mass_MeV
        @test r1.beta == r2.beta
        @test r1.particle.charge != r2.particle.charge

        cav = RFCavity("RF", 1.0, 1, 0.0)
        path = 20.0
        @test maxabs(rf_potential(cav, r1, path) - rf_potential(cav, r2, path)) > 1e-8
        m1 = element_lie_map(cav, r1, path)
        m2 = element_lie_map(cav, r2, path)
        @test poly_gap(m1, m2) > 1e-8
        @test m1 === element_lie_map(cav, r1, path)
        @test m1 !== m2
    end
end
