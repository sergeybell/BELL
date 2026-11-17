@testset "coordinate round-trip" begin
    lat = load_structure(joinpath(pkgdir(Bell), "examples", "lattices", "magnetic_2.bell"))
    ref = lat.ref
    C = path_length(lat)
    h = lat.harmonic
    # laboratory (x, x', y, y', φ, δ), every component nonzero
    u = [1.0e-3, -2.0e-3, 3.0e-3, 4.0e-3, 0.15, -1.0e-3]
    z = laboratory_to_canonical(u, ref, C, h)
    u_back = canonical_to_laboratory(z, ref, C, h)
    @test maximum(abs, u_back - u) < 1e-12
    # canonical (x, px, y, py, l, δ) → laboratory → canonical
    z0 = [1.0e-3, -2.0e-3, 3.0e-3, 4.0e-4, 0.02, 1.0e-3]
    z_back = laboratory_to_canonical(canonical_to_laboratory(z0, ref, C, h), ref, C, h)
    @test maximum(abs, z_back - z0) < 1e-12
end
