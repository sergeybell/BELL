# Reference numbers copied from a real run recorded in
# benchmarks/twiss/bell/results.csv. This test does not read that file.
# Linear optics: load_structure → ring_matrix → periodic_optics → optics_along.
# Chromaticity: chromaticity(factor_to_map(ring_lie_factor(magnetic_lattice(lat)))).
#
# That run used convert_optim(magnetic_2.opt; close_ring=true): the BM angle
# there is π/8. examples/lattices/magnetic_2.bell keeps the unscaled angle, so
# the total bend is short of 2π by ~1e-7. Qx, Qy and the betas still match at
# rtol 1e-8; alpha_c and Dx need the slightly wider bands below.

function _bell_twiss_refs()
    return Dict{String,Float64}(
        "C_m" => 117.388296,
        "Qx" => 2.46172014987,
        "Qy" => 2.46133732305,
        "alpha_c" => 0.230470282157,
        "betx0_m" => 19.3675224161,
        "bety0_m" => 4.054470469,
        "Dx0_m" => 0.0110941104046,
        "xi_x" => -1.46153752013,
        "xi_y" => -2.22548794389,
    )
end

@testset "optics regression magnetic_2" begin
    refs = _bell_twiss_refs()
    lat = load_structure(joinpath(pkgdir(Bell), "examples", "lattices", "magnetic_2.bell"))
    C = path_length(lat)
    R = ring_matrix(lat)
    po = periodic_optics(R, lat.ref, C)
    opt = optics_along(lat, po)
    @test abs(C - refs["C_m"]) < 1e-6
    @test isapprox(opt.Qx, refs["Qx"]; rtol = 1e-8)
    @test isapprox(opt.Qy, refs["Qy"]; rtol = 1e-8)
    @test isapprox(po.alpha_c, refs["alpha_c"]; rtol = 1e-7)
    @test isapprox(po.x.beta, refs["betx0_m"]; rtol = 1e-8)
    @test isapprox(po.y.beta, refs["bety0_m"]; rtol = 1e-8)
    @test isapprox(po.D4[1], refs["Dx0_m"]; rtol = 2e-5)
    xi_x, xi_y = chromaticity(factor_to_map(ring_lie_factor(magnetic_lattice(lat))))
    @test isapprox(xi_x, refs["xi_x"]; rtol = 1e-5)
    @test isapprox(xi_y, refs["xi_y"]; rtol = 1e-5)
end
