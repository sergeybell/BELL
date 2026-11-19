# element_maps.jl
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# Maps of one element. Drift, quadrupole, and sector bend come from
# lattices/fodo.bell (deuteron, 270 MeV). The sextupole and the two RF
# cavities are built here. A thin cavity is rf_kick. A thick cavity uses
# path_length in its s-Hamiltonian.
#
#   julia --project=examples examples/element_maps/run.jl

using Bell, Bell.DiffAlgebra, Bell.Integrators, Printf, LinearAlgebra, TaylorSeries

include(joinpath(@__DIR__, "..", "common", "print_util.jl"))

const ORDER = 3
const NSTEPS = 8
const METHOD = RK4()
const Z0 = [1.0e-3, -2.0e-4, 1.2e-3, 1.0e-4, 5.0e-4, 1.0e-3]
const STARTED = time()

function poly_gap(A, B, ord = ORDER)
    return maximum(maxabs(truncate_to(a, ord) - truncate_to(b, ord)) for (a, b) in zip(A, B))
end

function da_gap(Md, Mt, ord = ORDER)
    gap = 0.0
    for (ad, at) in zip(Md, Mt)
        bt = truncate_to(at, ord)
        for k in 0:ord, e in exponents_of_order(k)
            gap = max(gap, abs(coefficient(ad, e) - coeff(bt, e)))
        end
    end
    return gap
end

at_taylor(M, z) = [evaluate(m, z) for m in M]
at_da(M, z) = [Float64(M[i](z)) for i in eachindex(M)]

function print_matrix(R)
    for i in 1:6
        println("    ", join((@sprintf("% .6e", R[i, j]) for j in 1:6), " "))
    end
end

function report_magnet(e, ref, C)
    section("$e  ($(kind_string(e)))",
            "linear_matrix, Lie series of H, Taylor and DA integrations of S grad H,",
            "element_flow and exact_track_element on one numeric point.")
    R = linear_matrix(e, ref)
    @printf("  linear_matrix symplectic error = %.3e\n", symplectic_error(R))
    print_matrix(R)

    Mseries = lie_map(hamiltonian(e, ref), elem_length(e); order = ORDER)
    Mlie = element_lie_map(e, ref, C; order = ORDER)
    F = dragt_finn_factor(Mlie; name = e.name)
    ds = elem_length(e) / NSTEPS
    Mt = taylor_map_element(e, ref; ds = ds, method = METHOD, order = ORDER, path_length = C)
    Md = element_da_map(e, ref, ORDER; method = METHOD, nsteps = NSTEPS)
    z_flow = element_flow(e, ref, Z0; method = METHOD, nsteps = NSTEPS, path_length = C)
    z_ex = exact_track_element(Z0, e, ref, C; ds = ds, method = METHOD)

    @printf("  |lie_map(H) - element_lie_map|     = %.3e\n", poly_gap(Mseries, Mlie))
    @printf("  |linear_matrix - Lie linear part| = %.3e\n", maximum(abs, R - map_linear_part(Mlie)))
    @printf("  |factor_to_map - Lie map|         = %.3e\n", poly_gap(factor_to_map(F), Mlie))
    @printf("  |f3| = %.3e    |f4| = %.3e    symplectic error of Lie map = %.3e\n",
            maxabs(F.f3), maxabs(F.f4), symplectic_error(Mlie))
    println("  f3 = ", poly_string(F.f3; digits = 5))
    println("  f4 = ", poly_string(F.f4; digits = 5))
    @printf("  |Taylor - Lie| = %.3e    symplectic error Taylor = %.3e\n",
            poly_gap(Mt, Mlie), symplectic_error(Mt))
    @printf("  |DA - Lie|     = %.3e    symplectic error DA     = %.3e\n",
            da_gap(Md, Mlie), symplectic_error(Md))
    @printf("  |DA - Taylor|  = %.3e\n", da_gap(Md, Mt))
    @printf("  |element_flow (%d steps) - exact_track| = %.3e\n",
            NSTEPS, maximum(abs, z_flow - z_ex))
    @printf("  |Lie(z) - exact_track|    = %.3e\n", maximum(abs, at_taylor(Mlie, Z0) - z_ex))
    @printf("  |Taylor(z) - exact_track| = %.3e\n", maximum(abs, at_taylor(Mt, Z0) - z_ex))
    @printf("  |DA(z) - exact_track|     = %.3e\n", maximum(abs, at_da(Md, Z0) - z_ex))
    if e isa Drift
        z1 = element_flow(e, ref, Z0; method = METHOD, nsteps = 1, path_length = C)
        @printf("  |element_flow (1 step) - exact_track| = %.3e\n", maximum(abs, z1 - z_ex))
        note("exact_track_element of a Drift is one step. The $(NSTEPS)-step flow is a finer step of the same H.")
    end
    return nothing
end

function report_thin(e, ref, C)
    section("$e  (thin RF, L = 0)",
            "No continuous Hamiltonian. The map is rf_kick, and the Lie generator is rf_potential.",
            "element_flow rejects a thin cavity.")
    R = linear_matrix(e, ref)
    @printf("  linear_matrix symplectic error = %.3e    |R - I| = %.3e\n",
            symplectic_error(R), maximum(abs, R - I))
    G = rf_potential(e, ref, C)
    Mseries = lie_map(G - constant_term(G), 1.0; order = ORDER)
    Mlie = element_lie_map(e, ref, C; order = ORDER)
    F = dragt_finn_factor(Mlie; name = e.name)
    Mt = taylor_map_element(e, ref; ds = 1 / NSTEPS, method = METHOD, order = ORDER, path_length = C)
    Md = element_da_map(e, ref, ORDER; path_length = C)
    z_kick = rf_kick(Z0, e, ref, C)
    z_ex = exact_track_element(Z0, e, ref, C)

    @printf("  |lie_map(rf_potential) - element_lie_map| = %.3e\n", poly_gap(Mseries, Mlie))
    @printf("  |factor_to_map - Lie map|                 = %.3e\n", poly_gap(factor_to_map(F), Mlie))
    @printf("  |f3| = %.3e    |f4| = %.3e\n", maxabs(F.f3), maxabs(F.f4))
    println("  f3 = ", poly_string(F.f3; digits = 5))
    println("  f4 = ", poly_string(F.f4; digits = 5))
    @printf("  |Taylor - Lie| = %.3e    |DA - Lie| = %.3e\n", poly_gap(Mt, Mlie), da_gap(Md, Mlie))
    @printf("  |rf_kick - exact_track| = %.3e\n", maximum(abs, z_kick - z_ex))
    @printf("  |Lie(z) - rf_kick|      = %.3e\n", maximum(abs, at_taylor(Mlie, Z0) - z_kick))
    @printf("  |Taylor(z) - rf_kick|   = %.3e\n", maximum(abs, at_taylor(Mt, Z0) - z_kick))
    @printf("  |DA(z) - rf_kick|       = %.3e\n", maximum(abs, at_da(Md, Z0) - z_kick))
    print("  hamiltonian(thin, ref, C, z) : ")
    try
        hamiltonian(e, ref, C, Z0)
        println("returned a value")
    catch err
        println(sprint(showerror, err))
    end
    print("  element_flow(thin) : ")
    try
        element_flow(e, ref, Z0; path_length = C)
        println("returned a value")
    catch err
        println(sprint(showerror, err))
    end
    return nothing
end

function report_thick(e, ref, C)
    section("$e  (thick RF, L = $(e.L) m)",
            "H = energy_term - ps + rf_potential/L. path_length is the FODO circumference.",
            "linear_matrix of an RFCavity is the drift matrix of the same length.")
    R = linear_matrix(e, ref)
    @printf("  linear_matrix symplectic error = %.3e\n", symplectic_error(R))
    print_matrix(R)
    H = hamiltonian(e, ref, C)
    Mseries = lie_map(H, elem_length(e); order = ORDER)
    Mlie = element_lie_map(e, ref, C; order = ORDER)
    F = dragt_finn_factor(Mlie; name = e.name)
    ds = elem_length(e) / NSTEPS
    Mt = taylor_map_element(e, ref; ds = ds, method = METHOD, order = ORDER, path_length = C)
    Md = element_da_map(e, ref, ORDER; method = METHOD, nsteps = NSTEPS, path_length = C)
    z_flow = element_flow(e, ref, Z0; method = METHOD, nsteps = NSTEPS, path_length = C)
    z_ex = exact_track_element(Z0, e, ref, C; ds = ds, method = METHOD)

    @printf("  |lie_map(H) - element_lie_map|     = %.3e\n", poly_gap(Mseries, Mlie))
    @printf("  |linear_matrix - Lie linear part| = %.3e\n", maximum(abs, R - map_linear_part(Mlie)))
    @printf("  |factor_to_map - Lie map|         = %.3e\n", poly_gap(factor_to_map(F), Mlie))
    @printf("  |f3| = %.3e    |f4| = %.3e\n", maxabs(F.f3), maxabs(F.f4))
    println("  f3 = ", poly_string(F.f3; digits = 5))
    println("  f4 = ", poly_string(F.f4; digits = 5))
    @printf("  |Taylor - Lie| = %.3e    symplectic error Taylor = %.3e\n",
            poly_gap(Mt, Mlie), symplectic_error(Mt))
    @printf("  |DA - Lie|     = %.3e    symplectic error DA     = %.3e\n",
            da_gap(Md, Mlie), symplectic_error(Md))
    @printf("  |element_flow - exact_track| = %.3e\n", maximum(abs, z_flow - z_ex))
    @printf("  |Lie(z) - exact_track|       = %.3e\n", maximum(abs, at_taylor(Mlie, Z0) - z_ex))
    @printf("  |Taylor(z) - exact_track|    = %.3e\n", maximum(abs, at_taylor(Mt, Z0) - z_ex))
    @printf("  |DA(z) - exact_track|        = %.3e\n", maximum(abs, at_da(Md, Z0) - z_ex))
    return nothing
end


banner("Element maps: Lie series, Taylor, DA, and a numeric step")

lat = load_structure(joinpath(@__DIR__, "..", "lattices", "fodo.bell"))
ref = lat.ref
C = path_length(lat)
drift = first(e for e in lat.elements if e isa Drift)
qf = first(e for e in lat.elements if e isa Quadrupole && e.k1 > 0)
bend = first(e for e in lat.elements if e isa SBend)
sx = Sextupole("SF", 0.20, 4.0)
thin = RFCavity("RFthin", 0.05, lat.harmonic, 0.0)
thick = RFCavity("RFthick", 0.05, lat.harmonic, 0.0, 0.40)

section("1. Reference particle",
        "The FODO file is a deuteron at 270 MeV. C is the path length used by both cavities.")
println(describe_reference(ref))
@printf("  path length C = %.7f m    harmonic = %d\n", C, lat.harmonic)
@printf("  probe z = (x, px, y, py, l, delta) = %s\n", Z0)
@printf("  Taylor / DA order = %d    RK4 steps = %d\n", ORDER, NSTEPS)

report_magnet(drift, ref, C)
report_magnet(qf, ref, C)
report_magnet(sx, ref, C)
report_magnet(bend, ref, C)
report_thin(thin, ref, C)
report_thick(thick, ref, C)

@printf("\n  wall time %.2f s\n", time() - STARTED)
println("Done.")
