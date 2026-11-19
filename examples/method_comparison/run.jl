# method_comparison.jl
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# Same elements and the FODO cell and ring, orders 1-3.
#   Lie path:     element_lie_map / ring_lie_factor, then factor_to_map
#   DA from H:    element_da_map / one_turn_map_da
#   Taylor check: taylor_map_element / ring_taylor_map
#   numeric:      exact_track_element / exact_track_turn
# The generator table extracts R, f3, f4 from three order-3 maps.
# src/ has no DA-TaylorN converter. da_to_taylorn rebuilds the polynomial
# from coefficients and is local to this example.
#
#   julia --project=examples examples/method_comparison/run.jl

using Bell, Bell.DiffAlgebra, Bell.Integrators, Printf, TaylorSeries

include(joinpath(@__DIR__, "..", "common", "print_util.jl"))

const ORDERS = (1, 2, 3)
const NSTEPS = 8
const METHOD = RK4()
const STEP_DS = 0.25
const Z0 = [1.0e-3, -2.0e-4, 1.2e-3, 1.0e-4, 5.0e-4, 1.0e-3]
const STARTED = time()

# src/ has no DA-TaylorN converter. Rebuild the polynomial in the TaylorSeries
# space Bell opened for Z (degree JET_ORDER) by walking DA coefficients.
function da_to_taylorn(a)
    a.sp.nv == 6 || error("da_to_taylorn expects 6 variables, got $(a.sp.nv)")
    p = zero(Z[1])
    one_p = one(p)
    for m in eachindex(a.c)
        c = Float64(a.c[m])
        iszero(c) && continue
        ex = monomial_exponents(a.sp, m)
        term = one_p
        for i in 1:6
            e = Int(ex[i])
            e == 0 && continue
            term *= Z[i]^e
        end
        p += c * term
    end
    return p
end

to_taylor(M) = [truncate_to(da_to_taylorn(m), MAP_ORDER) for m in M]

function poly_gap(A, B, ord)
    return maximum(maxabs(truncate_to(a, ord) - truncate_to(b, ord)) for (a, b) in zip(A, B))
end

function check_converter()
    sp = map_space(3)
    vars = da_vars(sp)
    for i in 1:6
        maxabs(da_to_taylorn(vars[i]) - Z[i]) < 1e-15 ||
            error("da_to_taylorn failed on coordinate $i")
    end
    p = vars[1]^2 * vars[6] + 0.5 * vars[2]^3
    t = da_to_taylorn(p)
    abs(coeff(t, [2, 0, 0, 0, 0, 1]) - 1.0) < 1e-14 || error("da_to_taylorn failed on x^2*delta")
    abs(coeff(t, [0, 3, 0, 0, 0, 0]) - 0.5) < 1e-14 || error("da_to_taylorn failed on px^3")
    return nothing
end

function step_ds(e)
    L = elem_length(e)
    return L == 0 ? 1.0 : L / NSTEPS
end

function lie_of(e, ref, C, ord)
    return element_lie_map(e, ref, C; order = ord)
end

function taylor_of(e, ref, C, ord)
    ds = e isa RFCavity && e.L == 0 ? 1 / NSTEPS : step_ds(e)
    return taylor_map_element(e, ref; ds = ds, method = METHOD, order = ord, path_length = C)
end

function da_of(e, ref, C, ord)
    if e isa RFCavity && e.L == 0
        return element_da_map(e, ref, ord; path_length = C)
    end
    return element_da_map(e, ref, ord; method = METHOD, nsteps = NSTEPS, path_length = C)
end

at_taylor(M, z) = [evaluate(m, z) for m in M]
at_da(M, z) = [Float64(M[i](z)) for i in eachindex(M)]

function drop_constant(M)
    cmax = maximum(m -> abs(constant_term(m)), M)
    cmax <= 1e-14 && return M, cmax
    cmax < 1e-8 || error("constant part $cmax is too large to drop before dragt_finn_factor")
    return [m - constant_term(m) for m in M], cmax
end

function as_factor(M, name)
    prepared, cmax = drop_constant(M)
    return dragt_finn_factor(prepared; name = name), cmax
end

function print_factor_row(tag, pair, A, B)
    @printf("  %-12s %-10s %12.3e %12.3e %12.3e\n",
            tag, pair, maximum(abs, A.R - B.R), maxabs(A.f3 - B.f3), maxabs(A.f4 - B.f4))
    return nothing
end

function field_element(e, ref, N)
    sp = map_space(N)
    L = elem_length(e)
    z0 = collect(da_vars(sp))
    return integrate_fixed(element_rhs(e, ref), z0, 0.0, L, NSTEPS; method = METHOD).z
end

function field_one_turn(lattice, N)
    sp = map_space(N)
    ref = lattice.ref
    built = Dict()
    M = identity_map(sp)
    for e in lattice.elements
        Me = get!(built, e) do
            field_element(e, ref, N)
        end
        M = compose(Me, M)
    end
    return M
end

function composed_taylor(lattice, ord)
    M = collect(Z)
    C = path_length(lattice)
    built = Dict()
    for e in lattice.elements
        # One ds for the whole line, the same value ring_taylor_map receives.
        # A per-element L/NSTEPS cannot be passed through ring_taylor_map.
        Me = get!(built, e) do
            taylor_map_element(e, lattice.ref; ds = STEP_DS, method = METHOD,
                               order = ord, path_length = C)
        end
        M = compose_maps(Me, M; order = ord)
    end
    return M
end

function composed_lie(lattice, ord)
    M = collect(Z)
    C = path_length(lattice)
    for e in lattice.elements
        Me = lie_of(e, lattice.ref, C, ord)
        M = compose_maps(Me, M; order = ord)
    end
    return M
end

function report_element(e, ref, C)
    section("$e  ($(kind_string(e)))")
    t_num = @elapsed z_num = exact_track_element(Z0, e, ref, C; ds = step_ds(e), method = METHOD)
    @printf("  exact_track in %.4f s\n", t_num)
    if !(e isa RFCavity && e.L == 0)
        z_flow = element_flow(e, ref, Z0; method = METHOD, nsteps = NSTEPS, path_length = C)
        @printf("  |element_flow (%d steps) - exact_track| = %.3e\n", NSTEPS, maximum(abs, z_flow - z_num))
    else
        @printf("  element_flow is not used: a thin cavity is rf_kick, and exact_track calls it\n")
    end
    if e isa Drift
        z1 = element_flow(e, ref, Z0; method = METHOD, nsteps = 1, path_length = C)
        @printf("  |element_flow (1 step) - exact_track| = %.3e    (a drift track is one step)\n",
                maximum(abs, z1 - z_num))
    end
    @printf("  %3s %10s %10s %10s %10s %10s %10s %8s %8s %8s %10s %10s %10s\n",
            "ord", "Lie-Ty", "Lie-DA", "Ty-DA", "symL", "symD", "symT",
            "tL", "tD", "tT", "zL-zN", "zT-zN", "zD-zN")
    saved = nothing
    for ord in ORDERS
        tL = @elapsed Mlie = lie_of(e, ref, C, ord)
        tT = @elapsed Mt = taylor_of(e, ref, C, ord)
        tD = @elapsed Md = da_of(e, ref, C, ord)
        zL = at_taylor(Mlie, Z0)
        zT = at_taylor(Mt, Z0)
        zD = at_da(Md, Z0)
        @printf("  %3d %10.1e %10.1e %10.1e %10.1e %10.1e %10.1e %8.4f %8.4f %8.4f %10.1e %10.1e %10.1e\n",
                ord, poly_gap(Mlie, Mt, ord), poly_gap(Mlie, to_taylor(Md), ord), poly_gap(Mt, to_taylor(Md), ord),
                symplectic_error(Mlie), symplectic_error(Md), symplectic_error(Mt),
                tL, tD, tT, maximum(abs, zL - z_num), maximum(abs, zT - z_num), maximum(abs, zD - z_num))
        ord == 3 && (saved = (lie = Mlie, da = Md))
    end
    return saved
end

function report_line(label, lattice)
    section(label,
            "Lie: ring_lie_factor then factor_to_map at the printed order.",
            "DA: one_turn_map_da, RK4, $(NSTEPS) steps on every element.",
            "Taylor: taylor_map_element composed with one ds = $(STEP_DS) m (step count max(4, ceil(L/ds))).",
            "ring_taylor_map uses that same ds and always truncates at MAP_ORDER.",
            "Numeric: exact_track_turn with ds = $(STEP_DS) m. A drift is still one step.")
    tF = @elapsed F = ring_lie_factor(lattice)
    @printf("  ring_lie_factor in %.3f s    |R - ring_matrix| = %.3e    sym(R) = %.3e\n",
            tF, maximum(abs, F.R - ring_matrix(lattice)), symplectic_error(F.R))
    tN = @elapsed z_num = exact_track_turn(Z0, lattice; ds = STEP_DS, method = METHOD)
    @printf("  exact_track_turn in %.3f s\n", tN)
    @printf("  %3s %10s %10s %10s %10s %10s %10s %8s %8s %8s %10s %10s %10s\n",
            "ord", "Lie-Ty", "Lie-DA", "Ty-DA", "symL", "symD", "symT",
            "tL", "tD", "tT", "zL-zN", "zT-zN", "zD-zN")
    bundle = nothing
    for ord in ORDERS
        tL = @elapsed Mlie = factor_to_map(F; order = ord)
        tD = @elapsed Md = one_turn_map_da(lattice, ord; method = METHOD, nsteps = NSTEPS)
        tT = @elapsed Mt = composed_taylor(lattice, ord)
        @printf("  %3d %10.1e %10.1e %10.1e %10.1e %10.1e %10.1e %8.4f %8.4f %8.4f %10.1e %10.1e %10.1e\n",
                ord, poly_gap(Mlie, Mt, ord), poly_gap(Mlie, to_taylor(Md), ord), poly_gap(Mt, to_taylor(Md), ord),
                symplectic_error(Mlie), symplectic_error(Md), symplectic_error(Mt),
                tL, tD, tT,
                maximum(abs, at_taylor(Mlie, Z0) - z_num),
                maximum(abs, at_taylor(Mt, Z0) - z_num),
                maximum(abs, at_da(Md, Z0) - z_num))
        if ord == MAP_ORDER
            t_pub = @elapsed Mpub = ring_taylor_map(lattice; ds = STEP_DS)
            t_comp = @elapsed Mcomp = composed_lie(lattice, ord)
            @printf("  ring_taylor_map in %.3f s    |ring_taylor_map - composed taylor_map_element| = %.3e\n",
                    t_pub, poly_gap(Mpub, Mt, ord))
            @printf("  composed element_lie_map in %.3f s    |composed Lie - factor_to_map| = %.3e\n",
                    t_comp, poly_gap(Mcomp, Mlie, ord))
            bundle = (factor = F, lie = Mlie, composed = Mcomp, da = Md)
        end
    end
    return bundle
end

function generator_block(tag, lie_map, da_h, da_field)
    A, cA = as_factor(lie_map, tag * "-lie")
    B, cB = as_factor(to_taylor(da_h), tag * "-daH")
    C, cC = as_factor(to_taylor(da_field), tag * "-field")
    @printf("  %s    |f3|_A = %.3e    |f4|_A = %.3e    constants A/B/C = %.1e / %.1e / %.1e\n",
            tag, maxabs(A.f3), maxabs(A.f4), cA, cB, cC)
    @printf("  %-12s %-10s %12s %12s %12s\n", "object", "pair", "|dR|", "|df3|", "|df4|")
    print_factor_row(tag, "A-B", A, B)
    print_factor_row(tag, "A-C", A, C)
    print_factor_row(tag, "B-C", B, C)
    @printf("  symplectic error  A (Lie Taylor) = %.3e    B (DA of H) = %.3e    C (DA of fields) = %.3e\n",
            symplectic_error(lie_map), symplectic_error(da_h), symplectic_error(da_field))
    return (A = A, B = B, C = C)
end

function time_dual(q, ref)
    H(z) = hamiltonian(q, ref, z)
    rhs = element_rhs(q, ref)
    symplectic_gradient(H, Z0)
    rhs(Z0, 0.0)
    n = 4000
    t_dual = @elapsed for _ in 1:n
        symplectic_gradient(H, Z0)
    end
    t_field = @elapsed for _ in 1:n
        rhs(Z0, 0.0)
    end
    element_flow(q, ref, Z0; method = METHOD, nsteps = 40)
    integrate_fixed(rhs, collect(Z0), 0.0, elem_length(q), 40; method = METHOD)
    n2 = 40
    t_flow = @elapsed for _ in 1:n2
        element_flow(q, ref, Z0; method = METHOD, nsteps = 40)
    end
    t_rhs = @elapsed for _ in 1:n2
        integrate_fixed(rhs, collect(Z0), 0.0, elem_length(q), 40; method = METHOD)
    end
    @printf("  symplectic_gradient (Dual of H), %d calls: %.4f s\n", n, t_dual)
    @printf("  element_rhs (explicit field),    %d calls: %.4f s\n", n, t_field)
    @printf("  ratio Dual/field = %.3f\n", t_dual / t_field)
    @printf("  element_flow, %d x 40 RK4 steps: %.4f s\n", n2, t_flow)
    @printf("  integrate_fixed(element_rhs), same count: %.4f s\n", t_rhs)
    @printf("  ratio flow/field-stepper = %.3f\n", t_flow / t_rhs)
    return nothing
end


banner("Lie, DA, Taylor, and numeric maps at orders 1-3")

check_converter()
note("da_to_taylorn reproduced the coordinate functions and x^2*delta, px^3.")

lat = load_structure(joinpath(@__DIR__, "..", "lattices", "fodo.bell"))
ref = lat.ref
C = path_length(lat)
cell_elems = lat.elements[1:8]
cell = Lattice("fodo_cell", cell_elems, ref, lat.harmonic)
drift = first(e for e in cell_elems if e isa Drift)
qf = first(e for e in cell_elems if e isa Quadrupole && e.k1 > 0)
bend = first(e for e in cell_elems if e isa SBend)
sx = Sextupole("SF", 0.20, 4.0)
thin = RFCavity("RFthin", 0.05, lat.harmonic, 0.0)
thick = RFCavity("RFthick", 0.05, lat.harmonic, 0.0, 0.40)

section("1. Setup",
        "Deuteron at 270 MeV from fodo.bell. SF and the cavities are built here.",
        "RK4, $(NSTEPS) steps on each element map. Probe |z| is a few millimetres or 1e-3.",
        "symplectic_error of a TaylorN map always runs through degree MAP_ORDER - 1.",
        "An order-2 Taylor column is the cubic defect of that truncation, not a broken Lie series.",
        "A DA map is checked through its own order, so its order-2 column is a different cut.")
println(describe_reference(ref))
@printf("  C = %.7f m    cell elements = %d    ring elements = %d\n",
        C, length(cell_elems), length(lat.elements))
@printf("  MAP_ORDER = %d    JET_ORDER = %d    factors carry f3 and f4\n", MAP_ORDER, JET_ORDER)

function collect_elements(elements, ref, C, keep_name)
    kept = nothing
    for e in elements
        saved = report_element(e, ref, C)
        e.name == keep_name && (kept = saved)
    end
    return kept
end

sx_saved = collect_elements((drift, qf, sx, bend, thin, thick), ref, C, "SF")

cell_bundle = report_line("2. FODO cell (8 elements)", cell)
ring_bundle = report_line("3. FODO ring (16 cells)", lat)

section("4. Dual gradient versus the explicit field",
        "element_flow evaluates S grad H with Dual numbers at every stage.",
        "element_rhs is the explicit field kept for comparison. Times are after one warmup call.",
        "The quadrupole is $(qf.name).")
time_dual(qf, ref)

section("5. Dragt-Finn generators",
        "(A) Lie path map, then dragt_finn_factor.",
        "(B) DA map of H, da_to_taylorn, then dragt_finn_factor.",
        "(C) DA map of the field equations: integrate_fixed(element_rhs, da_vars, 0, L, $(NSTEPS); method = RK4).",
        "Same RK4 step count for (B) and (C). (A) is the Lie series of H, not that stepper.",
        "element_rhs has no RFCavity method, so the table is the sextupole, the cell, and the ring.")

t_field = @elapsed sx_field = field_element(sx, ref, MAP_ORDER)
@printf("  sextupole field map in %.3f s\n", t_field)
generator_block("SF", sx_saved.lie, sx_saved.da, sx_field)

t_cell_f = @elapsed cell_field = field_one_turn(cell, MAP_ORDER)
@printf("\n  cell field one-turn in %.3f s\n", t_cell_f)
generator_block("FODO cell", cell_bundle.lie, cell_bundle.da, cell_field)
@printf("  |composed Lie series - factor_to_map| before the factor split = %.3e\n",
        poly_gap(cell_bundle.composed, cell_bundle.lie, MAP_ORDER))

t_ring_f = @elapsed ring_field = field_one_turn(lat, MAP_ORDER)
@printf("\n  ring field one-turn in %.3f s\n", t_ring_f)
generator_block("FODO ring", ring_bundle.lie, ring_bundle.da, ring_field)

@printf("\n  wall time %.2f s\n", time() - STARTED)
println("Done.")
