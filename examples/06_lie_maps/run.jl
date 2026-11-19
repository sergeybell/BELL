# 06_lie_maps/run.jl
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# Dragt-Finn factors of lattices/fodo.bell, then lattices/magnetic_2.bell.
# build_factors, concatenate, product_factors, ring_lie_factor, factor_to_map,
# symplectic_error, track, and apply_factor. concatenate is checked against
# the composed Lie series, including the Poisson-bracket sign at degree 4.
# magnetic_2 compares that product with the one-turn Taylor map and with
# a short symplectic track.
#
#   julia --project=examples examples/06_lie_maps/run.jl

using Bell, Printf, LinearAlgebra, TaylorSeries

include(joinpath(@__DIR__, "..", "common", "print_util.jl"))

const Z0 = [1.0e-3, -2.0e-4, 1.2e-3, 1.0e-4, 5.0e-4, 1.0e-3]
const STARTED = time()

function poly_gap(A, B)
    return maximum(maxabs(truncate_to(a, MAP_ORDER) - truncate_to(b, MAP_ORDER)) for (a, b) in zip(A, B))
end

function bch_pair(ea, eb, ref)
    Ma = lie_map(hamiltonian(ea, ref), elem_length(ea))
    Mb = lie_map(hamiltonian(eb, ref), elem_length(eb))
    Fa = dragt_finn_factor(Ma; name = ea.name)
    Fb = dragt_finn_factor(Mb; name = eb.name)
    Mc = compose_maps(Mb, Ma)
    Fc = dragt_finn_factor(Mc; name = "compose")
    Fp = concatenate(Fa, Fb; name = ea.name * "+" * eb.name)
    g3 = pull_generator(Fb.f3, Fa.R)
    g4 = pull_generator(Fb.f4, Fa.R)
    f4_opposite = hom(Fa.f4 + g4 + poisson(g3, Fa.f3) / 2, 4)
    return (
        label = ea.name * "+" * eb.name,
        f3 = maxabs(Fp.f3 - Fc.f3),
        f4 = maxabs(Fp.f4 - Fc.f4),
        f4opp = maxabs(f4_opposite - Fc.f4),
        R = maximum(abs, Fp.R - Fc.R),
        map = poly_gap(factor_to_map(Fp), Mc),
    )
end

function pair_gap(ea, eb, ref, C)
    Ma = element_lie_map(ea, ref, C)
    Mb = element_lie_map(eb, ref, C)
    Fa = dragt_finn_factor(Ma; name = ea.name)
    Fb = dragt_finn_factor(Mb; name = eb.name)
    Fp = concatenate(Fa, Fb; name = ea.name * "+" * eb.name)
    Mc = compose_maps(Mb, Ma)
    Fc = dragt_finn_factor(Mc; name = "compose")
    return (
        label = ea.name * "+" * eb.name,
        R = maximum(abs, Fp.R - Fc.R),
        f3 = maxabs(Fp.f3 - Fc.f3),
        f4 = maxabs(Fp.f4 - Fc.f4),
        map = poly_gap(factor_to_map(Fp), Mc),
    )
end

function print_unique(lattice, factors)
    @printf("  %-8s %-12s %12s %12s %12s\n", "name", "kind", "|f3|", "|f4|", "sym(R)")
    seen = Set{String}()
    for (e, F) in zip(lattice.elements, factors)
        e.name in seen && continue
        push!(seen, e.name)
        @printf("  %-8s %-12s %12.3e %12.3e %12.3e\n",
                e.name, kind_string(e), maxabs(F.f3), maxabs(F.f4), symplectic_error(F.R))
    end
    return nothing
end


banner("Lie factors of the FODO ring and of magnetic_2")

lat = load_structure(joinpath(@__DIR__, "..", "lattices", "fodo.bell"))
ref = lat.ref
C = path_length(lat)
cell_elems = lat.elements[1:8]
cell = Lattice("fodo_cell", cell_elems, ref, lat.harmonic)
cell_names = [e.name for e in cell_elems]
all(i -> lat.elements[i].name == cell_names[mod1(i, 8)], eachindex(lat.elements)) ||
    error("fodo.bell sequence is not 16 copies of the 8-element cell")

section("1. FODO line",
        "Sixteen copies of one cell. Factors are one Dragt-Finn factor per element.")
println(describe_reference(ref))
@printf("  elements %d    cell %d    C = %.7f m\n", length(lat.elements), length(cell_elems), C)
@printf("  cell: %s\n", join(cell_names, " "))

t_factors = @elapsed factors = build_factors(lat)
@printf("  build_factors: %d factors in %.3f s\n", length(factors), t_factors)
print_unique(lat, factors)


section("2. concatenate against the composed Lie series",
        "concatenate(Fa, Fb) is Fa then Fb. compose_maps(Mb, Ma) applies Ma first.",
        "The sextupole is not in the file; it supplies a cubic generator.")

sx = Sextupole("SF", 0.20, 4.0)
qf = first(e for e in cell_elems if e isa Quadrupole && e.k1 > 0)
qd = first(e for e in cell_elems if e isa Quadrupole && e.k1 < 0)
bend = first(e for e in cell_elems if e isa SBend)
@printf("  %-12s %12s %12s %12s %12s\n", "pair", "|dR|", "|df3|", "|df4|", "|d map|")
for (ea, eb) in ((qf, bend), (qf, qd), (qf, sx), (sx, bend))
    r = pair_gap(ea, eb, ref, C)
    @printf("  %-12s %12.3e %12.3e %12.3e %12.3e\n", r.label, r.R, r.f3, r.f4, r.map)
end


section("3. product_factors and ring_lie_factor",
        "product_factors folds concatenate. ring_lie_factor is that product for the line.",
        "The cell product is the first eight factors of the ring.")

t_cell = @elapsed Fcell = ring_lie_factor(cell)
Feight = product_factors(factors[1:8]; name = "eight")
Fring = ring_lie_factor(lat)
Mcell = factor_to_map(Fcell)
Mring = factor_to_map(Fring)
Rring = ring_matrix(lat)

@printf("  ring_lie_factor(cell) in %.3f s\n", t_cell)
@printf("  |cell product - first 8 factors|   |dR| = %.3e  |df3| = %.3e  |df4| = %.3e\n",
        maximum(abs, Fcell.R - Feight.R), maxabs(Fcell.f3 - Feight.f3), maxabs(Fcell.f4 - Feight.f4))
@printf("  |cell R - ring_matrix(cell)|       = %.3e\n", maximum(abs, Fcell.R - ring_matrix(cell)))
@printf("  |ring R - ring_matrix|             = %.3e\n", maximum(abs, Fring.R - Rring))
@printf("  symplectic error of cell factor map = %.3e\n", symplectic_error(Mcell))
@printf("  symplectic error of ring factor map = %.3e\n", symplectic_error(Mring))
@printf("  symplectic error of ring matrix      = %.3e\n", symplectic_error(Rring))
@printf("  |f3| cell = %.3e    |f4| cell = %.3e\n", maxabs(Fcell.f3), maxabs(Fcell.f4))
@printf("  |f3| ring = %.3e    |f4| ring = %.3e\n", maxabs(Fring.f3), maxabs(Fring.f4))


section("4. apply_factor and track",
        "apply_factor is the exact flow of f4, then f3, then R.",
        "track is that step once per element. The numeric turn is exact_track_turn.")

Fq = dragt_finn_factor(element_lie_map(qf, ref, C); name = qf.name)
z_fac = apply_factor(Fq, Z0)
z_trk = track([Fq], Z0)
z_ex = exact_track_element(Z0, qf, ref, C)
@printf("  |apply_factor - track| on %s = %.3e\n", qf.name, maximum(abs, z_fac - z_trk))
@printf("  |apply_factor - exact_track_element| = %.3e\n", maximum(abs, z_fac - z_ex))
@printf("  symplectic error of apply_factor at z = %.3e\n",
        symplectic_error(z -> apply_factor(Fq, z), Z0))
@printf("  symplectic error of the polynomial factor map = %.3e\n", symplectic_error(factor_to_map(Fq)))

w = foldl((z, F) -> apply_factor(F, z), factors[1:8]; init = copy(Z0))
@printf("  |eight apply_factor - track of those factors| = %.3e\n",
        maximum(abs, w - track(factors[1:8], Z0)))

t_turn = @elapsed z_turn = track(factors, Z0)
z_num = exact_track_turn(Z0, lat)
@printf("  one FODO turn of track in %.3f s\n", t_turn)
@printf("  |track - exact_track_turn| = %.3e\n", maximum(abs, z_turn - z_num))
@printf("  delta in  = %.6e    delta out (Lie) = %.6e    delta out (numeric) = %.6e\n",
        Z0[6], z_turn[6], z_num[6])


section("5. magnetic_2 factors",
        "Same factor pipeline on the second lattice file.")

mag = load_structure(joinpath(@__DIR__, "..", "lattices", "magnetic_2.bell"))
@printf("  elements %d    unique %d    C = %.7f m\n",
        length(mag.elements), length(unique_elements(mag)), path_length(mag))
t_mag = @elapsed factors_m = build_factors(mag)
t_ring = @elapsed Fring_m = ring_lie_factor(mag)
@printf("  build_factors in %.3f s    ring_lie_factor in %.3f s\n", t_mag, t_ring)
print_unique(mag, factors_m)
Rm = ring_matrix(mag)
Mm = factor_to_map(Fring_m)
@printf("  |ring R - ring_matrix| = %.3e\n", maximum(abs, Fring_m.R - Rm))
@printf("  symplectic error of factor map = %.3e    of ring matrix = %.3e\n",
        symplectic_error(Mm), symplectic_error(Rm))
t_mturn = @elapsed z_m = track(factors_m, Z0)
z_mn = exact_track_turn(Z0, mag; ds = 0.5)
@printf("  one turn of track in %.3f s    |track - exact_track_turn| = %.3e\n",
        t_mturn, maximum(abs, z_m - z_mn))


section("6. Poisson bracket in concatenate",
        "lie_map(H, L) is the series of one element. concatenate(Fa, Fb) is Fa then Fb.",
        "The degree-4 generator uses (1/2) of the Poisson bracket of the pulled f3.",
        "The column |df4 opp| flips that sign and does not match the composed series.")

drift = first(e for e in cell_elems if e isa Drift)
@printf("  %-12s %12s %12s %12s %12s %12s\n",
        "pair", "|df3|", "|df4|", "|df4 opp|", "|dR|", "|d map|")
for (ea, eb) in ((drift, qf), (qf, qd), (qf, sx), (sx, drift), (bend, qf), (sx, sx))
    r = bch_pair(ea, eb, ref)
    @printf("  %-12s %12.3e %12.3e %12.3e %12.3e %12.3e\n",
            r.label, r.f3, r.f4, r.f4opp, r.R, r.map)
end

println()
Ma, Mb, Mc = (lie_map(hamiltonian(e, ref), elem_length(e)) for e in (qf, sx, qd))
Fa, Fb, Fc = (dragt_finn_factor(M; name = e.name) for (M, e) in zip((Ma, Mb, Mc), (qf, sx, qd)))
Mabc = compose_maps(Mc, compose_maps(Mb, Ma))
Fabc = product_factors((Fa, Fb, Fc); name = "QF+SF+QD")
Ffrom = dragt_finn_factor(Mabc; name = "compose")
@printf("  three-element fold QF+SF+QD: |df3| = %.3e  |df4| = %.3e  |d map| = %.3e\n",
        maxabs(Fabc.f3 - Ffrom.f3), maxabs(Fabc.f4 - Ffrom.f4), poly_gap(factor_to_map(Fabc), Mabc))


section("7. magnetic_2 elements, Taylor map, and a tracked action",
        "Each distinct element: linear part versus linear_matrix, and the cubic generator.",
        "The ring Lie product is compared with ring_taylor_map.",
        "Two hundred turns of track keep the Courant-Snyder action of a small amplitude.")

z_amp = [3e-3, 2e-4, 2e-3, 2e-4, 1e-3, 2e-3]
@printf("  %-6s %12s %12s %12s %12s %10s\n",
        "elem", "|lin-R|", "recon", "TM sym@0", "LF sym@amp", "nterms f3")
for e in unique_elements(mag)
    M = element_lie_map(e, mag.ref)
    F = dragt_finn_factor(M; name = e.name)
    lin = maximum(abs, map_linear_part(M) - linear_matrix(e, mag.ref))
    rec = poly_gap(factor_to_map(F), M)
    s0 = symplectic_error(jacobian_of_map(M, zeros(6)))
    sf = symplectic_error(jacobian_of_function(z -> apply_factor(F, z), z_amp))
    nt = count(!iszero, vcat([F.f3[3].coeffs]...))
    @printf("  %-6s %12.3e %12.3e %12.3e %12.3e %10d\n",
            e.name, lin, rec, s0, sf, nt)
end

for e in unique_elements(mag)
    e isa SBend || e isa Quadrupole || e isa Sextupole || continue
    e.name in ("BM", "QF", "SD") || continue
    println()
    println("  ", e.name, ": ", hamiltonian_formula(e))
    println("    H = ", poly_string(hamiltonian(e, mag.ref); digits = 5))
    F = dragt_finn_factor(element_lie_map(e, mag.ref))
    println("    f3 = ", poly_string(F.f3; digits = 5))
end

t_tm = @elapsed Mtm = ring_taylor_map(mag)
Mlie = factor_to_map(Fring_m)
d2 = maximum(maxabs(truncate_to(a - b, 2)) for (a, b) in zip(Mlie, Mtm))
println()
@printf("  Taylor map in %.2f s\n", t_tm)
@printf("  |R(Lie) - ring_matrix|           = %.3e\n", maximum(abs, Fring_m.R - Rm))
@printf("  |R(Taylor) - ring_matrix|        = %.3e\n", maximum(abs, map_linear_part(Mtm) - Rm))
@printf("  max |2nd-order Lie - Taylor|     = %.3e\n", d2)
@printf("  symplectic error of Taylor map at 0 = %.3e    at z_amp = %.3e\n",
        symplectic_error(jacobian_of_map(Mtm, zeros(6))),
        symplectic_error(jacobian_of_map(Mtm, z_amp)))

Cm = path_length(mag)
h_m = mag.harmonic
u0 = [2e-3, 0.0, 1e-3, 0.0, 0.0, 1e-3]
z0 = laboratory_to_canonical(u0, mag.ref, Cm, h_m)
po_m = periodic_optics(Rm, mag.ref, Cm)
co = closed_orbit(Mtm, z0[6])
Jx(w) = (x = w[1] - co[1]; p = w[2] - co[2];
         (po_m.x.beta * p^2 + 2po_m.x.alpha * x * p + (1 + po_m.x.alpha^2) / po_m.x.beta * x^2) / 2)
Jy(w) = (x = w[3] - co[3]; p = w[4] - co[4];
         (po_m.y.beta * p^2 + 2po_m.y.alpha * x * p + (1 + po_m.y.alpha^2) / po_m.y.beta * x^2) / 2)

n_action = 200
wact = copy(z0)
jx = Float64[]
jy = Float64[]
t_act = @elapsed for _ in 1:n_action
    global wact = track(factors_m, wact)
    push!(jx, Jx(wact))
    push!(jy, Jy(wact))
end
println()
@printf("  start laboratory (x, x', y, y', phi, delta) = %s\n",
        join([@sprintf("% .3e", v) for v in u0], " "))
@printf("  %d Lie turns in %.2f s\n", n_action, t_act)
@printf("  Jx mean %.4e    relative spread %.3e\n",
        sum(jx) / n_action, (maximum(jx) - minimum(jx)) / (sum(jx) / n_action))
@printf("  Jy mean %.4e    relative spread %.3e\n",
        sum(jy) / n_action, (maximum(jy) - minimum(jy)) / (sum(jy) / n_action))
@printf("  final delta = %.15e    start delta = %.15e\n", wact[6], z0[6])
note("The spread is the nonlinear map and the delta dependence of the Twiss, on a symplectic step.")

@printf("\n  wall time %.2f s\n", time() - STARTED)
println("Done.")
