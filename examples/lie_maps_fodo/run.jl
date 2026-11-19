# lie_maps_fodo.jl
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# Dragt-Finn factors of lattices/fodo.bell and lattices/magnetic_2.bell.
# build_factors, concatenate, product_factors, ring_lie_factor, factor_to_map,
# symplectic_error, track, and apply_factor. Example 17 remains the BCH check
# on its own pairs; one pair is repeated here on the FODO quadrupole.
#
#   julia --project=examples examples/lie_maps_fodo/run.jl

using Bell, Printf, LinearAlgebra, TaylorSeries

include(joinpath(@__DIR__, "..", "common", "print_util.jl"))

const Z0 = [1.0e-3, -2.0e-4, 1.2e-3, 1.0e-4, 5.0e-4, 1.0e-3]
const STARTED = time()

function poly_gap(A, B)
    return maximum(maxabs(truncate_to(a, MAP_ORDER) - truncate_to(b, MAP_ORDER)) for (a, b) in zip(A, B))
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


section("5. magnetic_2",
        "Same factor pipeline on the second lattice file. The one-turn Taylor map is left to example 11.")

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

@printf("\n  wall time %.2f s\n", time() - STARTED)
println("Done.")
