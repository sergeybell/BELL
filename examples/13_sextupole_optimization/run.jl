# 13_optional_sextupole_optimization.jl
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# Optional: place three sextupole families and solve two example targets.
# Other examples leave the lattice as given (sextupoles at zero).
#
#   julia --project=examples examples/13_sextupole_optimization/run.jl

using Bell, Printf, LinearAlgebra, TaylorSeries

include(joinpath(@__DIR__, "..", "common", "print_util.jl"))
include(joinpath(@__DIR__, "..", "common", "reference_values.jl"))

banner("Sextupole families and example optimisation")

of = convert_optim(OPTIM_FILE, STRUCTURE_FILE)
lat0 = load_structure(STRUCTURE_FILE)
ref = lat0.ref
C = circumference(lat0)
po = periodic_optics(ring_matrix(lat0), ref, C)
opt = optics_along(lat0, po; ds = 0.01)


section("1. Placement",
        "Maxima of D_x, beta_x, beta_y with D_x > 0.5 m.")

lat, fam_idx, rows, notes = place_sextupole_families(lat0, opt)

@printf("  %-4s %-7s %9s %9s %-5s %9s %8s %8s %8s\n",
        "fam", "max of", "s_max", "value", "slot", "s_slot", "beta_x", "beta_y", "D_x")
for r in rows
    @printf("  %-4s %-7s %9.3f %9.4f %-5s %9.3f %8.3f %8.3f %8.3f\n",
            r.family, r.quantity, r.s_max, r.value, r.slot, r.s_slot, r.betx, r.bety, r.D)
end

println()
foreach(n -> println("  note: ", n), notes)

println()
for f in SEXT_FAMILIES
    nmag = length(fam_idx[f])
    L = isempty(fam_idx[f]) ? 0 : lat.elements[fam_idx[f][1]].L
    println("  family $f: $nmag magnets, L = $L m")
end

@printf("  circumference unchanged: %.9f m\n", circumference(lat))


section("2. Response matrix (per unit k2 = 1 m^-3)",
        "Linear in k2: xi1, alpha1, eta1, D2.")

t0 = time()
base, cols = response(lat, fam_idx)
@printf("  (computed in %.1f s)\n", time() - t0)

println()
@printf("  %-8s %12s %12s %12s %12s %12s\n",
        "family", "dxi1_x", "dxi1_y", "dalpha1", "deta1", "dD2(0) [m]")
for f in SEXT_FAMILIES
    haskey(cols, f) || continue
    c = cols[f][1]
    @printf("  %-8s %12.5f %12.5f %12.5f %12.5f %12.5f\n",
            f, c[:xi1x], c[:xi1y], c[:alpha1], c[:eta1], c[:D2_0])
end


section("3. Solutions",
        "A: xi1_x = xi1_y = alpha1 = 0",
        "B: xi1_x = xi1_y = 0, min rms D2")

kA = solve_A(base, cols)
kB = solve_B(base, cols)

for (nm, k) in (("A", kA), ("B", kB))
    parts = String[]
    for f in SEXT_FAMILIES
        haskey(k, f) || continue
        k2L = k[f] * lat.elements[fam_idx[f][1]].L
        Bpp = k[f] * ref.brho_Tm
        push!(parts, "$f: k2 = $(round(k[f], sigdigits=6)) m^-3 (k2*L = $(round(k2L, sigdigits=5)) m^-2, B'' = $(round(Bpp, sigdigits=5)) T/m^2)")
    end
    println("  solution $nm: ", join(parts, ";  "))
end

mA = ring_metrics(set_family_strengths(lat, fam_idx, kA))
mB = ring_metrics(set_family_strengths(lat, fam_idx, kB))


section("4. Before / after",
        "All from the order-3 one-turn map; xi2 = d2Q/ddelta2.")

@printf("  %-14s %14s %14s %14s\n", "quantity", "natural", "A (xi,alpha1)", "B (xi,D2rms)")
for (lab, k) in [("Qx frac", :Qx), ("Qy frac", :Qy),
                 ("xi1_x", :xi1x), ("xi1_y", :xi1y),
                 ("xi2_x", :xi2x), ("xi2_y", :xi2y),
                 ("D(0) [m]", :D1_0), ("D2(0) [m]", :D2_0),
                 ("max|D2| [m]", :D2max), ("rms D2 [m]", :D2rms),
                 ("alpha0", :alpha0), ("alpha1", :alpha1),
                 ("eta0", :eta0), ("eta1", :eta1)]
    @printf("  %-14s %14.6f %14.6f %14.6f\n",
            lab, getfield(base, k), getfield(mA, k), getfield(mB, k))
end

println()
println("  (linear-model prediction check A: xi1_x, xi1_y, alpha1 should be ~0)")


section("5. Tracking check of solution A")

latA = set_family_strengths(lat, fam_idx, kA)
MA = mA.M

function exact_orbit(L, M, dl)
    u = closed_orbit(M, dl)
    J = [evaluate(differentiate(M[i], j), [u; 0.0; dl]) for i in 1:4, j in 1:4] - I
    for _ in 1:10
        z = exact_track_turn([u; 0.0; dl], L)
        du = -J \ (z[1:4] - u)
        u += du
        maximum(abs, du) < 1e-15 && break
    end
    return u
end

dls = [-1.5e-3, -1e-3, -0.5e-3, 0.0, 0.5e-3, 1e-3, 1.5e-3]
qx = Float64[]
qy = Float64[]

for dl in dls
    u = exact_orbit(latA, MA, dl)
    z = [u; 0.0; dl] .+ [1e-5, 0, 1e-5, 0, 0, 0]
    xx = Float64[]
    yy = Float64[]
    for n in 1:256
        z = exact_track_turn(z, latA)
        push!(xx, z[1] - u[1])
        push!(yy, z[3] - u[3])
    end
    push!(qx, measure_tune(xx))
    push!(qy, measure_tune(yy))
end

V2 = hcat(ones(7), dls, dls .^ 2, dls .^ 3)
fx = V2 \ qx
fy = V2 \ qy

@printf("  tracking fit A: xi1_x = %.5f  xi2_x = %.3f   xi1_y = %.5f  xi2_y = %.3f\n",
        fx[2], 2fx[3], fy[2], 2fy[3])
@printf("  map         A: xi1_x = %.5f  xi2_x = %.3f   xi1_y = %.5f  xi2_y = %.3f\n",
        mA.xi1x, mA.xi2x, mA.xi1y, mA.xi2y)

outdir = joinpath(@__DIR__, "output")
mkpath(outdir)

p1 = plot_chromatic(
    [(title = "solution A (xi1 = 0, alpha1 = 0)",
      d = dls, qx_track = qx, qy_track = qy,
      qx_poly = x -> mA.Qx + mA.xi1x * x + mA.xi2x * x^2 / 2,
      qy_poly = x -> mA.Qy + mA.xi1y * x + mA.xi2y * x^2 / 2)],
    joinpath(outdir, "chromaticity_corrected.png"))

sext_s = Float64[]
s = 0.0
for e in lat.elements
    global s
    (e isa Sextupole && e.name in SEXT_FAMILIES) && push!(sext_s, s + e.L / 2)
    s += elem_length(e)
end

p2 = plot_dispersion2(
    [base.along, mA.along, mB.along],
    ["natural", "A: xi=0, alpha1=0", "B: xi=0, min rms D2"],
    sext_s,
    joinpath(outdir, "sextupoles_dispersion2.png"))

println()
println("  written: ", abspath(p1))
println("  written: ", abspath(p2))

println()
println("Done.")
