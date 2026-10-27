# =============================================================================
# 04_optional_sextupole_optimization.jl - OPTIONAL: sextupole families + example correction
#   (the only script that calls Optimizers.jl; all other scripts use the lattice
#    as given, i.e. sextupoles at zero strength)
# =============================================================================
#
# WHAT
#   1. Closed optics of magnetic_2; local maxima of D_x, beta_x, beta_y in the
#      dispersive arcs (D_x > 0.5 m); nearest short slot (the existing 12 cm
#      sextupole slots or short drifts) gets a family:
#        SXD at D_x maxima, SXB at beta_x maxima, SYB at beta_y maxima.
#   2. Response of xi1_x, xi1_y, alpha1, eta1, D2(0) to each family's k2
#      ("what each family controls").
#   3. Solution A: xi1_x = xi1_y = 0 and alpha1 = 0.
#      Solution B: xi1_x = xi1_y = 0 and rms D2 over the ring minimal.
#   4. Before/after table; D, D2 plot; chromaticity tracking check for A.
#
# WHY
#   Chromatic decoherence ~ xi1, xi2; delta^2 path/orbit terms ~ alpha1, D2.
#   Sextupoles act chromatically through x = D delta, so the efficiency for
#   xi_x ~ beta_x D, for xi_y ~ beta_y D, and for D2 / alpha1 ~ D^2 (and D^3).
#
# HOW TO RUN:  julia --project=. scripts/04_sextupoles.jl
# WHAT YOU SEE: placement table, response matrix, before/after table,
#   output/sextupoles_dispersion2.png, output/chromaticity_corrected.png
# =============================================================================
include(joinpath(@__DIR__, "..", "src", "LieRing.jl"))
using .LieRing, Printf, LinearAlgebra, TaylorSeries
include(joinpath(@__DIR__, "reference_values.jl"))
section(t, lines...) = (println("\n--- ", t, " ", "-"^max(0, 66 - length(t))); foreach(l -> println("  > ", l), lines))
println("="^72); println(" LieRing OPTIONAL - sextupole families and example optimisation (targets A/B)"); println("="^72)
of = parse_optim(OPTIM_FILE); lat0 = of.lattice; ref = lat0.ref; C = circumference(lat0)
po = periodic_optics(ring_matrix(lat0), ref, C); opt = optics_along(lat0, po; ds = 0.01)

section("1. Placement", "Maxima of D_x, beta_x, beta_y with D_x > 0.5 m (sextupoles need D != 0 to act chromatically).")
lat, fam_idx, rows, notes = place_sextupole_families(lat0, opt)
@printf("  %-4s %-7s %9s %9s %-5s %9s %8s %8s %8s\n", "fam", "max of", "s_max", "value", "slot", "s_slot", "beta_x", "beta_y", "D_x")
for r in rows
    @printf("  %-4s %-7s %9.3f %9.4f %-5s %9.3f %8.3f %8.3f %8.3f\n", r.family, r.quantity, r.s_max, r.value, r.slot, r.s_slot, r.betx, r.bety, r.D)
end
foreach(n -> println("  note: ", n), notes)
for f in SEXT_FAMILIES
    println("  family $f: $(length(fam_idx[f])) magnets, L = $(isempty(fam_idx[f]) ? 0 : lat.elements[fam_idx[f][1]].L) m")
end
@printf("  circumference unchanged: %.9f m\n", circumference(lat))

section("2. Response matrix (per unit k2 = 1 m^-3 for the whole family)",
        "Linear in k2: xi1, alpha1, eta1, D2.  Shows what each family controls.")
t0 = time()
base, cols = response(lat, fam_idx)
@printf("  (computed in %.1f s)\n", time() - t0)
@printf("  %-8s %12s %12s %12s %12s %12s\n", "family", "dxi1_x", "dxi1_y", "dalpha1", "deta1", "dD2(0) [m]")
for f in SEXT_FAMILIES
    haskey(cols, f) || continue
    c = cols[f][1]
    @printf("  %-8s %12.5f %12.5f %12.5f %12.5f %12.5f\n", f, c[:xi1x], c[:xi1y], c[:alpha1], c[:eta1], c[:D2_0])
end

section("3. Solutions", "A: xi1_x = xi1_y = alpha1 = 0;  B: xi1_x = xi1_y = 0, min rms D2.")
kA = solve_A(base, cols); kB = solve_B(base, cols)
for (nm, k) in (("A", kA), ("B", kB))
    println("  solution $nm: ", join(["$f: k2 = $(round(k[f], sigdigits=6)) m^-3 (k2*L = $(round(k[f]*lat.elements[fam_idx[f][1]].L, sigdigits=5)) m^-2, B'' = $(round(k[f]*ref.brho_Tm, sigdigits=5)) T/m^2)" for f in SEXT_FAMILIES if haskey(k, f)], ";  "))
end
mA = ring_metrics(set_family_strengths(lat, fam_idx, kA))
mB = ring_metrics(set_family_strengths(lat, fam_idx, kB))

section("4. Before / after", "All from the order-3 one-turn map; xi2 = d2Q/ddelta2.")
@printf("  %-14s %14s %14s %14s\n", "quantity", "natural", "A (xi,alpha1)", "B (xi,D2rms)")
for (lab, k) in [("Qx frac", :Qx), ("Qy frac", :Qy), ("xi1_x", :xi1x), ("xi1_y", :xi1y), ("xi2_x", :xi2x), ("xi2_y", :xi2y),
                 ("D(0) [m]", :D1_0), ("D2(0) [m]", :D2_0), ("max|D2| [m]", :D2max), ("rms D2 [m]", :D2rms),
                 ("alpha0", :alpha0), ("alpha1", :alpha1), ("eta0", :eta0), ("eta1", :eta1)]
    @printf("  %-14s %14.6f %14.6f %14.6f\n", lab, getfield(base, k), getfield(mA, k), getfield(mB, k))
end
println("  (linear-model prediction check A: xi1_x, xi1_y, alpha1 should be ~0 -> see column A)")

section("5. Tracking check of solution A", "Exact-H tracking, tunes vs delta, as in script 03.")
latA = set_family_strengths(lat, fam_idx, kA); MA = mA.M
function exact_orbit(L, M, dl)
    u = closed_orbit(M, dl)
    J = [evaluate(differentiate(M[i], j), [u; 0.0; dl]) for i in 1:4, j in 1:4] - I
    for _ in 1:10
        z = exact_track_turn([u; 0.0; dl], L); du = -J \ (z[1:4] - u); u += du
        maximum(abs, du) < 1e-15 && break
    end
    u
end
dls = [-1.5e-3, -1e-3, -0.5e-3, 0.0, 0.5e-3, 1e-3, 1.5e-3]
qx = Float64[]; qy = Float64[]
for dl in dls
    u = exact_orbit(latA, MA, dl); z = [u; 0.0; dl] .+ [1e-5, 0, 1e-5, 0, 0, 0]; xx = Float64[]; yy = Float64[]
    for n in 1:256
        z = exact_track_turn(z, latA); push!(xx, z[1] - u[1]); push!(yy, z[3] - u[3])
    end
    push!(qx, measure_tune(xx)); push!(qy, measure_tune(yy))
end
V2 = hcat(ones(7), dls, dls .^ 2, dls .^ 3); fx = V2 \ qx; fy = V2 \ qy
@printf("  tracking fit A: xi1_x = %.5f  xi2_x = %.3f   xi1_y = %.5f  xi2_y = %.3f\n", fx[2], 2fx[3], fy[2], 2fy[3])
@printf("  map         A: xi1_x = %.5f  xi2_x = %.3f   xi1_y = %.5f  xi2_y = %.3f\n", mA.xi1x, mA.xi2x, mA.xi1y, mA.xi2y)
outdir = joinpath(@__DIR__, "..", "output"); mkpath(outdir)
p1 = plot_chromatic([(title = "solution A (xi1 = 0, alpha1 = 0)", d = dls, qx_track = qx, qy_track = qy,
                      qx_poly = x -> mA.Qx + mA.xi1x * x + mA.xi2x * x^2 / 2, qy_poly = x -> mA.Qy + mA.xi1y * x + mA.xi2y * x^2 / 2)],
                    joinpath(outdir, "chromaticity_corrected.png"))
sext_s = Float64[]; s = 0.0
for e in lat.elements
    global s
    (e isa Sextupole && e.name in SEXT_FAMILIES) && push!(sext_s, s + e.L / 2); s += elem_length(e)
end
p2 = plot_dispersion2([base.along, mA.along, mB.along], ["natural", "A: xi=0, alpha1=0", "B: xi=0, min rms D2"], sext_s,
                      joinpath(outdir, "sextupoles_dispersion2.png"))
println("  written: ", abspath(p1)); println("  written: ", abspath(p2))
println("\nDone.")
