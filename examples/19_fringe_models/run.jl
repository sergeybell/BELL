# 19_fringe_models/run.jl
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# Fringe-field models of Bell on the ring magnetic_2.  The lattice file is the
# same for every row; only the FringeModel passed at build time changes.
#
#   1. Catalogue of models: NoEdge, HardEdge (linear, nonlinear), SoftEdge
#      with every profile (Enge, Tanh, LinearRamp, Tabulated) and every
#      method (Splitting, Magnus, TaylorIntegrate).
#   2. Hard-edge limit: Enge edge with gap -> 0 tends to HardEdge(order = 2).
#   3. Fringe integral: the linear vertical focusing of a soft edge equals
#      the hard edge with FINT = int F(1-F) dz / gap.
#   4. Convergence in the number of slices n for each method and order.
#   5. Symplecticity of the edge maps.
#
#   julia --project=examples examples/19_fringe_models/run.jl

using Bell, Printf, Plots

include(joinpath(@__DIR__, "..", "common", "print_util.jl"))

ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")
outdir = joinpath(@__DIR__, "output"); mkpath(outdir)

lat = load_structure(joinpath(@__DIR__, "..", "lattices", "magnetic_2.bell"))
nbend = count(e -> e isa SBend, lat.elements)
bm = first(e for e in lat.elements if e isa SBend)

banner("Fringe-field models of Bell (lattice $(lat.name))")
@printf("  %d sector bends %s: L = %.4f m, angle = %.6f rad, rho = %.4f m, e1 = e2 = %g\n",
        nbend, bm.name, bm.L, bm.angle, 1 / curvature(bm), bm.e1)
note("Every table row names the FringeModel used (dipole / quadrupole edge).")

obs(m) = fringe_observables(lat, m)
row(label, o) = @printf("  %-46s %9.6f %9.6f %10.6f %10.6f %8.1e %8.1e %6.1f\n",
                        label, o.Qx, o.Qy, o.xix, o.xiy, o.symp_edge_lin, o.symp_edge_pt, o.seconds)
header() = (@printf("  %-46s %9s %9s %10s %10s %8s %8s %6s\n", "model", "Qx frac", "Qy frac",
                    "xi_x", "xi_y", "symp(R)", "symp(z0)", "t [s]"); println("  ", "-"^112))

# ---------------------------------------------------------------------------
G = 0.05                                  # full gap of the Enge edge [m]
I = fringe_integral(Enge(gap = G))        # int F(1-F) dz
profiles = [("Enge(gap=50 mm)", Enge(gap = G)),
            ("Tanh(lambda=2I)", Tanh(2I)),
            ("LinearRamp(width=6I)", LinearRamp(6I)),
            ("Tabulated(Enge 50 mm, 401 pts)",
             (Fe = profile_function(Enge(gap = G))[1]; zt = collect(range(-0.4, 0.4, length = 401));
              Tabulated(zt, Fe.(zt))))]
methods = [("Splitting(n=20,order=4)", Splitting(n = 20, order = 4)),
           ("Magnus(n=40,order=4)", Magnus(n = 40, order = 4)),
           ("TaylorIntegrate(tol=1e-10)", TaylorIntegrate(tol = 1e-10))]

section("1. Catalogue of edge models",
        "Same lattice, model chosen at build time.  Tanh and LinearRamp are scaled to the",
        @sprintf("same fringe integral I = int F(1-F) dz = %.6f m as the Enge edge (FINT = %.5f).", I, I / G),
        "symp(R): symplecticity error of the linear part of the edge maps;",
        "symp(z0): of the edge-map Jacobian at z0 = (1e-3, 1e-4, 1e-3, 1e-4, 0, 1e-3).")
header()
results = Dict{String,Any}()
for (lab, m) in [("NoEdge (default FringeModel())", FringeModel()),
                 ("HardEdge(order=1)", FringeModel(dipole = HardEdge(order = 1))),
                 ("HardEdge(order=1, FINT of Enge 50 mm)", FringeModel(dipole = HardEdge(order = 1, gap = G, fint = I / G))),
                 ("HardEdge(order=2)", FringeModel(dipole = HardEdge(order = 2))),
                 ("HardEdge(order=3)", FringeModel(dipole = HardEdge(order = 3))),
                 ("HardEdge(order=2, FINT of Enge 50 mm)", FringeModel(dipole = HardEdge(order = 2, gap = G, fint = I / G))),
                 ("HardEdge(2) + quad HardEdge(2)", FringeModel(dipole = HardEdge(order = 2), quad = HardEdge(order = 2)))]
    local o = obs(m); results[lab] = o; row(lab, o)
end
println("  ", "-"^112)
for (pl, p) in profiles, (ml, me) in methods
    lab = "Soft " * pl * " / " * ml
    local o = obs(FringeModel(dipole = SoftEdge(profile = p, method = me)))
    results[lab] = o
    @printf("  %s\n", lab); row("", o)
end
osq = obs(FringeModel(dipole = SoftEdge(profile = Enge(gap = G)), quad = HardEdge(order = 2)))
row("Soft Enge 50 mm + quad HardEdge(2)", osq)
note("Quadrupole hard edges are degree-4 generators: they change amplitude detuning,")
note("not tunes or first-order chromaticity.")

# ---------------------------------------------------------------------------
section("2. Hard-edge limit: Enge gap -> 0",
        "SoftEdge(Enge(gap), Splitting(n=20, order=4)) against HardEdge(order=2).")
hard = results["HardEdge(order=2)"]
gaps = [0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001]
@printf("  %10s %10s %12s %12s %12s %12s\n", "gap [m]", "FINT", "Qy frac", "xi_y", "dQy-hard", "dxi_y-hard")
gq, gx = Float64[], Float64[]
for g in gaps
    local o = obs(FringeModel(dipole = SoftEdge(profile = Enge(gap = g), method = Splitting(n = 20, order = 4))))
    push!(gq, o.Qy - hard.Qy); push!(gx, o.xiy - hard.xiy)
    @printf("  %10.4f %10.5f %12.8f %12.6f %12.2e %12.2e\n", g, fringe_integral(Enge(gap = g)) / g,
            o.Qy, o.xiy, o.Qy - hard.Qy, o.xiy - hard.xiy)
end
@printf("  HardEdge(order=2):                %12.8f %12.6f\n", hard.Qy, hard.xiy)
note("Both differences fall linearly with the gap: the soft edge tends to the hard edge.")

# ---------------------------------------------------------------------------
section("3. Linear focusing of a soft edge = hard edge with FINT",
        "psi = e - h g FINT (1 + sin^2 e)/cos e with FINT from the Enge profile.")
hf = results["HardEdge(order=1, FINT of Enge 50 mm)"]
se = results["Soft Enge(gap=50 mm) / Splitting(n=20,order=4)"]
@printf("  Qx: soft %.8f   hard+FINT %.8f   NoEdge %.8f\n", se.Qx, hf.Qx, results["NoEdge (default FringeModel())"].Qx)
@printf("  Qy: soft %.8f   hard+FINT %.8f   NoEdge %.8f\n", se.Qy, hf.Qy, results["NoEdge (default FringeModel())"].Qy)
note("A straight EFB seen from the curved orbit gives a gradient h0 F' theta whose integral")
note("h0^2 int F(1-F) dz is the FINT term; the horizontal weak-focusing deficit cancels.")

# ---------------------------------------------------------------------------
section("4. Convergence in the number of slices / sub-intervals n",
        "Enge 50 mm.  Reference: Splitting(n=128, order=4).  |xi_y - ref| and |Qy - ref|.")
refo = obs(FringeModel(dipole = SoftEdge(profile = Enge(gap = G), method = Splitting(n = 128, order = 4))))
conv = Dict{String,Tuple{Vector{Int},Vector{Float64}}}()
@printf("  %-22s %6s %14s %14s\n", "method", "n", "|dxi_y|", "|dQy|")
for (lab, mk, ns) in [("Splitting order 2", n -> Splitting(n = n, order = 2), [2, 4, 8, 16, 32, 64]),
                      ("Splitting order 4", n -> Splitting(n = n, order = 4), [2, 4, 8, 16, 32]),
                      ("Magnus order 2", n -> Magnus(n = n, order = 2), [2, 4, 8, 16, 32, 64]),
                      ("Magnus order 4", n -> Magnus(n = n, order = 4), [2, 4, 8, 16, 32])]
    errs = Float64[]
    for n in ns
        local o = obs(FringeModel(dipole = SoftEdge(profile = Enge(gap = G), method = mk(n))))
        push!(errs, max(abs(o.xiy - refo.xiy), 1e-16))
        @printf("  %-22s %6d %14.3e %14.3e\n", lab, n, abs(o.xiy - refo.xiy), abs(o.Qy - refo.Qy))
    end
    conv[lab] = (ns, errs)
end
note("Errors fall with n; at small n the slices do not resolve the profile and the")
note("decrease is irregular, then order 4 converges faster than order 2.")

# ---------------------------------------------------------------------------
section("5. Method agreement (Enge 50 mm)")
for (ml, _) in methods
    local o = results["Soft Enge(gap=50 mm) / " * ml]
    @printf("  %-28s xi_y = %.8f   diff to Splitting(128,4) = %.2e\n", ml, o.xiy, o.xiy - refo.xiy)
end
note("TaylorIntegrate (RK4 on the DA vector field) is not symplectic by construction;")
note("its symp(z0) column above shows the size of that defect.")

# ---------------------------------------------------------------------------
section("6. Plots")
zs = range(-0.3, 0.3, length = 600)
p1 = plot(xlabel = "z from EFB [m]", ylabel = "F(z)", title = "Edge profiles (equal fringe integral)", legend = :topleft)
for (pl, p) in profiles
    F = profile_function(p)[1]; plot!(p1, zs, F.(zs), label = pl)
end
plot!(p1, zs, z -> z >= 0 ? 1.0 : 0.0, label = "hard edge", ls = :dash, color = :black)
f1 = joinpath(outdir, "fringe_profiles.png"); savefig(p1, f1); println("  written: ", f1)

labs = ["NoEdge (default FringeModel())", "HardEdge(order=1)", "HardEdge(order=2)", "HardEdge(order=2, FINT of Enge 50 mm)",
        ["Soft " * pl * " / Splitting(n=20,order=4)" for (pl, _) in profiles]...]
short = ["NoEdge", "Hard 1", "Hard 2", "Hard 2+FINT", "Soft Enge", "Soft Tanh", "Soft Linear", "Soft Tab."]
p2 = scatter(1:length(labs), [results[l].xiy for l in labs], xticks = (1:length(labs), short), xrotation = 30,
             ylabel = "xi_y", title = "Vertical chromaticity by fringe model (Bell)", label = "xi_y", ms = 6)
f2 = joinpath(outdir, "fringe_xi_by_model.png"); savefig(p2, f2); println("  written: ", f2)

p3 = plot(gaps, abs.(gx), xscale = :log10, yscale = :log10, marker = :o, xlabel = "Enge gap [m]",
          ylabel = "|difference to HardEdge(2)|", title = "Hard-edge limit of the Enge edge", label = "|xi_y - xi_y hard|")
plot!(p3, gaps, abs.(gq), marker = :s, label = "|Qy - Qy hard|")
f3 = joinpath(outdir, "fringe_gap_limit.png"); savefig(p3, f3); println("  written: ", f3)

p4 = plot(xscale = :log10, yscale = :log10, xlabel = "n", ylabel = "|xi_y - reference|",
          title = "Convergence of soft-edge methods (Enge 50 mm)", legend = :bottomleft)
for (lab, (ns, errs)) in sort(collect(conv))
    plot!(p4, ns, errs, marker = :o, label = lab)
end
f4 = joinpath(outdir, "fringe_convergence.png"); savefig(p4, f4); println("  written: ", f4)

section("Summary",
        "NoEdge reproduces the plain ring; HardEdge(1) changes nothing for e = 0 without FINT;",
        "the nonlinear hard edge adds the chromatic term h y^2 x'/2 at each edge; soft edges",
        "converge to it as gap -> 0 and add the FINT focusing; all Lie-based maps are symplectic.")
