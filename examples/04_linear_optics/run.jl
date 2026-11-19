# 04_linear_optics/run.jl
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# Linear optics. lattices/fodo.bell comes first: element matrices, the one-turn
# matrix, periodic Twiss, phase advance, and dispersion. lattices/magnetic_2.bell
# follows: the same readings on a ring that also has sextupoles, and the
# normal-form residual of its one-turn matrix. Twiss figures go under output/.
#
#   julia --project=examples examples/04_linear_optics/run.jl

using Bell, Printf, LinearAlgebra, Plots
using Bell.NormalForm

include(joinpath(@__DIR__, "..", "common", "print_util.jl"))

banner("Linear optics of the FODO line")

lat = load_structure(joinpath(@__DIR__, "..", "lattices", "fodo.bell"))
ref = lat.ref
C = path_length(lat)
bend = total_bend_angle(lat)


section("1. Line",
        "Sixteen FODO cells. periodic_optics is the periodic solution of the transfer matrix.",
        "The sum of the bending angles is printed next to 2 pi.")

@printf("  elements              : %d\n", length(lat.elements))
@printf("  path length C         : %.7f m\n", C)
@printf("  total bending angle   : %.9f rad  (%.6f deg)    2 pi = %.9f rad\n",
        bend, rad2deg(bend), 2pi)
@printf("  bending angle / 2 pi  : %.6f\n", bend / 2pi)
println(describe_reference(ref))


section("2. Element matrices",
        "R = exp(L S A), A = Hessian of H at z = 0.",
        "symplectic error is max |R' S R - S|.")

@printf("  %-6s %-12s %8s %8s %8s %8s %8s %8s %10s\n",
        "name", "kind", "R11", "R12", "R21", "R33", "R16", "R56", "sympl.err")
for e in unique_elements(lat)
    R = linear_matrix(e, ref)
    @printf("  %-6s %-12s %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %10.1e\n",
            e.name, kind_string(e), R[1, 1], R[1, 2], R[2, 1], R[3, 3], R[1, 6], R[5, 6],
            symplectic_error(R))
end


section("3. One-turn matrix, stability, periodic optics",
        "A plane is stable when |Tr| < 2, with Tr the trace of its 2x2 block.",
        "periodic_optics reports mu in [0, 2 pi). The integer tune is accumulated later.")

R = ring_matrix(lat)
po = periodic_optics(R, ref, C)
trx = R[1, 1] + R[2, 2]
try_ = R[3, 3] + R[4, 4]

println("  one-turn matrix (x, px, y, py, l, delta):")
for i in 1:6
    println("    ", join([@sprintf("% .6e", R[i, j]) for j in 1:6], " "))
end
println()
@printf("  symplectic error of R        : %.3e\n", symplectic_error(R))
@printf("  x-y coupling                 : %.3e\n", po.coupling)
@printf("  Tr_x = % .8f    |Tr_x| < 2 : %s     cos(mu_x) = %.8f\n",
        trx, abs(trx) < 2 ? "yes" : "no", po.x.cosmu)
@printf("  Tr_y = % .8f    |Tr_y| < 2 : %s     cos(mu_y) = %.8f\n",
        try_, abs(try_) < 2 ? "yes" : "no", po.y.cosmu)
@printf("  stable                       : %s\n",
        po.x.stable && po.y.stable ? "yes" : "no")
println()
@printf("  beta_x  = %8.4f m    alpha_x = % .6e    mu_x = %.6f rad\n",
        po.x.beta, po.x.alpha, po.x.mu)
@printf("  beta_y  = %8.4f m    alpha_y = % .6e    mu_y = %.6f rad\n",
        po.y.beta, po.y.alpha, po.y.mu)
@printf("  fractional tunes from mu/2pi : Qx = %.8f    Qy = %.8f\n",
        po.x.mu / 2pi, po.y.mu / 2pi)
println()
@printf("  D_x  = % .6e m    Dpx = % .6e\n", po.D4[1], po.D4[2])
@printf("  D_y  = % .6e m    Dpy = % .6e\n", po.D4[3], po.D4[4])
println()
@printf("  (l|delta)_closed = % .6f m\n", po.l_delta_closed)
@printf("  alpha_c = %.8f    gamma_t = %.6f    1/gamma^2 = %.8f    slip eta = %.8f\n",
        po.alpha_c, po.gamma_t, 1 / ref.gamma^2, po.eta_slip)


section("4. One cell and the 16-cell line",
        "The file repeats one cell sixteen times. The cell matrix has the same Twiss",
        "at the start; its phase advance is one sixteenth of the accumulated tune.")

ncell = length(lat.elements) ÷ 16
length(lat.elements) == 16 * ncell || error("expected 16 identical cells")
for i in 1:length(lat.elements)
    isequal(lat.elements[i], lat.elements[mod1(i, ncell)]) ||
        error("element $i does not repeat the first cell")
end
cell = Lattice(lat.name * "_cell", lat.elements[1:ncell], ref, lat.harmonic)
Rc = ring_matrix(cell)
Cc = path_length(cell)
poc = periodic_optics(Rc, ref, Cc)
opt_cell = optics_along(cell, poc)

@printf("  cell length = %.4f m    cell bend = %.6f rad\n", Cc, total_bend_angle(cell))
@printf("  cell Qx = %.8f    cell Qy = %.8f\n", opt_cell.Qx, opt_cell.Qy)
@printf("  16 * cell tune          : Qx = %.8f    Qy = %.8f\n",
        16 * opt_cell.Qx, 16 * opt_cell.Qy)


section("5. Optics along the 16-cell line")

opt = optics_along(lat, po)
cl = opt.closure
@printf("  accumulated tunes       : Qx = %.8f    Qy = %.8f\n", opt.Qx, opt.Qy)
@printf("  difference from 16*cell : dQx = %.3e    dQy = %.3e\n",
        opt.Qx - 16 * opt_cell.Qx, opt.Qy - 16 * opt_cell.Qy)
println()
@printf("  closure  rel beta_x %.1e, alpha_x %.1e, rel beta_y %.1e, alpha_y %.1e\n",
        cl.betx, cl.alfx, cl.bety, cl.alfy)
@printf("           D %.1e, D' %.1e, fractional tune %.1e / %.1e\n",
        cl.D, cl.Dp, cl.Qx_frac, cl.Qy_frac)
println("  periodic Twiss closes: ", maximum(values(cl)) < 1e-8 ? "yes" : "no")
println()
@printf("  max beta_x = %.4f m    max beta_y = %.4f m\n",
        maximum(opt.data.betx), maximum(opt.data.bety))
@printf("  max D_x    = %.4f m    min D_x    = %.4f m\n",
        maximum(opt.data.D), minimum(opt.data.D))

function compaction_integral(along)
    acc = 0.0
    s = along.data.s
    D = along.data.D
    for k in 2:length(s)
        s_mid = (s[k] + s[k - 1]) / 2
        idx = findfirst(b -> b[1] <= s_mid < b[2], along.bounds)
        e = along.bounds[idx][3]
        acc += curvature(e) * (D[k] + D[k - 1]) / 2 * (s[k] - s[k - 1])
    end
    return acc / along.path_length
end

println()
@printf("  alpha_c from (1/C) int h D ds = %.8f\n", compaction_integral(opt))
@printf("  alpha_c from periodic_optics  = %.8f\n", po.alpha_c)


section("6. Twiss figure")

outdir = joinpath(@__DIR__, "output")
mkpath(outdir)
png = plot_twiss(opt, joinpath(outdir, "fodo_twiss.png");
                 title = "FODO, 16 cells: periodic linear optics")
println("  written: ", abspath(png))


section("7. magnetic_2",
        "Same lattice file, now with sextupoles. The calls are the ones used on the FODO line.")

mag_file = joinpath(@__DIR__, "..", "lattices", "magnetic_2.bell")
mag = load_structure(mag_file)
Cm = path_length(mag)
println("  file : ", abspath(mag_file))
println(describe_reference(mag.ref))
println()
for e in unique_elements(mag)
    n = count(isequal(e), mag.elements)
    @printf("  %-6s x%-3d %s\n", e.name, n, string(e))
end
println()
@printf("  elements              : %d\n", length(mag.elements))
@printf("  path length C         : %.7f m\n", Cm)
@printf("  total bending angle   : %.9f rad  (%.6f deg)\n",
        total_bend_angle(mag), rad2deg(total_bend_angle(mag)))
@printf("  dipoles: %d, quadrupoles: %d, sextupoles: %d\n",
        count(e -> e isa SBend, mag.elements),
        count(e -> e isa Quadrupole, mag.elements),
        count(e -> e isa Sextupole, mag.elements))


section("8. magnetic_2 element matrices",
        "R = exp(L S A), A = Hessian of H at z = 0. Sextupoles are included.")

@printf("  %-6s %-12s %8s %8s %8s %8s %8s %8s %10s\n",
        "name", "kind", "R11", "R12", "R21", "R33", "R16", "R56", "sympl.err")
for e in unique_elements(mag)
    local Relem = linear_matrix(e, mag.ref)
    @printf("  %-6s %-12s %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %10.1e\n",
            e.name, kind_string(e), Relem[1, 1], Relem[1, 2], Relem[2, 1], Relem[3, 3],
            Relem[1, 6], Relem[5, 6], symplectic_error(Relem))
end


section("9. magnetic_2 one-turn matrix",
        "rotation_residual is |A^-1 R A - rotation| from the linear normal form.")

Rm = ring_matrix(mag)
pom = periodic_optics(Rm, mag.ref, Cm)
trxm = Rm[1, 1] + Rm[2, 2]
trym = Rm[3, 3] + Rm[4, 4]
println("  one-turn matrix (x, px, y, py, l, delta):")
for i in 1:6
    println("    ", join([@sprintf("% .6e", Rm[i, j]) for j in 1:6], " "))
end
println()
@printf("  symplectic error of R        : %.3e\n", symplectic_error(Rm))
@printf("  normal-form residual         : %.3e\n", rotation_residual(Rm, pom))
@printf("  x-y coupling                 : %.3e\n", pom.coupling)
@printf("  Tr_x = % .8f    |Tr_x| < 2 : %s     cos(mu_x) = %.8f\n",
        trxm, abs(trxm) < 2 ? "yes" : "no", pom.x.cosmu)
@printf("  Tr_y = % .8f    |Tr_y| < 2 : %s     cos(mu_y) = %.8f\n",
        trym, abs(trym) < 2 ? "yes" : "no", pom.y.cosmu)
@printf("  stable                       : %s\n",
        pom.x.stable && pom.y.stable ? "yes" : "no")
println()
@printf("  beta_x  = %8.4f m    alpha_x = % .6e    mu_x = %.6f rad\n",
        pom.x.beta, pom.x.alpha, pom.x.mu)
@printf("  beta_y  = %8.4f m    alpha_y = % .6e    mu_y = %.6f rad\n",
        pom.y.beta, pom.y.alpha, pom.y.mu)
@printf("  fractional tunes from mu/2pi : Qx = %.8f    Qy = %.8f\n",
        pom.x.mu / 2pi, pom.y.mu / 2pi)
println()
@printf("  D_x  = % .6e m    Dpx = % .6e\n", pom.D4[1], pom.D4[2])
@printf("  D_y  = % .6e m    Dpy = % .6e\n", pom.D4[3], pom.D4[4])
println()
@printf("  (l|delta)_closed = % .6f m\n", pom.l_delta_closed)
@printf("  alpha_c = %.8f    gamma_t = %.6f    1/gamma^2 = %.8f    slip eta = %.8f\n",
        pom.alpha_c, pom.gamma_t, 1 / mag.ref.gamma^2, pom.eta_slip)


section("10. magnetic_2 optics along the ring")

optm = optics_along(mag, pom)
clm = optm.closure
@printf("  accumulated tunes       : Qx = %.8f    Qy = %.8f\n", optm.Qx, optm.Qy)
println()
@printf("  closure  rel beta_x %.1e, alpha_x %.1e, rel beta_y %.1e, alpha_y %.1e\n",
        clm.betx, clm.alfx, clm.bety, clm.alfy)
@printf("           D %.1e, D' %.1e, fractional tune %.1e / %.1e\n",
        clm.D, clm.Dp, clm.Qx_frac, clm.Qy_frac)
println("  periodic Twiss closes: ", maximum(values(clm)) < 1e-8 ? "yes" : "no")
println()
@printf("  max beta_x = %.4f m    max beta_y = %.4f m\n",
        maximum(optm.data.betx), maximum(optm.data.bety))
@printf("  max D_x    = %.4f m    min D_x    = %.4f m\n",
        maximum(optm.data.D), minimum(optm.data.D))
println()
@printf("  alpha_c from (1/C) int h D ds = %.8f\n", compaction_integral(optm))
@printf("  alpha_c from periodic_optics  = %.8f\n", pom.alpha_c)

png_m = plot_twiss(optm, joinpath(outdir, "magnetic_2_twiss.png");
                   title = "magnetic_2: periodic linear optics")
println()
println("  written: ", abspath(png_m))

println()
println("Done.")
