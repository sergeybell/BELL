# 10_linear_optics.jl
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# Closed linear optics of lattices/magnetic_2.bell.
# OptiM is converted first; the structure file is what Bell loads.
#
#   julia --project=examples examples/10_linear_optics/run.jl

using Bell, Printf, LinearAlgebra

include(joinpath(@__DIR__, "..", "common", "print_util.jl"))
include(joinpath(@__DIR__, "..", "common", "reference_values.jl"))

banner("Linear optics of magnetic_2")


section("1. Lattice file")

of = convert_optim(OPTIM_FILE, STRUCTURE_FILE)
lat = load_structure(STRUCTURE_FILE)

println("  file                 : ", of.path)
println("  header               : Energy = $(of.energy_MeV) MeV (kinetic), Mass = $(of.mass_MeV) MeV")
println("  sequence tokens      : $(length(of.sequence)), 'Number of periods' in file = $(of.nperiods_file), repetitions used = $(of.nrepeat_used)")

foreach(w -> println("  WARNING: ", w), of.warnings)


section("2. Reference particle")

println(describe_reference(lat.ref))

@printf("  (physical deuteron mass would be %.8f MeV -> B*rho = %.7f T*m; not used)\n",
        DEUTERON.mass_MeV, reference_from_kinetic(DEUTERON, of.energy_MeV).brho_Tm)


section("3. Elements")

for e in unique_elements(lat)
    n = count(isequal(e), lat.elements)
    @printf("  %-5s x%-3d %s\n", e.name, n, string(e))
end

C = circumference(lat)

println()
@printf("  number of elements    : %d\n", length(lat.elements))
@printf("  circumference C       : %.7f m\n", C)
@printf("  total bending angle   : %.9f deg\n", rad2deg(total_bend_angle(lat)))
@printf("  dipoles: %d, quads: %d, sextupoles: %d\n",
        count(e -> e isa SBend, lat.elements),
        count(e -> e isa Quadrupole, lat.elements),
        count(e -> e isa Sextupole, lat.elements))


section("4. Element linear maps",
        "R = exp(L S A), A = Hessian of H at z=0; symplectic error |R'SR - S|.")

for e in unique_elements(lat)
    R = linear_matrix(e, lat.ref)
    @printf("  %-5s  R11=% .6f R12=% .6f R21=% .6f R33=% .6f R16=% .6f R56=% .6f  sympl.err=%.1e\n",
            e.name, R[1,1], R[1,2], R[2,1], R[3,3], R[1,6], R[5,6], symplectic_error(R))
end


section("5. One-turn matrix and closed optics")

R = ring_matrix(lat)
po = periodic_optics(R, lat.ref, C)

println("  one-turn 6x6 matrix (canonical x, px, y, py, l, delta):")
for i in 1:6
    println("    ", join([@sprintf("% .6e", R[i, j]) for j in 1:6], " "))
end

println()
@printf("  symplectic error of one-turn matrix : %.2e\n", symplectic_error(R))
@printf("  x-y coupling (max off-diag block)    : %.2e\n", po.coupling)
@printf("  stability: cos(mu_x) = %.8f, cos(mu_y) = %.8f  -> %s\n",
        po.x.cosmu, po.y.cosmu,
        po.x.stable && po.y.stable ? "STABLE" : "UNSTABLE")

println()
@printf("  beta_x  = %.6f m   alpha_x = % .6e\n", po.x.beta, po.x.alpha)
@printf("  beta_y  = %.6f m   alpha_y = % .6e\n", po.y.beta, po.y.alpha)
@printf("  D_x     = % .6e m   D'_x = % .6e\n", po.D4[1], po.D4[2])
@printf("  D_y     = % .3e m   D'_y = % .3e\n", po.D4[3], po.D4[4])

println()
@printf("  (l|delta)_closed = % .6f m\n", po.l_delta_closed)
@printf("  alpha_c = %.8f   gamma_t = %.6f   1/gamma0^2 = %.8f   slip eta = %.8f\n",
        po.alpha_c, po.gamma_t, 1 / lat.ref.gamma^2, po.eta_slip)
@printf("  normal-form check: |A^-1 R A - rotation| = %.2e\n", rotation_residual(R, po))


section("6. Optics along the ring + CLOSED check")

opt = optics_along(lat, po)
cl = opt.closure

@printf("  rel. beta_x %.1e, alpha_x %.1e, rel. beta_y %.1e, alpha_y %.1e, D %.1e, D' %.1e, Qfrac %.1e/%.1e\n",
        cl.betx, cl.alfx, cl.bety, cl.alfy, cl.D, cl.Dp, cl.Qx_frac, cl.Qy_frac)

closed = maximum(values(cl)) < 1e-8
println("  CLOSED: ", closed ? "yes (all closure errors < 1e-8)" : "NO")

println()
@printf("  tunes Qx = %.8f   Qy = %.8f  (integer from accumulated phase advance)\n", opt.Qx, opt.Qy)
@printf("  max beta_x = %.4f m, max beta_y = %.4f m, max D = %.4f m, min D = %.4f m\n",
        maximum(opt.data.betx), maximum(opt.data.bety), maximum(opt.data.D), minimum(opt.data.D))

# Independent check: alpha_c = (1/C) int h D ds over bends.
acint = 0.0
for k in 2:length(opt.data.s)
    s_mid = (opt.data.s[k] + opt.data.s[k-1]) / 2
    idx = findfirst(b -> b[1] <= s_mid < b[2], opt.bounds)
    e = opt.bounds[idx][3]
    global acint += curvature(e) * (opt.data.D[k] + opt.data.D[k-1]) / 2 * (opt.data.s[k] - opt.data.s[k-1])
end

println()
@printf("  alpha_c from (1/C) int h*D ds (trapezoid)  = %.8f  (independent check)\n", acint / C)


section("7. Plot")

outdir = joinpath(@__DIR__, "output")
mkpath(outdir)
png = plot_twiss(opt, joinpath(outdir, "magnetic_2_twiss.png"); title = "magnetic_2: closed optics (Bell)")
println("  written: ", abspath(png))


section("8. Comparison with OptiM / COSY")

ini = of.initial
rows = [
    ("beta_x(s=0) [m]", po.x.beta, ini["BetaX[cm]"] / 100, "OptiM Initial"),
    ("beta_y(s=0) [m]", po.y.beta, ini["BetaY[cm]"] / 100, "OptiM Initial"),
    ("alpha_x(s=0)", po.x.alpha, ini["AlfaX"], "OptiM Initial"),
    ("alpha_y(s=0)", po.y.alpha, ini["AlfaY"], "OptiM Initial"),
    ("D_x(s=0) [m]", po.D4[1], ini["DispersX[cm]"] / 100, "OptiM Initial"),
    ("D'_x(s=0)", po.D4[2], ini["Dsp_PrimeX"], "OptiM Initial"),
    ("Qx fractional", mod(opt.Qx, 1), COSY_REF.Qx_frac[1], "COSY wp_natural.dat"),
    ("Qy fractional", mod(opt.Qy, 1), COSY_REF.Qy_frac[1], "COSY wp_natural.dat"),
    ("alpha_c", po.alpha_c, COSY_REF.alpha0[1], "COSY integrals_Icirc.dat"),
]

@printf("  %-18s %18s %18s %12s  %s\n", "quantity", "Bell", "reference", "diff", "source")
for (n, a, b, src) in rows
    @printf("  %-18s %18.10g %18.10g %12.3e  %s\n", n, a, b, a - b, src)
end

println()
println("Done.")
