# =============================================================================
# 09_lie_normal_form.jl - one-turn Lie product with thin RF, homological NF
# =============================================================================
include(joinpath(@__DIR__, "..", "src", "LieRing.jl"))
include(joinpath(@__DIR__, "reference_values.jl"))
using .LieRing, Printf, LinearAlgebra, TaylorSeries

println("="^72)
println(" One-turn Lie factor and normal form (homological operator on generators)")
println("="^72)

of = parse_optim(OPTIM_FILE)
ref = reference_from_kinetic(DEUTERON, 270.0)
lat = Lattice(of.lattice.name, of.lattice.elements, ref, 5)
C = circumference(lat)
cav = RFCavity("RF", 0.2, 5, Float64(pi))

t0 = time()
F = one_turn_factor_rf(lat, cav, ref, C)
@printf("  one-turn factor in %.2f s; |f3|=%.3e  |f4|=%.3e  |R-I| linear\n",
        time() - t0, maxabs(F.f3), maxabs(F.f4))
@printf("  |const of factor map| = %.3e (origin must stay)\n",
        maximum(abs(constant_term(m)) for m in factor_to_map(F)))

t0 = time()
nf = dragt_finn_nf6d(F)
@printf("  NF in %.2f s\n", time() - t0)
@printf("  Q = (%.6f, %.6f, %.6f)\n", (nf.mu ./ 2pi)...)
@printf("  |A^{-1} R A - rot| = %.3e   |A'SA - S| = %.3e\n", nf.Rn_residual, nf.symp_residual)
@printf("  |f3| before = %.3e   after F3 = %.3e\n", nf.f3_before, nf.f3_after)
@printf("  |f4 - h4| after F4 = %.3e\n", nf.f4_nonres_after)
@printf("  d (torus, from T) = %.8e  %.8e  %.8e\n", nf.d...)
@printf("  d (grid check)    = %.8e  %.8e  %.8e\n", nf.d_grid...)
@printf("  |d - d_grid|      = %.3e\n", maximum(abs, nf.d - nf.d_grid))

Mm = ring_taylor_map(lat)
pe = path_expansion(Mm, ref, C)
xi = chromaticity_orders(Mm)
d_note = [2pi / C * xi.x[2] / pe.eta0, 2pi / C * xi.y[2] / pe.eta0,
          -pe.eta1 * nf.sigma_delta_s^2 / pe.eta0]
@printf("  d (TOF formula)   = %.8e  %.8e  %.8e\n", d_note...)
@printf("  rel d − formula   = %.3e  %.3e  %.3e\n", ((nf.d ./ d_note) .- 1)...)

# f4' sign check: residual of cubics must vanish; opposite Poisson would not
Rrot = rotation_matrix_6(nf.mu)
f3n = pull_generator(F.f3, nf.A)
f3res = f3n + nf.F3 - pull_generator(nf.F3, Rrot)
@printf("  cubic homological residual |f3+F3-F3∘r| = %.3e\n", maxabs(f3res))

println("\nDone.")
