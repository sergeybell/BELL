# =============================================================================
# 06_coordinates.jl - user (x,x',y,y',phi,delta) <-> canonical (x,px,y,py,l,delta)
# =============================================================================
#
# WHAT
#   1. Shows that Constants.canonical_to_user / user_to_canonical are EXACT
#      nonlinear formulas (sqrt, 1/v(delta)) and that they are applied to the
#      TaylorN variables, i.e. as full polynomials truncated only at MAP_ORDER.
#   2. Round trip T o T^-1 - identity as polynomials: residual coefficients up to
#      MAP_ORDER must be 0 (machine precision); the first nonzero coefficients
#      appear only at order MAP_ORDER+1 (truncation, printed).
#   3. Point residual with the exact (non-polynomial) functions at large
#      amplitude, and the error of the TRUNCATED polynomial at the same point.
#   4. A few one-turn map coefficients in both systems (ring at 270 MeV).
#   5. Symplecticity: J^T S J - S of the one-turn map in canonical variables
#      (~1e-15) and in user variables (O(1): user variables are not canonical).
#
# HOW TO RUN: julia --project=. scripts/06_coordinates.jl
# WHAT YOU SEE: printed tables only.
# =============================================================================
include(joinpath(@__DIR__, "..", "src", "LieRing.jl"))
using .LieRing, Printf, LinearAlgebra, TaylorSeries
include(joinpath(@__DIR__, "reference_values.jl"))
section(t, lines...) = (println("\n--- ", t, " ", "-"^max(0, 66 - length(t))); foreach(l -> println("  > ", l), lines))
println("="^72); println(" LieRing - exact coordinate transformation user <-> canonical"); println("="^72)

of = parse_optim(OPTIM_FILE); lat = of.lattice; ref = lat.ref; C = circumference(lat); h = 5
section("1. Formulas used (src/Constants.jl)",
        "canonical -> user: x' = px/ps, y' = py/ps, ps = sqrt((1+delta)^2 - px^2 - py^2),",
        "                   phi = 2 pi h l / (C v(delta)/v0)",
        "user -> canonical: px = x'(1+delta)/sqrt(1+x'^2+y'^2), py likewise, l = phi C (v/v0) / (2 pi h)",
        "v(delta)/v0 = (1+delta)/sqrt(1 + beta0^2 ((1+delta)^2 - 1))  (exact, Constants.velocity_ratio)",
        "Both functions are generic Julia code: evaluated on TaylorN they give the full Taylor",
        "expansion (no linearisation); user_map truncates the result at MAP_ORDER = $(MAP_ORDER).")
T  = canonical_to_user(collect(Z), ref, C, h)
Ti = user_to_canonical(collect(Z), ref, C, h)
println("  x'(px,py,delta) to order 3 : ", poly_string(truncate_to(T[2], 3)))
println("  phi(l,delta)    to order 3 : ", poly_string(truncate_to(T[5], 3)))
nonlin = maximum(maxabs(T[i] - truncate_to(T[i], 1)) for i in 1:6)
@printf("  largest nonlinear coefficient in T: %.3e  (> 0 => conversion is nonlinear)\n", nonlin)

section("2. Polynomial round trip T(T^-1(u)) - u", "residual per component, split by order")
RT = [evaluate(T[i], Ti) for i in 1:6]
for i in 1:6
    r = RT[i] - Z[i]
    low = maxabs(truncate_to(r, MAP_ORDER)); hi = maxabs(r - truncate_to(r, MAP_ORDER))
    @printf("  component %d: max|coef| orders <= %d: %.2e   orders > %d (kept up to JET_ORDER): %.2e\n", i, MAP_ORDER, low, MAP_ORDER, hi)
end
section("3. Point residuals at large amplitude", "u = (5 mm, 20 mrad, 5 mm, 20 mrad, 1.0 rad, 1e-2)")
u = [5e-3, 0.02, 5e-3, 0.02, 1.0, 1e-2]
ue = canonical_to_user(user_to_canonical(u, ref, C, h), ref, C, h)
@printf("  exact functions round trip     : max|u - T(T^-1 u)| = %.2e\n", maximum(abs, ue - u))
zc = user_to_canonical(u, ref, C, h)
Tp = [truncate_to(c, MAP_ORDER) for c in T]
up = [evaluate(Tp[i], zc) for i in 1:6]
@printf("  polynomial T (order %d) vs exact T at that point: max diff = %.2e (truncation error, O(amp^%d))\n", MAP_ORDER, maximum(abs, up - u), MAP_ORDER + 1)
@printf("  compare: linearised T would give max diff = %.2e\n", maximum(abs, [evaluate(truncate_to(T[i], 1), zc) for i in 1:6] - u))

section("4. One-turn map coefficients in both systems (270 MeV, RF off)")
M = ring_taylor_map(lat); U = user_map(M, ref, C, h)
rows = [("x | px  vs x | x'", 1, [0,1,0,0,0,0]), ("x | delta", 1, [0,0,0,0,0,1]),
        ("px | delta  vs  x' | delta", 2, [0,0,0,0,0,1]), ("l | delta  vs  phi | delta", 5, [0,0,0,0,0,1]),
        ("l | delta^2  vs phi | delta^2", 5, [0,0,0,0,0,2]), ("l | px^2  vs  phi | x'^2", 5, [0,2,0,0,0,0]),
        ("px | px delta vs x'|x' delta", 2, [0,1,0,0,0,1]), ("x | px^2 delta vs x|x'^2 delta", 1, [0,2,0,0,0,1])]
@printf("  %-32s %16s %16s\n", "coefficient", "canonical", "user")
for (nm, i, e) in rows
    @printf("  %-32s %16.8e %16.8e\n", nm, coeff(M[i], e), coeff(U[i], e))
end
println("  (phi coefficients include the factor 2 pi h / C and the v(delta) terms)")
section("5. Symplecticity", "max |J^T S J - S| of the map Jacobian at z = (1e-3, 1e-3, 1e-3, 1e-3, 0, 1e-3)")
z0 = [1e-3, 1e-3, 1e-3, 1e-3, 0.0, 1e-3]
J(Mp) = [evaluate(differentiate(Mp[i], j), z0) for i in 1:6, j in 1:6]
@printf("  canonical map: %.2e  (zero up to truncation at order %d)\n", symplectic_error(J(M)), MAP_ORDER)
@printf("  user map     : %.2e  (not small: (x,x',y,y',phi,delta) is not a canonical set)\n", symplectic_error(J(U)))
println("\nDone.")
