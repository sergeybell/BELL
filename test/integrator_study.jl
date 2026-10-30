# test/integrator_study.jl - comparison of Bell.Integrators methods on the Henon-Heiles system
#   H = (px^2 + py^2)/2 + (x^2 + y^2)/2 + x^2 y - y^3/3,   split H = T(p) [drift] + V(q) [nonlinear kick]
# Run:  julia --project=. test/integrator_study.jl     -> prints table, writes test/output/integrator_study.md
using LinearAlgebra, Printf
using Bell, Bell.DiffAlgebra, Bell.Integrators

const NF = Ref(0)
H(z) = (z[3]^2 + z[4]^2) / 2 + (z[1]^2 + z[2]^2) / 2 + z[1]^2 * z[2] - z[2]^3 / 3
rhs(z, s) = (NF[] += 1; [z[3], z[4], -z[1] - 2z[1] * z[2], -z[2] - z[1]^2 + z[2]^2])
drift(z, h) = [z[1] + h * z[3], z[2] + h * z[4], z[3], z[4]]
kick(z, h) = (NF[] += 1; [z[1], z[2], z[3] - h * (z[1] + 2z[1] * z[2]), z[4] - h * (z[2] + z[1]^2 - z[2]^2)])
const FL = SplitFlows(drift, kick)
const z0 = [0.1, -0.2, 0.25, 0.1]            # E ~ 0.0622 < 1/6: bounded, mildly chaotic region
run(m, z, T, n) = m isa SplittingMethod ? integrate(FL, z, 0.0, T; method = m, nsteps = n).z :
                                          integrate(rhs, z, 0.0, T; method = m, nsteps = n).z
methods = [("Leapfrog", Leapfrog()), ("Yoshida4 (Forest-Ruth)", Yoshida(4)), ("Yoshida6", Yoshida(6)),
           ("Yoshida8", Yoshida(8)), ("RK4", RK4()), ("DP5", DormandPrince5()), ("DOP853", DormandPrince853()),
           ("GBS(4) ord 8", GraggBulirschStoer(4))]
T = 10.0
zref = integrate(rhs, z0, 0.0, T; method = GraggBulirschStoer(8), tol = 1e-16).z
S = [0 0 1 0; 0 0 0 1; -1 0 0 0; 0 -1 0 0]
sp = DASpace(4, 1)
rows = String[]
println(rpad("method", 24), "  err h=0.1  err h=0.05  ord(.4/.2) |JtSJ-S|  maxdE(T=1e4,h=0.1)  fev/step  time(T=1e4) s")
for (name, m) in methods
    e1 = norm(run(m, z0, T, 100) - zref); e2 = norm(run(m, z0, T, 200) - zref)
    p = log2(norm(run(m, z0, T, 25) - zref) / norm(run(m, z0, T, 50) - zref))   # order from h = 0.4 / 0.2
    zd = [da_var(sp, i; at = z0[i]) for i in 1:4]
    J = linear_part(run(m, zd, T, 100) .- 0.0)  # Jacobian of the T=10 map from the DA result
    serr = opnorm(J' * S * J - S)
    # long-term: 1e4 time units, h = 0.1, energy sampled every 10 units
    E0 = H(z0); z = copy(z0); dE = 0.0; NF[] = 0
    t = @elapsed for k in 1:1000
        z = run(m, z, 10.0, 100); dE = max(dE, abs(H(z) - E0))
    end
    fps = NF[] / 1e5
    line = @sprintf("%-24s  %9.2e  %9.2e  %5.2f  %9.2e  %9.2e  %6.1f  %7.2f", name, e1, e2, p, serr, dE, fps, t)
    println(line); push!(rows, line)
end
mkpath(joinpath(@__DIR__, "output"))
open(joinpath(@__DIR__, "output", "integrator_study.md"), "w") do io
    println(io, "| method | err h=0.1 | err h=0.05 | order (h 0.4→0.2) | ‖JᵀSJ−S‖ | max ΔE (T=1e4) | f/step | time, s |")
    println(io, "|---|---|---|---|---|---|---|---|")
    for l in rows
        f = split(l, r"\s{2,}"); println(io, "| ", join(f, " | "), " |")
    end
end