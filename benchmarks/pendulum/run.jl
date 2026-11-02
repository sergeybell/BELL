# Pendulum error against number of steps for every integrator in Bell.
#   julia --project=. benchmarks/pendulum/run.jl
using Printf
using Bell
using Bell.Integrators

pend(z, s) = [z[2], -sin(z[1])]
energy(z) = z[2]^2 / 2 - cos(z[1])
const theta = 0.3
const z0 = [theta, 0.0]
const E0 = energy(z0)
const period = 2pi * (1 + theta^2 / 16 + 11 * theta^4 / 3072)
const flows = SplitFlows((z, h) -> [z[1] + h * z[2], z[2]],
                         (z, h) -> [z[1], z[2] - h * sin(z[1])])

function report(name, n, z)
    dE = abs(energy(z) - E0)
    dist = hypot(z[1] - z0[1], z[2] - z0[2])
    @printf("  %-22s n=%4d  dE=%.3e  dist=%.3e\n", name, n, dE, dist)
end

println("mathematical pendulum, theta = $theta, period = $period")
for n in (25, 50, 100, 200)
    println("steps $n")
    report("RK4", n, integrate(pend, z0, 0.0, period; method = RK4(), nsteps = n).z)
    report("ImplicitMidpoint", n, integrate(pend, z0, 0.0, period; method = ImplicitMidpoint(), nsteps = n).z)
    report("ImplicitYoshida(4)", n, integrate(pend, z0, 0.0, period; method = ImplicitYoshida(4), nsteps = n).z)
    report("Yoshida(4) split", n, integrate(flows, z0, 0.0, period; method = Yoshida(4), nsteps = n).z)
    report("DormandPrince5", n, integrate(pend, z0, 0.0, period; method = DormandPrince5(), nsteps = n).z)
    report("DormandPrince853", n, integrate(pend, z0, 0.0, period; method = DormandPrince853(), nsteps = n).z)
    report("GBS(4)", n, integrate(pend, z0, 0.0, period; method = GraggBulirschStoer(4), nsteps = n).z)
end
