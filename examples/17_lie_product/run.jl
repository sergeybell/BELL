# 17_lie_product.jl
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# Product of Dragt-Finn factors through the Poisson bracket.
# concatenate(Fa, Fb) is checked against the Lie series of a pair of elements.
# Sign of the BCH term (1/2)[f3b o Ra, f3a] is fixed by this comparison.
#
#   julia --project=examples examples/17_lie_product/run.jl

using Bell, Printf, LinearAlgebra, TaylorSeries

include(joinpath(@__DIR__, "..", "common", "print_util.jl"))

err(a, b) = maxabs(a - b)
maperr(A, B) = maximum(err(a, b) for (a, b) in zip(A, B))

banner("Lie product of Dragt-Finn factors to degree 4")

ref = reference_from_kinetic(DEUTERON, 270.0)
d  = Drift("D", 2.0)
qf = Quadrupole("QF", 0.4,  0.3)
qd = Quadrupole("QD", 0.4, -0.3)
sx = Sextupole("SX", 0.2,  8.0)
sb = SBend("B", 1.5, 0.05 * 1.5, 0.0)

function factor_of(e)
    M = lie_map(hamiltonian(e, ref), elem_length(e))
    F = dragt_finn_factor(M; name = e.name)
    rec = maperr(factor_to_map(F), M)
    return M, F, rec
end


section("1. Single element: Lie series vs Dragt-Finn reconstruction")

@printf("  %-6s %12s %12s %12s\n", "elem", "recon", "|f3|", "|f4|")
for e in (d, qf, qd, sx, sb)
    M, F, rec = factor_of(e)
    @printf("  %-6s %12.3e %12.3e %12.3e\n", e.name, rec, maxabs(F.f3), maxabs(F.f4))
end


function check_pair(ea, eb; label = ea.name * "+" * eb.name)
    Ma, Fa, _ = factor_of(ea)
    Mb, Fb, _ = factor_of(eb)
    Mc = compose_maps(Mb, Ma)
    Fc = dragt_finn_factor(Mc; name = label)
    Fp = concatenate(Fa, Fb; name = label)
    # Opposite BCH sign, for diagnosis only.
    g3 = pull_generator(Fb.f3, Fa.R)
    g4 = pull_generator(Fb.f4, Fa.R)
    f4m = hom(Fa.f4 + g4 + poisson(g3, Fa.f3) / 2, 4)
    return (
        label = label,
        f3 = err(Fp.f3, Fc.f3),
        f4 = err(Fp.f4, Fc.f4),
        f4opp = err(f4m, Fc.f4),
        R = maximum(abs, Fp.R - Fc.R),
        map = maperr(factor_to_map(Fp), Mc),
    )
end


section("2. Pair of elements: concatenate vs composed Lie series")

@printf("  %-12s %12s %12s %12s %12s %12s\n",
        "pair", "|Δf3|", "|Δf4|", "|Δf4 opp|", "|ΔR|", "|Δmap|")
pairs = ((d, qf), (qf, qd), (qf, sx), (sx, d), (sb, qf), (qf, sb), (sx, sx))
for (ea, eb) in pairs
    r = check_pair(ea, eb)
    @printf("  %-12s %12.3e %12.3e %12.3e %12.3e %12.3e\n",
            r.label, r.f3, r.f4, r.f4opp, r.R, r.map)
end


section("3. Three elements (fold)")

Ma, Fa, _ = factor_of(qf)
Mb, Fb, _ = factor_of(sx)
Mc, Fc, _ = factor_of(qd)
Mabc = compose_maps(Mc, compose_maps(Mb, Ma))
Fabc = product_factors((Fa, Fb, Fc); name = "QF+SX+QD")
Ffrom = dragt_finn_factor(Mabc; name = "comp")

@printf("  |Δf3| = %.3e  |Δf4| = %.3e  |Δmap| = %.3e\n",
        err(Fabc.f3, Ffrom.f3), err(Fabc.f4, Ffrom.f4),
        maperr(factor_to_map(Fabc), Mabc))

println()
println("Done.")
