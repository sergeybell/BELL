# =============================================================================
# 06_poisson_bracket.jl - Poisson brackets and Lie operators
# =============================================================================
#
# WHAT THIS EXAMPLE TEACHES
#   How Hamiltonian dynamics can be done with algebra on polynomials: the
#   Poisson bracket, conserved quantities, and the Lie transformation
#   exp(:-tH:) which turns a Hamiltonian H into its time-t map (a Taylor map
#   as in example 05).  This is the language of Lie-algebraic beam dynamics
#   (Dragt, Forest, COSY INFINITY).
#
# PHASE SPACE AND POISSON BRACKET
# -------------------------------
# One degree of freedom: canonical coordinate q and momentum p.
# For two functions f(q,p), g(q,p) the Poisson bracket is
#
#     [f, g] = df/dq * dg/dp - df/dp * dg/dq
#
# Properties: [q, p] = 1, antisymmetry [f,g] = -[g,f], linearity, the Leibniz
# rule [f, g*h] = [f,g]*h + g*[f,h], and the Jacobi identity
# [f,[g,h]] + [g,[h,f]] + [h,[f,g]] = 0.
#
# Hamilton's equations can be written with the bracket:  dz/dt = [z, H]
# (z = q or p; check: [q,H] = dH/dp, [p,H] = -dH/dq).  A function F(q,p)
# is conserved,  dF/dt = [F, H] = 0,  iff its bracket with H vanishes.
#
# LIE OPERATOR AND LIE TRANSFORMATION
# -----------------------------------
# The Lie operator of f acts on g as        :f: g = [f, g]
# and its exponential is the series         exp(:f:) g = g + [f,g] + [f,[f,g]]/2! + ...
#
# SIGN CONVENTION used everywhere below:
#     dz/dt = [z, H] = [-H, z] = :-H: z     ==>    z(t) = exp(:-t H:) z
# i.e. exp(:-tH:) applied to the coordinate functions q, p gives the
# position and momentum after time t, as functions (polynomials) of the
# initial q, p.  So exp(:-tH:) IS the time-t flow (the Taylor map) of H.
#
# Order of products: exp(:f:) exp(:g:) q = exp(:f:) ( exp(:g:) q ).  This is
# the map "first exp(:f:), then exp(:g:)" (since exp(:f:) G(q,p) = G(exp(:f:)q,
# exp(:f:)p) - a Lie transformation acts on a function by acting on its arguments).
#
# Implementation: q and p are TaylorN polynomials (see example 05); the bracket
# is a one-liner with differentiate.
#
# SECTIONS
#   1. Basic brackets           - [q,p] = 1, antisymmetry, linearity, Leibniz, Jacobi
#   2. Conserved quantities     - [F, H] = 0 for the harmonic oscillator
#   3. Free particle            - Lie series terminates: q -> q + t p
#   4. Harmonic oscillator      - Lie series sums to cos/sin; error vs # of terms
#   5. Anharmonic oscillator    - H = p^2/2 + q^4/4, Lie series vs RK4 integration
#   6. Non-commutativity        - exp(:f:)exp(:g:) != exp(:g:)exp(:f:), BCH term
#   7. Symplecticity            - [Q, P] = det J = 1 for Hamiltonian maps
#
# HOW TO RUN (from the Julia folder)
#   julia --project=. examples/06_poisson_bracket.jl
#   (needs TaylorSeries, listed in Project.toml)
#
# EXPECTED OUTPUT / HOW TO READ IT
#   Polynomials in q, p printed in ASCII ("1.0*q + 2.0*p").  Identities print
#   0 (or 1 for [Q,P]); series results are printed next to the exact formula;
#   the error tables should decrease quickly with the number of terms K.

using TaylorSeries, Printf

const N = 20                                             # truncation order in (q, p)
const SPACE = JetSpace(order = N, variables = [:q, :p])
const q, p = variables(SPACE)

"Poisson bracket [f, g] (variable 1 = q, variable 2 = p)."
pb(f, g) = differentiate(f, 1) * differentiate(g, 2) - differentiate(f, 2) * differentiate(g, 1)

"exp(:f:) g, summed up to the term :f:^K g / K! (stops early if the terms vanish)."
function lie_exp(f, g; K = 30)
    result, term = g, g
    for n in 1:K
        term = pb(f, term) / n           # term = :f:^n g / n!
        iszero(term) && break
        result += term
    end
    return result
end

# --- ASCII printing of polynomials (TaylorSeries' own uses Unicode) ---------
function str(f; digits = 6)
    terms = String[]
    for d in 0:N, i in d:-1:0
        c = getcoeff(f, [i, d - i])
        iszero(c) && continue
        mono = join([k == 1 ? v : "$v^$k" for (v, k) in (("q", i), ("p", d - i)) if k > 0], "*")
        cs = string(round(c, sigdigits = digits))
        push!(terms, isempty(mono) ? cs : cs * "*" * mono)
    end
    isempty(terms) && return "0"
    return replace(join(terms, " + "), "+ -" => "- ")
end

# -----------------------------------------------------------------------------
# Printing helpers (ASCII only - the Windows console garbles Unicode).
#   banner(title, lines...)   big title + "What this example does"
#   section(title, lines...)  section heading + 1-3 lines of explanation
#   print_summary(lines...)   "Summary / what you learned" at the end
# -----------------------------------------------------------------------------
function banner(title, lines...)
    println("="^72)
    println(" ", title)
    println("="^72)
    println("What this example does:")
    foreach(l -> println("  ", l), lines)
end

function section(title, lines...)
    println()
    println("--- ", title, " ", "-"^max(0, 66 - length(title)))
    foreach(l -> println("  > ", l), lines)
    isempty(lines) || println()
end

function print_summary(lines...)
    println()
    println("="^72)
    println(" Summary / what you learned")
    println("="^72)
    foreach(l -> println("  * ", l), lines)
end

banner("Example 06: Poisson brackets and Lie operators (polynomials in q, p)",
       "Poisson bracket [f,g] = df/dq dg/dp - df/dp dg/dq on TaylorN polynomials.",
       "We check its algebraic rules, find conserved quantities ([F,H] = 0), and",
       "compute the time-t map of a Hamiltonian as the Lie series exp(:-tH:) q, p.",
       "Results are compared with exact formulas and a Runge-Kutta integration;",
       "finally we verify symplecticity [Q,P] = 1 of the resulting maps.")

# ---------------------------------------------------------------------------
section("1. Basic Poisson brackets",
        "[f,g] = df/dq*dg/dp - df/dp*dg/dq.  [q,p] = 1 defines canonical variables.",
        "The algebraic rules below must all print exactly 0.")
# ---------------------------------------------------------------------------
println("[q, p] = ", str(pb(q, p)), "    [p, q] = ", str(pb(p, q)),
        "    [q, q] = ", str(pb(q, q)), "    [p, p] = ", str(pb(p, p)))

f = q^2 * p
g = p^3 + q
h = q * p + q^3
println("f = ", str(f), ",   g = ", str(g), ",   h = ", str(h))
println("[f, g]                     = ", str(pb(f, g)))
println("antisymmetry  [f,g]+[g,f]  = ", str(pb(f, g) + pb(g, f)))
println("linearity     [2f+3g,h] - 2[f,h] - 3[g,h] = ", str(pb(2f + 3g, h) - 2pb(f, h) - 3pb(g, h)))
println("Leibniz       [f,g*h] - ([f,g]*h + g*[f,h]) = ", str(pb(f, g * h) - (pb(f, g) * h + g * pb(f, h))))
println("Jacobi        [f,[g,h]] + [g,[h,f]] + [h,[f,g]] = ",
        str(pb(f, pb(g, h)) + pb(g, pb(h, f)) + pb(h, pb(f, g))))

# ---------------------------------------------------------------------------
section("2. Conserved quantities",
        "dF/dt = [F, H]: F is conserved iff [F, H] = 0.  Harmonic oscillator",
        "H = (p^2+q^2)/2: H itself is conserved, q and p are not ([q,H] = p = dq/dt).")
# ---------------------------------------------------------------------------
H = (p^2 + q^2) / 2                      # harmonic oscillator
println("H = ", str(H))
println("[H, H]         = ", str(pb(H, H)), "    (energy is always conserved)")
println("[H, q^2 + p^2] = ", str(pb(H, q^2 + p^2)), "    (q^2 + p^2 = 2H is conserved)")
println("[q, H]         = ", str(pb(q, H)), "    (= dq/dt, q itself is not conserved)")
println("[p, H]         = ", str(pb(p, H)), "    (= dp/dt)")

# ---------------------------------------------------------------------------
section("3. Free particle H = p^2/2: the series terminates",
        "exp(:f:) g = g + [f,g] + [f,[f,g]]/2! + ...  with f = -tH.  Here [f,[f,q]] = 0,",
        "so the series stops and gives the exact drift q -> q + t p, p -> p.")
# ---------------------------------------------------------------------------
t = 2.0
f = -t * p^2 / 2                         # generator -tH
println("t = $t,  f = -t*H = ", str(f))
println(":f: q    = ", str(pb(f, q)))
println(":f:^2 q  = ", str(pb(f, pb(f, q))), "    -> all further terms are 0")
println("exp(:f:) q = ", str(lie_exp(f, q)), "    (exact: q + t*p)")
println("exp(:f:) p = ", str(lie_exp(f, p)), "    (exact: p, momentum conserved)")

# ---------------------------------------------------------------------------
section("4. Harmonic oscillator H = (p^2+q^2)/2: series -> cos, sin",
        "The Lie series no longer terminates; its terms sum to the rotation",
        "q -> q cos t + p sin t.  The table shows the error vs the number of terms K.")
# ---------------------------------------------------------------------------
# :-tH: q = t p,  :-tH: p = -t q, so the series for q gives
# q (1 - t^2/2! + t^4/4! - ...) + p (t - t^3/3! + ...) = q cos t + p sin t.
t = 1.0
f = -t * H
Q, P = lie_exp(f, q), lie_exp(f, p)
println("t = $t,  30 terms:")
println("  exp(:-tH:) q = ", str(Q, digits = 12))
println("  exact        = ", round(cos(t), sigdigits = 12), "*q + ", round(sin(t), sigdigits = 12), "*p")
println("  exp(:-tH:) p = ", str(P, digits = 12))
println("  exact        = -", round(sin(t), sigdigits = 12), "*q + ", round(cos(t), sigdigits = 12), "*p")

println("\nerror vs number of terms K in the series (t = $t):")
@printf("%4s  %12s\n", "K", "max error")
for K in 2:2:18
    QK, PK = lie_exp(f, q; K = K), lie_exp(f, p; K = K)
    err = maximum(abs.([getcoeff(QK, [1, 0]) - cos(t), getcoeff(QK, [0, 1]) - sin(t),
                        getcoeff(PK, [1, 0]) + sin(t), getcoeff(PK, [0, 1]) - cos(t)]))
    @printf("%4d  %12.3e\n", K, err)
end
println("(error ~ t^(K+1)/(K+1)!, the Taylor remainder of cos/sin)")

# ---------------------------------------------------------------------------
section("5. Anharmonic oscillator H = p^2/2 + q^4/4 vs numerical ODE",
        "A nonlinear system: the Lie series gives q(t), p(t) as polynomials in the",
        "initial q, p.  Evaluated at (0.8, 0.3), t = 0.5, it converges to RK4 as K grows.")
# ---------------------------------------------------------------------------
Ha = p^2 / 2 + q^4 / 4
println("H = ", str(Ha))
# The individual terms of the Lie series  Q = sum_k t^k/k! :-H:^k q
println("terms :-H:^k q  (multiply by t^k/k!):")
term = q
for k in 0:4
    println("  k = $k:  ", str(term))
    global term = pb(-Ha, term)
end

# Hand-written 4th-order Runge-Kutta for dq/dt = p, dp/dt = -q^3
rhs(z) = [z[2], -z[1]^3]
function rk4(z, t; steps = 10_000)
    dt = t / steps
    for _ in 1:steps
        k1 = rhs(z); k2 = rhs(z + dt / 2 * k1); k3 = rhs(z + dt / 2 * k2); k4 = rhs(z + dt * k3)
        z = z + dt / 6 * (k1 + 2k2 + 2k3 + k4)
    end
    return z
end

t, z0 = 0.5, [0.8, 0.3]
zref = rk4(z0, t)
println("\nt = $t, initial (q, p) = ", z0, ",  RK4 reference = ", zref)
@printf("%4s  %20s  %20s  %10s\n", "K", "q(t) from series", "p(t) from series", "max error")
# :-H:^k q has degree 2k+1, so K <= 9 terms are computed exactly at order N = 20.
for K in 1:9
    QK, PK = lie_exp(-t * Ha, q; K = K), lie_exp(-t * Ha, p; K = K)
    zK = [evaluate(QK, z0), evaluate(PK, z0)]
    @printf("%4d  %20.15f  %20.15f  %10.2e\n", K, zK..., maximum(abs.(zK - zref)))
end

# ---------------------------------------------------------------------------
section("6. Non-commutativity: exp(:f:)exp(:g:) vs exp(:g:)exp(:f:)",
        "The order of two Lie transformations matters.  By the BCH formula the",
        "difference is, to leading order, :[f,g]: q - the two printed lines must agree.")
# ---------------------------------------------------------------------------
e = 0.1
f, g = e * p^2 / 2, e * q^2 / 2
A = lie_exp(f, lie_exp(g, q))            # exp(:f:) exp(:g:) q
B = lie_exp(g, lie_exp(f, q))            # exp(:g:) exp(:f:) q
println("f = ", str(f), ",  g = ", str(g), "   (eps = $e)")
println("exp(:f:)exp(:g:) q = ", str(A))
println("exp(:g:)exp(:f:) q = ", str(B))
println("difference         = ", str(A - B))
# BCH: exp(:f:)exp(:g:) = exp(:f + g + [f,g]/2 + ...:), so to leading order
# the difference of the two orderings is :[f,g]: q  (order eps^2).
println("[f, g]             = ", str(pb(f, g)))
println(":[f,g]: q          = ", str(pb(pb(f, g), q)), "    (leading BCH term = the difference)")

# ---------------------------------------------------------------------------
section("7. Symplecticity: [Q, P] = 1 and det(Jacobian) = 1",
        "Any Hamiltonian map (q,p) -> (Q,P) keeps [Q,P] = 1 (= det J in 1D, area",
        "preservation), even if nonlinear.  A truncated Lie series violates it slightly.")
# ---------------------------------------------------------------------------
# For a Hamiltonian flow the new coordinates Q(q,p), P(q,p) satisfy [Q,P] = 1.
# In one degree of freedom [Q,P] = dQ/dq dP/dp - dQ/dp dP/dq = det(Jacobian).
jacdet(Q, P) = differentiate(Q, 1) * differentiate(P, 2) - differentiate(Q, 2) * differentiate(P, 1)

Q, P = lie_exp(-1.0 * H, q), lie_exp(-1.0 * H, p)
println("harmonic oscillator, t = 1:   [Q,P] = ", str(pb(Q, P), digits = 15))

# A nonlinear map built from two terminating Lie series:
#   kick  exp(:-tau*q^4/4:):  p -> p - tau*q^3
#   drift exp(:-tau*p^2/2:):  q -> q + tau*p
tau = 0.5
fk, fd = -tau * q^4 / 4, -tau * p^2 / 2
Q = lie_exp(fk, lie_exp(fd, q))          # first kick, then drift
P = lie_exp(fk, lie_exp(fd, p))
println("kick then drift, tau = $tau:")
println("  Q = ", str(Q))
println("  P = ", str(P))
println("  [Q,P] = ", str(pb(Q, P)), "   det J = ", str(jacdet(Q, P)), "   (exactly 1, although the map is nonlinear)")

# Truncated Lie series of the anharmonic oscillator: [Q,P] = 1 + O(t^(K+1)).
z0, t = [0.8, 0.3], 0.5
println("anharmonic oscillator, truncated series, t = $t, [Q,P] - 1 at (q,p) = ", z0, ":")
for K in 1:5        # K <= 5 keeps all products below the truncation order N = 20
    QK, PK = lie_exp(-t * Ha, q; K = K), lie_exp(-t * Ha, p; K = K)
    @printf("  K = %d:  %10.2e\n", K, evaluate(pb(QK, PK), z0) - 1)
end

print_summary("The Poisson bracket encodes Hamilton's equations: dz/dt = [z, H].",
              "exp(:-tH:) applied to q, p gives the time-t map as polynomials (a Taylor map).",
              "The Lie series is exact if it terminates, otherwise it converges like a Taylor series.",
              "Hamiltonian maps are symplectic: [Q, P] = 1; Lie transformations do not commute (BCH).")
