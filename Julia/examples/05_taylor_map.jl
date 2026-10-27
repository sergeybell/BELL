# =============================================================================
# 05_taylor_map.jl - what a Taylor map is and how to work with it
# =============================================================================
#
# WHAT THIS EXAMPLE TEACHES
#   The object at the heart of map-based beam dynamics (COSY INFINITY, DA):
#   a Taylor map = a vector of multivariate polynomials truncated at order N,
#   and how to build, inspect, evaluate, compose and invert such maps with the
#   package TaylorSeries.jl (new API: JetSpace + TaylorN).
#
# TAYLOR MAP
# ----------
# Take a (smooth) map of n variables, z_f = M(z_i), z = (z_1, ..., z_n).
# Expand every component of M in a Taylor series around z_i = 0 and keep all
# terms up to total degree N.  The result is a vector of n multivariate
# polynomials - this vector IS the Taylor map of order N:
#
#     z_f,k = sum over monomials   C(k | e_1 ... e_n) * z_1^e_1 * ... * z_n^e_n ,
#             e_1 + ... + e_n <= N
#
# So the map is stored as a table of numbers: one column per final coordinate
# k, one row per monomial (e_1 ... e_n).  The usual notation for an entry is
#
#     (final | monomial),   e.g. (x|xy) = coefficient of x*y in x_f.
#
# The degree-1 rows form the Jacobian matrix (linear part) of the map;
# higher rows are the nonlinear corrections.
#
# TRUNCATION
# ----------
# All arithmetic is done "modulo degree N+1": after every +, *, /, exp, ...
# the terms of degree > N are simply dropped.  Example (N = 3): x^2 * x^2 = 0.
# Because the coefficients of degree <= N never depend on the dropped terms,
# every result is EXACT up to degree N (this is "truncated power series
# algebra", a.k.a. differential algebra, DA).
#
# In TaylorSeries.jl such a truncated polynomial is a TaylorN object. A JetSpace
# fixes the number of variables and the order N; variables(space) returns the
# identity polynomials x, y, ... from which everything else is built.
#
# SECTIONS
#   1. Polynomials in x, y      - arithmetic, truncation, getcoeff, evaluate
#   2. Taylor expansions        - exp, sin, 1/(1+x), sqrt(1+x) applied to x
#   3. Derivatives / integrals  - differentiate, integrate (and truncation loss)
#   4. A map                    - vector of polynomials, coefficient table
#   5. Composition              - (M o R)(z) = M(R(z)), not commutative
#   6. Inverse map              - inverse_map and a fixed-point iteration
#   7. Jacobian                 - linear part and det J (area change)
#
# HOW TO RUN (from the Julia folder)
#   julia --project=. examples/05_taylor_map.jl
#   (needs TaylorSeries, listed in Project.toml)
#
# EXPECTED OUTPUT / HOW TO READ IT
#   Polynomials are printed in ASCII as "1.0 + 2.0*x - 3.0*y + 0.5*y^2".
#   Maps are printed as a coefficient table: one row per monomial, e.g. the
#   row "(x*y)" and column "y_f" holds the coefficient (y|xy).  Checks print
#   differences that should be exactly 0 (or ~|z|^(N+1) where a truncation
#   effect is demonstrated on purpose).

using TaylorSeries, Printf

const N = 3                                             # truncation order
const SPACE = JetSpace(order = N, variables = [:x, :y])
const x, y = variables(SPACE)

# ---------------------------------------------------------------------------
# Small helpers for readable ASCII output (TaylorSeries' own printing uses
# Unicode superscripts, which some consoles garble).
# ---------------------------------------------------------------------------

"All exponent vectors (e_1, ..., e_nv) with e_1 + ... + e_nv == d (descending)."
function monomials(nv, d)
    nv == 1 && return [[d]]
    return [[i; rest] for i in d:-1:0 for rest in monomials(nv - 1, d - i)]
end

"Monomial name, e.g. [2,1] -> \"x^2*y\", [0,0] -> \"1\"."
function monomial_name(e, names = get_variable_names(SPACE))
    parts = [e[i] == 1 ? names[i] : "$(names[i])^$(e[i])" for i in eachindex(e) if e[i] > 0]
    return isempty(parts) ? "1" : join(parts, "*")
end

"Polynomial as ASCII string, e.g. \"1 + 2*x - 0.5*x*y\"."
function str(f; digits = 6)
    nv, terms = get_numvars(f), String[]
    for d in 0:N, e in monomials(nv, d)
        c = getcoeff(f, e)
        iszero(c) && continue
        cs = string(round(c, sigdigits = digits))
        push!(terms, e == zeros(Int, nv) ? cs : cs * "*" * monomial_name(e))
    end
    isempty(terms) && return "0"
    return replace(join(terms, " + "), "+ -" => "- ")
end

"Coefficient table (final | monomial) of a map M = vector of TaylorN."
function print_map(M, outnames)
    @printf("%-10s", "monomial")
    for o in outnames
        @printf("%14s", o * "_f")
    end
    println()
    for d in 0:N, e in monomials(get_numvars(M[1]), d)
        c = [getcoeff(m, e) for m in M]
        all(iszero, c) && continue               # skip empty rows, as COSY does
        @printf("%-10s", "(" * monomial_name(e) * ")")
        for ci in c
            @printf("%14.6f", ci)
        end
        println()
    end
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

banner("Example 05: Taylor maps with TaylorSeries.jl (order N = $N, variables x, y)",
       "A Taylor map is a vector of polynomials in the initial coordinates, cut",
       "at total degree N; all arithmetic drops terms above N, and everything",
       "kept is exact.  We build polynomials, expand functions, differentiate,",
       "and treat a nonlinear 2D map: coefficient table, composition, inverse",
       "and Jacobian - the basic operations of map-based beam dynamics.")

# ---------------------------------------------------------------------------
section("1. Polynomials in x, y (order N = $N)",
        "x and y are the identity polynomials; +, * work as usual but every term",
        "of degree > $N is dropped.  getcoeff(f, [i,j]) is the coefficient of x^i y^j,",
        "evaluate(f, [x0, y0]) the value at a point.")
# ---------------------------------------------------------------------------
f = 1 + 2x - 3y + x * y + 0.5y^2
g = x - y^2
println("f        = ", str(f))
println("g        = ", str(g))
println("f + g    = ", str(f + g))
println("f * g    = ", str(f * g), "    (degree-4 terms dropped)")
println("coefficient of x*y in f : getcoeff(f, [1,1]) = ", getcoeff(f, [1, 1]))
println("coefficient of y^2 in f : getcoeff(f, [0,2]) = ", getcoeff(f, [0, 2]))
println("constant term of f      : constant_term(f)   = ", constant_term(f))
println("f at (x,y) = (0.1, 0.2) : evaluate(f, [0.1, 0.2]) = ", evaluate(f, [0.1, 0.2]))
println("   direct formula       : ", 1 + 2 * 0.1 - 3 * 0.2 + 0.1 * 0.2 + 0.5 * 0.2^2)

println("\ntruncation at order $N (look for the dropped terms):")
println("x^2       = ", str(x^2))
println("x^2 * x   = ", str(x^2 * x))
println("x^2 * x^2 = ", str(x^2 * x^2), "    (degree 4 > N, so it is dropped)")
println("(1+x)^5   = ", str((1 + x)^5), "    (only 1 + 5x + 10x^2 + 10x^3 survive)")

# ---------------------------------------------------------------------------
section("2. Taylor expansion of functions around 0",
        "Applying exp, sin, ... to the polynomial x gives its Taylor series up to",
        "order $N.  Compared with the textbook series: difference must be 0.")
# ---------------------------------------------------------------------------
# Applying a function to the identity polynomial x gives its Taylor series.
expansions = [
    ("exp(x)",      exp(x),        1 + x + x^2 / 2 + x^3 / 6,          "1 + x + x^2/2 + x^3/6"),
    ("sin(x+y)",    sin(x + y),    (x + y) - (x + y)^3 / 6,            "(x+y) - (x+y)^3/6"),
    ("1/(1+x)",     1 / (1 + x),   1 - x + x^2 - x^3,                  "1 - x + x^2 - x^3"),
    ("sqrt(1+x)",   sqrt(1 + x),   1 + x / 2 - x^2 / 8 + x^3 / 16,     "1 + x/2 - x^2/8 + x^3/16"),
]
for (name, series, known, formula) in expansions
    println(rpad(name, 10), " = ", str(series))
    println(rpad("", 10), "   known: ", formula, ",  max coeff. difference = ",
            norm(series - known, Inf))
end

# ---------------------------------------------------------------------------
section("3. Derivatives and integrals",
        "differentiate(h, i) = dh/dz_i, integrate(h, i) = antiderivative in z_i.",
        "Integration raises the degree, so terms that would exceed N are lost.")
# ---------------------------------------------------------------------------
h = x^3 + 2x^2 * y - y^2 + 4x
println("h          = ", str(h))
println("dh/dx      = ", str(differentiate(h, 1)), "    differentiate(h, 1)")
println("dh/dy      = ", str(differentiate(h, 2)), "    differentiate(h, 2)")
println("int h dx   = ", str(integrate(h, 1)), "    integrate(h, 1); degree > $N dropped")
println("d/dx int h dx = ", str(differentiate(integrate(h, 1), 1)),
        "    (= h minus its degree-$N part: x^3 + 2x^2*y were lost in the integral)")

# ---------------------------------------------------------------------------
section("4. A map = vector of polynomials",
        "M(x, y) = (x + y^2, y + x*y).  Table: row = monomial, column = final",
        "coordinate; entry = (final|monomial).  Degree-1 rows = linear (matrix) part.")
# ---------------------------------------------------------------------------
# A simple nonlinear 2D map:  M(x, y) = (x + y^2,  y + x*y)
M = [x + y^2, y + x * y]
println("M(x,y) = (", str(M[1]), ",  ", str(M[2]), ")")
println("\ncoefficient table (final | monomial):")
print_map(M, ["x", "y"])
println("e.g. (x|y^2) = ", getcoeff(M[1], [0, 2]), ",  (y|x*y) = ", getcoeff(M[2], [1, 1]))
z0 = [0.1, 0.2]
println("\nM applied to the point z0 = ", z0, " :  ", [evaluate(m, z0) for m in M])

# ---------------------------------------------------------------------------
section("5. Composition of maps",
        "(M o R)(z) = M(R(z)): substitute R's polynomials into M (first R, then M).",
        "Composition of a beamline = composition of element maps; order matters.")
# ---------------------------------------------------------------------------
# A second map: rotation by angle th plus a quadratic term.
th = 0.3
c, s = cos(th), sin(th)
R = [c * x + s * y, -s * x + c * y + x^2]
println("R(x,y) = (", str(R[1]), ",  ", str(R[2]), ")")

# (A o B)(z) = A(B(z)): substitute the polynomials of B into the variables of A.
compose(A, B) = [evaluate(a, B) for a in A]

MR = compose(M, R)          # first R, then M
RM = compose(R, M)          # first M, then R
println("\nM o R  (first R, then M):")
print_map(MR, ["x", "y"])
println("\nR o M  (first M, then R):")
print_map(RM, ["x", "y"])
println("\nmax |coeff(M o R) - coeff(R o M)| = ", round(maximum(norm.(MR - RM, Inf)), sigdigits = 4),
        "   -> composition is NOT commutative")

# Check against applying the two maps one after the other to a point.
step  = [evaluate(r, z0) for r in R]
two   = [evaluate(m, step) for m in M]
once  = [evaluate(m, z0) for m in MR]
println("\npoint z0 = ", z0)
println("  M(R(z0)) step by step : ", two)
println("  (M o R)(z0)           : ", once)
println("  difference            : ", maximum(abs.(two - once)),
        "   (~|z0|^4: terms of degree > $N were dropped in M o R)")

# ---------------------------------------------------------------------------
section("6. Inverse map to order $N",
        "Minv with M o Minv = identity up to order $N: from inverse_map and from",
        "the fixed-point iteration Minv = id - Nl(Minv).  Both must agree.")
# ---------------------------------------------------------------------------
Minv = inverse_map(M)                   # built into TaylorSeries
println("inverse_map(M):")
print_map(Minv, ["x", "y"])

# The same by fixed-point iteration. Write M(z) = z + Nl(z) with Nl the
# nonlinear part (here the linear part is the identity). Then
#     Minv(w) = w - Nl(Minv(w)),
# and each iteration makes one more order correct.
Nl  = M - [x, y]
Mi  = [x, y]
for k in 1:N
    global Mi = [x, y] - compose(Nl, Mi)
end
println("fixed-point iteration agrees: max diff = ", maximum(norm.(Mi - Minv, Inf)))

println("\ncheck  M o Minv  = identity:  ", str.(compose(M, Minv)))
println("check  Minv o M  = identity:  ", str.(compose(Minv, M)))
println("(both exact up to order $N; the identity is all that survives truncation)")

# ---------------------------------------------------------------------------
section("7. Jacobian (linear part) and its determinant",
        "J_ij = dM_i/dz_j (polynomials).  det J = local area change: 1 everywhere",
        "for a symplectic (Hamiltonian) map in 1 degree of freedom; here it is not.")
# ---------------------------------------------------------------------------
# J_ij = d M_i / d z_j  - each entry is itself a polynomial.
J = [differentiate(M[i], j) for i in 1:2, j in 1:2]
println("J(x,y) = [ ", str(J[1, 1]), "    ", str(J[1, 2]), " ;")
println("           ", str(J[2, 1]), "    ", str(J[2, 2]), " ]")
detJ = J[1, 1] * J[2, 2] - J[1, 2] * J[2, 1]
println("det J(x,y) = ", str(detJ), "    (a polynomial: the map changes areas)")
J0 = constant_term.(J)
println("linear part at the origin (degree-1 rows of the table) = ", J0, ",  det = ", J0[1, 1] * J0[2, 2] - J0[1, 2] * J0[2, 1])

JR = [differentiate(R[i], j) for i in 1:2, j in 1:2]
detR = JR[1, 1] * JR[2, 2] - JR[1, 2] * JR[2, 1]
println("for R: det J = ", str(detR), "    (rotation part has det 1, the x^2 term adds -2s*x)")

print_summary("A Taylor map of order N = table of coefficients (final|monomial), degree <= N.",
              "TaylorN arithmetic truncates at N; results up to order N are exact (DA).",
              "Maps are composed by substitution (M o R != R o M) and can be inverted to order N.",
              "The degree-1 part is the transfer matrix; det J measures area (phase-space) change.")
