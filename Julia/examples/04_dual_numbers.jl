# =============================================================================
# 04_dual_numbers.jl - automatic differentiation with hand-made dual numbers
# =============================================================================
#
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# THE PROBLEM
#   Given a function that is computed by a PROGRAM (a chain of +, *, /, sin,
#   exp, ... - possibly thousands of lines, e.g. a particle tracking code),
#   get its value AND its exact derivative at a point x0, without deriving any
#   formula by hand and without finite differences (f(x0+h)-f(x0))/h, which
#   are only approximate.
#
# THE IDEA: DUAL NUMBERS
#   A dual number is  a + b*eps  with the rule  eps^2 = 0  (eps is a formal
#   symbol, not a small real number; in the code it is called ε, typed
#   \varepsilon<TAB>).  Taylor's formula then terminates after the linear term:
#
#       f(x0 + eps) = f(x0) + f'(x0) * eps          (exactly, since eps^2 = 0)
#
#   So: feed  x0 + eps  into the unchanged program and read the eps
#   coefficient of the result - that is f'(x0), exact to machine precision.
#   To make this work we only teach Julia the arithmetic of duals:
#       (a + b eps) + (c + d eps) = (a + c) + (b + d) eps
#       (a + b eps) * (c + d eps) = a c + (a d + b c) eps        (product rule)
#       (a + b eps) / (c + d eps) = a/c + (b c - a d)/c^2 eps    (quotient rule)
#       g(a + b eps)              = g(a) + g'(a) b eps           (chain rule)
#   for g = sin, cos, exp, sqrt.  Julia's multiple dispatch then runs any
#   existing function (written for ordinary numbers) on duals automatically.
#
# THE BEAM-OPTICS APPLICATION
#   A linear beamline maps (x, x') at the entrance to (x, x') at the exit.
#   Push the "particle" (x0 + eps, x0') through the element functions: the eps
#   coefficients of the final x and x' are the matrix elements (x|x) and
#   (x'|x) - the first column of the transfer matrix.  Seeding eps in x0'
#   instead gives the second column (x|x'), (x'|x').  We check the result
#   against the usual matrix product M = D(L2) * Lens(f) * D(L1).
#   This is first-order, one-variable differential algebra (DA); COSY INFINITY
#   does the same with all Taylor coefficients up to order n in many
#   variables (see examples 05 and 06).
#
# SECTIONS
#   0. Definition of the Dual type and its arithmetic (product, quotient,
#      chain rule) - this is the whole "AD library", about 20 lines
#   1. Basic checks: x^2 and x^3 + 2x at x0 = 3, sin(x) e^x / x at x0 = 1
#   2. Tracking a dual particle through drift - thin lens - drift
#   3. Comparison with the transfer-matrix product
#
# HOW TO RUN (from the Julia folder)
#   julia --project=. examples/04_dual_numbers.jl
#   (no packages needed, "julia examples/04_dual_numbers.jl" also works)
#
# EXPECTED OUTPUT / HOW TO READ IT
#   Duals are printed as  "value + derivative eps", e.g.  x^2 at x0 = 3 gives
#   "9.0 + 6.0eps": value 9, derivative 2*3 = 6.  x^3 + 2x gives
#   "33.0 + 29.0eps" (27 + 6 and 3*9 + 2).  For the beamline, the eps parts
#   form the matrix [-0.5 2.5; -0.5 0.5], identical to the matrix product;
#   (x|x) = 1 - L2/f = -0.5.

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

# =============================================================================
# 0. The Dual type: value a and eps-coefficient b (the derivative)
# =============================================================================
struct Dual <: Number
    a::Float64   # value
    b::Float64   # derivative (coefficient of ε)
end
Dual(a::Real) = Dual(a, 0.0)            # a constant has zero derivative
const ε = Dual(0.0, 1.0)                # the formal symbol eps, eps^2 = 0

value(d::Dual) = d.a
deriv(d::Dual) = d.b
# printed as "a + b eps" (ASCII "eps", so the Windows console shows it correctly)
Base.show(io::IO, d::Dual) = print(io, d.a, (d.b < 0 ? " - " : " + "), abs(d.b), "eps")

# Arithmetic: sum rule, product rule, quotient rule
Base.:+(x::Dual, y::Dual) = Dual(x.a + y.a, x.b + y.b)
Base.:-(x::Dual, y::Dual) = Dual(x.a - y.a, x.b - y.b)
Base.:-(x::Dual)          = Dual(-x.a, -x.b)
Base.:*(x::Dual, y::Dual) = Dual(x.a * y.a, x.a * y.b + x.b * y.a)
Base.:/(x::Dual, y::Dual) = Dual(x.a / y.a, (x.b * y.a - x.a * y.b) / y.a^2)

# Mixing with ordinary numbers: Julia's promotion system converts Real -> Dual,
# so 2 * x, x + 1, 1 / x all work automatically with the methods above.
Base.convert(::Type{Dual}, x::Real) = Dual(x)
Base.promote_rule(::Type{Dual}, ::Type{<:Real}) = Dual

# Elementary functions: chain rule  g(a + b eps) = g(a) + g'(a) b eps
Base.sin(x::Dual)  = Dual(sin(x.a),  cos(x.a) * x.b)
Base.cos(x::Dual)  = Dual(cos(x.a), -sin(x.a) * x.b)
Base.exp(x::Dual)  = Dual(exp(x.a),  exp(x.a) * x.b)
Base.sqrt(x::Dual) = Dual(sqrt(x.a), x.b / (2 * sqrt(x.a)))

banner("Example 04: dual numbers = exact derivatives of a program",
       "Problem: a function is given as a program; we want its value AND its",
       "exact derivative at x0 without deriving formulas or using finite",
       "differences.  Trick: feed x0 + eps (with eps^2 = 0) into the program;",
       "the result is f(x0) + f'(x0) eps - just read off the eps coefficient.",
       "Then the same trick gives the transfer matrix of a simple beamline.")

# =============================================================================
# 1. Basic checks: compare the eps coefficient with the derivative by hand
# =============================================================================
section("1. Value and derivative of simple functions",
        "Printed as 'value + derivative eps'.  We evaluate at x0 + eps and compare",
        "the eps coefficient with the derivative computed by hand.")
x = 3 + ε
println("x            = ", x, "        (x0 = 3 seeded with eps: dx/dx = 1)")
println("x^2          = ", x^2, "        (expect 9 + 6eps:  3^2 = 9,  2*3 = 6)")
f(x) = x^3 + 2x
println("x^3 + 2x     = ", f(x), "       (expect 33 + 29eps: 27+6, 3*9+2)")
g(x) = sin(x) * exp(x) / x
y = g(1.0 + ε)
dg_exact = (cos(1) * exp(1) + sin(1) * exp(1)) - sin(1) * exp(1)   # d/dx[sin e^x / x] at 1
println("sin(x)exp(x)/x at x0 = 1: ", y)
println("  analytic derivative = ", dg_exact, ",  difference = ", deriv(y) - dg_exact,
        "   (exact to machine precision)")

# =============================================================================
# 2. Particle tracking with duals: drift L1 -> thin lens (focal f) -> drift L2
#    The element maps are ordinary Julia functions of (x, x') written for real
#    numbers; they work on duals without any change.
# =============================================================================
drift(x, xp, L) = (x + L * xp, xp)
thin_lens(x, xp, f) = (x, xp - x / f)

function line(x, xp; L1 = 1.0, fl = 2.0, L2 = 3.0)
    x, xp = drift(x, xp, L1)
    x, xp = thin_lens(x, xp, fl)
    x, xp = drift(x, xp, L2)
    return x, xp
end

section("2. Tracking a dual particle: drift(1) + thin lens(f=2) + drift(3)",
        "Start on axis (x0 = x0' = 0) and seed eps in x0: the eps parts of the final",
        "(x, x') are (x|x) and (x'|x).  Seeding eps in x0' gives (x|x') and (x'|x').")
x0, xp0 = 0.0, 0.0                       # reference particle on axis
# Seed ε in x: the ε-coefficients of the result are (x|x) and (x'|x)
xf, xpf = line(x0 + ε, Dual(xp0))
println("x0 = 0 + 1eps  ->  x = ", xf, ",  x' = ", xpf)
# Seed ε in x': the ε-coefficients are (x|x') and (x'|x')
xf2, xpf2 = line(Dual(x0), xp0 + ε)
println("x0'= 0 + 1eps  ->  x = ", xf2, ",  x' = ", xpf2)

# =============================================================================
# 3. The same transfer matrix from the usual matrix product
#    M = D(L2) * Lens(f) * D(L1)   (the first element acts first, i.e. rightmost)
# =============================================================================
section("3. Check against the transfer-matrix product",
        "Collect the eps coefficients into a 2x2 matrix and compare it with",
        "M = D(3) * Lens(2) * D(1).  They must be identical.")
M_dual = [deriv(xf) deriv(xf2); deriv(xpf) deriv(xpf2)]
D(L) = [1 L; 0 1]
Lens(f) = [1 0; -1/f 1]
M_mat = D(3.0) * Lens(2.0) * D(1.0)
println("map from duals   (x|x) (x|x'); (x'|x) (x'|x') = ", M_dual)
println("matrix product   D(3)*Lens(2)*D(1)             = ", M_mat)
println("identical? ", M_dual ≈ M_mat)
println("(x|x) = 1 - L2/f = ", 1 - 3.0 / 2.0)

print_summary("A dual number a + b eps with eps^2 = 0 carries a value and a derivative.",
              "Run an unchanged program on x0 + eps: the eps coefficient is f'(x0), exact.",
              "Pushing x0 + eps through a beamline gives the transfer matrix (x|x), (x'|x), ...",
              "This is 1st-order DA; examples 05/06 extend it to all orders and variables.")
