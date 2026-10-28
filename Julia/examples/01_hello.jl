# =============================================================================
# 01_hello.jl - first steps in Julia
# =============================================================================
#
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# WHAT THIS EXAMPLE TEACHES
#   The very basics of the language, enough to read the later physics examples:
#   printing, variables and their types, arithmetic (including the difference
#   between "/" and integer division), powers, math functions and constants,
#   defining functions, if/elseif/else, for loops, arrays and string
#   interpolation ("$x" inside a string).
#
# BACKGROUND (no physics here - just the language)
#   * Julia picks the type of a value automatically: 7 is an Int64, 2.5 is a
#     Float64, "Julia" is a String.  typeof(x) tells you the type.
#   * "/" always gives a floating-point result (7/2 = 3.5); div(7,2) = 3 is the
#     integer quotient and 7 % 2 = 1 the remainder.
#   * Angles are in radians: sin(pi/6) = 1/2.  Floating-point rounding means you
#     may see 0.49999999999999994 instead of 0.5 - that is normal.
#   * Arrays are indexed from 1 (not 0 as in C/Python); a[end] is the last one.
#   * A small check formula used below:  sum_{i=1}^{n} i^2 = n(n+1)(2n+1)/6.
#
# SECTIONS
#   1. Hello World            - println
#   2. Variables and types    - typeof
#   3. Arithmetic             - + - * /, div, %, ^, sqrt
#   4. Constants & functions  - pi, sin, cos, exp, round
#   5. Functions              - long form (function ... end) and one-line form
#   6. if / elseif / else
#   7. for loops              - sum of squares, factorial
#   8. Arrays                 - indexing, length, sum, sort, push!
#   9. String interpolation   - "$a + $b = $(a + b)"
#
# HOW TO RUN (from the Julia folder)
#   julia --project=. examples/01_hello.jl
#   (this script needs no packages, so plain "julia examples/01_hello.jl" works too)
#
# EXPECTED OUTPUT
#   A banner, then one block per section: a heading, a short explanation, and
#   the printed values (e.g. "7 / 2 = 3.5", "div(7, 2) = 3", "sum of squares
#   1..10 = 385").  Each line shows the Julia expression on the left and its
#   value on the right, so you can compare with what you expected.
#   At the end a short summary of what you learned.

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

banner("Example 01: first steps in Julia",
       "Shows the basic building blocks of the language: printing, variables and",
       "their types, arithmetic, math functions, user-defined functions, if/else,",
       "for loops, arrays and string interpolation.  Every printed line has the",
       "form  'expression = value'  so you can check the result yourself.")

# ---------------------------------------------------------------
# 1. Hello World: println prints its arguments followed by a newline.
# ---------------------------------------------------------------
section("1. Hello World",
        "println(...) prints all its arguments and then a newline.")
println("Hello, World!")

# ---------------------------------------------------------------
# 2. Variables: no declarations needed, Julia infers the type from the value.
# ---------------------------------------------------------------
section("2. Variables and their types",
        "No declarations needed: the type is chosen from the value.",
        "Int64 = integer, Float64 = double-precision real, String = text.")
a = 7              # an integer (Int64)
b = 2.5            # a floating-point number (Float64)
name = "Julia"     # a string (String)
println("a = ", a, "  has type ", typeof(a))
println("b = ", b, "  has type ", typeof(b))
println("name = ", name, "  has type ", typeof(name))

# ---------------------------------------------------------------
# 3. Arithmetic.  Mixing Int and Float64 gives Float64.  Note the three kinds
#    of division: "/" (always real), div (integer quotient), % (remainder).
# ---------------------------------------------------------------
section("3. Basic arithmetic",
        "Int + Float64 gives a Float64.  '/' always returns a real number,",
        "div(a, b) is the integer quotient and a % b the remainder.")
println("a + b = ", a + b)
println("a - b = ", a - b)
println("a * b = ", a * b)

# Division: "/" always gives a Float64, "div" gives an integer
println("7 / 2     = ", 7 / 2)        # 3.5
println("div(7, 2) = ", div(7, 2))    # 3   (integer division)
println("7 % 2     = ", 7 % 2)        # 1   (remainder)

# Powers and square root
println("2^10      = ", 2^10)
println("b^2       = ", b^2)
println("sqrt(16)  = ", sqrt(16))     # result is a Float64: 4.0

# ---------------------------------------------------------------
# 4. Built-in constants and math functions (angles in radians).
#    Results like 0.49999999999999994 are ordinary floating-point rounding.
# ---------------------------------------------------------------
section("4. Constants and math functions",
        "Angles are in radians.  Values such as 0.49999999999999994 instead of",
        "0.5 are normal floating-point rounding (about 1e-16), not a bug.")
println("pi        = ", float(pi))   # "pi" is a special constant; float() shows its value
println("sin(pi/6) = ", sin(pi / 6))  # 0.5 (up to rounding)
println("cos(pi/3) = ", cos(pi / 3))  # 0.5 (up to rounding)
println("cos(0)    = ", cos(0))
println("exp(1)    = ", exp(1))
println("round(pi, digits=3) = ", round(pi, digits = 3))

# ---------------------------------------------------------------
# 5. Functions: the long form (function ... end) for several lines, the short
#    form f(x) = ... for one-liners.  Both define the same kind of object.
# ---------------------------------------------------------------
section("5. Functions",
        "Long form:  function circle_area(r) ... end     (area = pi r^2)",
        "Short form: square(x) = x^2")
# Long form
function circle_area(r)
    return pi * r^2
end

# Short (one-line) form
square(x) = x^2

println("circle_area(2.0) = ", circle_area(2.0))
println("square(9)        = ", square(9))

# ---------------------------------------------------------------
# 6. if / elseif / else: choose a branch depending on a condition.
# ---------------------------------------------------------------
section("6. if / elseif / else",
        "The sign of x = -3 is tested; exactly one branch is executed.")
x = -3
if x > 0
    println("x = $x is positive")
elseif x < 0
    println("x = $x is negative")
else
    println("x = $x is zero")
end

# ---------------------------------------------------------------
# 7. for loops: sum of squares (checked against the closed formula) and a
#    factorial written by hand (checked against the built-in factorial).
# ---------------------------------------------------------------
section("7. for loops",
        "Sum of squares 1^2 + ... + 10^2, compared with n(n+1)(2n+1)/6, and a",
        "hand-written factorial compared with the built-in factorial(5).")
# Sum of squares 1^2 + 2^2 + ... + 10^2
total = 0
for i in 1:10              # 1:10 is the range 1, 2, ..., 10 (both ends included)
    global total += i^2    # "global" is needed only at the top level of a script
end
println("sum of squares 1..10 = ", total, "   (formula n(n+1)(2n+1)/6 = ", div(10 * 11 * 21, 6), ")")

# Factorial written as a function with a loop (inside a function no "global" needed)
function my_factorial(n)
    result = 1
    for i in 2:n
        result *= i        # same as result = result * i
    end
    return result
end
println("my_factorial(5) = ", my_factorial(5), "   (built-in factorial(5) = ", factorial(5), ")")

# ---------------------------------------------------------------
# 8. A small array (indices start at 1!).  Functions ending in "!" modify
#    their argument (push! appends in place).
# ---------------------------------------------------------------
section("8. Arrays",
        "Indices start at 1, numbers[end] is the last element.  push! appends",
        "in place (a '!' in a function name means: it modifies its argument).")
numbers = [3, 1, 4, 1, 5, 9]
println("numbers        = ", numbers)
println("first element  = ", numbers[1])
println("last element   = ", numbers[end])
println("length         = ", length(numbers))
println("sum            = ", sum(numbers))
println("maximum        = ", maximum(numbers))
println("sorted         = ", sort(numbers))

push!(numbers, 2)          # append an element ("!" means: the function modifies its argument)
println("after push!    = ", numbers)

for (i, value) in enumerate(numbers)    # loop over index and value together
    println("  numbers[$i] = $value")
end

# ---------------------------------------------------------------
# 9. String interpolation: $variable or $(expression) inside a string
#    (use $(name) when the next character could be part of a name, e.g. "!")
# ---------------------------------------------------------------
section("9. String interpolation",
        "Inside a string, \$x inserts the value of x and \$(expr) the value of",
        "an expression, e.g. \"\$a + \$b = \$(a + b)\".")
r = 1.5
println("Hello from $(name)! A circle of radius $r has area $(round(circle_area(r), digits = 4)).")
println("$a + $b = $(a + b)")

print_summary("Julia infers types (Int64, Float64, String); '/' is real division, div/% integer.",
              "Functions: 'function f(x) ... end' or 'f(x) = ...'; loops 'for i in 1:n ... end'.",
              "Arrays start at index 1; '!' functions (push!) modify their argument.",
              "\"\$x\" and \"\$(expr)\" insert values into strings - handy for printing results.")
