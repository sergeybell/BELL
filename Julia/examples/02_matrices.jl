# =============================================================================
# 02_matrices.jl - vectors and matrices in Julia
# =============================================================================
#
# WHAT THIS EXAMPLE TEACHES
#   How to create and use vectors and matrices with the standard library
#   LinearAlgebra - the tools needed later for transfer matrices in beam optics.
#
# BACKGROUND
#   * [1.0, 2.0, 3.0] (commas) is a column vector; [1.0 2.0; 3.0 4.0] (spaces
#     separate columns, ";" separates rows) is a 2x2 matrix.  Indices are
#     (row, column) and start at 1.
#   * A * B is the matrix product, (A B)_ij = sum_k A_ik B_kj; in general
#     A B != B A.  A .* B (with a dot) multiplies element by element.
#   * det(A) = a11 a22 - a12 a21 for 2x2, tr(A) = sum of the diagonal,
#     inv(A) is the inverse: A inv(A) = I.  I is the identity matrix.
#   * A linear system A x = b is solved with x = A \ b (more accurate and
#     faster than inv(A) * b).
#   * Eigenvalues/eigenvectors:  S v = lambda v.  A symmetric matrix has real
#     eigenvalues.  A 2D rotation R(theta) = [cos -sin; sin cos] has det = 1,
#     R' R = I and the complex eigenvalues exp(+-i theta).
#
# SECTIONS
#   1. Vectors                 - +, scalar *, dot product, norm
#   2. Matrices                - literals, zeros, ones, comprehension
#   3. Indexing                - A[i, j], whole rows/columns with ":"
#   4. Products                - A*B vs B*A vs A.*B, A^2 vs A.^2
#   5. det, tr, inv, transpose
#   6. Identity I              - A + I, isapprox (A*inv(A) ~ I)
#   7. Solving A x = b         - backslash operator
#   8. Eigenvalues             - eigvals, eigen
#   9. 2D rotation matrix      - an orthogonal matrix with det = 1
#
# HOW TO RUN (from the Julia folder)
#   julia --project=. examples/02_matrices.jl
#   (only standard libraries are used, "julia examples/02_matrices.jl" also works)
#
# EXPECTED OUTPUT
#   For each section a heading, a short explanation, then the results.
#   Matrices are printed row by row.  Tiny numbers like 8.88e-16 are
#   floating-point rounding (should be 0).  Checks print "true" when the
#   identity holds, e.g. "A * inv(A) approx. equal to I ? true".

using LinearAlgebra   # standard library: det, tr, inv, eigvals, I, norm, ...

# -----------------------------------------------------------------------------
# Printing helpers (ASCII only - the Windows console garbles Unicode).
#   banner(title, lines...)   big title + "What this example does"
#   section(title, lines...)  section heading + 1-3 lines of explanation
#   print_summary(lines...)   "Summary / what you learned" at the end
#   showmat(M)                print a matrix row by row (without Julia's
#                             "2x2 Matrix{Float64}" header, which uses a
#                             Unicode multiplication sign)
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

showmat(M) = (Base.print_matrix(IOContext(stdout, :compact => true, :limit => true), M); println())

banner("Example 02: vectors and matrices",
       "Creates vectors and matrices and applies the basic linear-algebra",
       "operations: products (matrix vs element-wise), transpose, determinant,",
       "trace, inverse, solving A x = b, eigenvalues and a 2D rotation matrix.",
       "Each result is followed by a check you can verify by hand.")

# ---------------------------------------------------------------
# 1. Vectors: arithmetic is component-wise, dot and norm come from LinearAlgebra.
# ---------------------------------------------------------------
section("1. Vectors",
        "v = [1,2,3], w = [4,5,6].  dot(v, w) = 1*4 + 2*5 + 3*6 = 32,",
        "norm(v) = sqrt(1 + 4 + 9) = 3.7417.")
v = [1.0, 2.0, 3.0]          # commas -> column vector
w = [4.0, 5.0, 6.0]
println("v = ", v)
println("w = ", w)
println("v[1] = ", v[1], ",  v[end] = ", v[end], ",  length(v) = ", length(v))
println("v + w      = ", v + w)
println("2 * v      = ", 2 * v)
println("dot(v, w)  = ", dot(v, w))      # scalar product 1*4 + 2*5 + 3*6 = 32
println("norm(v)    = ", norm(v))        # length sqrt(1 + 4 + 9)

# ---------------------------------------------------------------
# 2. Matrices: spaces separate columns, semicolons (or new lines) separate rows
# ---------------------------------------------------------------
section("2. Matrices",
        "Literal [1 2; 3 4]: spaces separate columns, ';' separates rows.",
        "zeros(m,n), ones(m,n) and comprehensions [f(i,j) for i in .., j in ..].")
A = [1.0 2.0;
     3.0 4.0]
B = [0.0 1.0;
     1.0 0.0]
println("A ="); showmat(A)
println("B ="); showmat(B)
println("size(A) = ", size(A))           # (rows, columns)

# Other ways to create matrices
println("zeros(2, 3) ="); showmat(zeros(2, 3))
println("ones(2, 2) ="); showmat(ones(2, 2))
println("3x3 matrix from a formula M[i, j] = 10i + j:")
showmat([10i + j for i in 1:3, j in 1:3])

# ---------------------------------------------------------------
# 3. Indexing (row, column), 1-based
# ---------------------------------------------------------------
section("3. Indexing",
        "A[i, j] = row i, column j (starting at 1); ':' selects a whole row/column.")
println("A[1, 2]  = ", A[1, 2])       # row 1, column 2
println("A[2, :]  = ", A[2, :])         # whole second row
println("A[:, 1]  = ", A[:, 1])         # whole first column

# ---------------------------------------------------------------
# 4. Matrix product vs element-wise product.  The dot (.*, .^) means
#    "apply element by element"; without the dot it is linear algebra.
# ---------------------------------------------------------------
section("4. Matrix product vs element-wise product",
        "A * B is the matrix product (order matters: A*B != B*A);",
        "A .* B and A .^ 2 act element by element; A^2 = A * A.")
println("A * B  (matrix product, rows times columns) ="); showmat(A * B)
println("B * A  (order matters: A*B != B*A in general) ="); showmat(B * A)
println("A .* B (element-wise product, note the dot) ="); showmat(A .* B)
println("A * v[1:2] (matrix times vector) = ", A * v[1:2])
println("A .^ 2 (each element squared) ="); showmat(A .^ 2)
println("A^2 = A * A ="); showmat(A^2)

# ---------------------------------------------------------------
# 5. Transpose, determinant, trace, inverse
# ---------------------------------------------------------------
section("5. Transpose, determinant, trace, inverse",
        "det(A) = 1*4 - 2*3 = -2, tr(A) = 1 + 4 = 5.  A * inv(A) should be the",
        "identity; an entry like 8.9e-16 is rounding error (really 0).")
println("A' (transpose) ="); showmat(A')
println("det(A) = ", det(A))            # 1*4 - 2*3 = -2
println("tr(A)  = ", tr(A))             # trace = 1 + 4 = 5
println("inv(A) ="); showmat(inv(A))
println("A * inv(A) ="); showmat(A * inv(A))   # identity (up to rounding)

# ---------------------------------------------------------------
# 6. Identity matrix
# ---------------------------------------------------------------
section("6. Identity matrix I",
        "I adapts its size automatically.  'approximately equal' (isapprox,",
        "written with the approx operator in the code) ignores rounding errors.")
# "I" is a special identity object that adapts its size automatically
println("A + I (adds 1 to the diagonal) ="); showmat(A + I)
# ≈ (type \approx<TAB>) means "approximately equal", same as isapprox(x, y)
println("A * inv(A) approx. equal to I ? ", A * inv(A) ≈ I)
I3 = Matrix(1.0I, 3, 3)                         # explicit 3x3 identity matrix
println("3x3 identity ="); showmat(I3)

# ---------------------------------------------------------------
# 7. Solving a linear system A x = b
# ---------------------------------------------------------------
#   1 x1 + 2 x2 = 5
#   3 x1 + 4 x2 = 6
section("7. Solving a linear system A x = b",
        "System  x1 + 2 x2 = 5,  3 x1 + 4 x2 = 6.  Exact solution x = (-4, 4.5);",
        "x = A \\ b solves it without forming inv(A).  Check: A * x == b.")
b = [5.0, 6.0]
x = A \ b                      # backslash solves A x = b (better than inv(A) * b)
println("solution of A x = b:  x = ", x)
println("check A * x = ", A * x, "  (should equal b = ", b, ")")

# ---------------------------------------------------------------
# 8. Eigenvalues and eigenvectors:  S v = lambda v
# ---------------------------------------------------------------
section("8. Eigenvalues and eigenvectors",
        "S = [2 1; 1 2] is symmetric -> real eigenvalues 1 and 3 with",
        "eigenvectors (-1, 1)/sqrt(2) and (1, 1)/sqrt(2) (the columns below).")
S = [2.0 1.0;
     1.0 2.0]                  # symmetric matrix -> real eigenvalues
println("S ="); showmat(S)
println("eigvals(S) = ", eigvals(S))         # 1 and 3
E = eigen(S)                                 # values and vectors together
println("eigenvectors (columns) ="); showmat(E.vectors)
println("check S * v1 = lambda1 * v1: ", S * E.vectors[:, 1] ≈ E.values[1] * E.vectors[:, 1])

# ---------------------------------------------------------------
# 9. Example: 2D rotation matrix
# ---------------------------------------------------------------
# R(theta) rotates a vector by the angle theta (counter-clockwise).
section("9. 2D rotation matrix",
        "R(theta) = [cos -sin; sin cos] with theta = 90 deg turns (1,0) into (0,1).",
        "Rotations keep lengths and areas: det R = 1, R' R = I, eigenvalues exp(+-i theta).")
theta = pi / 2                               # 90 degrees
R = [cos(theta) -sin(theta);
     sin(theta)  cos(theta)]
p = [1.0, 0.0]
println("rotation by 90 deg of p = ", p, "  ->  R * p = ", round.(R * p, digits = 12))
println("det(R) = ", det(R), "   (rotations keep areas: det = 1)")
println("R' * R approx. equal to I ? ", R' * R ≈ I, "   (rotation matrices are orthogonal: inv(R) = R')")
println("eigvals(R) = ", eigvals(R), "   (complex: exp(+-i*theta))")

print_summary("Vectors [a, b, c], matrices [a b; c d]; indices (row, col) start at 1.",
              "A * B is the matrix product (not commutative), A .* B is element-wise.",
              "det, tr, inv, eigvals/eigen come from LinearAlgebra; solve A x = b with A \\ b.",
              "Compare floating-point results with isapprox (approx), not with ==.")
