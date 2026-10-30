# =============================================================================
# lie_algebra.jl - truncated polynomial space, Poisson bracket, Dragt-Finn product
# =============================================================================
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# WHAT
#   Polynomial space in the 6 canonical variables (x, px, y, py, l, delta).
#   Hamiltonians are kept to total degree JET_ORDER = 4; maps obtained from
#   the Lie series are cut at MAP_ORDER = 3.
#   poisson, lie_exp, lie_map, LieFactor, dragt_finn_factor, concatenate,
#   product_factors, factor_to_map.
#
# WHY
#   The ring is a product of Lie transformations, not a substitution of
#   truncated maps.  Each factor is symplectic; the truncated polynomial
#   expansion of the product is symplectic only up to the dropped degree.
#
# =============================================================================
module LieAlgebra

using TaylorSeries, LinearAlgebra

export SPACE, Z, JET_ORDER, MAP_ORDER, truncate_to, hom, maxabs, map_linear_part,
       poly_string, map_string, coeff, exponents_of_order,
       poisson, lie_exp, lie_map, LieFactor, dragt_finn_factor, concatenate,
       factor_to_map, compose_maps, SYMPLECTIC_S, identity_factor,
       pull_generator, product_factors, rotation_matrix_6



"""
    JET_ORDER :: Int = 4

Maximum total degree of every `TaylorN` polynomial in Bell. Hamiltonians are
expanded to this degree, so maps of degree [`MAP_ORDER`](@ref) = 3 are exact
(a degree-4 H gives the cubic terms of the map). Changing it changes the
whole polynomial space ([`SPACE`](@ref)).
"""
const JET_ORDER = 4   # phase 2: H to order 4 -> maps exact to order 3
"""
    MAP_ORDER :: Int = 3

Degree at which all transfer maps are truncated (`truncate_to(·, MAP_ORDER)`).
Dragt-Finn factors therefore have generators f3 (degree 3) and f4 (degree 4).
"""
const MAP_ORDER = 3
"""
    SPACE

TaylorSeries `JetSpace` with 6 variables (x, px, y, py, l, d) and order
[`JET_ORDER`](@ref). Defined once at load time; all polynomials share it.
"""
const SPACE = JetSpace(order = JET_ORDER, variables = [:x, :px, :y, :py, :l, :d])
"""
    Z :: Vector{TaylorN{Float64}}

The six canonical variables as polynomials: `Z[1] = x`, `Z[2] = px`, `Z[3] = y`,
`Z[4] = py`, `Z[5] = l`, `Z[6] = δ`. `collect(Z)` is the identity map;
a map is a 6-vector of polynomials in `Z`.
"""
const Z = variables(SPACE)

"""
    truncate_to(p::TaylorN, n) -> TaylorN

Keep only the homogeneous parts of degree 0…n of `p` (still in the order-
[`JET_ORDER`](@ref) space). Used after every composition to stay at MAP_ORDER.
"""
truncate_to(p::TaylorN, n::Int) = sum(TaylorN(p[k], JET_ORDER) for k in 0:n)

"""
    hom(p::TaylorN, k) -> TaylorN

Homogeneous part of total degree `k` of `p`.
"""
hom(p::TaylorN, k::Int) = TaylorN(p[k], JET_ORDER)

"""
    maxabs(p) -> Float64

Largest absolute coefficient of a polynomial (or `abs(x)` for a number);
the norm used in convergence tests and residual reports.
"""
maxabs(p::TaylorN) = maximum(k -> maximum(abs, p[k].coeffs; init = 0.0), 0:JET_ORDER)
maxabs(x::Real) = abs(x)

"""
    coeff(p, e::Vector{Int}) -> Float64

Coefficient of the monomial x^e1 px^e2 y^e3 py^e4 l^e5 δ^e6 in `p`.
Example: `coeff(M[1], [0,0,0,0,0,1])` is the dispersion (x|δ).
"""
coeff(p::TaylorN, e::Vector{Int}) = getcoeff(p, e)

"""
    map_linear_part(M) -> Matrix{Float64} (6×6)

Linear part of a polynomial map: R[i, j] = ∂M_i/∂z_j at z = 0.
"""
map_linear_part(M::AbstractVector{<:TaylorN}) =
    [constant_term(differentiate(M[i], j)) for i in 1:6, j in 1:6]

"""
    exponents_of_order(k) -> Vector{Vector{Int}}

All exponent vectors (length 6, non-negative, sum `k`) in a fixed
lexicographic order. This fixes the monomial basis used by the homological
equation of the normal form (`homological_operator`).
For k = 3 there are 56 monomials, for k = 4 there are 126.
"""
function exponents_of_order(k::Int)
    out = Vector{Vector{Int}}()
    function rec(prefix, left, pos)
        if pos == 6
            push!(out, vcat(prefix, left)); return
        end
        for a in left:-1:0
            rec(vcat(prefix, a), left - a, pos + 1)
        end
    end
    rec(Int[], k, 1)
    return out
end

"""
    poly_string(p; names = VAR_NAMES, tol = …) -> String

ASCII rendering of a polynomial, terms with |coefficient| below `tol` skipped.
"""
function poly_string(p::TaylorN; names = ["x", "px", "y", "py", "l", "delta"],
                     digits = 6, tol = 1e-14, maxorder = JET_ORDER)
    terms = String[]
    for k in 0:maxorder, e in exponents_of_order(k)
        c = getcoeff(p, e)
        abs(c) < tol && continue
        mono = join([e[i] == 1 ? names[i] : "$(names[i])^$(e[i])" for i in 1:6 if e[i] > 0], "*")
        cs = string(round(c, sigdigits = digits))
        push!(terms, isempty(mono) ? cs : cs * "*" * mono)
    end
    isempty(terms) && return "0"
    return replace(join(terms, " + "), "+ -" => "- ")
end

"""
    map_string(M; names = VAR_NAMES, kw...) -> String

All six components of a map, one line each, via [`poly_string`](@ref).
"""
function map_string(M; names = ["x", "px", "y", "py", "l", "delta"], kw...)
    join(["  $(rpad(names[i], 6)) = " * poly_string(M[i]; names = names, kw...) for i in 1:6], "\n")
end



const QIDX = (1, 3, 5)
const PIDX = (2, 4, 6)

"""
    SYMPLECTIC_S :: Matrix{Float64} (6×6)

Symplectic unit matrix in the ordering (x, px, y, py, l, δ): block-diagonal
with blocks [0 1; −1 0]. Hamilton's equations: dz/ds = S ∇H.
A matrix M is symplectic iff Mᵀ S M = S (see `symplectic_error`).
"""
const SYMPLECTIC_S = let S = zeros(6, 6)
    for k in 1:3
        S[2k-1, 2k] = 1.0; S[2k, 2k-1] = -1.0
    end
    S
end

"""
    poisson(f, g) -> TaylorN

Poisson bracket [f, g] = Σ_k (∂f/∂q_k ∂g/∂p_k − ∂f/∂p_k ∂g/∂q_k),
q = (x, y, l), p = (px, py, δ). The Lie operator is :f: g = [f, g].
Degree: deg[f, g] = deg f + deg g − 2.
"""
function poisson(f::TaylorN, g::TaylorN)
    s = zero(f)
    for (q, p) in zip(QIDX, PIDX)
        s += differentiate(f, q) * differentiate(g, p) - differentiate(f, p) * differentiate(g, q)
    end
    return s
end

"""
    lie_exp(f, g; nmax = 400, tol = 1e-18) -> TaylorN

Lie transformation exp(:f:) g = g + [f, g] + [f, [f, g]]/2! + … .

# Algorithm
Terms are generated recursively, `term_n = [f, term_{n−1}]/n`, and summed until
`maxabs(term) ≤ tol·max(1, maxabs(result))`; a warning is issued after `nmax`
terms. Since the space is truncated at JET_ORDER, the series terminates
exactly when f has degree ≥ 3; for quadratic f (linear part) it converges
like the exponential series.
"""
function lie_exp(f::TaylorN, g::TaylorN; nmax::Int = 400, tol = 1e-18)
    result, term = g, g
    for n in 1:nmax
        term = poisson(f, term) / n
        result += term
        maxabs(term) <= tol * max(1.0, maxabs(result)) && return result
    end
    @warn "lie_exp did not converge in $nmax terms"
    return result
end

"""
    lie_map(H, L) -> Vector{TaylorN}

Transfer map of an element with s-independent Hamiltonian `H` over length `L` [m]:
z ↦ exp(−L :H:) z (component-wise), truncated to [`MAP_ORDER`](@ref).
Since dz/ds = [z, H] = −:H: z, this is the exact flow of H to that order.
"""
lie_map(H::TaylorN, L::Real) = [truncate_to(lie_exp(-L * H, Z[i]), MAP_ORDER) for i in 1:6]

"""
    compose_maps(A, B) -> Vector{TaylorN}

Composition A ∘ B (B is applied first, then A): component i is A_i(B(z)),
truncated to [`MAP_ORDER`](@ref). Truncation makes the result symplectic only
up to the dropped degree; prefer [`concatenate`](@ref) for long products.
"""
compose_maps(A, B) = [truncate_to(evaluate(a, collect(B)), MAP_ORDER) for a in A]

"""
    LieFactor(name, R, f3, f4)

Dragt-Finn factorised map M = R ∘ N3 ∘ P4, with N3 = exp(:f3:) z and P4 = exp(:f4:) z.
Applied to coordinates: first P4, then N3, then the matrix R.

# Fields
- `name::String` – label (element or product name).
- `R::Matrix{Float64}` – 6×6 linear (symplectic) part.
- `f3::TaylorN{Float64}` – homogeneous cubic generator.
- `f4::TaylorN{Float64}` – homogeneous quartic generator.
- `field::Vector{TaylorN}`, `field4::Vector{TaylorN}` – vector fields
  dz/dt = [f, z] = −S ∇f of f3 and f4 (computed by the constructor). Their
  time-1 flows, integrated symplectically by `Tracking`, are N3 and P4.

Each factor is exactly symplectic as a map (exponentials of Hamiltonian vector
fields); only its generators are truncated.
"""
struct LieFactor
    name::String
    R::Matrix{Float64}
    f3::TaylorN{Float64}
    f4::TaylorN{Float64}
    field::Vector{TaylorN{Float64}}
    field4::Vector{TaylorN{Float64}}
end

_vfield(f) = (g = [differentiate(f, j) for j in 1:6]; [-sum(SYMPLECTIC_S[i, j] * g[j] for j in 1:6) for i in 1:6])

LieFactor(name, R, f3, f4 = zero(Z[1])) = LieFactor(name, Matrix{Float64}(R), f3, f4, _vfield(f3), _vfield(f4))

"""
    identity_factor(name = "identity") -> LieFactor

Factor with R = I and f3 = f4 = 0 (the start of a product).
"""
identity_factor(name = "identity") = LieFactor(name, Matrix{Float64}(I, 6, 6), zero(Z[1]))

_euler(g, k) = (Z[1] * g[2] - Z[2] * g[1] + Z[3] * g[4] - Z[4] * g[3] + Z[5] * g[6] - Z[6] * g[5]) / k

"""
    dragt_finn_factor(M; name = "") -> LieFactor

Extract the Dragt-Finn factorisation of a polynomial map `M` without constant
term (the origin must be a fixed point).

# Algorithm
1. R = linear part of M; N = R⁻¹ M (so N = identity + degree ≥ 2).
2. The degree-2 part N₂ of N is the vector field of f3. f3 is recovered by
   the Euler-homogeneity formula f = (Σ_k q_k g_{p_k} − p_k g_{q_k})/deg,
   where g is the degree-2 part of N (helper `_euler`).
3. If MAP_ORDER ≥ 3, compute P = exp(−:f3:) ∘ N. Its degree-3 part gives f4 in
   the same way.

# Throws
An error if `M` has a constant (closed-orbit) part. Shift the map first
(`center_map_at`).
"""
function dragt_finn_factor(M::AbstractVector{<:TaylorN}; name = "")
    maximum(abs(constant_term(m)) for m in M) > 1e-14 && error("map has a constant (closed-orbit) part")
    R = map_linear_part(M)
    Ri = inv(R)
    N = [sum(Ri[i, j] * M[j] for j in 1:6) for i in 1:6]
    f3 = _euler([TaylorN(N[i][2], JET_ORDER) for i in 1:6], 3)
    f4 = zero(Z[1])
    if MAP_ORDER >= 3
        N3inv = [truncate_to(lie_exp(-f3, Z[i]), MAP_ORDER) for i in 1:6]
        P = compose_maps(N3inv, N)
        f4 = _euler([TaylorN(P[i][3], JET_ORDER) for i in 1:6], 4)
    end
    return LieFactor(name, R, f3, f4)
end

"""
    linear_on_Z(R) -> Vector{TaylorN}

Polynomials (R z)_i = Σ_j R[i, j] Z[j]: the linear map `R` as a polynomial map.
"""
linear_on_Z(R) = [sum(R[i, j] * Z[j] for j in 1:6) for i in 1:6]

"""
    pull_generator(f, R) -> TaylorN

Generator pulled through a linear map: (f ∘ R)(z) = f(R z). Used to move a
later generator to the front of a product, and to express generators in
normal-form coordinates (R = A).
"""
pull_generator(f::TaylorN, R) = evaluate(f, linear_on_Z(R))

"""
    rotation_matrix_6(mu) -> Matrix{Float64} (6×6)

Block-diagonal rotation with phase advances `mu = (μx, μy, μs)` [rad], blocks
[cos μ  sin μ; −sin μ  cos μ]: the linear normal form of a stable one-turn matrix.
"""
function rotation_matrix_6(mu)
    R = zeros(6, 6)
    for k in 1:3
        c, s = cos(mu[k]), sin(mu[k])
        R[2k-1:2k, 2k-1:2k] = [c s; -s c]
    end
    return R
end

"""
    concatenate(Fa, Fb; name) -> LieFactor

Factor of the product "Fa first, then Fb", correct to degree 4.

# Algorithm
The second generators are pulled through Ra (g3 = f3b ∘ Ra, g4 = f4b ∘ Ra), then
combined by Baker-Campbell-Hausdorff:

- R  = Rb · Ra
- f3 = f3a + g3
- f4 = f4a + g4 − ½ [g3, f3a]

The sign follows from the order of application (f3a first). Nested brackets are
of degree ≥ 5 and are dropped. The result is again an exactly symplectic map.
This is how the one-turn factor of the ring is built ([`product_factors`](@ref)),
instead of composing truncated polynomials.
"""
function concatenate(Fa::LieFactor, Fb::LieFactor; name = Fa.name * "+" * Fb.name)
    g3 = pull_generator(Fb.f3, Fa.R)
    g4 = pull_generator(Fb.f4, Fa.R)
    f3 = hom(Fa.f3 + g3, 3)
    f4 = hom(Fa.f4 + g4 - poisson(g3, Fa.f3) / 2, 4)
    return LieFactor(name, Fb.R * Fa.R, f3, f4)
end

"""
    product_factors(factors; name = "product") -> LieFactor

Left fold of [`concatenate`](@ref) starting from [`identity_factor`](@ref);
the first factor in the list is applied first. Throws if an entry is not a `LieFactor`.
"""
function product_factors(factors; name = "product")
    acc = identity_factor(name)
    for F in factors
        F isa LieFactor || error("product_factors expects LieFactor, got $(typeof(F))")
        acc = concatenate(acc, F; name)
    end
    return acc
end

"""
    factor_to_map(F) -> Vector{TaylorN}

Polynomial map of a factor: M = R ∘ exp(:f3:) ∘ exp(:f4:), truncated to
[`MAP_ORDER`](@ref). This is the "reading" of the one-turn factor from which
dispersion, path length and chromaticity are extracted.
"""
function factor_to_map(F::LieFactor)
    N3 = [truncate_to(lie_exp(F.f3, Z[j]), MAP_ORDER) for j in 1:6]
    P4 = [truncate_to(lie_exp(F.f4, Z[j]), MAP_ORDER) for j in 1:6]
    N = compose_maps(N3, P4)
    return [truncate_to(sum(F.R[i, j] * N[j] for j in 1:6), MAP_ORDER) for i in 1:6]
end

end # module
