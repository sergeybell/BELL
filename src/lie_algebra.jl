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



const JET_ORDER = 4   # phase 2: H to order 4 -> maps exact to order 3
const MAP_ORDER = 3
const SPACE = JetSpace(order = JET_ORDER, variables = [:x, :px, :y, :py, :l, :d])
"The six canonical variables as TaylorN polynomials: Z[1] = x, ..., Z[6] = delta."
const Z = variables(SPACE)

"Keep only the homogeneous parts of order 0..n of a TaylorN polynomial."
truncate_to(p::TaylorN, n::Int) = sum(TaylorN(p[k], JET_ORDER) for k in 0:n)

"Homogeneous part of total degree k."
hom(p::TaylorN, k::Int) = TaylorN(p[k], JET_ORDER)

"Largest absolute coefficient of a polynomial (used for convergence / error checks)."
maxabs(p::TaylorN) = maximum(k -> maximum(abs, p[k].coeffs; init = 0.0), 0:JET_ORDER)
maxabs(x::Real) = abs(x)

"Coefficient of the monomial z^e, e = exponent vector of length 6."
coeff(p::TaylorN, e::Vector{Int}) = getcoeff(p, e)

"Linear part (6x6 matrix R[i,j] = d M_i / d z_j at z = 0) of a polynomial map."
map_linear_part(M::AbstractVector{<:TaylorN}) =
    [constant_term(differentiate(M[i], j)) for i in 1:6, j in 1:6]

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

"ASCII string of a polynomial; names = variable names, terms below `tol` are skipped."
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

"Print all components of a map, one line each."
function map_string(M; names = ["x", "px", "y", "py", "l", "delta"], kw...)
    join(["  $(rpad(names[i], 6)) = " * poly_string(M[i]; names = names, kw...) for i in 1:6], "\n")
end



const QIDX = (1, 3, 5)
const PIDX = (2, 4, 6)

const SYMPLECTIC_S = let S = zeros(6, 6)
    for k in 1:3
        S[2k-1, 2k] = 1.0; S[2k, 2k-1] = -1.0
    end
    S
end

"Poisson bracket [f, g] = sum_k df/dq_k dg/dp_k - df/dp_k dg/dq_k."
function poisson(f::TaylorN, g::TaylorN)
    s = zero(f)
    for (q, p) in zip(QIDX, PIDX)
        s += differentiate(f, q) * differentiate(g, p) - differentiate(f, p) * differentiate(g, q)
    end
    return s
end

"exp(:f:) g summed as a Lie series; stops when a term is below tol (relative)."
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

"Order-2 transfer map z -> exp(-L :H:) z of an element with s-independent H."
lie_map(H::TaylorN, L::Real) = [truncate_to(lie_exp(-L * H, Z[i]), MAP_ORDER) for i in 1:6]

"Composition A o B (apply B first, then A), truncated to MAP_ORDER."
compose_maps(A, B) = [truncate_to(evaluate(a, collect(B)), MAP_ORDER) for a in A]

"""
Dragt-Finn factor  M = R o N3 o P4  (functions: first P4, then N3, then R),
N3 = exp(:f3:) z,  P4 = exp(:f4:) z.  `field3`/`field4` are the vector fields
dz/dt = [f, z] = -S grad f  whose time-1 flows are N3 and P4 (Tracking.jl).
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

identity_factor(name = "identity") = LieFactor(name, Matrix{Float64}(I, 6, 6), zero(Z[1]))

_euler(g, k) = (Z[1] * g[2] - Z[2] * g[1] + Z[3] * g[4] - Z[4] * g[3] + Z[5] * g[6] - Z[6] * g[5]) / k

"Extract (R, f3, f4) from a polynomial map M (f4 only if MAP_ORDER >= 3)."
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

"Linear image of the identity: (R z)_i as TaylorN polynomials."
linear_on_Z(R) = [sum(R[i, j] * Z[j] for j in 1:6) for i in 1:6]

"Generator pulled through a linear map: (f ∘ R)(z) = f(R z)."
pull_generator(f::TaylorN, R) = evaluate(f, linear_on_Z(R))

"6×6 block-diagonal rotation with phase advances mu (planes x, y, s)."
function rotation_matrix_6(mu)
    R = zeros(6, 6)
    for k in 1:3
        c, s = cos(mu[k]), sin(mu[k])
        R[2k-1:2k, 2k-1:2k] = [c s; -s c]
    end
    return R
end

"""
Factor of the product 'Fa first, then Fb' to degree 4.

The maps are R ∘ exp(:f3:) ∘ exp(:f4:), i.e. the coordinate map
z ↦ R (exp(:f3:) (exp(:f4:) z)).  After pulling the second generator
through Ra, the cubics combine by BCH of those coordinate maps:

    f3 = f3a + f3b ∘ Ra
    f4 = f4a + f4b ∘ Ra − (1/2)[f3b ∘ Ra, f3a]

The minus is required by the order of application (f3a first, then
f3b ∘ Ra); it is the same as +(1/2)[f3a, f3b ∘ Ra].  Nested brackets
are degree 5 and are dropped.
"""
function concatenate(Fa::LieFactor, Fb::LieFactor; name = Fa.name * "+" * Fb.name)
    g3 = pull_generator(Fb.f3, Fa.R)
    g4 = pull_generator(Fb.f4, Fa.R)
    f3 = hom(Fa.f3 + g3, 3)
    f4 = hom(Fa.f4 + g4 - poisson(g3, Fa.f3) / 2, 4)
    return LieFactor(name, Fb.R * Fa.R, f3, f4)
end

"Fold concatenate: first factor in the list is applied first."
function product_factors(factors; name = "product")
    acc = identity_factor(name)
    for F in factors
        F isa LieFactor || error("product_factors expects LieFactor, got $(typeof(F))")
        acc = concatenate(acc, F; name)
    end
    return acc
end

"Polynomial map of a factor: M = R o exp(:f3:) o exp(:f4:), truncated to MAP_ORDER."
function factor_to_map(F::LieFactor)
    N3 = [truncate_to(lie_exp(F.f3, Z[j]), MAP_ORDER) for j in 1:6]
    P4 = [truncate_to(lie_exp(F.f4, Z[j]), MAP_ORDER) for j in 1:6]
    N = compose_maps(N3, P4)
    return [truncate_to(sum(F.R[i, j] * N[j] for j in 1:6), MAP_ORDER) for i in 1:6]
end

end # module
