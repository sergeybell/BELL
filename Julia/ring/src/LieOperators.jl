# =============================================================================
# LieOperators.jl - Poisson bracket, Lie series, Dragt-Finn factors
# =============================================================================
#
# WHAT
#   poisson(f, g)          [f, g] in the canonical variables
#   lie_exp(f, g)          exp(:f:) g = g + [f,g] + [f,[f,g]]/2! + ...
#                          (series summed until the terms vanish)
#   lie_map(H, L)          transfer map of an s-independent element of length L:
#                          z(L) = exp(-L :H:) z   (since dz/ds = [z,H] = -:H: z),
#                          truncated to MAP_ORDER = 2
#   LieFactor              factorised (Dragt-Finn) form  M = R o N3 o P4,
#                          N3 = exp(:f3:) z, P4 = exp(:f4:) z  (f4 since phase 2)
#                          i.e.  M(z) = R * N(z),  N = exp(:f3:) z
#   dragt_finn_factor(M)   extract R (6x6), cubic f3 and quartic f4 from an
#                          order-3 map (P4 = N3^-1 o R^-1 o M = z + g3 + ...,
#                          f4 = (1/4) sum_k (q_k g3_{p_k} - p_k g3_{q_k}))
#   concatenate(Fa, Fb)    factor of "Fa first, then Fb":
#                          R = Rb*Ra,  f3 = f3a + f3b o Ra   (exact to order 2;
#                          f4 of a product (BCH term) not implemented -> f4 = 0)
#   factor_to_map(F)       back to an order-2 polynomial map
#   SYMPLECTIC_S           the 6x6 symplectic form S (blocks [0 1; -1 0])
#
# WHY
#   The Lie-algebraic description keeps the map symplectic by construction:
#   R = exp(L S A) is a symplectic matrix and exp(:f3:) is a canonical
#   transformation for ANY f3.  Tracking (Tracking.jl) applies each factor
#   exactly-symplectically, factor by factor.
#
# HOW f3 IS EXTRACTED (own derivation)
#   N = R^-1 M = z + g(z) + O(3), g quadratic.  exp(:f3:) q = q - df3/dp + ...,
#   exp(:f3:) p = p + df3/dq + ...  ->  df3/dq_k = g_{p_k}, df3/dp_k = -g_{q_k};
#   Euler's theorem for a homogeneous cubic gives
#       f3 = (1/3) sum_k ( q_k g_{p_k} - p_k g_{q_k} ).
#   The reconstruction error (M vs factor_to_map) is printed by scripts/02 and
#   is zero only if M is symplectic to order 2 - an independent check.
#
# WHAT YOU SEE
#   Nothing printed here.
# =============================================================================
module LieOperators

using TaylorSeries, LinearAlgebra
using ..Jets

export poisson, lie_exp, lie_map, LieFactor, dragt_finn_factor, concatenate,
       factor_to_map, compose_maps, SYMPLECTIC_S, identity_factor

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

"Factor of the product 'Fa first, then Fb' to order 2 (cubic generators add after transport)."
function concatenate(Fa::LieFactor, Fb::LieFactor; name = Fa.name * "+" * Fb.name)
    Raz = [sum(Fa.R[i, j] * Z[j] for j in 1:6) for i in 1:6]
    f3 = Fa.f3 + evaluate(Fb.f3, Raz)
    return LieFactor(name, Fb.R * Fa.R, f3)
end

"Polynomial map of a factor: M = R o N3 o P4, truncated to MAP_ORDER."
function factor_to_map(F::LieFactor)
    N3 = [truncate_to(lie_exp(F.f3, Z[j]), MAP_ORDER) for j in 1:6]
    P4 = [truncate_to(lie_exp(F.f4, Z[j]), MAP_ORDER) for j in 1:6]
    N = compose_maps(N3, P4)
    return [truncate_to(sum(F.R[i, j] * N[j] for j in 1:6), MAP_ORDER) for i in 1:6]
end

end # module