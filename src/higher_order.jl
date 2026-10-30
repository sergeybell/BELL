# =============================================================================
# MapAnalysis.jl - read coefficients of an already-built polynomial map
# =============================================================================
#
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# WHAT
#   Readers of a one-turn polynomial M.  They do not integrate the lattice.
#   closed_orbit, tunes_at_delta     off-momentum fixed point and tunes
#   dispersion_orders, path_expansion  D1, D2, alpha0/1, eta0/1
#   chromaticity_orders              Q(delta) by central difference
#   dispersion_along                 D(s), D2(s) via Lie element polynomials
#   observables(lat)                 one Lie product; Twiss, eta, xi from it
#
# WHY
#   Linear optics and second order are two readings of one map, not stacked
#   layers.  The RK4 Taylor map is a check, not the source of these numbers.
#
# WHAT YOU SEE
#   Nothing printed.
# =============================================================================
module HigherOrder

using TaylorSeries, LinearAlgebra
using ..Coordinates, ..Elements, ..Structure, ..LieAlgebra, ..Linear, ..Tracking

export closed_orbit, tunes_at_delta, dispersion_orders, path_expansion,
       chromaticity_orders, dispersion_along, velocity_expansion, observables

_evalpt(M, z) = [evaluate(M[i], z) for i in 1:6]
_jac(M, z) = [evaluate(differentiate(M[i], j), z) for i in 1:6, j in 1:6]

"4D closed orbit (x,px,y,py) at fixed delta; l is irrelevant (no RF)."
function closed_orbit(M, d::Real; iters = 30)
    u = zeros(4)
    for _ in 1:iters
        z = [u; 0.0; d]
        F = _evalpt(M, z)[1:4] .- u
        J = _jac(M, z)[1:4, 1:4] - I
        du = -J \ F
        u .+= du
        maximum(abs, du) < 1e-16 && break
    end
    return u
end

function _frac_tune(m)
    c = (m[1, 1] + m[2, 2]) / 2
    s = sign(m[1, 2]) * sqrt(max(0.0, 1 - c^2))
    return mod(atan(s, c), 2pi) / 2pi
end

"Fractional tunes of the Jacobian on the closed orbit at momentum deviation delta."
function tunes_at_delta(M, d::Real)
    u = closed_orbit(M, d)
    J = _jac(M, [u; 0.0; d])
    return (_frac_tune(J[1:2, 1:2]), _frac_tune(J[3:4, 3:4]), u)
end

"delta^k coefficient helper: quadratic part of M_i evaluated at vector w."
_q(M, i, w) = evaluate(TaylorN(M[i][2], JET_ORDER), w)

function dispersion_orders(M)
    R = map_linear_part(M)
    A = I - R[1:4, 1:4]
    D1 = A \ R[1:4, 6]
    w = [D1; 0.0; 1.0]
    D2 = A \ [_q(M, i, w) for i in 1:4]
    return (D1 = D1, D2 = D2, w = w, R = R)
end

"v/v0 - 1 = c1 d + c2 d^2 + c3 d^3 (exact Taylor of velocity_ratio)."
function velocity_expansion(ref::RefParticle)
    g2 = ref.gamma^2; b2 = ref.beta^2
    c1 = 1 / g2
    c2 = -1.5 * b2 / g2
    c3 = 0.5 * b2 * (5b2 - 1) / g2
    return c1, c2, c3
end

"""
path_expansion to order 2 (alpha0, alpha1) from the truncated map.
alpha_i defined by C(delta)/C0 = 1 + alpha0 d + alpha1 d^2 + alpha2 d^3
with Delta l = -C(d) + C0 (v/v0) on the closed orbit.
"""
function path_expansion(M, ref::RefParticle, C::Real)
    d = dispersion_orders(M)
    R = d.R
    l1 = dot(R[5, 1:4], d.D1) + R[5, 6]
    l2 = dot(R[5, 1:4], d.D2) + _q(M, 5, d.w)
    c1, c2, c3 = velocity_expansion(ref)
    alpha0 = -l1 / C + c1
    alpha1 = -l2 / C + c2
    eta0 = -l1 / C
    eta1 = -l2 / C + (l1 / C) * c1
    eta1_formula = alpha1 - alpha0 * c1 - c2 + c1^2
    return (l1 = l1, l2 = l2, alpha0 = alpha0, alpha1 = alpha1, eta0 = eta0, eta1 = eta1,
            eta1_formula = eta1_formula, gamma_t = 1 / sqrt(abs(alpha0)), c1 = c1, c2 = c2, c3 = c3)
end

function chromaticity_orders(M; dd = 1e-4)
    q0 = tunes_at_delta(M, 0.0); qp = tunes_at_delta(M, dd); qm = tunes_at_delta(M, -dd)
    f(k) = (q0[k], (qp[k] - qm[k]) / (2dd), (qp[k] - 2 * q0[k] + qm[k]) / dd^2)
    return (x = f(1), y = f(2))     # (Q0_frac, xi1, xi2 = d2Q/ddelta2)
end

"D, D', D2, D2' at each element exit; orbit pushed by Lie element polynomials."
function dispersion_along(lat::Lattice, D1, D2)
    d = Z[6]
    u = [D1[i] * d + D2[i] * d^2 for i in 1:4]
    z = [u; zero(d); d]
    s = 0.0
    out = (s = [0.0], D = [D1[1]], Dp = [D1[2]], D2 = [D2[1]], D2p = [D2[2]])
    for e in lat.elements
        Me = element_lie_map(e, lat.ref)
        z = [truncate_to(evaluate(m, z), 2) for m in Me]
        s += elem_length(e)
        push!(out.s, s)
        push!(out.D, coeff(z[1], [0, 0, 0, 0, 0, 1])); push!(out.Dp, coeff(z[2], [0, 0, 0, 0, 0, 1]))
        push!(out.D2, coeff(z[1], [0, 0, 0, 0, 0, 2])); push!(out.D2p, coeff(z[2], [0, 0, 0, 0, 0, 2]))
    end
    return out
end

"""
    observables(lat) -> NamedTuple

One-turn Lie factor of the lattice, its polynomial, and the readings of that
polynomial: periodic Twiss from F.R, eta from path_expansion, xi from
chromaticity_orders.  The RK4 Taylor map is not used.
"""
function observables(lat::Lattice)
    C = circumference(lat)
    F = ring_lie_factor(lat)
    M = factor_to_map(F)
    po = periodic_optics(F.R, lat.ref, C)
    pe = path_expansion(M, lat.ref, C)
    xi = chromaticity_orders(M)
    return (F = F, M = M, po = po, pe = pe, xi = xi, C = C)
end

end # module
