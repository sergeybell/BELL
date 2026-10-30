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

"""
    closed_orbit(M, δ; iters = 30) -> Vector{Float64} (4)

Transverse closed orbit u = (x, px, y, py) of the one-turn map `M` at fixed
momentum deviation δ: solves M_{1:4}(u, 0, δ) = u by Newton iteration with the
exact Jacobian of the polynomial (stops at |Δu| < 1e-16). l is irrelevant
without RF.
"""
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

"""
    tunes_at_delta(M, δ) -> (Qx_frac, Qy_frac, u)

Fractional tunes from the 2×2 blocks of the Jacobian of `M` on the closed orbit
`u` at momentum δ (cos 2πQ = tr/2, sign of the sine from the (1,2) element).
"""
function tunes_at_delta(M, d::Real)
    u = closed_orbit(M, d)
    J = _jac(M, [u; 0.0; d])
    return (_frac_tune(J[1:2, 1:2]), _frac_tune(J[3:4, 3:4]), u)
end

"delta^k coefficient helper: quadratic part of M_i evaluated at vector w."
_q(M, i, w) = evaluate(TaylorN(M[i][2], JET_ORDER), w)

"""
    dispersion_orders(M) -> (D1, D2, w, R)

First- and second-order periodic dispersion at the start of the map.
The closed orbit is expanded as u(δ) = D1 δ + D2 δ² + O(δ³).

# Algorithm
- D1 = (I − R₄)⁻¹ R_{1:4,6};
- D2 = (I − R₄)⁻¹ q(w), with q the quadratic part of M_{1:4} evaluated at
  w = (D1, 0, 1).

Returns 4-vectors `D1`, `D2` (x, px, y, py) [m, 1], plus `w` and the linear part `R`.
"""
function dispersion_orders(M)
    R = map_linear_part(M)
    A = I - R[1:4, 1:4]
    D1 = A \ R[1:4, 6]
    w = [D1; 0.0; 1.0]
    D2 = A \ [_q(M, i, w) for i in 1:4]
    return (D1 = D1, D2 = D2, w = w, R = R)
end

@doc raw"""
    velocity_expansion(ref) -> (c1, c2, c3)

Exact Taylor coefficients ``v/v_0 - 1 = c_1\delta + c_2\delta^2 + c_3\delta^3``:
``c_1 = 1/\gamma_0^2``, ``c_2 = -3\beta_0^2/(2\gamma_0^2)``, ``c_3 = \beta_0^2(5\beta_0^2 - 1)/(2\gamma_0^2)``.
"""
function velocity_expansion(ref::RefParticle)
    g2 = ref.gamma^2; b2 = ref.beta^2
    c1 = 1 / g2
    c2 = -1.5 * b2 / g2
    c3 = 0.5 * b2 * (5b2 - 1) / g2
    return c1, c2, c3
end

@doc raw"""
    path_expansion(M, ref, C) -> NamedTuple

Momentum compaction and phase slip to second order, read from the one-turn map.

# Algorithm
On the dispersion orbit, ``l(\delta) = l_1\delta + l_2\delta^2`` with
``l_1 = R_{5,1:4}\cdot D_1 + R_{56}`` and ``l_2 = R_{5,1:4}\cdot D_2 + q_5(w)``. Since
``\Delta l = -\Delta C + C_0(v/v_0 - 1)``:

- ``\alpha_0 = -l_1/C + c_1``, ``\alpha_1 = -l_2/C + c_2`` (``c_i`` from [`velocity_expansion`](@ref));
- ``\eta_0 = -l_1/C``, ``\eta_1 = -l_2/C + (l_1/C)c_1``.

Definitions: ``C(\delta)/C_0 = 1 + \alpha_0\delta + \alpha_1\delta^2``, ``\Delta T/T_0 = \eta_0\delta + \eta_1\delta^2``.

# Returns
`(l1, l2, alpha0, alpha1, eta0, eta1, eta1_formula, gamma_t, c1, c2, c3)`;
`eta1_formula` is the same ``\eta_1`` from ``\alpha_0, \alpha_1`` (consistency check).
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

@doc raw"""
    chromaticity_orders(M; dd = 1e-4) -> (x = (Q0, ξ1, ξ2), y = (Q0, ξ1, ξ2))

Fractional tune and first and second chromaticity per plane, measured on the
polynomial map (not from an integral formula). For ``\delta \in \{-dd, 0, +dd\}``:
closed orbit by Newton ([`closed_orbit`](@ref)), Jacobian of ``M`` on that orbit, tune
from its 2×2 blocks ([`tunes_at_delta`](@ref)). Then

- ``\xi_1 = \dfrac{Q(+dd) - Q(-dd)}{2\,dd}``;
- ``\xi_2 = \dfrac{Q(+dd) - 2Q(0) + Q(-dd)}{dd^2}``, so ``Q(\delta) \approx Q_0 + \xi_1\delta + \xi_2\delta^2/2``.

The tunes are computed on the off-momentum closed orbit, so sextupole feed-down
is included. ``\xi_1`` is exact for the order-3 map; ``\xi_2`` includes all terms present in it.
"""
function chromaticity_orders(M; dd = 1e-4)
    q0 = tunes_at_delta(M, 0.0); qp = tunes_at_delta(M, dd); qm = tunes_at_delta(M, -dd)
    f(k) = (q0[k], (qp[k] - qm[k]) / (2dd), (qp[k] - 2 * q0[k] + qm[k]) / dd^2)
    return (x = f(1), y = f(2))     # (Q0_frac, xi1, xi2 = d2Q/ddelta2)
end

"""
    dispersion_along(lat, D1, D2) -> NamedTuple (s, D, Dp, D2, D2p)

D, D', D2, D2' at each element exit. The orbit u = D1 δ + D2 δ² (as a polynomial
in δ) is pushed through the Lie polynomials of the elements
([`element_lie_map`](@ref)), and the δ and δ² coefficients of x and px are read.
"""
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
