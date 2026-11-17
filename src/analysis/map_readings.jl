# =============================================================================
# map_readings.jl - read coefficients of an already-built map
# =============================================================================
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# Linear optics of a matrix lives in linear.jl. This file reads coefficients
# of a polynomial or DA map up to the order that map already has. It does
# not integrate the lattice.
# =============================================================================
module MapReadings

using TaylorSeries, LinearAlgebra
using ..Coordinates, ..Elements, ..Structure, ..LieAlgebra, ..Linear, ..Tracking, ..DATracking, ..DiffAlgebra

export closed_orbit, tunes_at_delta, dispersion_orders, path_expansion,
       chromaticity_orders, dispersion_along, velocity_expansion, observables,
       chromaticity, momentum_compaction, chromaticity_series, find_optics_maxima

_evalpt(M, z) = [evaluate(M[i], z) for i in 1:6]
_jac(M, z) = [evaluate(differentiate(M[i], j), z) for i in 1:6, j in 1:6]

"""
    closed_orbit(M, δ; iters = 30) -> Vector{Float64} (4)

Transverse closed orbit u = (x, px, y, py) of a one-turn polynomial `M` at
fixed δ. Newton iteration on M_{1:4}(u, 0, δ) = u with the exact Jacobian of
the polynomial, at most `iters` steps, stop at |Δu| < 1e-16. This is the
fixed point of an already-built map. Building the map is the expensive step.
`periodic_optics` is the linear Twiss of a matrix; this function is the
off-momentum fixed point of a nonlinear map. Integer tunes are accumulated
in [`optics_along`](@ref).
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
`u` at momentum δ. The integer part of the tune is the accumulated phase in
[`optics_along`](@ref).
"""
function tunes_at_delta(M, d::Real)
    u = closed_orbit(M, d)
    J = _jac(M, [u; 0.0; d])
    return (_frac_tune(J[1:2, 1:2]), _frac_tune(J[3:4, 3:4]), u)
end

_q(M, i, w) = evaluate(TaylorN(M[i][2], JET_ORDER), w)

"""
    dispersion_orders(M) -> (D1, D2, w, R)

First- and second-order periodic dispersion at the start of the map.
The closed orbit is expanded as u(δ) = D1 δ + D2 δ² + O(δ³).
Linear dispersion of a matrix is [`periodic_optics`](@ref). This reading
adds the δ² term from the polynomial.
"""
function dispersion_orders(M)
    R = map_linear_part(M)
    A = I - R[1:4, 1:4]
    D1 = A \ R[1:4, 6]
    w = [D1; 0.0; 1.0]
    D2 = A \ [_q(M, i, w) for i in 1:4]
    return (D1 = D1, D2 = D2, w = w, R = R)
end

"""
    velocity_expansion(ref) -> (c1, c2, c3)

Taylor coefficients of v/v0 − 1 through δ³. Truncated at this order by design.
"""
function velocity_expansion(ref::RefParticle)
    g2 = ref.gamma^2; b2 = ref.beta^2
    c1 = 1 / g2
    c2 = -1.5 * b2 / g2
    c3 = 0.5 * b2 * (5b2 - 1) / g2
    return c1, c2, c3
end

"""
    path_expansion(M, ref, path_length) -> NamedTuple

Momentum compaction and phase slip to second order, read from the one-turn map.
"""
function path_expansion(M, ref::RefParticle, path_length::Real)
    d = dispersion_orders(M)
    R = d.R
    l1 = dot(R[5, 1:4], d.D1) + R[5, 6]
    l2 = dot(R[5, 1:4], d.D2) + _q(M, 5, d.w)
    c1, c2, c3 = velocity_expansion(ref)
    alpha0 = -l1 / path_length + c1
    alpha1 = -l2 / path_length + c2
    eta0 = -l1 / path_length
    eta1 = -l2 / path_length + (l1 / path_length) * c1
    eta1_formula = alpha1 - alpha0 * c1 - c2 + c1^2
    return (l1 = l1, l2 = l2, alpha0 = alpha0, alpha1 = alpha1, eta0 = eta0, eta1 = eta1,
            eta1_formula = eta1_formula, gamma_t = 1 / sqrt(abs(alpha0)), c1 = c1, c2 = c2, c3 = c3)
end

"""
    chromaticity_orders(M; dd = 1e-4) -> (x = (Q0, ξ1, ξ2), y = (Q0, ξ1, ξ2))

Fractional tune and first and second chromaticity by finite difference on the
closed orbit of the polynomial `M`. The series of arbitrary order is
[`chromaticity_series`](@ref).
"""
function chromaticity_orders(M; dd = 1e-4)
    q0 = tunes_at_delta(M, 0.0); qp = tunes_at_delta(M, dd); qm = tunes_at_delta(M, -dd)
    f(k) = (q0[k], (qp[k] - qm[k]) / (2dd), (qp[k] - 2 * q0[k] + qm[k]) / dd^2)
    return (x = f(1), y = f(2))
end

"""
    chromaticity(M; dd = 1e-5)

Linear chromaticities ξ = dQ/dδ from [`chromaticity_orders`](@ref).
"""
function chromaticity(M; dd = 1e-5)
    ch = chromaticity_orders(M; dd = dd)
    return (ch.x[2], ch.y[2])
end

"""
    momentum_compaction(M, ref, path_length)

α0, α1 from [`path_expansion`](@ref).
"""
function momentum_compaction(M, ref::RefParticle, path_length::Real)
    p = path_expansion(M, ref, path_length)
    return (p.alpha0, p.alpha1)
end

"""
    dispersion_along(lattice, D1, D2) -> NamedTuple

D, D', D2, D2' at each element exit.
"""
function dispersion_along(lattice::Lattice, D1, D2)
    d = Z[6]
    u = [D1[i] * d + D2[i] * d^2 for i in 1:4]
    z = [u; zero(d); d]
    s = 0.0
    along = (s = [0.0], D = [D1[1]], Dp = [D1[2]], D2 = [D2[1]], D2p = [D2[2]])
    for e in lattice.elements
        Me = e isa RFCavity ? element_lie_map(Drift(e.name, e.L), lattice.ref) : element_lie_map(e, lattice.ref)
        z = [truncate_to(evaluate(m, z), 2) for m in Me]
        s += elem_length(e)
        push!(along.s, s)
        push!(along.D, coeff(z[1], [0, 0, 0, 0, 0, 1])); push!(along.Dp, coeff(z[2], [0, 0, 0, 0, 0, 1]))
        push!(along.D2, coeff(z[1], [0, 0, 0, 0, 0, 2])); push!(along.D2p, coeff(z[2], [0, 0, 0, 0, 0, 2]))
    end
    return along
end

"""
    observables(lattice) -> NamedTuple

One-turn Lie factor, its polynomial, and the readings of that polynomial.
"""
function observables(lattice::Lattice)
    Lpath = path_length(lattice)
    F = ring_lie_factor(magnetic_lattice(lattice))
    M = factor_to_map(F)
    po = periodic_optics(F.R, lattice.ref, Lpath)
    pe = path_expansion(M, lattice.ref, Lpath)
    xi = chromaticity_orders(M)
    return (F = F, M = M, po = po, pe = pe, xi = xi, path_length = Lpath)
end

"""
    find_optics_maxima(opt; frac = 0.8, dmin = 0.5) -> Vector{NamedTuple}

Local maxima of D_x, β_x, β_y in the table from [`optics_along`](@ref).
A maximum counts if it exceeds `frac` × global maximum; maxima closer than
`dmin` [m] are merged. This reading does not know sextupole family names.
"""
function find_optics_maxima(opt; frac = 0.8, dmin = 0.5)
    d = opt.data
    found = NamedTuple[]
    for (q, arr) in (("D_x", d.D), ("beta_x", d.betx), ("beta_y", d.bety))
        mask = d.D .> dmin
        g = maximum(arr[mask])
        last_s = -1e9
        for k in 2:length(arr)-1
            (mask[k] && arr[k] >= arr[k-1] && arr[k] >= arr[k+1] && arr[k] >= frac * g) || continue
            d.s[k] - last_s < 0.5 && continue
            last_s = d.s[k]
            push!(found, (quantity = q, s = d.s[k], value = arr[k], betx = d.betx[k], bety = d.bety[k], D = d.D[k]))
        end
    end
    return found
end

@doc raw"""
    chromaticity_series(M; iters = N + 2) -> NamedTuple
    chromaticity_series(lattice, N; cache = true) -> NamedTuple

Exact chromaticity series of a one-turn DA map. The Lattice method drops thin
`RFCavity` elements and uses [`one_turn_map_da`](@ref) (field ODEs on DA).
"""
function chromaticity_series(M::AbstractVector{<:DA}; iters::Integer = M[1].sp.no + 2)
    sp = M[1].sp; N = sp.no
    z = da_vars(sp); d = z[6]; O = zero(d)
    A = inv(linear_part(M)[1:4, 1:4] - I)
    u = [O, O, O, O]
    for _ in 1:iters
        Mu = compose(M[1:4], [u; O; d])
        r = Mu .- u
        u = [u[i] - sum(A[i, j] * r[j] for j in 1:4) for i in 1:4]
    end
    Ms = compose(M[1:4], [u[1] + z[1], u[2] + z[2], u[3] + z[3], u[4] + z[4], z[5], d])
    on = [O, O, O, O, O, d]
    R = [compose([deriv(Ms[i], j)], on)[1] for i in 1:4, j in 1:4]
    function tune(a, b)
        c = (R[a, a] + R[b, b]) / 2
        c0 = constant_part(c); abs(c0) < 1 || error("chromaticity_series: plane $a unstable, |tr R/2| = $(abs(c0))")
        mu0 = acos(c0); constant_part(R[a, b]) < 0 && (mu0 = 2pi - mu0)
        mu = da_const(sp, mu0)
        for _ in 1:N + 2
            mu = mu + (cos(mu) - c) / sin(mu)
        end
        q = [coefficient(mu, [0, 0, 0, 0, 0, n]) / (2pi) for n in 0:N-1]
        q, [factorial(n) * q[n+1] for n in 1:N-1]
    end
    qx, xx = tune(1, 2); qy, xy = tune(3, 4)
    (x = qx, y = qy, xi_x = xx, xi_y = xy, orbit = u, nmax = N - 1)
end
function chromaticity_series(lattice::Lattice, N::Integer; cache::Bool = true)
    chromaticity_series(one_turn_map_da(magnetic_lattice(lattice), N; cache = cache))
end

end # module
