# =============================================================================
# TaylorMap.jl - truncated Taylor map of the ring (COSY-style analogue)
# =============================================================================
#
# WHAT
#   An INDEPENDENT route to the ring map, used to cross-check the Lie code:
#   taylor_map_element(e, ref)  integrates Hamilton's equations
#                               dz/ds = S grad H(z)
#                               with classical RK4 in truncated polynomial
#                               (TaylorN) arithmetic, starting from the
#                               identity map - i.e. the map is obtained by
#                               "differential-algebra" integration, not by a
#                               Lie series.
#   ring_taylor_map(lat)        composition M_ring = M_n o ... o M_1,
#                               truncated to MAP_ORDER = 2 after every step.
#   user_map(M, ref, C, h)      the same map in the user coordinates
#                               (x, x', y, y', phi, delta):  U = T o M o T^-1.
#   closed_orbit(M, delta)      4D fixed point at fixed delta (Newton).
#   chromaticity(M)             dQ/ddelta from the tunes of the Jacobian on the
#                               off-momentum closed orbit (central difference).
#   momentum_compaction(M, ...) alpha_0, alpha_1 of  dC/C = a0 d + a1 d^2.
#
# WHY
#   COSY INFINITY produces exactly this kind of object (a truncated Taylor
#   map).  Building it in our own coordinates lets us compare (a) with our
#   Lie maps (must agree to order 2 up to integration error) and (b) with
#   COSY OUTPUT numbers (tunes, chromaticity, alpha_c) at the end.
#   Note: a truncated Taylor map is NOT exactly symplectic (see Symplectic.jl).
#
# WHAT YOU SEE
#   Nothing printed here; scripts/02 prints the comparisons.
# =============================================================================
module TaylorMap

using TaylorSeries, LinearAlgebra
using ..Constants, ..Elements, ..Hamiltonians, ..Jets, ..LieOperators

export taylor_map_element, element_maps, ring_taylor_map, user_map, closed_orbit, tunes_at_delta,
       chromaticity, momentum_compaction

_V(H) = (g = [differentiate(H, j) for j in 1:6]; [sum(SYMPLECTIC_S[i, j] * g[j] for j in 1:6) for i in 1:6])
_ev(V, z) = [evaluate(V[i], z) for i in 1:6]

"Order-2 Taylor map of one element by RK4 in TaylorN arithmetic (step <= ds)."
function taylor_map_element(e::AbstractElement, ref::RefParticle; ds = 0.01)
    e isa Marker && return collect(Z)
    e isa RFCavity && error("RF not included in the phase-1 Taylor map")
    L = elem_length(e)
    V = _V(hamiltonian(e, ref))
    n = max(4, ceil(Int, L / ds)); h = L / n
    z = collect(Z)
    for _ in 1:n
        k1 = _ev(V, z)
        k2 = _ev(V, z .+ (h / 2) .* k1)
        k3 = _ev(V, z .+ (h / 2) .* k2)
        k4 = _ev(V, z .+ h .* k3)
        z = [truncate_to(z[i] + h / 6 * (k1[i] + 2k2[i] + 2k3[i] + k4[i]), JET_ORDER) for i in 1:6]
    end
    return [truncate_to(zi, MAP_ORDER) for zi in z]
end

const _TMCACHE = Dict{Any,Any}()

"Taylor maps of all elements in order (cached per distinct element value and energy)."
element_maps(lat::Lattice; ds = 0.01) =
    [get!(() -> taylor_map_element(e, lat.ref; ds = ds), _TMCACHE, (e, lat.ref.p0_MeV, ds)) for e in lat.elements]

"Ring map by composing the element Taylor maps M_ring = M_n o ... o M_1."
function ring_taylor_map(lat::Lattice; ds = 0.01)
    M = collect(Z)
    for Me in element_maps(lat; ds = ds)
        M = compose_maps(Me, M)
    end
    return M
end
"Map in user coordinates: U = T o M o T^-1, T = canonical_to_user."
function user_map(M, ref::RefParticle, C::Real, h::Real)
    Tinv = [truncate_to(c, MAP_ORDER) for c in user_to_canonical(collect(Z), ref, C, h)]
    T    = [truncate_to(c, MAP_ORDER) for c in canonical_to_user(collect(Z), ref, C, h)]
    return compose_maps(T, compose_maps(M, Tinv))
end

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

"Chromaticities xi = dQ/ddelta (central difference, exact to first order for an order-2 map)."
function chromaticity(M; dd = 1e-5)
    qp = tunes_at_delta(M, dd); qm = tunes_at_delta(M, -dd)
    return ((qp[1] - qm[1]) / (2dd), (qp[2] - qm[2]) / (2dd))
end

"""
dC/C(delta) on the closed orbit: dC/C = -dl/C + (v/v0 - 1); fit a0, a1 from
delta = -dd, 0, +dd.  Returns (alpha0, alpha1).
"""
function momentum_compaction(M, ref::RefParticle, C::Real; dd = 1e-4)
    f(d) = (u = closed_orbit(M, d); z = [u; 0.0; d]; -_evalpt(M, z)[5] / C + (velocity_ratio(d, ref) - 1))
    fp, fc, fm = f(dd), f(0.0), f(-dd)   # (not "f0": 2f0 would parse as the Float32 literal 2.0)
    return ((fp - fm) / (2dd), (fp - 2 * fc + fm) / (2dd^2))
end

end # module
