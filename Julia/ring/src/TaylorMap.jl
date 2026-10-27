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
#                               truncated to MAP_ORDER after every step.
#   user_map(M, ref, C, h)      the same map in the user coordinates
#                               (x, x', y, y', phi, delta):  U = T o M o T^-1.
#   chromaticity, momentum_compaction  wrappers around MapAnalysis readings
#                               of whatever polynomial is passed in.
#
# WHY
#   Independent check of the Lie product: same Hamiltonians, integrated as
#   polynomials and composed by substitution.  The one-turn object of the
#   ring is the product of Lie factors, not this map.
#
# WHAT YOU SEE
#   Nothing printed here; scripts/02 prints the comparisons.
# =============================================================================
module TaylorMap

using TaylorSeries, LinearAlgebra
using ..Constants, ..Elements, ..Hamiltonians, ..Jets, ..LieOperators, ..MapAnalysis

export taylor_map_element, element_maps, ring_taylor_map, user_map,
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

"xi = dQ/ddelta from chromaticity_orders of the given polynomial."
function chromaticity(M; dd = 1e-5)
    ch = chromaticity_orders(M; dd = dd)
    return (ch.x[2], ch.y[2])
end

"alpha0, alpha1 from path_expansion of the given polynomial."
function momentum_compaction(M, ref::RefParticle, C::Real; dd = 1e-4)
    p = path_expansion(M, ref, C)
    return (p.alpha0, p.alpha1)
end

end # module
