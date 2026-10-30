# =============================================================================
# differential_algebra.jl - TaylorN truncated-series cross-check of the Lie maps
#   (module name DifferentialAlgebra kept for compatibility; this is NOT the general
#    DA engine - that is Bell.DiffAlgebra in src/da/DiffAlgebra.jl, with ODE
#    integrators in Bell.Integrators, src/integrators/Integrators.jl)
# =============================================================================
#
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
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
@doc raw"""
    DifferentialAlgebra

Truncated-series (TaylorSeries.jl `TaylorN`) cross-check route for the Lie maps:
element maps by RK4 in truncated polynomial arithmetic, composed into the ring map
and compared with the Lie product.  Scope is limited to that cross-check; the
general differential algebra `{}_nD_v` (arbitrary `n`, `v` at run time) is
Bell.DiffAlgebra, and generic ODE integrators are Bell.Integrators.
"""
module DifferentialAlgebra

using TaylorSeries, LinearAlgebra
using ..Coordinates, ..Elements, ..Structure, ..Hamiltonians, ..LieAlgebra, ..HigherOrder

export taylor_map_element, element_maps, ring_taylor_map, user_map,
       chromaticity, momentum_compaction

_V(H) = (g = [differentiate(H, j) for j in 1:6]; [sum(SYMPLECTIC_S[i, j] * g[j] for j in 1:6) for i in 1:6])
_ev(V, z) = [evaluate(V[i], z) for i in 1:6]

"""
    taylor_map_element(e, ref; ds = 0.01) -> Vector{TaylorN}

Taylor map of one element obtained by integrating Hamilton's equations
dz/ds = S ∇H in polynomial (differential-algebra) arithmetic.

# Algorithm
Start from the identity map z = Z. Do n = max(4, ⌈L/ds⌉) classical RK4 steps of
length h = L/n with the polynomial vector field S ∇H. Truncate to JET_ORDER after
every step and to MAP_ORDER at the end.

RK4 is of order 4 in h and not symplectic. This map is an independent check of the
Lie maps, not the production route. Markers give the identity; thin RF raises an error.
"""
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

"""
    element_maps(lat; ds = 0.01) -> Vector

[`taylor_map_element`](@ref) for every element of `lat`, cached per (element, p0, ds).
"""
element_maps(lat::Lattice; ds = 0.01) =
    [get!(() -> taylor_map_element(e, lat.ref; ds = ds), _TMCACHE, (e, lat.ref.p0_MeV, ds)) for e in lat.elements]

"""
    ring_taylor_map(lat; ds = 0.01) -> Vector{TaylorN}

One-turn Taylor map M_n ∘ ⋯ ∘ M_1 by [`compose_maps`](@ref) of the RK4 element maps.
Used as a cross-check of the Lie product, and in `one_turn_map_rf`.
"""
function ring_taylor_map(lat::Lattice; ds = 0.01)
    M = collect(Z)
    for Me in element_maps(lat; ds = ds)
        M = compose_maps(Me, M)
    end
    return M
end
"""
    user_map(M, ref, C, h) -> Vector{TaylorN}

The map `M` written in user coordinates (x, x', y, y', φ, δ):
U = T ∘ M ∘ T⁻¹, where T = [`canonical_to_user`](@ref) (evaluated as polynomials).
"""
function user_map(M, ref::RefParticle, C::Real, h::Real)
    Tinv = [truncate_to(c, MAP_ORDER) for c in user_to_canonical(collect(Z), ref, C, h)]
    T    = [truncate_to(c, MAP_ORDER) for c in canonical_to_user(collect(Z), ref, C, h)]
    return compose_maps(T, compose_maps(M, Tinv))
end

"""
    chromaticity(M; dd = 1e-5)

Linear chromaticities ξ = dQ/dδ from [`chromaticity_orders`](@ref) of the polynomial `M`.
"""
function chromaticity(M; dd = 1e-5)
    ch = chromaticity_orders(M; dd = dd)
    return (ch.x[2], ch.y[2])
end

"""
    momentum_compaction(M, ref, C; dd = 1e-4)

α0, α1 from [`path_expansion`](@ref) of the polynomial `M`.
"""
function momentum_compaction(M, ref::RefParticle, C::Real; dd = 1e-4)
    p = path_expansion(M, ref, C)
    return (p.alpha0, p.alpha1)
end

end # module
