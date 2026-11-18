# =============================================================================
# fringe_maps.jl - edge maps built from the fringe Hamiltonian
# =============================================================================
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# HardEdge is a thin Lie generator. SoftEdge integrates the s-dependent
# fringe Hamiltonian: Splitting and Magnus stay products of Lie maps;
# TaylorIntegrate calls element_flow, which evaluates S ∇H(s) with the dual
# symplectic gradient at each stage. A pole-face angle is a thin wedge at
# local s = 0, the effective field boundary, between the two halves.
# =============================================================================
module FringeMaps

using TaylorSeries
using ..Coordinates, ..Elements, ..Structure, ..LieAlgebra, ..Hamiltonians, ..Fringe,
      ..Symplectic, ..Tracking, ..MapReadings, ..Integrators, ..ElementFlow
import ..Tracking: build_factors, ring_lie_factor

export edge_map, fringe_observables

# Parent, side, model key, full reference particle.
const _ECacheKey = Tuple{Any,Symbol,Any,RefParticle}
const _ECACHE = Dict{_ECacheKey,Vector{TaylorN{Float64}}}()

_as_generator(H::TaylorN) = H - constant_term(H)

# Cache key of an edge model: repr() for analytic models, plus the table contents for a
# Tabulated profile (its repr is only "Tabulated(N points)", so two different tables of the
# same length used to share one cache entry).
_model_key(m::EdgeModel) = repr(m)
_model_key(m::SoftEdge) = _soft_model_key(m, m.profile)
_soft_model_key(m::SoftEdge, ::FringeProfile) = repr(m)
_soft_model_key(m::SoftEdge, p::Tabulated) = (repr(m), Tuple(p.z), Tuple(p.B))

# ----------------------------------------------------------------------------- hard edges
function _hard(e::SBend, side, m::HardEdge, ref)
    x, px, y, py, l, d = Z
    h = curvature(e)
    ang = side === :in ? e.e1 : e.e2
    g = m.gap === nothing ? e.gap : m.gap
    fi = m.fint === nothing ? e.fint : m.fint
    psi = ang - h * g * fi * (1 + sin(ang)^2) / cos(ang)
    te, tp = tan(ang), tan(psi)
    K = -h * te * x^2 / 2
    if m.order <= 1
        K += h * tp * y^2 / 2
    else
        ps = sqrt((1 + d)^2 - px^2 - py^2)
        t = px / ps                                 # tan of the slope angle x'
        T = side === :in ? (tp + t) / (1 - tp * t) : (tp - t) / (1 + tp * t)
        K += h * y^2 / 2 * T
        m.order == 2 && (K = truncate_to(K, 3))
    end
    K -= constant_term(K)
    return lie_map(K, 1.0)
end

function _hard(e::Quadrupole, side, m::HardEdge, ref)
    m.order <= 1 && return collect(Z)               # no linear hard-edge effect
    x, px, y, py, l, d = Z
    s = side === :in ? 1.0 : -1.0
    K = s * e.k1 / (12 * (1 + d)) * (x^3 * px + 3x * y^2 * px - 3x^2 * y * py - y^3 * py)
    return lie_map(K, 1.0)
end

_hard(e::AbstractElement, side, m::HardEdge, ref) =
    error("HardEdge is implemented for SBend and Quadrupole, not for $(kind_string(e)) $(e.name)")

# ----------------------------------------------------------------------------- soft edges
_slice_count(E::Fringe._Edge, n::Int, s0::Real, s1::Real) =
    max(1, round(Int, n * abs(float(s1) - float(s0)) / (2 * E.W)))

function _soft(E::Fringe._Edge, m::Splitting, ref, k1, s0::Real, s1::Real)
    span = float(s1) - float(s0)
    span == 0 && return collect(Z)
    n = _slice_count(E, m.n, s0, s1)
    ws = splitting_weights(m.order >= 4 ? Yoshida(4) : Leapfrog())
    hs = span / n
    M = collect(Z)
    a = float(s0)
    for _ in 1:n
        for w in ws
            b = a + w * hs
            M = compose_maps(lie_map(_as_generator(Fringe._slice_hamiltonian(E, a, b, ref, Z; k1 = k1)), b - a), M)
            a = b
        end
    end
    return M
end

function _soft(E::Fringe._Edge, m::Magnus, ref, k1, s0::Real, s1::Real)
    span = float(s1) - float(s0)
    span == 0 && return collect(Z)
    n = _slice_count(E, m.n, s0, s1)
    hs = span / n
    wwin = max(E.ep, abs(hs) * sqrt(3) / 6)
    Hs(s) = _as_generator(fringe_hamiltonian(E, ref, s, Z; k1 = k1, w = wwin))
    M = collect(Z)
    for k in 0:n-1
        mid = float(s0) + (k + 0.5) * hs
        if m.order <= 2
            Om = hs * Hs(mid)
        else
            H1 = Hs(mid - hs * sqrt(3) / 6)
            H2 = Hs(mid + hs * sqrt(3) / 6)
            Om = hs / 2 * (H1 + H2) - sqrt(3) / 12 * hs^2 * poisson(H1, H2)
        end
        M = compose_maps(lie_map(Om, 1.0), M)
    end
    return M
end

function _flow_map(E::Fringe._Edge, ref, k1, s0::Real, s1::Real, n::Int)
    n >= 1 || return collect(Z)
    hs = (float(s1) - float(s0)) / n
    w = max(E.ep, abs(hs) / 4)
    H(z, s) = fringe_hamiltonian(E, ref, s, z; k1 = k1, w = w)
    zf = element_flow(H, collect(Z); s0 = s0, s1 = s1, method = RK4(), nsteps = n)
    return [truncate_to(zi, MAP_ORDER) for zi in zf]
end

function _soft(E::Fringe._Edge, m::TaylorIntegrate, ref, k1, s0::Real, s1::Real)
    frac = abs(float(s1) - float(s0)) / (2 * E.W)
    n = max(1, round(Int, 16 * frac))
    nmax = max(n, round(Int, m.nmax * frac))
    M = _flow_map(E, ref, k1, s0, s1, n)
    while 2n <= nmax
        M2 = _flow_map(E, ref, k1, s0, s1, 2n)
        diff = maximum(maxabs(M2[i] - M[i]) for i in 1:6)
        M, n = M2, 2n
        diff < m.tol && break
    end
    return M
end

# Local s = 0 is the effective field boundary. The thin wedge is applied there,
# between the upstream and downstream halves. Placing it before or after the
# whole fringe puts the kick on the field-free end, a distance W off the EFB.
function _soft_region(E::Fringe._Edge, method::FringeMethod, ref, k1, ang)
    ang == 0 && return _soft(E, method, ref, k1, -E.W, E.W)
    left = _soft(E, method, ref, k1, -E.W, 0.0)
    right = _soft(E, method, ref, k1, 0.0, E.W)
    Wg = lie_map(-E.h0 * tan(ang) * (Z[1]^2 - Z[3]^2) / 2, 1.0)
    return compose_maps(right, compose_maps(Wg, left))
end

function _soft_edge(e::SBend, side, m::SoftEdge, ref)
    E = Fringe.edge_field(e, m.profile, side)
    ang = side === :in ? e.e1 : e.e2
    S = _soft_region(E, m.method, ref, e.k1, ang)
    W = E.W
    Bm = lie_map(hamiltonian(e, ref), -W)
    Dm = lie_map(hamiltonian(Drift("_fringe_drift", 1.0), ref), -W)
    return side === :in ? compose_maps(Bm, compose_maps(S, Dm)) : compose_maps(Dm, compose_maps(S, Bm))
end
_soft_edge(e::AbstractElement, side, m::SoftEdge, ref) =
    error("SoftEdge is implemented for SBend only, not for $(kind_string(e)) $(e.name)")

_edge_polynomial(e, side, ::NoEdge, ref) = collect(Z)
_edge_polynomial(e, side, m::HardEdge, ref) = _hard(e, side, m, ref)
_edge_polynomial(e, side, m::SoftEdge, ref) = _soft_edge(e, side, m, ref)

"""
    edge_map(E::EdgeMap, ref) -> Vector{TaylorN}

Polynomial map (order MAP_ORDER) of a thin edge item: identity for `NoEdge`, the
hard-edge Lie map for `HardEdge`, or the composed soft-edge map for `SoftEdge`.
Cached per (parent magnet, side, model, reference particle).
"""
function edge_map(E::EdgeMap, ref::RefParticle)
    return get!(_ECACHE, (E.parent, E.side, _model_key(E.model), ref)::_ECacheKey) do
        _edge_polynomial(E.parent, E.side, E.model, ref)
    end
end

# ----------------------------------------------------------------------------- ring
_sequence_factor(it::EdgeMap, lattice, cache) =
    get!(() -> dragt_finn_factor(edge_map(it, lattice.ref); name = it.name), cache,
         (it.parent, it.side, _model_key(it.model)))
_sequence_factor(it::AbstractElement, lattice, cache) =
    get!(() -> dragt_finn_factor(element_lie_map(it, lattice.ref, path_length(lattice)); name = it.name),
         cache, it)

"One factor per item of fringe_sequence (edges are separate factors)."
function build_factors(lattice::Lattice, model::FringeModel)
    cache = Dict{Any,LieFactor}()
    factors = LieFactor[]
    for it in fringe_sequence(lattice, model)
        push!(factors, _sequence_factor(it, lattice, cache))
    end
    return factors
end

"One-turn Lie factor of the ring with the given fringe model."
ring_lie_factor(lattice::Lattice, model::FringeModel; name = lattice.name) =
    product_factors(build_factors(lattice, model); name)

_edge_symplecticity(::AbstractElement, ref, z0) = (0.0, 0.0)
function _edge_symplecticity(it::EdgeMap, ref, z0)
    E = edge_map(it, ref)
    return (symplectic_error(map_linear_part(E)), symplectic_error(jacobian_of_map(E, z0)))
end

"""
    fringe_observables(lattice, model) -> NamedTuple

Fractional tunes, first/second-order chromaticity of the one-turn map and the
symplecticity of the edge maps (linear part, and Jacobian at a test point z0).
"""
function fringe_observables(lattice::Lattice, model::FringeModel = FringeModel(); dd = 1e-4,
                            z0 = [1e-3, 1e-4, 1e-3, 1e-4, 0.0, 1e-3])
    t0 = time()
    F = ring_lie_factor(lattice, model)
    M = factor_to_map(F)
    ch = chromaticity_orders(M; dd)
    se_lin, se_pt = 0.0, 0.0
    for it in fringe_sequence(lattice, model)
        sl, sp = _edge_symplecticity(it, lattice.ref, z0)
        se_lin = max(se_lin, sl)
        se_pt = max(se_pt, sp)
    end
    return (Qx = ch.x[1], Qy = ch.y[1], xix = ch.x[2], xiy = ch.y[2], xi2x = ch.x[3], xi2y = ch.y[3],
            symp_edge_lin = se_lin, symp_edge_pt = se_pt, symp_ring = symplectic_error(F.R),
            seconds = time() - t0, model = model)
end

end # module
