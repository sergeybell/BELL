# =============================================================================
# lie_flow.jl - element maps to arbitrary order by the Lie flow operator
# =============================================================================
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# WHAT
#   A second, independent map route next to the Dragt-Finn / Lie route of
#   lie_algebra.jl (which is NOT changed):
#     hamiltonian_da(e, ref, sp)   the FULL s-Hamiltonian of e (exact square root,
#                                  curvature, energy term) as a Bell.DiffAlgebra DA
#     lie_flow_map(e, ref, N)      M = exp(-L :H:) z, summed to machine precision,
#                                  exact through order N (H built at order N+1)
#     element_flow_maps(lat, N)    the same for every element (cached)
#     one_turn_map_flow(lat, N)    M_ring = M_n o ... o M_1 at order N
#   Drift, Quadrupole, Sextupole, SBend body (no edge / fringe terms), Marker,
#   Cavity (as a drift, RF off).  Thin RFCavity: error.
#
# WHY
#   The Lie route expands H to degree 4 (maps of order 3); here H is never
#   expanded by hand and the order N is a run-time parameter.
# =============================================================================
module LieFlow

using LinearAlgebra
using ..DiffAlgebra
using ..Coordinates, ..Elements, ..Structure

export hamiltonian_da, lie_flow_map, element_flow_maps, one_turn_map_flow, flow_space

const _SPACES = Dict{Tuple{Int,Int},DASpace}()
"""
    flow_space(N) -> DASpace

The (cached, shared) 6-variable DA space of order `N` used by the flow maps, so that maps of
different elements can be composed. Variables: ``z = (x, p_x, y, p_y, l, \\delta)``.
"""
flow_space(N::Integer) = get!(() -> DASpace(6, N), _SPACES, (6, Int(N)))

@doc raw"""
    hamiltonian_da(e, ref, sp) -> DA

Full s-Hamiltonian of `e` in the DA space `sp` (6 variables ``x, p_x, y, p_y, l, \delta``), with
``p_s = \sqrt{(1+\delta)^2 - p_x^2 - p_y^2}`` and ``e(\delta)`` = `energy_term` — the same
expressions as `hamiltonian`, but evaluated in ``{}_ND_6`` to the order of `sp` (no degree-4 cut):

| element | ``H`` |
|:--|:--|
| `Drift`, `Cavity` | ``e(\delta) - p_s`` |
| `Quadrupole` | ``e(\delta) - p_s + \tfrac{k_1}{2}(x^2-y^2)`` |
| `Sextupole` | ``e(\delta) - p_s + \tfrac{k_2}{6}(x^3-3xy^2)`` |
| `SBend` (body) | ``e(\delta) - (1+hx)p_s + hx + \tfrac{h^2x^2}{2} + \tfrac{k_1}{2}(x^2-y^2)`` |

The constant is removed (it does not act on the map).
"""
function hamiltonian_da(e::AbstractElement, ref::RefParticle, sp::DASpace)
    x, px, y, py, l, d = da_vars(sp)
    ps = sqrt((1 + d)^2 - px * px - py * py); en = energy_term(d, ref)
    H = if e isa Drift || e isa Cavity
        en - ps
    elseif e isa Quadrupole
        en - ps + e.k1 / 2 * (x * x - y * y)
    elseif e isa Sextupole
        en - ps + e.k2 / 6 * (x * x * x - 3x * y * y)
    elseif e isa SBend
        h = curvature(e)
        en - (1 + h * x) * ps + h * x + h^2 / 2 * x * x + e.k1 / 2 * (x * x - y * y)
    elseif e isa Marker
        zero(x)
    else
        error("hamiltonian_da: no s-independent Hamiltonian for $(typeof(e))")
    end
    H - constant_part(H)
end

@doc raw"""
    lie_flow_map(e, ref, N; tol = eps(), nmax = 60, nsquare = :auto) -> Vector{DA}

Transfer map of the s-independent element `e` to order ``N`` by the Lie flow operator
```math
\mathcal M_e = e^{-L{:}H{:}}\,z = \sum_{n\ge0}\frac{(-L)^n}{n!}{:}H{:}^n z,\qquad
{:}H{:}\,g = [H, g],
```
with the full Hamiltonian ([`hamiltonian_da`](@ref)) built in ``{}_{N+1}D_6`` and the series
its vector field ``F = S\nabla H`` kept through order ``N``, and the series summed in ``{}_ND_6`` to
machine precision (`DiffAlgebra.vector_field_flow`, scaling-and-squaring when ``L\,\|\partial F/\partial z\|``
is large); the result is in [`flow_space`](@ref)`(N)`
and is exact through order ``N`` up to roundoff. `SBend`: body only, edge angles ``e_1, e_2``
and fringe fields are not included. `Marker`: identity. Thin `RFCavity`: error.
"""
function lie_flow_map(e::AbstractElement, ref::RefParticle, N::Integer; tol = eps(), nmax = 60, nsquare = :auto)
    e isa RFCavity && error("lie_flow_map: thin RFCavity is a kick, not an s-independent element")
    sp = flow_space(N)
    e isa Marker && return identity_map(sp)
    L = elem_length(e); L == 0 && return identity_map(sp)
    H = hamiltonian_da(e, ref, flow_space(N + 1))                    # H through N+1 ...
    F = [change_order(f, sp) for f in hamiltonian_vector_field(H)]   # ... F = S grad H through N
    vector_field_flow(F, L; tol = tol, nmax = nmax, nsquare = nsquare) # series summed at order N
end

const _FLOWCACHE = Dict{Any,Vector{DA{Float64}}}()

"""
    element_flow_maps(lat, N; cache = true) -> Vector{Vector{DA}}

[`lie_flow_map`](@ref) of every element of `lat`, cached per (element, p0, N). SBend edge
angles are ignored (a warning is printed once per call if any is nonzero).
"""
function element_flow_maps(lat::Lattice, N::Integer; cache::Bool = true)
    any(e -> e isa SBend && (e.e1 != 0 || e.e2 != 0), lat.elements) &&
        @warn "element_flow_maps: SBend edge angles are not included in the flow maps (body only)"
    [cache ? get!(() -> lie_flow_map(e, lat.ref, N), _FLOWCACHE, (e, lat.ref.p0_MeV, Int(N))) :
             lie_flow_map(e, lat.ref, N) for e in lat.elements]
end

@doc raw"""
    one_turn_map_flow(lat, N; cache = true) -> Vector{DA}

One-turn map ``\mathcal M = \mathcal M_n\circ\dots\circ\mathcal M_1`` at order ``N`` from
[`element_flow_maps`](@ref) (`compose`, truncated to ``N`` in ``{}_ND_6``). Thin RF cavities are
not allowed (s-independent route). Independent of the Dragt-Finn product of `lie_algebra.jl`.
"""
function one_turn_map_flow(lat::Lattice, N::Integer; cache::Bool = true)
    M = identity_map(flow_space(N))
    for Me in element_flow_maps(lat, N; cache = cache)
        M = compose(Me, M)
    end
    M
end

end # module