# =============================================================================
# lie_flow.jl - element maps to arbitrary order by the Lie flow operator
# =============================================================================
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# WHAT
#   Element maps to order N by the Lie flow operator. N is a run-time argument.
#   The TaylorSeries path in lie_algebra.jl stays at generators f3, f4.
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
``p_s = \sqrt{(1+\delta)^2 - p_x^2 - p_y^2}`` and [`energy_term`](@ref) — the same
expressions as `hamiltonian`, evaluated to the order of `sp`:

| element | ``H`` |
|:--|:--|
| `Drift`, `LongCavity` | ``\mathrm{energy\_term}(\delta) - p_s`` |
| `Quadrupole` | ``\mathrm{energy\_term}(\delta) - p_s + \tfrac{k_1}{2}(x^2-y^2)`` |
| `Sextupole` | ``\mathrm{energy\_term}(\delta) - p_s + \tfrac{k_2}{6}(x^3-3xy^2)`` |
| `SBend` (body) | ``\mathrm{energy\_term}(\delta) - (1+hx)p_s + hx + \tfrac{h^2x^2}{2} + \tfrac{k_1}{2}(x^2-y^2)`` |

The constant is removed (it does not act on the map).
"""
function hamiltonian_da(e::AbstractElement, ref::RefParticle, sp::DASpace)
    x, px, y, py, l, d = da_vars(sp)
    ps = sqrt((1 + d)^2 - px * px - py * py); en = energy_term(d, ref)
    H = if e isa Drift || e isa LongCavity
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
``\mathcal M_e = e^{-L{:}H{:}} z``. H is built at order N+1; the series of F = S ∇H is
summed at order N. Truncation of a Lie series and truncation of a Taylor polynomial
are different operations. `SBend`: body only. Thin `RFCavity`: error.
"""
function lie_flow_map(e::AbstractElement, ref::RefParticle, N::Integer; tol = eps(), nmax = 60, nsquare = :auto)
    e isa RFCavity && error("lie_flow_map: thin RFCavity is a kick, not an s-independent element")
    sp = flow_space(N)
    e isa Marker && return identity_map(sp)
    L = elem_length(e); L == 0 && return identity_map(sp)
    H = hamiltonian_da(e, ref, flow_space(N + 1))
    F = [change_order(f, sp) for f in hamiltonian_vector_field(H)]
    vector_field_flow(F, L; tol = tol, nmax = nmax, nsquare = nsquare)
end

const _FLOWCACHE = Dict{Any,Vector{DA{Float64}}}()

"""
    element_flow_maps(lattice, N; cache = true) -> Vector{Vector{DA}}

[`lie_flow_map`](@ref) of every element of `lattice`, cached per (element, p0, N).
"""
function element_flow_maps(lattice::Lattice, N::Integer; cache::Bool = true)
    any(e -> e isa SBend && (e.e1 != 0 || e.e2 != 0), lattice.elements) &&
        @warn "element_flow_maps: SBend edge angles are not included in the flow maps (body only)"
    [cache ? get!(() -> lie_flow_map(e, lattice.ref, N), _FLOWCACHE, (e, lattice.ref.p0_MeV, Int(N))) :
             lie_flow_map(e, lattice.ref, N) for e in lattice.elements]
end

"""
    one_turn_map_flow(lattice, N; cache = true) -> Vector{DA}

One-turn map at order N from composition of element flow maps.
"""
function one_turn_map_flow(lattice::Lattice, N::Integer; cache::Bool = true)
    M = identity_map(flow_space(N))
    for Me in element_flow_maps(lattice, N; cache = cache)
        M = compose(Me, M)
    end
    M
end

end # module
