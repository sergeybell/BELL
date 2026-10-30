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

export hamiltonian_da, lie_flow_map, element_flow_maps, one_turn_map_flow, flow_space, chromaticity_series

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

@doc raw"""
    chromaticity_series(M; iters = N + 2) -> NamedTuple
    chromaticity_series(lat, N; cache = true) -> NamedTuple

Exact chromaticity series of a one-turn DA map ``\mathcal M`` (variables
``x, p_x, y, p_y, l, \delta``, order ``N``, e.g. [`one_turn_map_flow`](@ref)) - no finite
differences, no fit:

1. parametric closed orbit ``u(\delta)``, ``\mathcal M_{1..4}(u(\delta), 0, \delta) = u(\delta)``, as a DA
   series in ``\delta`` (Newton iteration with the constant matrix ``(R_0 - I)^{-1}``, one order
   per step);
2. the map shifted to the orbit, ``\tilde{\mathcal M}(z) = \mathcal M(u + z_\perp, l, \delta) - u``, and its
   linear part on the orbit ``R_{ij}(\delta) = \partial\tilde{\mathcal M}_i/\partial z_j|_{z_\perp = 0}`` as DA in
   ``\delta`` (exact through ``\delta^{N-1}``);
3. ``\cos\mu(\delta) = \tfrac12\operatorname{tr}R_{\rm plane}(\delta)``, solved for ``\mu`` by Newton in the DA,
   started from the branch ``\mu_0 \in (0,\pi)`` if ``R_{12}(0) > 0``, ``\mu_0 \in (\pi, 2\pi)`` otherwise;
   ``Q = \mu/2\pi`` (fractional tune).

Returns `(x = q, y = q, xi_x, xi_y, orbit, nmax)`: `q[n+1]` ``= q_n`` with
``Q(\delta) = \sum_n q_n\delta^n``, `xi[n]` ``= \xi_n = n!\,q_n = d^nQ/d\delta^n`` for
``n = 1\dots N-1`` (``q_N`` is not exact and is not returned), `orbit` = ``u(\delta)`` (4 DA).
Planes are treated as uncoupled (2x2 blocks). The `Lattice` method drops thin `RFCavity`
elements (``\delta`` is a parameter) and uses [`one_turn_map_flow`](@ref) (no SBend edges/fringes).
"""
function chromaticity_series(M::AbstractVector{<:DA}; iters::Integer = order(M[1]) + 2)
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
function chromaticity_series(lat::Lattice, N::Integer; cache::Bool = true)
    l2 = Lattice(lat.name, AbstractElement[e for e in lat.elements if !(e isa RFCavity)], lat.ref, lat.harmonic)
    chromaticity_series(one_turn_map_flow(l2, N; cache = cache))
end

end # module