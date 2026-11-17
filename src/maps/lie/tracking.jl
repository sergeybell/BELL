# =============================================================================
# tracking.jl - symplectic tracking with factorised Lie maps
# =============================================================================
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# The ring object is a product of Lie factors. Tracking applies each factor
# (flow of f4, then f3, then R). RFCavity maps come from its Hamiltonian like any element (v0.17).
# A Taylor-series map of the same Hamiltonian is obtained by integrating
# hamiltonian_field with a chosen ODE method; it stores no copy of H.
# =============================================================================
module Tracking

using TaylorSeries, LinearAlgebra
using ..Coordinates, ..Elements, ..Structure, ..Hamiltonians, ..LieAlgebra, ..Integrators

export build_factors, apply_factor, track, element_lie_map, ring_lie_factor,
       taylor_map_element, element_maps, ring_taylor_map, laboratory_map

const _MCACHE = Dict{Any,Any}()

"""
    element_lie_map(e, ref; order = MAP_ORDER) -> Vector{TaylorN}

Map exp(−L :H:) z of one element ([`lie_map`](@ref) of its [`hamiltonian`](@ref)),
truncated to `order`; identity for markers. Cached per (element, p0, order).
"""
element_lie_map(e::AbstractElement, ref::RefParticle, path_length::Real = 0.0; order::Int = MAP_ORDER) =
    get!(_MCACHE, (e, ref.p0_MeV, e isa RFCavity ? Float64(path_length) : 0.0, order)) do
        if e isa Marker
            collect(Z)
        elseif e isa RFCavity && e.L == 0
            G = rf_potential(e, ref, path_length)
            lie_map(G - constant_term(G), 1.0; order)       # thin limit: exp(-:G:) = rf_kick
        else
            lie_map(hamiltonian(e, ref, path_length), elem_length(e); order)
        end
    end

"""
    ring_lie_factor(lattice; name = lattice.name) -> LieFactor

One-turn Dragt-Finn factor of the magnetic ring: [`product_factors`](@ref) of
[`build_factors`](@ref)(lattice), RF cavities included (use `magnetic_lattice` for the magnetic ring).
"""
function ring_lie_factor(lattice::Lattice; name = lattice.name)
    return product_factors(LieFactor[F for F in build_factors(lattice)]; name)
end

"""
    build_factors(lattice) -> Vector

One [`LieFactor`](@ref) per element. RF cavities too: thin = exp(-:G:) of [`rf_potential`](@ref) (= [`rf_kick`](@ref)),
thick = exp(-L:H:).
"""
function build_factors(lattice::Lattice)
    cache = Dict{Any,LieFactor}()
    factors = LieFactor[]
    Lpath = path_length(lattice)
    for e in lattice.elements
        push!(factors, get!(() -> dragt_finn_factor(element_lie_map(e, lattice.ref, Lpath); name = e.name), cache, e))
    end
    return factors
end

"""
    apply_nonlinear(F, z) -> Vector

Apply the non-linear part of a factor, N3 ∘ P4, each by the implicit midpoint
rule of [`implicit_midpoint`](@ref Bell.Integrators.implicit_midpoint).
"""
function apply_nonlinear(F::LieFactor, z::AbstractVector)
    ev4 = (u, s) -> [evaluate(F.field4[i], u) for i in 1:6]
    ev3 = (u, s) -> [evaluate(F.field[i], u) for i in 1:6]
    w = iszero(F.f4) ? collect(z) : integrate_fixed(ev4, collect(z), 0.0, 1.0, 2; method = ImplicitMidpoint()).z
    return integrate_fixed(ev3, w, 0.0, 1.0, 2; method = ImplicitMidpoint()).z
end

"""
    apply_factor(F::LieFactor, z) -> Vector

Track one point through a factor: w = N3(P4(z)), then R·w.
"""
apply_factor(F::LieFactor, z) = (w = apply_nonlinear(F, z); [sum(F.R[i, j] * w[j] for j in 1:6) for i in 1:6])

"""
    track(factors, z; turns = 1, ref = nothing, path_length = 0.0) -> Vector

Track a canonical point `z` for `turns` turns. Lie factors use [`apply_factor`](@ref);
Every element (RF included) is a Lie factor; `ref`/`path_length` are kept for compatibility and ignored.
"""
function track(factors, z; turns::Int = 1, ref = nothing, path_length = 0.0)
    w = collect(z)
    for n in 1:turns, F in factors
        w = apply_factor(F, w)
    end
    return w
end

const _TMCACHE = Dict{Any,Any}()

"""
    taylor_map_element(e, ref; ds = 0.01, method = RK4(), order = MAP_ORDER) -> Vector{TaylorN}

Taylor map of one element: integrate ``dz/ds = S\\nabla H`` of [`hamiltonian`](@ref)
with `method` in polynomial arithmetic, starting from the identity. Truncated
to `order`. This is a check of the Lie map, not the production route.
"""
function taylor_map_element(e::AbstractElement, ref::RefParticle; ds = 0.01,
                            method::ODEMethod = RK4(), order::Int = MAP_ORDER, path_length::Real = 0.0)
    e isa Marker && return collect(Z)
    thin = e isa RFCavity && e.L == 0
    L = thin ? 1.0 : elem_length(e)
    V = hamiltonian_field(thin ? rf_potential(e, ref, path_length) : hamiltonian(e, ref, path_length))
    n = max(4, ceil(Int, L / ds))
    f = (z, s) -> [evaluate(V[i], z) for i in 1:6]
    r = integrate_fixed(f, collect(Z), 0.0, L, n; method)
    return [truncate_to(zi, order) for zi in r.z]
end

"""
    element_maps(lattice; ds = 0.01) -> Vector

[`taylor_map_element`](@ref) for every element of `lattice`, cached per (element, p0, ds).
"""
element_maps(lattice::Lattice; ds = 0.01) =
    [get!(() -> taylor_map_element(e, lattice.ref; ds = ds, path_length = path_length(lattice)), _TMCACHE, (e, lattice.ref.p0_MeV, path_length(lattice), ds)) for e in lattice.elements]

"""
    ring_taylor_map(lattice; ds = 0.01) -> Vector{TaylorN}

One-turn Taylor map by composition of the element maps.
"""
function ring_taylor_map(lattice::Lattice; ds = 0.01)
    M = collect(Z)
    for Me in element_maps(lattice; ds = ds)
        M = compose_maps(Me, M)
    end
    return M
end

"""
    laboratory_map(M, ref, path_length, h) -> Vector{TaylorN}

The map `M` written in laboratory coordinates:
U = T ∘ M ∘ T⁻¹, where T = [`canonical_to_laboratory`](@ref).
"""
function laboratory_map(M, ref::RefParticle, path_length::Real, h::Real)
    Tinv = [truncate_to(c, MAP_ORDER) for c in laboratory_to_canonical(collect(Z), ref, path_length, h)]
    T    = [truncate_to(c, MAP_ORDER) for c in canonical_to_laboratory(collect(Z), ref, path_length, h)]
    return compose_maps(T, compose_maps(M, Tinv))
end

end # module
