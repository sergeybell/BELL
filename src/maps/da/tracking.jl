# =============================================================================
# maps/da/tracking.jl - magnet maps from the s-Hamiltonian on differential algebra
# =============================================================================
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# The model of an element is its Hamiltonian H, and dz/ds = S ∇H.
# Differential algebra is an algebra used to compute the map, not a second
# physical model. A magnet is integrated from H by ElementFlow, on the identity
# jet, to a chosen order. A Taylor map is symplectic only up to its truncation
# order. The field equations remain for comparison and for tests. A thin RF
# cavity is the kick. A thick cavity integrates H and needs the path length.
# The Lie computation uses the same H as Lie factors.
# =============================================================================
module DATracking

using LinearAlgebra
using ..DiffAlgebra, ..Integrators
using ..Coordinates, ..Elements, ..Structure, ..FieldTrack, ..ElementFlow

export map_space, element_da_map, element_da_maps, one_turn_map_da

const _SPACES = Dict{Tuple{Int,Int},DASpace}()

"""
    map_space(N) -> DASpace

Cached shared 6-variable DA space of order `N` for element maps
``z = (x, p_x, y, p_y, l, \\delta)``.
"""
map_space(N::Integer) = get!(() -> DASpace(6, N), _SPACES, (6, Int(N)))

@doc raw"""
    element_da_map(e, ref, N; method = DormandPrince853(), tol = 1e-14, nsteps = nothing) -> Vector{DA}

Transfer map of element `e` to order ``N``. A magnet and a thick
[`RFCavity`](@ref) integrate ``dz/ds = S\\nabla H`` ([`element_flow`](@ref))
on the DA identity jet. A [`Marker`](@ref) is the identity. A thin cavity
is [`rf_kick`](@ref) on the DA variables. An [`RFCavity`](@ref) needs
`path_length`. `SBend`: body only (no edge angles).
"""
function element_da_map(e::Marker, ref::RefParticle, N::Integer;
                        method::ODEMethod = DormandPrince853(), tol = 1e-14,
                        nsteps::Union{Nothing,Int} = nothing, path_length = nothing, ds = 0.01)
    return identity_map(map_space(N))
end

function element_da_map(e::RFCavity, ref::RefParticle, N::Integer;
                        method::ODEMethod = DormandPrince853(), tol = 1e-14,
                        nsteps::Union{Nothing,Int} = nothing, path_length = nothing, ds = 0.01)
    path_length === nothing && error("element_da_map(::RFCavity): pass path_length")
    z0 = collect(da_vars(map_space(N)))
    e.L == 0 && return rf_kick(z0, e, ref, path_length)
    plen = Float64(path_length)
    if nsteps === nothing
        return element_flow(e, ref, z0; s1 = e.L, method = method, tol = tol, path_length = plen)
    else
        return element_flow(e, ref, z0; s1 = e.L, method = method, tol = tol,
                            nsteps = nsteps, path_length = plen)
    end
end

function element_da_map(e::AbstractElement, ref::RefParticle, N::Integer;
                        method::ODEMethod = DormandPrince853(), tol = 1e-14,
                        nsteps::Union{Nothing,Int} = nothing, path_length = nothing, ds = 0.01)
    sp = map_space(N)
    L = elem_length(e)
    L == 0 && return identity_map(sp)
    z0 = collect(da_vars(sp))
    plen = path_length === nothing ? 0.0 : Float64(path_length)
    if nsteps === nothing
        return element_flow(e, ref, z0; s1 = L, method = method, tol = tol, path_length = plen)
    else
        return element_flow(e, ref, z0; s1 = L, method = method, tol = tol,
                            nsteps = nsteps, path_length = plen)
    end
end

# Key: element, reference particle, order, method value, tol, nsteps.
# `path_length` and `ds` do not change a non-RF magnet map (RF is not cached).
const _DACacheKey = Tuple{Any,RefParticle,Int,Any,Float64,Union{Nothing,Int}}
const _DACACHE = Dict{_DACacheKey,Vector{DA{Float64}}}()

function _da_cache_key(e, ref::RefParticle, N::Integer, kw)
    method = get(kw, :method, DormandPrince853())
    tol = Float64(get(kw, :tol, 1e-14))
    nsteps = get(kw, :nsteps, nothing)
    nkey = nsteps isa Integer ? Int(nsteps) : nsteps
    return (e, ref, Int(N), method, tol, nkey)::_DACacheKey
end

"""
    element_da_maps(lattice, N; cache = true, kw...) -> Vector{Vector{DA}}

[`element_da_map`](@ref) of every element. An [`RFCavity`](@ref) is built with
`path_length(lattice)` and is not cached. With `cache = true`, every other
element is stored under `(element, reference particle, N, method, tol, nsteps)`. `method` is
the integrator value, so two methods of one type with different parameters do
not share an entry. `path_length` and `ds` do not change a non-RF magnet map
and are not part of the key. A warning is given when an [`SBend`](@ref) has a
nonzero edge angle: the map is the body only.
"""
function element_da_maps(lattice::Lattice, N::Integer; cache::Bool = true, kw...)
    any(e -> e isa SBend && (e.e1 != 0 || e.e2 != 0), lattice.elements) &&
        @warn "element_da_maps: SBend edge angles are not included (body only)"
    out = Vector{DA{Float64}}[]
    for e in lattice.elements
        if e isa RFCavity
            rfkw = merge(kw, (; path_length = path_length(lattice)))
            push!(out, element_da_map(e, lattice.ref, N; rfkw...))
            continue
        end
        M = if cache
            get!(() -> element_da_map(e, lattice.ref, N; kw...), _DACACHE,
                 _da_cache_key(e, lattice.ref, N, kw))
        else
            element_da_map(e, lattice.ref, N; kw...)
        end
        push!(out, M)
    end
    out
end

"""
    one_turn_map_da(lattice, N; cache = true, kw...) -> Vector{DA}

One-turn DA map at order N from composition of element maps.
RF cavities are included (thin `rf_kick`, thick flow of `H`); use `magnetic_lattice` for the magnetic ring.
"""
function one_turn_map_da(lattice::Lattice, N::Integer; cache::Bool = true, kw...)
    M = identity_map(map_space(N))
    for Me in element_da_maps(lattice, N; cache = cache, kw...)
        M = compose(Me, M)
    end
    M
end

end # module
