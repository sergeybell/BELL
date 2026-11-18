# =============================================================================
# maps/da/tracking.jl - maps from the field equations on differential algebra
# =============================================================================
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# The model of an element is its Hamiltonian H, and dz/ds = S ∇H.
# Differential algebra is an algebra used to compute the map, not a second
# physical model. This computation integrates the field equations of motion,
# the alternative formulation kept for comparison, on series initial
# conditions to a chosen order. It does not use H. A Taylor map is symplectic
# only up to its truncation order. The Lie computation uses H and is not
# unified with this one.
# =============================================================================
module DATracking

using LinearAlgebra
using ..DiffAlgebra, ..Integrators
using ..Coordinates, ..Elements, ..Structure, ..Fields, ..FieldTrack

export map_space, element_da_map, element_da_maps, one_turn_map_da

const _SPACES = Dict{Tuple{Int,Int},DASpace}()

"""
    map_space(N) -> DASpace

Cached shared 6-variable DA space of order `N` for field maps
``z = (x, p_x, y, p_y, l, \\delta)``.
"""
map_space(N::Integer) = get!(() -> DASpace(6, N), _SPACES, (6, Int(N)))

@doc raw"""
    element_da_map(e, ref, N; method = DormandPrince853(), tol = 1e-14, nsteps = nothing) -> Vector{DA}

Transfer map of element `e` to order ``N`` by integrating the field equations of
motion on DA initial conditions ``z(0) =`` identity jet. A [`Marker`](@ref) is the
identity. An [`RFCavity`](@ref) needs `path_length`: thin uses [`rf_kick`](@ref),
thick is drift/kick/drift. `SBend`: body only (no edge angles).
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
    return exact_track_element(collect(da_vars(map_space(N))), e, ref, path_length; ds = ds)
end

function element_da_map(e::AbstractElement, ref::RefParticle, N::Integer;
                        method::ODEMethod = DormandPrince853(), tol = 1e-14,
                        nsteps::Union{Nothing,Int} = nothing, path_length = nothing, ds = 0.01)
    sp = map_space(N)
    L = elem_length(e)
    L == 0 && return identity_map(sp)
    z0 = collect(da_vars(sp))
    f = element_rhs(e, ref)
    if nsteps === nothing
        return integrate(f, z0, 0.0, L; method = method, tol = tol).z
    else
        return integrate_fixed(f, z0, 0.0, L, nsteps; method = method).z
    end
end

# Key: element, p0, order, method value, tol, nsteps.
# `path_length` and `ds` do not change a non-RF field map (RF is not cached).
const _DACacheKey = Tuple{Any,Float64,Int,Any,Float64,Union{Nothing,Int}}
const _DACACHE = Dict{_DACacheKey,Vector{DA{Float64}}}()

function _da_cache_key(e, p0::Real, N::Integer, kw)
    method = get(kw, :method, DormandPrince853())
    tol = Float64(get(kw, :tol, 1e-14))
    nsteps = get(kw, :nsteps, nothing)
    nkey = nsteps isa Integer ? Int(nsteps) : nsteps
    return (e, Float64(p0), Int(N), method, tol, nkey)::_DACacheKey
end

"""
    element_da_maps(lattice, N; cache = true, kw...) -> Vector{Vector{DA}}

[`element_da_map`](@ref) of every element. An [`RFCavity`](@ref) is built with
`path_length(lattice)` and is not cached. With `cache = true`, every other
element is stored under `(element, p0, N, method, tol, nsteps)`. `method` is
the integrator value, so two methods of one type with different parameters do
not share an entry. `path_length` and `ds` do not change a non-RF field map
and are not part of the key. A warning is given when an [`SBend`](@ref) has a
nonzero edge angle: the map is the body only.
"""
function element_da_maps(lattice::Lattice, N::Integer; cache::Bool = true, kw...)
    any(e -> e isa SBend && (e.e1 != 0 || e.e2 != 0), lattice.elements) &&
        @warn "element_da_maps: SBend edge angles are not included (body only)"
    out = Vector{DA{Float64}}[]
    for e in lattice.elements
        if e isa RFCavity
            push!(out, element_da_map(e, lattice.ref, N; path_length = path_length(lattice)))
            continue
        end
        M = if cache
            get!(() -> element_da_map(e, lattice.ref, N; kw...), _DACACHE,
                 _da_cache_key(e, lattice.ref.p0_MeV, N, kw))
        else
            element_da_map(e, lattice.ref, N; kw...)
        end
        push!(out, M)
    end
    out
end

"""
    one_turn_map_da(lattice, N; cache = true, kw...) -> Vector{DA}

One-turn DA map at order N from composition of element field maps.
RF cavities are included (exact `rf_kick` / thick splitting); use `magnetic_lattice` for the magnetic ring.
"""
function one_turn_map_da(lattice::Lattice, N::Integer; cache::Bool = true, kw...)
    M = identity_map(map_space(N))
    for Me in element_da_maps(lattice, N; cache = cache, kw...)
        M = compose(Me, M)
    end
    M
end

end # module
