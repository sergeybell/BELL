# =============================================================================
# da/tracking.jl - element maps by integrating field ODEs on DA
# =============================================================================
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# Same equations as Fields / FieldTrack, integrated on DA variables to order N.
# Does not build or call Hamiltonians. Lie-path maps stay in Tracking / LieAlgebra.
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
motion on DA initial conditions ``z(0) =`` identity jet. Thin [`RFCavity`](@ref):
use [`rf_kick`](@ref) on DA variables. `SBend`: body only (no edge angles).
"""
function element_da_map(e::AbstractElement, ref::RefParticle, N::Integer;
                        method::ODEMethod = DormandPrince853(), tol = 1e-14,
                        nsteps::Union{Nothing,Int} = nothing, path_length = nothing, ds = 0.01)
    if e isa RFCavity
        path_length === nothing && error("element_da_map(::RFCavity): pass path_length")
        return exact_track_element(collect(da_vars(map_space(N))), e, ref, path_length; ds = ds)
    end
    sp = map_space(N)
    e isa Marker && return identity_map(sp)
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

const _DACACHE = Dict{Any,Vector{DA{Float64}}}()

"""
    element_da_maps(lattice, N; cache = true, kw...) -> Vector{Vector{DA}}

[`element_da_map`](@ref) of every non-RF element; thin RF slots are skipped
with a warning if present (compose RF separately via [`rf_kick`](@ref)).
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
                 (e, lattice.ref.p0_MeV, Int(N), typeof(get(kw, :method, DormandPrince853()))))
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
