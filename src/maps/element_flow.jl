# =============================================================================
# element_flow.jl - one integration of a magnet from its s-Hamiltonian
# =============================================================================
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# dz/ds = S ∇H for a Drift, Quadrupole, Sextupole, or SBend. Numbers and
# TaylorN evaluate S ∇H at each stage. A DA jet builds that vector field once
# and evaluates the polynomial on the current state. A Marker and L == 0 are
# the identity. An RF cavity is not this flow. The field equations stay in
# Fields for comparison and for tests.
# =============================================================================
module ElementFlow

using ..Coordinates, ..Elements, ..Hamiltonians, ..Integrators, ..DiffAlgebra
using ..Coordinates: SYMPLECTIC_S

export element_flow

# Monomial tables for the polynomial vector field. One order above the map:
# deriv drops the top degree, and the map needs that degree of S ∇H.
const _SPACES = Dict{Tuple{Int,Int},DASpace}()

function _field_space(nv::Integer, no::Integer)
    return get!(() -> DASpace(Int(nv), Int(no)), _SPACES, (Int(nv), Int(no)))
end

_da_state(z) = !isempty(z) && all(zi -> zi isa DA, z)

function _da_rhs(e, ref, z0)
    sp = z0[1].sp
    all(zi -> zi.sp == sp, z0) ||
        throw(ArgumentError("element_flow: DA components must share one space"))
    spV = _field_space(sp.nv, sp.no + 1)
    V = vector_field_on_identity(z -> hamiltonian(e, ref, z), spV; S = SYMPLECTIC_S)
    return (z, s) -> compose(V, z)
end

function _state_rhs(e, ref)
    return (z, s) -> symplectic_gradient(w -> hamiltonian(e, ref, w), z; S = SYMPLECTIC_S)
end

function _fixed_steps(s0, s1, nsteps, ds)
    if nsteps !== nothing
        n = Int(nsteps)
        n >= 1 || throw(ArgumentError("element_flow: nsteps must be >= 1, got $n"))
        return n
    end
    ds === nothing && return nothing
    step = abs(Float64(ds))
    step > 0 || throw(ArgumentError("element_flow: ds must be nonzero"))
    return max(1, ceil(Int, abs(Float64(s1) - Float64(s0)) / step))
end

@doc raw"""
    element_flow(e, ref, z0; s0 = 0, s1 = elem_length(e),
                 method = DormandPrince853(), tol = 1e-14,
                 nsteps = nothing, ds = nothing, path_length = 0) -> Vector

Integrate one magnet from its s-Hamiltonian, ``dz/ds = S\nabla H``, with
[`SYMPLECTIC_S`](@ref) and [`hamiltonian`](@ref)`(e, ref, z)`.

`z0` is a canonical 6-vector of numbers, `TaylorN`, or DA. A [`Marker`](@ref),
an element of length 0, or a call with `s1 == s0`, returns a copy of `z0`.

On numbers and on `TaylorN` each stage evaluates
[`symplectic_gradient`](@ref). On DA the polynomial vector field is built once
with [`vector_field_on_identity`](@ref), in a space one order above `z0`,
because [`deriv`](@ref) drops the top degree. The integrator evaluates that
polynomial on the current jet.

`nsteps` is a fixed step count. Otherwise `ds` is the step length. If both are
omitted, an adaptive method uses `tol` as the local error tolerance. A
fixed-step method then needs `nsteps` or `ds`.

`path_length` is accepted so callers can share keywords with the trackers.
A magnetic Hamiltonian does not use it. An [`RFCavity`](@ref) is rejected.
"""
function element_flow(e::AbstractElement, ref::RefParticle, z0;
                      s0::Real = 0.0, s1::Real = elem_length(e),
                      method::ODEMethod = DormandPrince853(), tol::Real = 1e-14,
                      nsteps::Union{Nothing,Integer} = nothing,
                      ds::Union{Nothing,Real} = nothing, path_length::Real = 0.0)
    z = collect(z0)
    (elem_length(e) == 0 || s1 == s0) && return z
    rhs = _da_state(z) ? _da_rhs(e, ref, z) : _state_rhs(e, ref)
    n = _fixed_steps(s0, s1, nsteps, ds)
    if n === nothing
        is_adaptive(method) ||
            throw(ArgumentError("$(typeof(method)) is fixed-step: pass nsteps or ds"))
        return integrate(rhs, z, s0, s1; method = method, tol = tol).z
    end
    return integrate_fixed(rhs, z, s0, s1, n; method = method).z
end

function element_flow(::Marker, ::RefParticle, z0;
                      s0::Real = 0.0, s1::Real = 0.0,
                      method::ODEMethod = DormandPrince853(), tol::Real = 1e-14,
                      nsteps::Union{Nothing,Integer} = nothing,
                      ds::Union{Nothing,Real} = nothing, path_length::Real = 0.0)
    return collect(z0)
end

function element_flow(::RFCavity, ::RefParticle, ::Any;
                      s0::Real = 0.0, s1::Real = 0.0,
                      method::ODEMethod = DormandPrince853(), tol::Real = 1e-14,
                      nsteps::Union{Nothing,Integer} = nothing,
                      ds::Union{Nothing,Real} = nothing, path_length::Real = 0.0)
    throw(ArgumentError("element_flow does not integrate an RFCavity; use exact_track_element"))
end

end # module
