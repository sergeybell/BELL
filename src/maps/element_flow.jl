# =============================================================================
# element_flow.jl - one integration of a magnet from its s-Hamiltonian
# =============================================================================
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# dz/ds = S ∇H for a Drift, Quadrupole, Sextupole, SBend, or a thick RFCavity.
# Numbers and TaylorN evaluate S ∇H at each stage. A DA jet builds that vector
# field once and evaluates the polynomial on the current state. A Marker and a
# magnet with L == 0 are the identity. A thin RF cavity is rejected. The field
# equations stay in Fields for comparison and for tests.
#
# element_flow(H, z0; s0, s1) integrates an explicit scalar H(z, s). Each stage
# evaluates the dual symplectic gradient of w -> H(w, s). A soft edge uses that
# path because the fringe Hamiltonian depends on s.
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

function _flow_rhs(H, z)
    if _da_state(z)
        sp = z[1].sp
        all(zi -> zi.sp == sp, z) ||
            throw(ArgumentError("element_flow: DA components must share one space"))
        spV = _field_space(sp.nv, sp.no + 1)
        V = vector_field_on_identity(H, spV; S = SYMPLECTIC_S)
        return (w, s) -> compose(V, w)
    end
    return (w, s) -> symplectic_gradient(H, w; S = SYMPLECTIC_S)
end

function _integrate_rhs(z, s0, s1, method, tol, nsteps, ds, rhs)
    n = _fixed_steps(s0, s1, nsteps, ds)
    if n === nothing
        is_adaptive(method) ||
            throw(ArgumentError("$(typeof(method)) is fixed-step: pass nsteps or ds"))
        return integrate(rhs, z, s0, s1; method = method, tol = tol).z
    end
    return integrate_fixed(rhs, z, s0, s1, n; method = method).z
end

function _integrate_flow(z, s0, s1, method, tol, nsteps, ds, H)
    return _integrate_rhs(z, s0, s1, method, tol, nsteps, ds, _flow_rhs(H, z))
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
    element_flow(H, z0; s0 = 0, s1, method = DormandPrince853(), tol = 1e-14,
                 nsteps = nothing, ds = nothing) -> Vector

Integrate one element from its s-Hamiltonian, ``dz/ds = S\nabla H``, with
[`SYMPLECTIC_S`](@ref). A magnet uses [`hamiltonian`](@ref)`(e, ref, z)`.
A thick [`RFCavity`](@ref) uses [`hamiltonian`](@ref)`(e, ref, path_length, z)`
and needs a positive `path_length`. A thin cavity is rejected: the map is
[`rf_kick`](@ref) and the generator is [`rf_potential`](@ref).

`element_flow(H, z0; s0, s1)` integrates a scalar `H(z, s)`. `s1` is required.
Each stage evaluates [`symplectic_gradient`](@ref) of `w -> H(w, s)`, including
on `TaylorN` and on DA. The magnet methods stay on the autonomous path: one
DA vector field, built once.

`z0` is a canonical 6-vector of numbers, `TaylorN`, or DA. A [`Marker`](@ref),
a magnet of length 0, or a call with `s1 == s0`, returns a copy of `z0`.

On numbers and on `TaylorN` each stage evaluates
[`symplectic_gradient`](@ref). On DA the polynomial vector field is built once
with [`vector_field_on_identity`](@ref), in a space one order above `z0`,
because [`deriv`](@ref) drops the top degree. The integrator evaluates that
polynomial on the current jet.

`nsteps` is a fixed step count. Otherwise `ds` is the step length. If both are
omitted, an adaptive method uses `tol` as the local error tolerance. A
fixed-step method then needs `nsteps` or `ds`.

`path_length` is accepted so callers can share keywords with the trackers.
A magnetic Hamiltonian does not use it. A thick cavity does.
"""
function element_flow(e::AbstractElement, ref::RefParticle, z0;
                      s0::Real = 0.0, s1::Real = elem_length(e),
                      method::ODEMethod = DormandPrince853(), tol::Real = 1e-14,
                      nsteps::Union{Nothing,Integer} = nothing,
                      ds::Union{Nothing,Real} = nothing, path_length::Real = 0.0)
    z = collect(z0)
    (elem_length(e) == 0 || s1 == s0) && return z
    return _integrate_flow(z, s0, s1, method, tol, nsteps, ds, w -> hamiltonian(e, ref, w))
end

function element_flow(::Marker, ::RefParticle, z0;
                      s0::Real = 0.0, s1::Real = 0.0,
                      method::ODEMethod = DormandPrince853(), tol::Real = 1e-14,
                      nsteps::Union{Nothing,Integer} = nothing,
                      ds::Union{Nothing,Real} = nothing, path_length::Real = 0.0)
    return collect(z0)
end

function element_flow(e::RFCavity, ref::RefParticle, z0;
                      s0::Real = 0.0, s1::Real = elem_length(e),
                      method::ODEMethod = DormandPrince853(), tol::Real = 1e-14,
                      nsteps::Union{Nothing,Integer} = nothing,
                      ds::Union{Nothing,Real} = nothing, path_length::Real = 0.0)
    e.L > 0 || throw(ArgumentError(
        "element_flow does not integrate a thin RFCavity; use rf_kick or rf_potential"))
    z = collect(z0)
    s1 == s0 && return z
    path_length > 0 || throw(ArgumentError(
        "element_flow(::RFCavity) needs a positive path_length"))
    return _integrate_flow(z, s0, s1, method, tol, nsteps, ds,
                           w -> hamiltonian(e, ref, path_length, w))
end

function element_flow(H::Function, z0;
                      s0::Real = 0.0, s1::Real,
                      method::ODEMethod = DormandPrince853(), tol::Real = 1e-14,
                      nsteps::Union{Nothing,Integer} = nothing,
                      ds::Union{Nothing,Real} = nothing)
    z = collect(z0)
    s1 == s0 && return z
    rhs = (w, s) -> symplectic_gradient(u -> H(u, s), w; S = SYMPLECTIC_S)
    return _integrate_rhs(z, s0, s1, method, tol, nsteps, ds, rhs)
end

end # module
