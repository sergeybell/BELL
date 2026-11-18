# =============================================================================
# symplectic.jl - one check: max |J^T S J - S|
# =============================================================================
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# symplectic_error accepts anything that reduces to a Jacobian or a series map:
# a numeric matrix, a tracking step z -> z, a TaylorSeries polynomial, a DA map.
# The defect is always max |J^T S J - S|. For a series, that maximum runs over
# every retained coefficient. DA arithmetic stays in DiffAlgebra; this file
# is the call site. The TaylorN method, and the Jacobians that probe a map
# with TaylorN, are added by LieAlgebra so this module does not import it.
# SYMPLECTIC_S is the coordinate pairing in Coordinates.
# =============================================================================
module Symplectic

using LinearAlgebra
using ..Coordinates: SYMPLECTIC_S
using ..DiffAlgebra: DA, da_symplectic_error

export symplectic_error, jacobian_of_map, jacobian_of_function

"""
    symplectic_error(J::AbstractMatrix) -> Float64
    symplectic_error(M::AbstractVector{<:TaylorN}) -> Float64
    symplectic_error(M::AbstractVector{<:DA}; kwargs...) -> Float64
    symplectic_error(f, z0) -> Float64

Single symplecticity check, ``\\max|J^T S J - S|``, for every object Bell tracks:

- a numeric Jacobian;
- a `TaylorN` map (maximum over coefficients through degree `MAP_ORDER - 1`;
  that method is defined in `LieAlgebra`);
- a DA map, computed by [`da_symplectic_error`](@ref Bell.DiffAlgebra.da_symplectic_error);
- a tracking step `f` at the point `z0` (Jacobian from [`jacobian_of_function`](@ref)).

Zero (up to round-off) for an exactly symplectic linear map. A truncated series
is symplectic only through its retained degree.
"""
function symplectic_error end

symplectic_error(J::AbstractMatrix) = maximum(abs, J' * SYMPLECTIC_S * J - SYMPLECTIC_S)

symplectic_error(M::AbstractVector{<:DA}; kwargs...) = da_symplectic_error(M; kwargs...)

symplectic_error(f, z0::AbstractVector) = symplectic_error(jacobian_of_function(f, z0))

"""
    jacobian_of_map(M, z0) -> Matrix{Float64}

Jacobian ∂M_i/∂z_j of a polynomial map evaluated at the point `z0` (exact
differentiation of the polynomial). The method is defined in `LieAlgebra`.
"""
function jacobian_of_map end

"""
    jacobian_of_function(f, z0) -> Matrix{Float64}

Exact Jacobian of a generic tracking function `f` (any code that works on
`TaylorN`, e.g. `z -> track(factors, z)` or `exact_track_turn`) at `z0`.
The method is defined in `LieAlgebra`.

# Algorithm
Calls `f(z0 .+ Z)` once in polynomial arithmetic; the linear coefficients of the
output are the Jacobian. This is automatic differentiation; no finite differences.
"""
function jacobian_of_function end

end # module
