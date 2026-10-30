# =============================================================================
# Symplectic.jl - symplecticity checks  || J^T S J - S ||
# =============================================================================
#
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# WHAT
#   symplectic_error(J)              max |J^T S J - S| for a 6x6 Jacobian J
#   jacobian_of_map(M, z0)           Jacobian of a polynomial map at point z0
#   jacobian_of_function(f, z0)      Jacobian of ANY generic function f(z)
#                                    (e.g. a tracking routine) at z0, computed
#                                    exactly by feeding z0 + Z (TaylorN) in -
#                                    automatic differentiation, no finite steps.
#
# WHY
#   A truncated Taylor map is symplectic only up to its truncation order
#   (error grows like |z|^order at amplitude), whereas the factorised Lie map
#   applied element by element must be symplectic to machine precision at any
#   amplitude.  These functions make that difference visible (scripts/02).
#
# WHAT YOU SEE
#   Nothing printed here.
# =============================================================================
module Symplectic

using TaylorSeries, LinearAlgebra
using ..LieAlgebra

export symplectic_error, jacobian_of_map, jacobian_of_function

"""
    symplectic_error(J) -> Float64

max |Jᵀ S J − S| over all entries: 0 for an exactly symplectic matrix. Applied
to Jacobians of maps or tracking functions to measure the symplecticity defect.
"""
symplectic_error(J::AbstractMatrix) = maximum(abs, J' * SYMPLECTIC_S * J - SYMPLECTIC_S)

"""
    jacobian_of_map(M, z0) -> Matrix{Float64}

Jacobian ∂M_i/∂z_j of a polynomial map evaluated at the point `z0` (exact
differentiation of the polynomial).
"""
jacobian_of_map(M, z0::AbstractVector) =
    [evaluate(differentiate(M[i], j), collect(float.(z0))) for i in 1:6, j in 1:6]

"""
    jacobian_of_function(f, z0) -> Matrix{Float64}

Exact Jacobian of a generic tracking function `f` (any code that works on
`TaylorN`, e.g. `z -> track(factors, z)` or `exact_track_turn`) at `z0`.

# Algorithm
Calls `f(z0 .+ Z)` once in polynomial arithmetic; the linear coefficients of the
output are the Jacobian. This is automatic differentiation; no finite differences.
"""
function jacobian_of_function(f, z0::AbstractVector)
    zt = [z0[i] + Z[i] for i in 1:6]
    out = f(zt)
    return [constant_term(differentiate(out[i], j)) for i in 1:6, j in 1:6]
end

end # module
