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

symplectic_error(J::AbstractMatrix) = maximum(abs, J' * SYMPLECTIC_S * J - SYMPLECTIC_S)

jacobian_of_map(M, z0::AbstractVector) =
    [evaluate(differentiate(M[i], j), collect(float.(z0))) for i in 1:6, j in 1:6]

function jacobian_of_function(f, z0::AbstractVector)
    zt = [z0[i] + Z[i] for i in 1:6]
    out = f(zt)
    return [constant_term(differentiate(out[i], j)) for i in 1:6, j in 1:6]
end

end # module
