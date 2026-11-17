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
# is the call site.
# =============================================================================
module Symplectic

using TaylorSeries, LinearAlgebra
using ..LieAlgebra
using ..DiffAlgebra: DA, da_symplectic_error

export symplectic_error, jacobian_of_map, jacobian_of_function

"""
    symplectic_error(object) -> Float64
    symplectic_error(f, z0) -> Float64

Symplecticity defect ``\\max|J^T S J - S|``. The object is the thing being
checked: a numeric Jacobian, a polynomial map, a DA map, or a tracking step
`f` at point `z0`. Zero (up to round-off) for an exactly symplectic linear map.
A truncated Taylor map is symplectic only through its retained degree.
"""
symplectic_error(J::AbstractMatrix) = maximum(abs, J' * SYMPLECTIC_S * J - SYMPLECTIC_S)

function symplectic_error(M::AbstractVector{<:TaylorN})
    J = [differentiate(M[i], j) for i in 1:6, j in 1:6]
    err = 0.0
    for a in 1:6, b in 1:6
        s = zero(M[1])
        for i in 1:6, j in 1:6
            SYMPLECTIC_S[i, j] == 0 && continue
            s += SYMPLECTIC_S[i, j] * (J[i, a] * J[j, b])
        end
        err = max(err, maxabs(truncate_to(s - SYMPLECTIC_S[a, b], MAP_ORDER - 1)))
    end
    return err
end

symplectic_error(M::AbstractVector{<:DA}; kwargs...) = da_symplectic_error(M; kwargs...)

symplectic_error(f, z0::AbstractVector) = symplectic_error(jacobian_of_function(f, z0))

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
