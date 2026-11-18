# =============================================================================
# lie_series.jl - canonical pairs and the DA symplecticity defect
#   (included into module DiffAlgebra; pure mathematics, no accelerator concepts)
# =============================================================================
#
# WHAT
#   canonical_pairs(nv)          internal default pairs (1,2), (3,4), ...
#   da_symplectic_error(M)       max coefficient of J^T S J - S, J = dM/dz (as DA)
#
# The series flow (poisson_bracket, lie_series, lie_flow, vector_field_flow,
# change_order) and hamiltonian_vector_field had no callers. Magnet maps
# integrate S ∇H. The field equations remain for comparison. The Lie path
# integrates TaylorN Hamiltonians in LieAlgebra.
# =============================================================================

# Default pairs (q_k, p_k) = (2k-1, 2k). Not part of the public API.
canonical_pairs(nv::Integer) = [(2k - 1, 2k) for k in 1:(nv ÷ 2)]

@doc raw"""
    da_symplectic_error(M; pairs) -> Float64

``\max \big|\,[J^T S J - S]_\alpha\big|`` over all matrix elements and all coefficients
``\alpha`` (orders ``0 \dots N-1``), where ``J_{ij} = \partial M_i/\partial z_j`` as DA and ``S`` is
the symplectic form of `pairs`. The default pairs are ``(2k-1, 2k)``. Parameter
rows and columns are excluded. Zero (up to roundoff) for a map that is symplectic
through order ``N`` — this tests all orders, not only the linear part.

This is the DA method of `symplectic_error` (`Bell.Symplectic`), which is the single check
for a numeric Jacobian, a `TaylorN` map, a DA map, or a tracking step.
"""
function da_symplectic_error(M::AbstractVector{<:DA}; pairs = canonical_pairs(M[1].sp.nv))
    idx = vcat([[q, p] for (q, p) in pairs]...); n = length(idx)
    S = zeros(n, n); for k in 1:length(pairs); S[2k-1, 2k] = 1; S[2k, 2k-1] = -1; end
    J = [deriv(M[idx[a]], idx[b]) for a in 1:n, b in 1:n]
    err = 0.0
    for a in 1:n, b in 1:n
        s = zero(M[1])
        for i in 1:n, j in 1:n
            S[i, j] == 0 && continue
            s += S[i, j] * (J[i, a] * J[j, b])
        end
        # degree N of J^T S J would need degree N+1 of M: compare through N-1 only
        err = max(err, coeffnorm(truncate_da(s - S[a, b], M[1].sp.no - 1)))
    end
    err
end
