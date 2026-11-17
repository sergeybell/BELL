# =============================================================================
# lie_series.jl - canonical pairs and the DA symplecticity defect
#   (included into module DiffAlgebra; pure mathematics, no accelerator concepts)
# =============================================================================
#
# WHAT
#   canonical_pairs(nv)          default canonical pairs (1,2), (3,4), ...
#   hamiltonian_vector_field(H)  F = S grad H, i.e. dz/ds = [z, H] = F(z)
#   da_symplectic_error(M)       max coefficient of J^T S J - S, J = dM/dz (as DA)
#
# The series flow (poisson_bracket, lie_series, lie_flow, vector_field_flow,
# change_order) had no callers in the package. The field / DA maps integrate
# element_rhs; the Lie path integrates TaylorN Hamiltonians in LieAlgebra.
# =============================================================================

@doc raw"""
    canonical_pairs(nv) -> Vector{Tuple{Int,Int}}

Default canonical pairs ``(q_k, p_k) = (2k-1, 2k)``, ``k = 1..\lfloor nv/2\rfloor``; e.g. for
``z = (x, p_x, y, p_y, l, \delta)``: ``(x,p_x), (y,p_y), (l,\delta)``. Variables not in a pair are
parameters (their flow is the identity).
"""
canonical_pairs(nv::Integer) = [(2k - 1, 2k) for k in 1:nv÷2]

@doc raw"""
    hamiltonian_vector_field(H; pairs) -> Vector{DA}

``F = S\,\nabla H``: ``F_{q_k} = \partial H/\partial p_k``, ``F_{p_k} = -\partial H/\partial q_k``,
zero for parameters. Hamilton's equations are ``z' = [z, H] = F(z)`` and, for any ``g``,
``-{:}H{:}\,g = [g, H] = F\cdot\nabla g``.
"""
function hamiltonian_vector_field(H::DA; pairs = canonical_pairs(H.sp.nv))
    F = [zero(H) for _ in 1:H.sp.nv]
    for (q, p) in pairs
        F[q] = deriv(H, p); F[p] = -deriv(H, q)
    end
    F
end

@doc raw"""
    da_symplectic_error(M; pairs = canonical_pairs(nv)) -> Float64

``\max \big|\,[J^T S J - S]_\alpha\big|`` over all matrix elements and all coefficients
``\alpha`` (orders ``0 \dots N-1``), where ``J_{ij} = \partial M_i/\partial z_j`` as DA and ``S`` is
the symplectic form of `pairs` (parameter rows/columns are excluded). Zero (up to roundoff)
for a map that is symplectic through order ``N`` — this tests all orders, not only the linear part.

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
