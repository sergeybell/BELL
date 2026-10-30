# =============================================================================
# lie_series.jl - Poisson brackets, Lie operators and Lie-series flows on DA
#   (included into module DiffAlgebra; pure mathematics, no accelerator concepts)
# =============================================================================
#
# WHAT
#   canonical_pairs(nv)          default canonical pairs (1,2), (3,4), ...
#   poisson_bracket(f, g)        [f, g] = sum_k df/dq_k dg/dp_k - df/dp_k dg/dq_k
#   hamiltonian_vector_field(H)  F = S grad H, i.e. dz/ds = [z, H] = F(z)
#   lie_series(H, g, t)          exp(t :H:) g = sum t^n/n! :H:^n g      (:H: g = [H, g])
#   lie_flow(H, L)               the flow exp(-L :H:) z of H over L, all components,
#                                summed to machine precision, with scaling-and-squaring
#   change_order(a, sp)          re-embed a DA in a space of another order (truncate / pad)
#   da_symplectic_error(M)          max coefficient of J^T S J - S, J = dM/dz (as DA)
# =============================================================================

@doc raw"""
    canonical_pairs(nv) -> Vector{Tuple{Int,Int}}

Default canonical pairs ``(q_k, p_k) = (2k-1, 2k)``, ``k = 1..\lfloor nv/2\rfloor``; e.g. for
``z = (x, p_x, y, p_y, l, \delta)``: ``(x,p_x), (y,p_y), (l,\delta)``. Variables not in a pair are
parameters (their flow is the identity).
"""
canonical_pairs(nv::Integer) = [(2k - 1, 2k) for k in 1:nv÷2]

@doc raw"""
    poisson_bracket(f, g; pairs = canonical_pairs(nv)) -> DA

```math
[f, g] = \sum_k \left(\frac{\partial f}{\partial q_k}\frac{\partial g}{\partial p_k}
       - \frac{\partial f}{\partial p_k}\frac{\partial g}{\partial q_k}\right).
```
The Lie operator is ``{:}f{:}\,g = [f, g]``. In ``{}_ND_v`` each derivative loses the top order,
so ``[f,g]`` is exact through order ``N-1`` in general, and through ``N`` when ``f`` (or ``g``)
has no terms of degree ``< 2`` (see [`lie_flow`](@ref)).
"""
function poisson_bracket(f::DA, g::DA; pairs = canonical_pairs(f.sp.nv))
    s = zero(f)
    for (q, p) in pairs
        s += deriv(f, q) * deriv(g, p) - deriv(f, p) * deriv(g, q)
    end
    s
end

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

# D g = F . grad g  (only the nonzero components of F)
_lie_D(F, idx, g) = (s = zero(g); for i in idx; s += F[i] * deriv(g, i); end; s)

@doc raw"""
    lie_series(H, g, t; pairs, tol = eps(), nmax = 200) -> DA

```math
e^{t{:}H{:}}\,g = \sum_{n\ge 0}\frac{t^n}{n!}{:}H{:}^n g, \qquad {:}H{:}\,g = [H, g],
```
terms ``u_n = \tfrac{t}{n}{:}H{:}\,u_{n-1}`` summed until two consecutive terms satisfy
``\|u_n\| \le \text{tol}\,\|\text{sum}\|`` (coefficient norm) or a term vanishes identically.
With ``t = -L`` this is the value of ``g`` after the flow of ``H`` over the length ``L``.
No scaling is applied here; see [`lie_flow`](@ref).
"""
function lie_series(H::DA, g::DA, t::Real; pairs = canonical_pairs(H.sp.nv), tol = eps(), nmax = 200)
    F = hamiltonian_vector_field(H; pairs = pairs)
    r, _ = _lie_sum(F, _nzcomp(F), g, -t, tol, nmax)
    r
end
_nzcomp(F) = [i for i in eachindex(F) if !iszero(F[i].c)]
function _lie_sum(F, idx, g, L, tol, nmax)
    acc = g; u = g; small = 0
    for n in 1:nmax
        u = (L / n) * _lie_D(F, idx, u)
        acc += u
        nu = coeffnorm(u)
        nu == 0 && return acc, n
        small = nu <= tol * coeffnorm(acc) ? small + 1 : 0
        small >= 2 && return acc, n
    end
    return acc, -1                      # not converged
end

@doc raw"""
    lie_flow(H, L; pairs, tol = eps(), nmax = 60, nsquare = :auto) -> Vector{DA}

Transfer map of the s-independent Hamiltonian ``H`` over the length ``L``:
```math
\mathcal M = e^{-L{:}H{:}}\,z = \sum_{n\ge0}\frac{(-L)^n}{n!}{:}H{:}^n z
           = \sum_{n\ge0}\frac{L^n}{n!}(F\cdot\nabla)^n z, \qquad F = S\nabla H,
```
for every component of ``z`` (parameters are left unchanged), summed to machine precision.

**Exactness.** If ``H`` has no linear terms (the reference orbit is a solution, ``F(0) = 0``),
``F\cdot\nabla`` does not lower the degree, so every term is exact in ``{}_ND_v`` except the
top degree ``N``, whose coefficients would need the degree-``N+1`` part of ``H``. **Build ``H``
in a space of order ``N+1`` and keep the result through order ``N``** ([`change_order`](@ref));
then the map is exact through ``N`` up to roundoff. A warning is issued if ``F(0) \ne 0``.

**Scaling and squaring.** The series converges like ``e^{L\rho}``, ``\rho = \|\partial F/\partial z(0)\|``.
With `nsquare = :auto`, ``k = \max(0, \lceil\log_2(L\rho/0.5)\rceil)`` and, if the series has
not converged in `nmax` terms, ``k`` is increased; the map for ``L/2^k`` is then squared
``k`` times, ``\mathcal M_L = (\mathcal M_{L/2^k})^{\circ 2^k}`` (exact: the flow is a one-parameter group).
Pass an integer to force ``k``.
"""
function lie_flow(H::DA, L::Real; pairs = canonical_pairs(H.sp.nv), tol = eps(), nmax::Integer = 60,
                  nsquare = :auto)
    sp = H.sp; F = hamiltonian_vector_field(H; pairs = pairs); idx = _nzcomp(F)
    maximum(abs(F[i].c[1]) for i in eachindex(F)) > 1e-13 * max(1.0, maximum(coeffnorm, F)) &&
        @warn "lie_flow: H has linear terms (F(0) ≠ 0); the top order of the map is not exact"
    A = [F[i].c[1+j] for i in 1:sp.nv, j in 1:sp.nv]
    rho = opnorm(A, Inf)
    k = nsquare === :auto ? max(0, ceil(Int, log2(max(abs(L) * rho, 1e-300) / 0.5))) : Int(nsquare)
    X = identity_map(sp)
    while true
        h = L / 2^k; M = similar(X); ok = true
        for i in 1:sp.nv
            M[i], n = _lie_sum(F, idx, X[i], h, tol, nmax)
            n < 0 && (ok = false; break)
        end
        if ok || nsquare !== :auto || k > 40
            ok || @warn "lie_flow: series not converged (nsquare = $k)"
            for _ in 1:k
                M = compose(M, M)
            end
            return M
        end
        k += 1
    end
end

@doc raw"""
    change_order(a, sp) -> DA

The DA `a` re-embedded in the space `sp` (same number of variables, any order): monomials of
degree ``> \text{order}(sp)`` are dropped, missing ones are zero.
"""
function change_order(a::DA, sp::DASpace)
    sp.nv == a.sp.nv || throw(ArgumentError("change_order: different number of variables"))
    c = zeros(eltype(a.c), nmonomials(sp))
    for m in 1:nmonomials(sp)
        sp.deg[m] > a.sp.no && break
        c[m] = coefficient(a, view(sp.exps, :, m))
    end
    DA(sp, c)
end

@doc raw"""
    da_symplectic_error(M; pairs = canonical_pairs(nv)) -> Float64

``\max \big|\,[J^T S J - S]_\alpha\big|`` over all matrix elements and all coefficients
``\alpha`` (orders ``0 \dots N-1``), where ``J_{ij} = \partial M_i/\partial z_j`` as DA and ``S`` is
the symplectic form of `pairs` (parameter rows/columns are excluded). Zero (up to roundoff)
for a map that is symplectic through order ``N`` — this tests all orders, not only the linear part.
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