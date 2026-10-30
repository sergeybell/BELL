# =============================================================================
# da/DiffAlgebra.jl - truncated power series algebra (differential algebra, DA)
# =============================================================================
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# WHAT
#   A pure-mathematics module: the algebra  _N D_v  of polynomials in v
#   variables truncated at total degree N (Berz, "Modern Map Methods in
#   Particle Beam Physics", ch. 2).  Both v and N are runtime parameters of a
#   DASpace.  Nothing here knows about accelerators.
#
#   DASpace(nv, no)      monomial table (graded order), index lookup
#   DA{T}                element of the algebra (dense coefficient vector)
#   + - * / ^            ring operations, truncated at degree N
#   exp log sqrt sin cos tan atan sinh cosh tanh inv, real powers
#                        by the Taylor series of f about the constant part
#   deriv, antideriv     d/dx_i and int dx_i
#   a(x), compose        evaluation at a point (numbers or DA), composition
#   inverse_map          inverse of a map with invertible linear part
#   coefficient, constant_part, linear_part, truncate_da
# =============================================================================
@doc raw"""
    DiffAlgebra

Truncated power-series algebra `{}_nD_v`: `v` variables, order `n` (both run-time),
arithmetic, intrinsics, `\partial`/`\partial^{-1}`, evaluation, composition and map inversion.
See the Developer's Guide page *Differential algebra*.
"""
module DiffAlgebra

using LinearAlgebra, Printf

export DASpace, DA, da_var, da_vars, da_const, coefficient, constant_part, linear_part,
       truncate_da, deriv, antideriv, compose, inverse_map, identity_map, nmonomials,
       monomial_exponents, coeffnorm, order, nvars,
       canonical_pairs, poisson_bracket, hamiltonian_vector_field, lie_series, lie_flow, vector_field_flow,
       change_order, da_symplectic_error

# ----------------------------------------------------------------------------- space
@doc raw"""
    DASpace(nv, no)

Truncated power-series algebra ``{}_N D_v`` with ``v`` = `nv` variables and maximal total
degree ``N`` = `no`. An element is
```math
a(x) = \sum_{|\alpha| \le N} a_\alpha\, x^\alpha, \qquad x^\alpha = x_1^{\alpha_1}\cdots x_v^{\alpha_v},
```
with ``\binom{v+N}{N}`` coefficients. Products are truncated: monomials with
``|\alpha| > N`` are dropped, so ``{}_N D_v`` is the quotient of the polynomial ring by
the ideal of degree ``> N``; it is the algebra of ``N``-jets of smooth functions
(the DA of Berz). Monomials are stored in graded order: degree 0, then 1, ..., ``N``,
lexicographic inside a degree.

Index lookup uses the packed key ``k(\alpha) = \sum_i \alpha_i (N+1)^{i-1}``, which is
additive (``k(\alpha+\beta) = k(\alpha)+k(\beta)`` when ``|\alpha+\beta| \le N``); keys go to
a dense table when ``(N+1)^v \le 2^{24}``, otherwise to a `Dict`.
"""
struct DASpace
    nv::Int
    no::Int
    exps::Matrix{Int}          # nv x n
    deg::Vector{Int}
    key::Vector{Int}
    pw::Vector{Int}            # (no+1)^(i-1)
    dense::Vector{Int32}       # key+1 -> index (or empty)
    dict::Dict{Int,Int}
    degstart::Vector{Int}      # degstart[d+1] = first index of degree d; degstart[no+2] = n+1
    parent::Vector{Int}        # monomial m = parent[m] * x_{pvar[m]}
    pvar::Vector{Int}
end

function _exps_of_degree(nv, d)
    nv == 1 && return [[d]]
    out = Vector{Vector{Int}}()
    for a in d:-1:0, r in _exps_of_degree(nv - 1, d - a)
        push!(out, [a; r])
    end
    out
end

function DASpace(nv::Integer, no::Integer)
    nv >= 1 || throw(ArgumentError("nv must be >= 1"))
    no >= 0 || throw(ArgumentError("no must be >= 0"))
    E = Vector{Vector{Int}}(); ds = Int[]
    for d in 0:no
        push!(ds, length(E) + 1); append!(E, _exps_of_degree(nv, d))
    end
    push!(ds, length(E) + 1)
    n = length(E); exps = reduce(hcat, E); exps = reshape(exps, nv, n)
    pw = [(no + 1)^(i - 1) for i in 1:nv]
    key = [sum(exps[:, m] .* pw) for m in 1:n]
    deg = [sum(exps[:, m]) for m in 1:n]
    big = (no + 1.0)^nv > 2.0^24
    dense = big ? Int32[] : zeros(Int32, (no + 1)^nv)
    dict = Dict{Int,Int}()
    for m in 1:n
        big ? (dict[key[m]] = m) : (dense[key[m]+1] = m)
    end
    look(k) = big ? dict[k] : Int(dense[k+1])
    parent = zeros(Int, n); pvar = zeros(Int, n)
    for m in 2:n
        j = findfirst(>(0), exps[:, m]); pvar[m] = j; parent[m] = look(key[m] - pw[j])
    end
    DASpace(nv, no, exps, deg, key, pw, dense, dict, ds, parent, pvar)
end

@inline _index(sp::DASpace, k::Int) = isempty(sp.dense) ? sp.dict[k] : Int(@inbounds sp.dense[k+1])

"Number of monomials ``\\binom{v+N}{N}`` of the space."
nmonomials(sp::DASpace) = size(sp.exps, 2)
"Exponent vector ``\\alpha`` of monomial number `m`."
monomial_exponents(sp::DASpace, m::Integer) = sp.exps[:, m]
Base.show(io::IO, sp::DASpace) = print(io, "DASpace(nv = ", sp.nv, ", no = ", sp.no, ", ", nmonomials(sp), " monomials)")
Base.:(==)(a::DASpace, b::DASpace) = a.nv == b.nv && a.no == b.no

# ----------------------------------------------------------------------------- element
@doc raw"""
    DA{T}(sp, c)          DA(sp, c::Number)

Element of ``{}_N D_v``: coefficient vector `c` (length [`nmonomials`](@ref), graded
order) in the space `sp`. `DA <: Number`, so generic numerical code (e.g. the ODE
integrators) runs on it unchanged. `DA(sp, x::Number)` is the constant ``x``;
[`da_var`](@ref) gives the coordinate functions ``x_i``.

The constant part ``a_0 = a(0)`` is the "real" value; the rest is the nilpotent part
``\bar a = a - a_0`` (``\bar a^{N+1} = 0``).
"""
struct DA{T<:Number} <: Number
    sp::DASpace
    c::Vector{T}
end
DA(sp::DASpace, x::Number) = (c = zeros(typeof(float(x)), nmonomials(sp)); c[1] = x; DA(sp, c))
"Constant DA ``x`` in `sp`."
da_const(sp::DASpace, x::Number) = DA(sp, x)
@doc raw"""
    da_var(sp, i; at = 0.0) -> DA

Coordinate function ``x_i`` (plus the constant `at`): ``a_0 = `` `at`, ``\partial a/\partial x_i = 1``.
"""
function da_var(sp::DASpace, i::Integer; at::Real = 0.0)
    1 <= i <= sp.nv || throw(ArgumentError("variable $i out of 1:$(sp.nv)"))
    c = zeros(Float64, nmonomials(sp)); c[1] = at
    sp.no >= 1 && (c[1+i] = 1.0)   # degree-1 monomials are x_1..x_v in this order
    DA(sp, c)
end
"All coordinate functions ``(x_1,\\dots,x_v)``, optionally about the point `at`."
da_vars(sp::DASpace; at = zeros(sp.nv)) = [da_var(sp, i; at = at[i]) for i in 1:sp.nv]
"Identity map ``x \\mapsto x`` of `sp` (vector of [`da_var`](@ref))."
identity_map(sp::DASpace) = da_vars(sp)
"Number of variables of the space of `a`."
nvars(a::DA) = a.sp.nv
"Truncation order ``N`` of the space of `a`."
order(a::DA) = a.sp.no

_check(a::DA, b::DA) = (a.sp === b.sp || a.sp == b.sp || throw(ArgumentError("DA of different spaces")); a.sp)

Base.zero(a::DA{T}) where {T} = DA(a.sp, zeros(T, length(a.c)))
Base.one(a::DA{T}) where {T} = (c = zeros(T, length(a.c)); c[1] = one(T); DA(a.sp, c))
Base.copy(a::DA) = DA(a.sp, copy(a.c))
Base.float(a::DA) = a
Base.convert(::Type{DA{T}}, a::DA) where {T} = DA(a.sp, convert(Vector{T}, a.c))
Base.promote_rule(::Type{DA{T}}, ::Type{S}) where {T,S<:Real} = DA{promote_type(T, S)}
Base.promote_rule(::Type{DA{T}}, ::Type{DA{S}}) where {T,S} = DA{promote_type(T, S)}
Base.:(==)(a::DA, b::DA) = a.sp == b.sp && a.c == b.c
Base.isapprox(a::DA, b::DA; atol = 0.0, rtol = sqrt(eps())) = isapprox(a.c, b.c; atol, rtol)

"Coefficient ``a_\\alpha`` for the exponent vector `alpha`."
function coefficient(a::DA, alpha::AbstractVector{<:Integer})
    length(alpha) == a.sp.nv || throw(ArgumentError("exponent length"))
    sum(alpha) > a.sp.no && return zero(eltype(a.c))
    a.c[_index(a.sp, sum(alpha .* a.sp.pw))]
end
"Constant part ``a_0 = a(0)``."
constant_part(a::DA) = a.c[1]
constant_part(x::Number) = x
"Coefficient norm ``\\max_\\alpha |a_\\alpha|`` (used for step control of DA integration)."
coeffnorm(a::DA) = maximum(abs, a.c)
coeffnorm(x::Number) = abs(x)
@doc raw"""
    linear_part(M::AbstractVector{<:DA}) -> Matrix

Jacobian at the origin, ``L_{ij} = \partial M_i/\partial x_j(0)``.
"""
linear_part(M::AbstractVector{<:DA}) = [M[i].c[1+j] for i in eachindex(M), j in 1:M[1].sp.nv]
"Drop all monomials of degree > `n` (``n \\le N``)."
function truncate_da(a::DA, n::Integer)
    c = copy(a.c); c[a.sp.degstart[min(n, a.sp.no)+2]:end] .= 0; DA(a.sp, c)
end

# ----------------------------------------------------------------------------- arithmetic
Base.:+(a::DA, b::DA) = (_check(a, b); DA(a.sp, a.c .+ b.c))
Base.:-(a::DA, b::DA) = (_check(a, b); DA(a.sp, a.c .- b.c))
Base.:-(a::DA) = DA(a.sp, -a.c)
Base.:+(a::DA, x::Number) = (c = copy(a.c) .+ zero(x); c[1] += x; DA(a.sp, c))
Base.:+(x::Number, a::DA) = a + x
Base.:-(a::DA, x::Number) = a + (-x)
Base.:-(x::Number, a::DA) = (-a) + x
Base.:*(a::DA, x::Number) = DA(a.sp, a.c .* x)
Base.:*(x::Number, a::DA) = a * x
Base.:/(a::DA, x::Number) = DA(a.sp, a.c ./ x)
Base.:/(x::Number, a::DA) = x * inv(a)
Base.:/(a::DA, b::DA) = a * inv(b)

@doc raw"""
    a * b

Truncated Cauchy product ``(ab)_\gamma = \sum_{\alpha+\beta=\gamma,\ |\gamma|\le N} a_\alpha b_\beta``.
Only pairs with ``|\alpha| + |\beta| \le N`` are visited (graded storage).
"""
function Base.:*(a::DA, b::DA)
    sp = _check(a, b)
    T = promote_type(eltype(a.c), eltype(b.c))
    c = zeros(T, length(a.c)); ds = sp.degstart; key = sp.key; N = sp.no
    @inbounds for i in 1:length(a.c)
        ai = a.c[i]; iszero(ai) && continue
        jmax = ds[N-sp.deg[i]+2] - 1; ki = key[i]
        for j in 1:jmax
            bj = b.c[j]; iszero(bj) && continue
            c[_index(sp, ki + key[j])] += ai * bj
        end
    end
    DA(sp, c)
end

function Base.:^(a::DA, n::Integer)
    n < 0 && return inv(a)^(-n)
    r = one(a); p = a
    while n > 0
        isodd(n) && (r = r * p); n >>= 1; n > 0 && (p = p * p)
    end
    r
end

# ----------------------------------------------------------------------------- intrinsics
@doc raw"""
    _series(a, f) -> DA

``f(a) = \sum_{k=0}^{N} f_k\, \bar a^{\,k}``, ``\bar a = a - a_0``, with ``f_k = f^{(k)}(a_0)/k!``
(Horner). Exact in ``{}_N D_v`` because ``\bar a^{N+1} = 0``. This is how every intrinsic
function is extended to DA.
"""
function _series(a::DA, f::AbstractVector)
    h = a - a.c[1]; N = a.sp.no
    r = DA(a.sp, f[N+1] + zero(eltype(a.c)))
    for k in N-1:-1:0
        r = r * h + f[k+1]
    end
    r
end

Base.exp(a::DA) = (e = exp(a.c[1]); _series(a, [e / factorial(big(k)) for k in 0:a.sp.no] .|> Float64))
function Base.log(a::DA)
    a0 = a.c[1]; a0 > 0 || throw(DomainError(a0, "log of DA with non-positive constant part"))
    _series(a, [k == 0 ? log(a0) : (-1)^(k + 1) / (k * a0^k) for k in 0:a.sp.no])
end
@doc raw"""
    a^p  (real p)

``(a_0+\bar a)^p = a_0^p \sum_k \binom{p}{k} (\bar a/a_0)^k``, requires ``a_0 > 0``
(``a_0 \ne 0`` for integer `p`, which uses repeated squaring instead).
"""
function Base.:^(a::DA, p::Real)
    isinteger(p) && return a^Int(p)
    a0 = a.c[1]; a0 > 0 || throw(DomainError(a0, "non-integer power of DA needs a positive constant part"))
    f = zeros(a.sp.no + 1); b = 1.0
    for k in 0:a.sp.no
        f[k+1] = b * a0^(p - k); b *= (p - k) / (k + 1)
    end
    _series(a, f)
end
Base.sqrt(a::DA) = a^0.5
Base.cbrt(a::DA) = a^(1 / 3)
@doc raw"""
    inv(a)

``1/a = \sum_k (-1)^k \bar a^k / a_0^{k+1}``; requires ``a_0 \ne 0`` (otherwise ``a`` is a
zero divisor and has no inverse in ``{}_N D_v``).
"""
function Base.inv(a::DA)
    a0 = a.c[1]; iszero(a0) && throw(DomainError(a0, "DA with zero constant part is not invertible"))
    _series(a, [(-1)^k / a0^(k + 1) for k in 0:a.sp.no])
end
Base.sin(a::DA) = (a0 = a.c[1]; _series(a, [sin(a0 + k * pi / 2) / factorial(big(k)) for k in 0:a.sp.no] .|> Float64))
Base.cos(a::DA) = (a0 = a.c[1]; _series(a, [cos(a0 + k * pi / 2) / factorial(big(k)) for k in 0:a.sp.no] .|> Float64))
Base.tan(a::DA) = sin(a) / cos(a)
Base.sinh(a::DA) = (e = exp(a); (e - inv(e)) / 2)
Base.cosh(a::DA) = (e = exp(a); (e + inv(e)) / 2)
Base.tanh(a::DA) = sinh(a) / cosh(a)
@doc raw"""
    atan(a)

``\arctan(a_0+\bar a) = \arctan a_0 + \arctan t``, ``t = \bar a/(1 + a_0 a)`` (addition
theorem), and ``\arctan t = \sum_k (-1)^k t^{2k+1}/(2k+1)`` for the nilpotent ``t``.
"""
function Base.atan(a::DA)
    a0 = a.c[1]; t = (a - a0) / (1 + a0 * a)
    f = zeros(a.sp.no + 1); for k in 0:2:a.sp.no-1; f[k+2] = (-1)^(k ÷ 2) / (k + 1); end
    _series(t, f) + atan(a0)
end
Base.abs(a::DA) = a.c[1] < 0 ? -a : a
Base.isless(a::DA, b::DA) = isless(a.c[1], b.c[1])
Base.isless(a::DA, x::Real) = isless(a.c[1], x)
Base.isless(x::Real, a::DA) = isless(x, a.c[1])

# ----------------------------------------------------------------------------- calculus
@doc raw"""
    deriv(a, i) -> DA

``\partial a/\partial x_i``. The result is exact to degree ``N-1`` (the degree-``N`` part of
``a`` only determines degree ``N-1`` of the derivative); its degree-``N`` coefficients are 0.
"""
function deriv(a::DA, i::Integer)
    sp = a.sp; c = zeros(eltype(a.c), length(a.c))
    @inbounds for m in 2:length(a.c)
        e = sp.exps[i, m]; (e == 0 || iszero(a.c[m])) && continue
        c[_index(sp, sp.key[m] - sp.pw[i])] += e * a.c[m]
    end
    DA(sp, c)
end
@doc raw"""
    antideriv(a, i) -> DA

``\int_0^{x_i} a\, dx_i`` (zero on the hyperplane ``x_i = 0``); terms that would reach degree
``N+1`` are dropped, so ``\partial_i \int a\,dx_i = a`` up to degree ``N-1``.
"""
function antideriv(a::DA, i::Integer)
    sp = a.sp; c = zeros(eltype(a.c), length(a.c))
    @inbounds for m in 1:length(a.c)
        (sp.deg[m] >= sp.no || iszero(a.c[m])) && continue
        c[_index(sp, sp.key[m] + sp.pw[i])] += a.c[m] / (sp.exps[i, m] + 1)
    end
    DA(sp, c)
end

# ----------------------------------------------------------------------------- evaluation
"Values of all monomials ``x^\\alpha`` at `x` (numbers or DA), graded order, one product each."
function _monomial_values(sp::DASpace, x::AbstractVector)
    length(x) == sp.nv || throw(ArgumentError("need $(sp.nv) arguments"))
    T = promote_type(map(typeof, x)...)
    v = Vector{T}(undef, nmonomials(sp)); v[1] = one(x[1])
    for m in 2:length(v)
        v[m] = v[sp.parent[m]] * x[sp.pvar[m]]
    end
    v
end
@doc raw"""
    a(x)

Evaluate ``a(x) = \sum_\alpha a_\alpha x^\alpha`` at `x` (vector of numbers, or of DA of another
space: this is composition).
"""
(a::DA)(x::AbstractVector) = (v = _monomial_values(a.sp, x); sum(a.c[m] * v[m] for m in eachindex(v) if !iszero(a.c[m]); init = zero(v[1]) * 0))
@doc raw"""
    compose(M, N) -> Vector{DA}

Composition ``(M \circ N)(x) = M(N(x))`` of DA maps: `M` has `nv` components in space ``S``,
`N` is a vector of ``v_S`` DA in any space. If ``N`` has a nonzero constant part the result is
the exact re-expansion only when ``M`` is a polynomial (it always is in ``{}_N D_v``); the usual
case is ``N(0) = 0``.
"""
function compose(M::AbstractVector{<:DA}, N::AbstractVector)
    N isa AbstractVector{DA{Float64}} && eltype(M) == DA{Float64} && return _compose_da(M, N)
    v = _monomial_values(M[1].sp, N)
    [sum(Mi.c[m] * v[m] for m in eachindex(v) if !iszero(Mi.c[m]); init = zero(v[1]) * 0) for Mi in M]
end
# c = a * b (coefficient arrays of space sp), loops over the nonzero entries bnz of b only
function _mul_into!(c::Vector{Float64}, sp::DASpace, a::Vector{Float64}, b::Vector{Float64}, bnz)
    fill!(c, 0.0); ds = sp.degstart; key = sp.key; N = sp.no
    @inbounds for i in eachindex(a)
        ai = a[i]; iszero(ai) && continue
        jmax = ds[N-sp.deg[i]+2] - 1; ki = key[i]
        for j in bnz
            j > jmax && break
            c[_index(sp, ki + key[j])] += ai * b[j]
        end
    end
    c
end
# M o N for DA arguments: monomial table of N built in place (one product per monomial, each
# against the sparse N[pvar]); if N(0) = 0, monomials of degree > order(N) vanish and are skipped;
# the result is accumulated in place.  Same values as the generic path.
function _compose_da(M::AbstractVector{DA{Float64}}, N::AbstractVector{DA{Float64}})
    spM = M[1].sp; spN = N[1].sp; length(N) == spM.nv || throw(ArgumentError("need $(spM.nv) arguments"))
    len = nmonomials(spN); nm = nmonomials(spM)
    used = falses(nm); for Mi in M, m in 1:nm; iszero(Mi.c[m]) || (used[m] = true); end
    n0 = all(iszero(n.c[1]) for n in N)
    # a monomial value is needed if it is used or is the parent of a needed one
    need = copy(used); for m in nm:-1:2; need[m] && (need[spM.parent[m]] = true); end
    Nnz = [findall(!iszero, n.c) for n in N]
    V = Vector{Vector{Float64}}(undef, nm); V[1] = zeros(len); V[1][1] = 1.0
    for m in 2:nm
        (!need[m] || (n0 && spM.deg[m] > spN.no)) && (V[m] = Float64[]; continue)
        par = V[spM.parent[m]]
        isempty(par) && (V[m] = Float64[]; continue)
        V[m] = _mul_into!(zeros(len), spN, par, N[spM.pvar[m]].c, Nnz[spM.pvar[m]])
    end
    out = [zeros(len) for _ in M]
    for (o, Mi) in zip(out, M), m in 1:nm
        c = Mi.c[m]; (iszero(c) || isempty(V[m])) && continue
        v = V[m]; @inbounds @simd for i in eachindex(v); o[i] += c * v[i]; end
    end
    [DA(spN, o) for o in out]
end
@doc raw"""
    inverse_map(M) -> Vector{DA}

Inverse of a map ``M: x \mapsto M(x)``, ``M(0) = 0`` (an origin-preserving map; subtract the
constant part first), with invertible linear part ``L``. Writing ``M = L + \mathcal N``
(``\mathcal N`` = terms of degree ``\ge 2``), the inverse is the fixed point of
```math
M^{-1} = L^{-1}\big(\mathrm{id} - \mathcal N \circ M^{-1}\big),
```
and each iteration fixes one more order, so ``N`` iterations give ``M^{-1}`` exactly in
``{}_N D_v`` (Berz, Modern Map Methods, sec. 2.3.1).
"""
function inverse_map(M::AbstractVector{<:DA})
    sp = M[1].sp; length(M) == sp.nv || throw(ArgumentError("map must have nv components"))
    all(iszero(Mi.c[1]) for Mi in M) || throw(ArgumentError("inverse_map needs M(0) = 0"))
    L = linear_part(M); Li = inv(L); I_ = identity_map(sp)
    NL = [Mi - sum(L[i, j] * I_[j] for j in 1:sp.nv) for (i, Mi) in enumerate(M)]
    X = [sum(Li[i, j] * I_[j] for j in 1:sp.nv) for i in 1:sp.nv]
    for _ in 1:sp.no
        R = I_ .- compose(NL, X)
        X = [sum(Li[i, j] * R[j] for j in 1:sp.nv) for i in 1:sp.nv]
    end
    X
end

# ----------------------------------------------------------------------------- printing
function Base.show(io::IO, a::DA)
    nz = count(!iszero, a.c)
    print(io, "DA(nv=", a.sp.nv, ", no=", a.sp.no, ", ", nz, " nonzero; a0 = ", a.c[1], ")")
end
function Base.show(io::IO, ::MIME"text/plain", a::DA)
    println(io, "DA in ", a.sp, ":")
    println(io, "    I   COEFFICIENT              ORDER  EXPONENTS")
    k = 0
    for m in eachindex(a.c)
        iszero(a.c[m]) && continue
        k += 1
        @printf(io, "%5d  %23.16e  %3d    %s\n", k, a.c[m], a.sp.deg[m], join(a.sp.exps[:, m], " "))
    end
    k == 0 && println(io, "  (zero)")
end

include("lie_series.jl")

end # module
