# =============================================================================
# dual.jl - forward-mode dual numbers and the phase-space gradient
#   (included into module DiffAlgebra; pure mathematics, no accelerator types)
# =============================================================================
#
# WHAT
#   Dual(value, partial)     value + ε·partial, ε² = 0
#   phase_space_gradient     ∂f/∂z_j at the current z, one probe per coordinate
#   symplectic_gradient      S ∇f; S is an argument (default: block pairing)
#   map_vector_field         S ∇f on a DA coordinate jet, ∇ taken by deriv
#   vector_field_on_identity the same jet on the identity of a DASpace
#
# The dual partial differentiates the expression. deriv differentiates a jet
# that is already truncated at degree N, so its degree-N part is zero. The two
# agree through order N-1.
# =============================================================================

@doc raw"""
    Dual(value, partial)
    Dual(value)

Forward-mode number ``\mathrm{value} + \varepsilon\,\mathrm{partial}`` with ``\varepsilon^2 = 0``.
`value` and `partial` share one numeric type (a float, a `TaylorN`, or a DA).
`Dual(value)` has a zero partial.

Arithmetic is the sum, product, and quotient rules, integer powers, and the
chain rule for `sqrt`, `inv`, `exp`, `log`, `sin`, and `cos`.
"""
struct Dual{T<:Number} <: Number
    value::T
    partial::T
end

function Dual(value::Number, partial::Number)
    T = promote_type(typeof(value), typeof(partial))
    return Dual{T}(convert(T, value), convert(T, partial))
end
Dual(value::Number) = Dual(value, zero(value))

Base.zero(a::Dual) = Dual(zero(a.value), zero(a.partial))
Base.one(a::Dual) = Dual(one(a.value), zero(a.partial))
Base.:(==)(a::Dual, b::Dual) = a.value == b.value && a.partial == b.partial
function Base.isapprox(a::Dual, b::Dual; kwargs...)
    return isapprox(a.value, b.value; kwargs...) && isapprox(a.partial, b.partial; kwargs...)
end
Base.show(io::IO, a::Dual) = print(io, "Dual(", a.value, ", ", a.partial, ")")

Base.:+(a::Dual, b::Dual) = Dual(a.value + b.value, a.partial + b.partial)
Base.:+(a::Dual, x::Number) = Dual(a.value + x, a.partial)
Base.:+(x::Number, a::Dual) = Dual(x + a.value, a.partial)
Base.:-(a::Dual) = Dual(-a.value, -a.partial)
Base.:-(a::Dual, b::Dual) = Dual(a.value - b.value, a.partial - b.partial)
Base.:-(a::Dual, x::Number) = Dual(a.value - x, a.partial)
Base.:-(x::Number, a::Dual) = Dual(x - a.value, -a.partial)
function Base.:*(a::Dual, b::Dual)
    return Dual(a.value * b.value, a.partial * b.value + a.value * b.partial)
end
Base.:*(a::Dual, x::Number) = Dual(a.value * x, a.partial * x)
Base.:*(x::Number, a::Dual) = Dual(x * a.value, x * a.partial)
function Base.:/(a::Dual, b::Dual)
    v = a.value / b.value
    d = (a.partial * b.value - a.value * b.partial) / (b.value * b.value)
    return Dual(v, d)
end
Base.:/(a::Dual, x::Number) = Dual(a.value / x, a.partial / x)
Base.:/(x::Number, a::Dual) = x * inv(a)

function Base.inv(a::Dual)
    v = inv(a.value)
    return Dual(v, -a.partial * v * v)
end

function Base.:^(a::Dual, n::Integer)
    n < 0 && return inv(a)^(-n)
    r = one(a)
    p = a
    m = n
    while m > 0
        isodd(m) && (r *= p)
        m >>= 1
        m > 0 && (p *= p)
    end
    return r
end

function Base.sqrt(a::Dual)
    v = sqrt(a.value)
    return Dual(v, a.partial / (2 * v))
end

function Base.exp(a::Dual)
    v = exp(a.value)
    return Dual(v, a.partial * v)
end

Base.log(a::Dual) = Dual(log(a.value), a.partial / a.value)
Base.sin(a::Dual) = Dual(sin(a.value), a.partial * cos(a.value))
Base.cos(a::Dual) = Dual(cos(a.value), -a.partial * sin(a.value))

@doc raw"""
    phase_space_gradient(f, z::AbstractVector) -> Vector

``\partial f/\partial z_j`` at the current point `z`, one forward-mode probe per
coordinate. `f` returns a scalar. `z` may hold numbers, `TaylorN`, or DA; each
partial has that same scalar type.

On a truncated series the dual rule differentiates the expression.
[`deriv`](@ref) differentiates the jet after truncation, so its degree-``N``
part is zero. The two agree through order ``N-1``, and at order ``N`` when `f`
is a polynomial of degree at most ``N``. The same rule holds for a shifted
coordinate jet (`da_vars` with a nonzero `at`): the constant part of each
partial is ``\partial f/\partial z_j`` at that point.
"""
function phase_space_gradient(f, z::AbstractVector)
    isempty(z) && throw(ArgumentError("z is empty"))
    return [_partial(f, z, j) for j in eachindex(z)]
end

function _probe(z, j)
    return [Dual(z[i], i == j ? one(z[i]) : zero(z[i])) for i in eachindex(z)]
end

_partial(f, z, j) = _derivative_part(f(_probe(z, j)))
_derivative_part(y::Dual) = y.partial
_derivative_part(y::Number) = zero(y)
_derivative_part(y) = throw(ArgumentError("f(z) must be a scalar, got $(typeof(y))"))

"""Block pairing with blocks `[0 1; -1 0]`. `n` must be positive and even."""
function _canonical_pairing(n::Integer)
    n > 0 && iseven(n) ||
        throw(ArgumentError("default pairing needs a positive even length, got $n"))
    S = zeros(n, n)
    for k in 1:(n ÷ 2)
        S[2k - 1, 2k] = 1.0
        S[2k, 2k - 1] = -1.0
    end
    return S
end

function _apply_pairing(S, g::AbstractVector)
    n = length(g)
    size(S, 1) == n && size(S, 2) == n ||
        throw(ArgumentError("S must be $(n)×$(n), got $(size(S))"))
    return [_pairing_row(S, g, i) for i in 1:n]
end

function _pairing_row(S, g, i)
    s = S[i, 1] * g[1]
    for j in 2:length(g)
        s += S[i, j] * g[j]
    end
    return s
end

@doc raw"""
    symplectic_gradient(f, z::AbstractVector; S = pairing(length(z))) -> Vector

``S \nabla f`` at the current `z`. The gradient is [`phase_space_gradient`](@ref).

`S` is an `n×n` matrix, ``n = `` `length(z)`. The default is the block-diagonal
pairing whose blocks are ``\begin{bmatrix} 0 & 1 \\ -1 & 0 \end{bmatrix}``;
`n` must then be even. Pass `S` to use another matrix.
"""
function symplectic_gradient(f, z::AbstractVector; S = _canonical_pairing(length(z)))
    return _apply_pairing(S, phase_space_gradient(f, z))
end

function _check_coordinate_jet(z::AbstractVector{<:DA})
    isempty(z) && throw(ArgumentError("z is empty"))
    sp = z[1].sp
    length(z) == sp.nv || throw(ArgumentError("expected $(sp.nv) components, got $(length(z))"))
    for (i, zi) in enumerate(z)
        zi.sp === sp || throw(ArgumentError("components must share one DASpace"))
        for m in 2:length(zi.c)
            expected = (sp.no >= 1 && m == 1 + i) ? one(eltype(zi.c)) : zero(eltype(zi.c))
            zi.c[m] == expected || throw(ArgumentError(
                "z must be a coordinate jet (da_vars or identity_map); " *
                "component $i is not a coordinate",
            ))
        end
    end
    return sp
end

function _as_da(y::DA, sp::DASpace)
    y.sp == sp || throw(ArgumentError("f(z) must be a DA in the space of z"))
    return y
end
_as_da(y::Real, sp::DASpace) = da_const(sp, y)
_as_da(y, ::DASpace) = throw(ArgumentError("f(z) must be a scalar DA, got $(typeof(y))"))

@doc raw"""
    map_vector_field(f, z::AbstractVector{<:DA}; S = pairing(length(z))) -> Vector{DA}

``S \nabla f`` on a coordinate jet `z` ([`identity_map`](@ref) or [`da_vars`](@ref),
including a nonzero expansion point). `f(z)` is the jet of the scalar and
``\nabla`` is [`deriv`](@ref), so degree ``N`` of each component is zero.

`z` is that coordinate jet, not an arbitrary map. [`phase_space_gradient`](@ref)
on the same jet differentiates the expression; the two agree through order ``N-1``.
`S` is the same argument as in [`symplectic_gradient`](@ref).
"""
function map_vector_field(f, z::AbstractVector{<:DA}; S = _canonical_pairing(length(z)))
    sp = _check_coordinate_jet(z)
    h = _as_da(f(z), sp)
    g = [deriv(h, j) for j in 1:length(z)]
    return _apply_pairing(S, g)
end

@doc raw"""
    vector_field_on_identity(f, sp::DASpace; S = pairing(sp.nv)) -> Vector{DA}

[`map_vector_field`](@ref) of `f` on [`identity_map`](@ref)`(sp)`.
"""
function vector_field_on_identity(f, sp::DASpace; S = _canonical_pairing(sp.nv))
    return map_vector_field(f, identity_map(sp); S = S)
end
