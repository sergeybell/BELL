# Implicit midpoint rule, and Yoshida composition of that rule.
# Distinct from Yoshida splitting of two exact sub-flows (yoshida.jl).

@doc raw"""
    implicit_midpoint(f, z, s, hs; maxiter = 40, tol = 1e-16) -> z'

One implicit midpoint step ``z' = z + h\,f((z+z')/2,\, s + h/2)``.
Fixed-point iteration; for real numbers it stops at `|Δ| ≤ tol`.
Symplectic for any, including s-dependent, Hamiltonian vector field; order 2.
"""
function implicit_midpoint(f, z, s, hs; maxiter::Integer = 40, tol::Real = 1e-16)
    zn = z .+ hs .* f(z, s + hs / 2)
    for _ in 1:maxiter
        znew = z .+ hs .* f((z .+ zn) ./ 2, s + hs / 2)
        if eltype(znew) <: Real
            done = maximum(abs.(znew .- zn)) <= tol
            zn = znew
            done && break
        else
            zn = znew
        end
    end
    return zn
end

@doc raw"""
    ImplicitMidpoint()

Implicit midpoint rule, order 2. Symplectic for a Hamiltonian vector field.
Fixed step only.
"""
struct ImplicitMidpoint <: ODEMethod end
method_order(::ImplicitMidpoint) = 2

function _step(::ImplicitMidpoint, f, z, s, h)
    return (implicit_midpoint(f, z, s, h), nothing, 1)
end

@doc raw"""
    ImplicitYoshida(order = 4)

Yoshida triple-jump composition of the implicit midpoint rule. Order 4 is
```math
S_4(h) = S_2(w_1 h)\,S_2(w_0 h)\,S_2(w_1 h),\qquad
w_1 = 1/(2-2^{1/3}),\quad w_0 = 1-2w_1.
```
This is the method used by exact tracking of an unexpanded Hamiltonian.
It is not [`Yoshida`](@ref), which composes two exact sub-flows.
"""
struct ImplicitYoshida <: ODEMethod
    order::Int
    function ImplicitYoshida(o::Integer = 4)
        (o >= 2 && iseven(o)) || throw(ArgumentError("ImplicitYoshida order must be even >= 2"))
        new(o)
    end
end
method_order(m::ImplicitYoshida) = m.order

"Composition weights of [`ImplicitYoshida`](@ref) (sum = 1)."
function implicit_yoshida_weights(order::Integer)
    w = [1.0]
    for n in 1:(order ÷ 2 - 1)
        w1 = 1 / (2 - 2^(1 / (2n + 1))); w0 = 1 - 2w1
        w = vcat(w1 .* w, w0 .* w, w1 .* w)
    end
    w
end
implicit_yoshida_weights(m::ImplicitYoshida) = implicit_yoshida_weights(m.order)

function _step(m::ImplicitYoshida, f, z, s, h)
    s_cur = s
    zn = z
    w = implicit_yoshida_weights(m)
    for wi in w
        hs = wi * h
        zn = implicit_midpoint(f, zn, s_cur, hs)
        s_cur += hs
    end
    return (zn, nothing, length(w))
end
