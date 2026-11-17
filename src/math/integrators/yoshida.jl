# Splitting methods for two exact sub-flows: Leapfrog and Yoshida composition.
# This is not the implicit-midpoint Yoshida of implicit_midpoint.jl.

@doc raw"""
    SplitFlows(flowA, flowB)

A split vector field ``f = f_A + f_B`` (e.g. ``H = T(p) + V(q)``) given by its two **exact**
sub-flows: `flowA(z, h)` and `flowB(z, h)` return ``\varphi^A_h(z)`` and ``\varphi^B_h(z)``
(drift and kick for a separable Hamiltonian). Pass it to [`integrate`](@ref) with a
[`SplittingMethod`](@ref). If both sub-flows are symplectic, so is every composition.
"""
struct SplitFlows{FA,FB}
    flowA::FA
    flowB::FB
end

"Supertype of composition (splitting) methods."
abstract type SplittingMethod <: ODEMethod end

@doc raw"""
    Leapfrog()

Strang splitting / Störmer-Verlet, order 2:
``\Phi_h = \varphi^A_{h/2}\circ\varphi^B_h\circ\varphi^A_{h/2}``.
"""
struct Leapfrog <: SplittingMethod end

@doc raw"""
    Yoshida(order)   order = 2, 4, 6, 8, ...

Triple-jump composition (Yoshida 1990, Suzuki): from the symmetric method ``S_{2n}``,
```math
S_{2n+2}(h) = S_{2n}(w_1h)\,S_{2n}(w_0h)\,S_{2n}(w_1h),\qquad
w_1 = \frac{1}{2-2^{1/(2n+1)}},\quad w_0 = 1-2w_1,
```
starting from [`Leapfrog`](@ref) (``S_2``). `Yoshida(4)` is the Forest-Ruth integrator
(3 kicks); order ``2n`` needs ``3^{n-1}`` leapfrog stages.
This composition acts on two exact sub-flows. The Yoshida composition of the
implicit midpoint rule is [`ImplicitYoshida`](@ref).
"""
struct Yoshida <: SplittingMethod
    order::Int
    function Yoshida(o::Integer)
        (o >= 2 && iseven(o)) || throw(ArgumentError("Yoshida order must be even >= 2"))
        new(o)
    end
end
"Forest-Ruth 4th-order integrator (= `Yoshida(4)`)."
ForestRuth() = Yoshida(4)

method_order(::Leapfrog) = 2
method_order(m::Yoshida) = m.order

"Leapfrog sub-step weights of a composition method (sum = 1)."
function splitting_weights(m::Yoshida)
    w = [1.0]
    for n in 1:(m.order ÷ 2 - 1)
        w1 = 1 / (2 - 2^(1 / (2n + 1))); w0 = 1 - 2w1
        w = vcat(w1 .* w, w0 .* w, w1 .* w)
    end
    w
end
splitting_weights(::Leapfrog) = [1.0]

function _split_step(p::SplitFlows, w, z, h)
    for wi in w
        z = p.flowA(z, wi * h / 2); z = p.flowB(z, wi * h); z = p.flowA(z, wi * h / 2)
    end
    z
end

@doc raw"""
    integrate(p::SplitFlows, z0, s0, s1; method = Yoshida(4), h = nothing, nsteps = nothing)

Fixed-step composition integration of a split problem (see [`SplitFlows`](@ref)); `method` is
[`Leapfrog`](@ref) or [`Yoshida`](@ref). Adjacent half drifts are not merged (simplicity over
the saving of one A-flow per stage).
"""
function integrate(p::SplitFlows, z0::AbstractVector, s0::Real, s1::Real; method::SplittingMethod = Yoshida(4),
                   h = nothing, nsteps = nothing, kwargs...)
    n = nsteps !== nothing ? Int(nsteps) : h !== nothing ? max(1, ceil(Int, abs(s1 - s0) / h)) :
        error("splitting methods are fixed-step: give h or nsteps")
    hh = (s1 - s0) / n; w = splitting_weights(method); z = collect(z0)
    for _ in 1:n
        z = _split_step(p, w, z, hh)
    end
    return (z = z, nsteps = n, nrejected = 0, nfev = 2 * length(w) * n)
end
