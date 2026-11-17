# =============================================================================
# integrators/Integrators.jl - general ODE integrators (pure mathematics)
# =============================================================================
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# One file per method. The module only holds the common driver integrate.
# =============================================================================
@doc raw"""
    Integrators

Selectable ODE integrators, generic over `Float64` and [DA](@ref Bell.DiffAlgebra.DA):
RK4, Dormand-Prince 5(4), DOP853, Gragg-Bulirsch-Stoer, implicit midpoint,
Yoshida composition of the implicit midpoint, and splitting methods
(Leapfrog, Yoshida 4/6/8) for user-supplied exact sub-flows.  Entry point integrate.
"""
module Integrators

using LinearAlgebra
using ..DiffAlgebra: DA, coeffnorm

export ODEMethod, RK4, DormandPrince5, DormandPrince853, GraggBulirschStoer,
       ImplicitMidpoint, ImplicitYoshida, implicit_midpoint, implicit_yoshida_weights,
       integrate, integrate_fixed, errnorm, method_order,
       SplitFlows, SplittingMethod, Leapfrog, Yoshida, ForestRuth, splitting_weights

"Supertype of the ODE methods of [`integrate`](@ref)."
abstract type ODEMethod end

@doc raw"""
    errnorm(x)

Magnitude used by the step control: ``|x|`` for numbers; for a DA the coefficient norm
``\max_\alpha |x_\alpha|`` (all Taylor coefficients are controlled, not only the constant part).
Pass `norm = x -> abs(constant_part(x))` to [`integrate`](@ref) to control only the
reference trajectory.
"""
errnorm(x::Number) = abs(x)
errnorm(x::DA) = coeffnorm(x)

_axpy(z, h, k) = [z[i] + h * k[i] for i in eachindex(z)]
_lin(z, h, cs, ks) = [z[i] + h * sum(cs[j] * ks[j][i] for j in eachindex(cs) if !iszero(cs[j]); init = zero(z[i])) for i in eachindex(z)]

_err(::ODEMethod, e, z, zn, tol, rtol, norm) =
    maximum(norm(e[i]) / (tol + rtol * max(norm(z[i]), norm(zn[i]))) for i in eachindex(z))

include("rk4.jl")
include("dormand_prince.jl")
include("gragg_bulirsch_stoer.jl")
include("implicit_midpoint.jl")
include("yoshida.jl")

@doc raw"""
    integrate(f, z0, s0, s1; method = DormandPrince5(), tol = 1e-12, rtol = tol, h = nothing,
              nsteps = nothing, h0 = nothing, hmax = abs(s1 - s0), maxsteps = 10^6,
              norm = errnorm) -> NamedTuple (z, nsteps, nrejected, nfev)

Integrate ``dz/ds = f(z, s)`` from ``s_0`` to ``s_1``; `z0` is a vector of any number type
(`Float64`, [`DA`](@ref), ...), `f(z, s)` returns a vector of the same length.

- Fixed step: give `h` or `nsteps` (or use a method with no error estimate:
  `RK4()`, `ImplicitMidpoint()`, `ImplicitYoshida()`).
- Adaptive (DormandPrince5, GraggBulirschStoer): local error per step
  ``\mathrm{err} = \max_i \|e_i\| / (\mathrm{tol} + \mathrm{rtol}\,\|z_i\|)`` with
  ``\|\cdot\|`` = `norm` (for DA: coefficient norm, see [`errnorm`](@ref)); the step is
  accepted if ``\mathrm{err} \le 1`` and the next step is
  ``h \cdot \min(5, \max(0.2,\ 0.9\,\mathrm{err}^{-1/q}))``, ``q`` = order of the estimate
  (4 for DP5, ``2k-1`` for GBS).
"""
function integrate(f, z0::AbstractVector, s0::Real, s1::Real; method::ODEMethod = DormandPrince5(),
                   tol::Real = 1e-12, rtol::Real = tol, h = nothing, nsteps = nothing, h0 = nothing,
                   hmax = abs(s1 - s0), maxsteps::Integer = 10^6, norm = errnorm)
    z = collect(z0)
    s1 == s0 && return (z = z, nsteps = 0, nrejected = 0, nfev = 0)
    fixed = method isa RK4 || method isa ImplicitMidpoint || method isa ImplicitYoshida ||
            h !== nothing || nsteps !== nothing
    if fixed
        n = nsteps !== nothing ? Int(nsteps) : h !== nothing ? max(1, ceil(Int, abs(s1 - s0) / h)) :
            error("$(typeof(method)) is fixed-step: give h or nsteps")
        return integrate_fixed(f, z, s0, s1, n; method)
    end
    dir = sign(s1 - s0); s = float(s0)
    hc = h0 === nothing ? dir * min(abs(s1 - s0), method isa GraggBulirschStoer ? abs(s1 - s0) / 4 : 1e-2 * abs(s1 - s0)) : dir * abs(h0)
    q = _errexp(method)
    nacc = 0; nrej = 0; nf = 0
    while dir * (s1 - s) > 0
        nacc + nrej > maxsteps && error("integrate: maxsteps exceeded at s = $s")
        abs(hc) > abs(s1 - s) && (hc = s1 - s)
        zn, e, k = _step(method, f, z, s, hc); nf += k
        err = _err(method, e, z, zn, tol, rtol, norm)
        if err <= 1
            s = (abs(s1 - s - hc) < 1e-14 * abs(s1 - s0)) ? float(s1) : s + hc; z = zn; nacc += 1
        else
            nrej += 1
        end
        fac = err == 0 ? 5.0 : min(5.0, max(0.2, 0.9 * err^(-1 / (q + 1))))
        hc = dir * min(abs(hmax), abs(hc) * fac)
    end
    return (z = z, nsteps = nacc, nrejected = nrej, nfev = nf)
end

"""
    integrate_fixed(f, z0, s0, s1, n; method = RK4()) -> NamedTuple (z, nsteps, nrejected, nfev)

`n` equal steps of `method` (any [`ODEMethod`](@ref); embedded estimates are ignored).
"""
function integrate_fixed(f, z0::AbstractVector, s0::Real, s1::Real, n::Integer; method::ODEMethod = RK4())
    z = collect(z0); h = (s1 - s0) / n; nf = 0
    for i in 0:n-1
        z, _, k = _step(method, f, z, s0 + i * h, h); nf += k
    end
    return (z = z, nsteps = n, nrejected = 0, nfev = nf)
end

end # module
