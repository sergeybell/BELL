# =============================================================================
# integrators/Integrators.jl - general ODE integrators (pure mathematics)
# =============================================================================
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# WHAT
#   integrate(f, z0, s0, s1; method, tol)   dz/ds = f(z, s)
#   Methods (all explicit, generic over the number type: Float64, DA, ...):
#     RK4()                    classical Runge-Kutta, order 4, fixed step
#     DormandPrince5()         embedded 5(4) pair (FSAL), adaptive or fixed
#     DormandPrince853()       order 8 with 5th/3rd-order error estimators (DOP853)
#     Leapfrog(), Yoshida(2n)  symplectic compositions for SplitFlows (exact sub-flows)
#     GraggBulirschStoer(k)    extrapolated modified midpoint, order 2k,
#                              embedded error estimate, adaptive or fixed
#   The symplectic integrators of Bell (symplectic.jl, exact_integrator.jl)
#   are not touched; this module is independent of the accelerator code.
# =============================================================================
@doc raw"""
    Integrators

Selectable ODE integrators, generic over `Float64` and [DA](@ref Bell.DiffAlgebra.DA):
RK4, Dormand-Prince 5(4), DOP853, Gragg-Bulirsch-Stoer, and splitting methods
(Leapfrog, Yoshida 4/6/8) for user-supplied exact sub-flows.  Entry point integrate.
"""
module Integrators

using LinearAlgebra
using ..DiffAlgebra: DA, coeffnorm

export ODEMethod, RK4, DormandPrince5, DormandPrince853, GraggBulirschStoer, integrate, integrate_fixed, errnorm,
       method_order, SplitFlows, SplittingMethod, Leapfrog, Yoshida, ForestRuth, splitting_weights

"Supertype of the ODE methods of [`integrate`](@ref)."
abstract type ODEMethod end

@doc raw"""
    RK4()

Classical 4-stage Runge-Kutta method, order 4. Fixed step only.
```math
z_{n+1} = z_n + \tfrac{h}{6}(k_1 + 2k_2 + 2k_3 + k_4).
```
"""
struct RK4 <: ODEMethod end

@doc raw"""
    DormandPrince5()

Dormand-Prince embedded Runge-Kutta pair 5(4) (Dormand & Prince 1980), 7 stages with the
first-same-as-last property; the 5th-order solution is propagated (local extrapolation), the
difference to the 4th-order one is the local error estimate.
"""
struct DormandPrince5 <: ODEMethod end

@doc raw"""
    GraggBulirschStoer(k = 8)

Gragg-Bulirsch-Stoer extrapolation: modified midpoint rule with ``n_j = 2j`` substeps
(``j = 1..k``), whose error has an expansion in ``h^2`` only (Gragg), extrapolated to
``h \to 0`` by the Aitken-Neville scheme in ``h^2``:
```math
T_{j,m+1} = T_{j,m} + \frac{T_{j,m} - T_{j-1,m}}{(n_j/n_{j-m})^2 - 1}.
```
``T_{k,k}`` has order ``2k``; ``\|T_{k,k} - T_{k,k-1}\|`` is the embedded error estimate.
An explicit method of arbitrary even order (``k = 8`` gives order 16).
"""
struct GraggBulirschStoer <: ODEMethod
    k::Int
end
GraggBulirschStoer() = GraggBulirschStoer(8)

"Classical order of the method (GBS: ``2k``)."
method_order(::RK4) = 4
method_order(::DormandPrince5) = 5
method_order(m::GraggBulirschStoer) = 2m.k

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

# ----------------------------------------------------------------------------- one steps
function _step(::RK4, f, z, s, h)
    k1 = f(z, s); k2 = f(_axpy(z, h / 2, k1), s + h / 2)
    k3 = f(_axpy(z, h / 2, k2), s + h / 2); k4 = f(_axpy(z, h, k3), s + h)
    return ([z[i] + h / 6 * (k1[i] + 2k2[i] + 2k3[i] + k4[i]) for i in eachindex(z)], nothing, 4)
end

const DP_C = (0.0, 1 / 5, 3 / 10, 4 / 5, 8 / 9, 1.0, 1.0)
const DP_A = ((), (1 / 5,), (3 / 40, 9 / 40), (44 / 45, -56 / 15, 32 / 9),
              (19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729),
              (9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656),
              (35 / 384, 0.0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84))
const DP_B5 = (35 / 384, 0.0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84, 0.0)
const DP_B4 = (5179 / 57600, 0.0, 7571 / 16695, 393 / 640, -92097 / 339200, 187 / 2100, 1 / 40)

function _step(::DormandPrince5, f, z, s, h)
    ks = Any[f(z, s)]
    for st in 2:7
        push!(ks, f(_lin(z, h, DP_A[st], ks[1:st-1]), s + DP_C[st] * h))
    end
    z5 = _lin(z, h, DP_B5, ks)
    e = [h * sum((DP_B5[j] - DP_B4[j]) * ks[j][i] for j in 1:7) for i in eachindex(z)]
    return (z5, e, 7)
end

function _midpoint(f, z, s, H, n)
    h = H / n
    zm = z; z1 = _axpy(z, h, f(z, s))
    for m in 1:n-1
        zn = _axpy(zm, 2h, f(z1, s + m * h)); zm = z1; z1 = zn
    end
    [(zm[i] + z1[i] + h * f(z1, s + H)[i]) / 2 for i in eachindex(z)]  # Gragg smoothing
end

function _step(m::GraggBulirschStoer, f, z, s, H)
    k = m.k; ns = [2j for j in 1:k]
    T = Vector{Any}(undef, k); nf = 0
    for j in 1:k
        Tj = Vector{Any}(undef, j)
        Tj[1] = _midpoint(f, z, s, H, ns[j]); nf += ns[j] + 1
        for q in 1:j-1
            r = (ns[j] / ns[j-q])^2 - 1
            Tj[q+1] = [Tj[q][i] + (Tj[q][i] - T[j-1][q][i]) / r for i in eachindex(z)]
        end
        T[j] = Tj
    end
    zk = T[k][k]
    e = k > 1 ? [zk[i] - T[k][k-1][i] for i in eachindex(z)] : nothing
    return ([x for x in zk], e, nf)
end

# ----------------------------------------------------------------------------- DOP853
@doc raw"""
    DormandPrince853()

Explicit Runge-Kutta method of order 8 by Dormand & Prince (Hairer, Nørsett & Wanner,
*Solving ODE I*, sec. II.10, code DOP853): 12 stages, error estimator combining embedded 5th-
and 3rd-order solutions,
```math
\mathrm{err} = \frac{\mathrm{err}_5^2}{\sqrt{\mathrm{err}_5^2 + 0.01\,\mathrm{err}_3^2}},
```
step factor ``\mathrm{err}^{-1/8}``. Coefficients as published (also used by SciPy's `DOP853`).
"""
struct DormandPrince853 <: ODEMethod end
method_order(::DormandPrince853) = 8

const D8_C = [0.0, 0.526001519587677318785587544488e-01, 0.789002279381515978178381316732e-01,
    0.118350341907227396726757197510, 0.281649658092772603273242802490, 0.333333333333333333333333333333,
    0.25, 0.307692307692307692307692307692, 0.651282051282051282051282051282, 0.6,
    0.857142857142857142857142857142, 1.0, 1.0]
const D8_A = let A = zeros(13, 13)
    s(i, j, v) = (A[i+1, j+1] = v)
    s(1,0,5.26001519587677318785587544488e-2)
    s(2,0,1.97250569845378994544595329183e-2); s(2,1,5.91751709536136983633785987549e-2)
    s(3,0,2.95875854768068491816892993775e-2); s(3,2,8.87627564304205475450678981324e-2)
    s(4,0,2.41365134159266685502369798665e-1); s(4,2,-8.84549479328286085344864962717e-1); s(4,3,9.24834003261792003115737966543e-1)
    s(5,0,3.7037037037037037037037037037e-2); s(5,3,1.70828608729473871279604482173e-1); s(5,4,1.25467687566822425016691814123e-1)
    s(6,0,3.7109375e-2); s(6,3,1.70252211019544039314978060272e-1); s(6,4,6.02165389804559606850219397283e-2); s(6,5,-1.7578125e-2)
    s(7,0,3.70920001185047927108779319836e-2); s(7,3,1.70383925712239993810214054705e-1); s(7,4,1.07262030446373284651809199168e-1)
    s(7,5,-1.53194377486244017527936158236e-2); s(7,6,8.27378916381402288758473766002e-3)
    s(8,0,6.24110958716075717114429577812e-1); s(8,3,-3.36089262944694129406857109825); s(8,4,-8.68219346841726006818189891453e-1)
    s(8,5,2.75920996994467083049415600797e1); s(8,6,2.01540675504778934086186788979e1); s(8,7,-4.34898841810699588477366255144e1)
    s(9,0,4.77662536438264365890433908527e-1); s(9,3,-2.48811461997166764192642586468); s(9,4,-5.90290826836842996371446475743e-1)
    s(9,5,2.12300514481811942347288949897e1); s(9,6,1.52792336328824235832596922938e1); s(9,7,-3.32882109689848629194453265587e1)
    s(9,8,-2.03312017085086261358222928593e-2)
    s(10,0,-9.3714243008598732571704021658e-1); s(10,3,5.18637242884406370830023853209); s(10,4,1.09143734899672957818500254654)
    s(10,5,-8.14978701074692612513997267357); s(10,6,-1.85200656599969598641566180701e1); s(10,7,2.27394870993505042818970056734e1)
    s(10,8,2.49360555267965238987089396762); s(10,9,-3.0467644718982195003823669022)
    s(11,0,2.27331014751653820792359768449); s(11,3,-1.05344954667372501984066689879e1); s(11,4,-2.00087205822486249909675718444)
    s(11,5,-1.79589318631187989172765950534e1); s(11,6,2.79488845294199600508499808837e1); s(11,7,-2.85899827713502369474065508674)
    s(11,8,-8.87285693353062954433549289258); s(11,9,1.23605671757943030647266201528e1); s(11,10,6.43392746015763530355970484046e-1)
    s(12,0,5.42937341165687622380535766363e-2); s(12,5,4.45031289275240888144113950566); s(12,6,1.89151789931450038304281599044)
    s(12,7,-5.8012039600105847814672114227); s(12,8,3.1116436695781989440891606237e-1); s(12,9,-1.52160949662516078556178806805e-1)
    s(12,10,2.01365400804030348374776537501e-1); s(12,11,4.47106157277725905176885569043e-2)
    A
end
const D8_B = D8_A[13, 1:12]
const D8_E3 = let e = [D8_B; 0.0]
    e[1] -= 0.244094488188976377952755905512; e[9] -= 0.733846688281611857341361741547
    e[12] -= 0.220588235294117647058823529412e-1; e
end
const D8_E5 = let e = zeros(13)
    e[1] = 0.1312004499419488073250102996e-1; e[6] = -0.1225156446376204440720569753e+1
    e[7] = -0.4957589496572501915214079952; e[8] = 0.1664377182454986536961530415e+1
    e[9] = -0.3503288487499736816886487290; e[10] = 0.3341791187130174790297318841
    e[11] = 0.8192320648511571246570742613e-1; e[12] = -0.2235530786388629525884427845e-1; e
end

function _step(::DormandPrince853, f, z, s, h)
    ks = Any[f(z, s)]
    for st in 2:12
        push!(ks, f(_lin(z, h, view(D8_A, st, 1:st-1), ks[1:st-1]), s + D8_C[st] * h))
    end
    zn = _lin(z, h, D8_B, ks)
    push!(ks, f(zn, s + h))
    e5 = [h * sum(D8_E5[j] * ks[j][i] for j in 1:13 if !iszero(D8_E5[j])) for i in eachindex(z)]
    e3 = [h * sum(D8_E3[j] * ks[j][i] for j in 1:13 if !iszero(D8_E3[j])) for i in eachindex(z)]
    return (zn, (e5, e3), 13)
end

_err(::ODEMethod, e, z, zn, tol, rtol, norm) =
    maximum(norm(e[i]) / (tol + rtol * max(norm(z[i]), norm(zn[i]))) for i in eachindex(z))
function _err(::DormandPrince853, e, z, zn, tol, rtol, norm)
    e5 = _err(RK4(), e[1], z, zn, tol, rtol, norm); e3 = _err(RK4(), e[2], z, zn, tol, rtol, norm)
    d = sqrt(e5^2 + 0.01 * e3^2); d == 0 ? 0.0 : e5^2 / d
end
_errexp(m::GraggBulirschStoer) = 2m.k - 1
_errexp(::DormandPrince5) = 4
_errexp(::DormandPrince853) = 7
# ----------------------------------------------------------------------------- drivers
@doc raw"""
    integrate(f, z0, s0, s1; method = DormandPrince5(), tol = 1e-12, rtol = tol, h = nothing,
              nsteps = nothing, h0 = nothing, hmax = abs(s1 - s0), maxsteps = 10^6,
              norm = errnorm) -> NamedTuple (z, nsteps, nrejected, nfev)

Integrate ``dz/ds = f(z, s)`` from ``s_0`` to ``s_1``; `z0` is a vector of any number type
(`Float64`, [`DA`](@ref), ...), `f(z, s)` returns a vector of the same length.

- Fixed step: give `h` or `nsteps` (or use `RK4()`, which has no error estimate).
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
    if method isa RK4 || h !== nothing || nsteps !== nothing
        n = nsteps !== nothing ? Int(nsteps) : h !== nothing ? max(1, ceil(Int, abs(s1 - s0) / h)) :
            error("RK4 is fixed-step: give h or nsteps")
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

# ----------------------------------------------------------------------------- splitting
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

end # module