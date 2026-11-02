# =============================================================================
# fringe.jl - fringe-field (edge) models of magnets, chosen at call time
# =============================================================================
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# The lattice stays physical: an SBend carries L, angle, e1, e2 and optional
# gap / FINT. No element carries a fringe flag. The model is chosen when the
# ring is built. FringeModel() is NoEdge everywhere and leaves the ring as is.
#
# HardEdge is a thin lens at each end. The body Hamiltonian is not changed.
# Linear (order 1): dpx = h tan(e) x, dpy = -h tan(psi) y. Higher order keeps
# a Lie generator K from the step-field B_s.
#
# SoftEdge is a different layout. Half-width W is read from the profile (Enge:
# first z where F is below 1e-10 and 1-F is below 1e-10, step gap/20, cap 40
# gaps). The edge map is Body(-W) o Soft o Drift(-W): a slice of solid body at
# depth W is removed and replaced by the s-dependent Hamiltonian. W depends on
# the profile and the gap, not on magnet strength (the profile is normalised;
# the strength sits in h0). Subtracting the body on depth W is the account of
# field leaking into the magnet; if notes later say otherwise, the formulas
# stay until that is reviewed.
#
# enge_default is an unverified default coefficient set stored in this
# repository. There is no citation for it; it is not a standard.
#
# RK4 of the soft-edge field and Yoshida weights of splitting come from
# Integrators. Simpson of a profile integral is quadrature, not an ODE method.
# =============================================================================
module Fringe

using TaylorSeries, LinearAlgebra
using ..Coordinates, ..Elements, ..Structure, ..LieAlgebra, ..Hamiltonians, ..Linear,
      ..Symplectic, ..Tracking, ..MapReadings, ..Integrators
import ..Tracking: build_factors, ring_lie_factor

export EdgeModel, NoEdge, HardEdge, SoftEdge,
       FringeProfile, Enge, Tanh, LinearRamp, Tabulated, enge_default,
       FringeMethod, Splitting, Magnus, TaylorIntegrate,
       FringeModel, EdgeMap, edge_model, fringe_sequence, edge_map,
       profile_function, fringe_integral, fringe_observables, fringe_hamiltonian

# ----------------------------------------------------------------------------- models
"""
    EdgeModel

Supertype of edge (fringe-field) models: [`NoEdge`](@ref), [`HardEdge`](@ref),
[`SoftEdge`](@ref). A model is chosen per magnet family in a [`FringeModel`](@ref)
when the ring is built; elements themselves carry no fringe flag.
"""
abstract type EdgeModel end
"""
    FringeProfile

Supertype of longitudinal field profiles F(z) (0 outside, 1 inside the magnet):
[`Enge`](@ref), [`Tanh`](@ref), [`LinearRamp`](@ref), [`Tabulated`](@ref).
They are normalised and centred on the effective field boundary by [`profile_function`](@ref).
"""
abstract type FringeProfile end
"""
    FringeMethod

Supertype of integrators for the s-dependent soft-edge Hamiltonian:
[`Splitting`](@ref), [`Magnus`](@ref), [`TaylorIntegrate`](@ref).
"""
abstract type FringeMethod end

"No fringe field: the edge is the identity."
struct NoEdge <: EdgeModel end

@doc raw"""
    HardEdge(; order = 2, gap = nothing, fint = nothing)

Thin hard-edge model at the magnet ends; the map is ``e^{-{:}K{:}}z`` (symplectic).

# Fields
- `order::Int` – 1: linear dipole-edge kick, ``K = -\tfrac{h\tan\varepsilon}{2}x^2 + \tfrac{h\tan\psi}{2}y^2``
  (``p_x \mathrel{+}= h\tan\varepsilon\,x``, ``p_y \mathrel{-}= h\tan\psi\,y``) with
  ``\psi = \varepsilon - hg\,\mathrm{FINT}(1+\sin^2\varepsilon)/\cos\varepsilon``;
  2 or 3: nonlinear edge, the vertical term uses ``\tan(\psi \pm x')`` with the exact
  ``x' = p_x/p_s``; ``K`` is kept to degree 3 (order 2) or 4 (order 3).
- `gap`, `fint` – override the magnet's full gap [m] and FINT (`nothing` = use the magnet's).

For a `Quadrupole` (order ≥ 2) the cubic hard-edge generator
``K = \pm\frac{k_1}{12(1+\delta)}(x^3p_x + 3xy^2p_x - 3x^2yp_y - y^3p_y)`` is used; order 1 is the identity.
Other element types raise an error.
"""
struct HardEdge <: EdgeModel
    order::Int
    gap::Union{Nothing,Float64}
    fint::Union{Nothing,Float64}
end
HardEdge(; order::Int = 2, gap = nothing, fint = nothing) =
    HardEdge(order, gap === nothing ? nothing : Float64(gap), fint === nothing ? nothing : Float64(fint))

"""
    enge_default :: Vector{Float64}

Unverified default Enge coefficients of a dipole edge (z in units of the full
gap). Stored in this repository with no citation; not a published standard.
"""
const enge_default = [0.478959, 1.911289, -1.185953, 1.630554, -1.082657, 0.318111]

@doc raw"""
    Enge(coeffs = enge_default; gap = 0.05)

Enge profile ``F(z) = 1/\big(1 + \exp P(-z/g)\big)``, ``P(t) = \sum_{k\ge 1} c_k t^{k-1}``, ``g`` = full gap [m].
`coeffs::Vector{Float64}` are the Enge coefficients (dimensionless, ``z`` in units of the gap).
Default `coeffs` is [`enge_default`](@ref), an unverified set without a citation.
"""
struct Enge <: FringeProfile
    coeffs::Vector{Float64}
    gap::Float64
end
Enge(coeffs::AbstractVector = enge_default; gap = 0.05) = Enge(collect(Float64, coeffs), Float64(gap))

@doc raw"""
    Tanh(lambda)

Profile ``F(z) = \tfrac12\big(1 + \tanh(z/\lambda)\big)``; `lambda` [m] is the fall-off length.
"""
struct Tanh <: FringeProfile
    lambda::Float64
end

"Trapezoid edge: F rises linearly over `width`."
struct LinearRamp <: FringeProfile
    width::Float64
end

"Tabulated profile B(z); z along the beam, any units of B (normalised to max |B|)."
struct Tabulated <: FringeProfile
    z::Vector{Float64}
    B::Vector{Float64}
end

"""
    Splitting(; n = 20, order = 4)

Soft-edge integration by operator splitting. The region [−W, W] is cut into `n`
slices with piecewise-constant field (slice averages of F', F''). Each slice is
the exact Lie map exp(−Δs :H_k:). With `order = 4` every slice is a Yoshida
triple jump (weights W1, W0, W1). A product of Lie maps, so symplectic to the
truncation order.
"""
struct Splitting <: FringeMethod
    n::Int
    order::Int
end
Splitting(; n::Int = 20, order::Int = 4) = Splitting(n, order)

@doc raw"""
    Magnus(; n = 10, order = 4)

Soft-edge integration by the Magnus expansion on `n` sub-intervals of length ``\Delta s``.
- Order 2: ``\Omega = \Delta s\,H(s_m)``.
- Order 4: two Gauss nodes ``H_{1,2} = H(s_m \mp \Delta s\sqrt3/6)``,
  ``\Omega = \Delta s(H_1 + H_2)/2 - (\sqrt3/12)\Delta s^2[H_1, H_2]``.

Each interval map is ``e^{-{:}\Omega{:}}``.
"""
struct Magnus <: FringeMethod
    n::Int
    order::Int
end
Magnus(; n::Int = 10, order::Int = 4) = Magnus(n, order)

"""
    TaylorIntegrate(; tol = 1e-10, nmax = 512)

Soft-edge integration of the polynomial vector field with RK4 in TaylorN. Starts
with 16 steps and halves the step until the map changes by less than `tol`
(max coefficient) or `nmax` steps are reached. Not symplectic by construction;
it is the reference for the other methods.
"""
struct TaylorIntegrate <: FringeMethod
    tol::Float64
    nmax::Int
end
TaylorIntegrate(; tol = 1e-10, nmax::Int = 512) = TaylorIntegrate(Float64(tol), nmax)

"""
    SoftEdge(; profile = Enge(), method = Splitting())

Soft-edge dipole model: the field b(s) = h0 F(s) is integrated through [−W, W]
around the effective field boundary (EFB) using `profile::FringeProfile`, with
`method::FringeMethod`.

# Layout
The soft region replaces "drift W, then body W" of the hard-edge layout, so the
lattice length is unchanged:
- EdgeIn = Body(−W) ∘ Soft ∘ Drift(−W);
- EdgeOut = Drift(−W) ∘ Soft_out ∘ Body(−W).

A pole-face angle adds a linear wedge at the EFB. The magnet must satisfy 2W < L.
Implemented for `SBend` only.
"""
struct SoftEdge <: EdgeModel
    profile::FringeProfile
    method::FringeMethod
end
SoftEdge(; profile::FringeProfile = Enge(), method::FringeMethod = Splitting()) = SoftEdge(profile, method)

"""
    FringeModel(; dipole = NoEdge(), quad = NoEdge(), sext = NoEdge(), override = Dict())

Choice of [`EdgeModel`](@ref) per magnet family, with per-name overrides
(`override::Dict{String,EdgeModel}`, keys are element names).
`FringeModel()` reproduces the plain ring exactly.

```julia
model = FringeModel(dipole = SoftEdge(profile = Enge(gap = 0.05)),
                    quad = HardEdge(order = 2),
                    override = Dict("BM3" => HardEdge()))
F = ring_lie_factor(lattice, model)
```
"""
struct FringeModel
    dipole::EdgeModel
    quad::EdgeModel
    sext::EdgeModel
    override::Dict{String,EdgeModel}
end
FringeModel(; dipole::EdgeModel = NoEdge(), quad::EdgeModel = NoEdge(), sext::EdgeModel = NoEdge(),
            override = Dict{String,EdgeModel}()) =
    FringeModel(dipole, quad, sext, Dict{String,EdgeModel}(String(k) => v for (k, v) in override))

"""
    edge_model(m::FringeModel, e) -> EdgeModel

Model that applies to element `e`: the override for `e.name` if present, else the
family entry (`dipole` for SBend, `quad` for Quadrupole, `sext` for Sextupole),
else `NoEdge()`.
"""
function edge_model(m::FringeModel, e::AbstractElement)
    haskey(m.override, e.name) && return m.override[e.name]
    e isa SBend && return m.dipole
    e isa Quadrupole && return m.quad
    e isa Sextupole && return m.sext
    return NoEdge()
end

Base.show(io::IO, ::NoEdge) = print(io, "NoEdge()")
Base.show(io::IO, m::HardEdge) = print(io, "HardEdge(order=", m.order,
    m.gap === nothing ? "" : ", gap=$(m.gap)", m.fint === nothing ? "" : ", fint=$(m.fint)", ")")
Base.show(io::IO, p::Enge) = print(io, "Enge(gap=", p.gap, p.coeffs == enge_default ? "" : ", coeffs=$(p.coeffs)", ")")
Base.show(io::IO, p::Tanh) = print(io, "Tanh(", p.lambda, ")")
Base.show(io::IO, p::LinearRamp) = print(io, "LinearRamp(", p.width, ")")
Base.show(io::IO, p::Tabulated) = print(io, "Tabulated(", length(p.z), " points)")
Base.show(io::IO, m::Splitting) = print(io, "Splitting(n=", m.n, ", order=", m.order, ")")
Base.show(io::IO, m::Magnus) = print(io, "Magnus(n=", m.n, ", order=", m.order, ")")
Base.show(io::IO, m::TaylorIntegrate) = print(io, "TaylorIntegrate(tol=", m.tol, ")")
Base.show(io::IO, m::SoftEdge) = print(io, "SoftEdge(", m.profile, ", ", m.method, ")")

# ----------------------------------------------------------------------------- profiles
_raw(p::Enge, z) = (u = -z / p.gap; P = sum(p.coeffs[k] * u^(k - 1) for k in eachindex(p.coeffs));
                   P > 700 ? 0.0 : 1 / (1 + exp(P)))
_raw(p::Tanh, z) = (1 + tanh(z / p.lambda)) / 2
_raw(p::LinearRamp, z) = clamp(z / p.width + 0.5, 0.0, 1.0)
function _raw(p::Tabulated, z)
    zz, B = p.z, p.B ./ maximum(abs, p.B)
    if B[1] > B[end]                       # falling table: mirror to rising
        zz = -reverse(zz); B = reverse(B)
    end
    z <= zz[1] && return 0.0 + B[1] * 0
    z >= zz[end] && return 1.0
    k = searchsortedlast(zz, z)
    t = (z - zz[k]) / (zz[k+1] - zz[k])
    return B[k] + t * (B[k+1] - B[k])
end

function _halfwidth(p::Enge)
    z = 0.0
    while z < 40p.gap && (_raw(p, -z) > 1e-10 || 1 - _raw(p, z) > 1e-10)
        z += p.gap / 20
    end
    return z
end
_halfwidth(p::Tanh) = 12p.lambda
_halfwidth(p::LinearRamp) = p.width / 2
_halfwidth(p::Tabulated) = maximum(abs, p.z)

function _simpson(f, a, b; n = 4000)
    h = (b - a) / n
    s = f(a) + f(b)
    for k in 1:n-1
        s += (isodd(k) ? 4 : 2) * f(a + k * h)
    end
    return s * h / 3
end

"""
profile_function(p) -> (F, W, shift)
Normalised profile F(z), z = 0 at the effective field boundary (equal area:
int (F - step) dz = 0), F = 0 below -W and 1 above W.
"""
function profile_function(p::FringeProfile)
    W0 = _halfwidth(p)
    c = 0.0
    for _ in 1:40
        L = W0 + abs(c)
        G = _simpson(z -> _raw(p, z + c), -L, 0.0) + _simpson(z -> _raw(p, z + c) - 1, 0.0, L)
        c -= G
        abs(G) < 1e-13 && break
    end
    W = W0 + abs(c)
    return (z -> _raw(p, z + c), W, c)
end

"Fringe integral I = int F (1 - F) dz  (= g * FINT for a dipole of full gap g)."
function fringe_integral(p::FringeProfile)
    F, W, _ = profile_function(p)
    return _simpson(z -> F(z) * (1 - F(z)), -W, W; n = 8000)
end

# ----------------------------------------------------------------------------- sequence
"Thin edge of a parent magnet; side = :in or :out."
struct EdgeMap
    name::String
    parent::AbstractElement
    side::Symbol
    model::EdgeModel
end

const _LAST_FRINGE = Ref("")
function _announce_fringe(model::FringeModel)
    (model.dipole isa NoEdge && model.quad isa NoEdge && model.sext isa NoEdge &&
     isempty(model.override)) && return
    key = repr((model.dipole, model.quad, model.sext, sort(collect(keys(model.override)))))
    _LAST_FRINGE[] == key && return
    _LAST_FRINGE[] = key
    @info "fringe model enabled" dipole = model.dipole quad = model.quad sext = model.sext override = model.override
end

"Element list with EdgeMap items around every magnet whose model is not NoEdge."
function fringe_sequence(lattice::Lattice, model::FringeModel)
    _announce_fringe(model)
    expanded = Any[]
    for e in lattice.elements
        m = edge_model(model, e)
        if m isa NoEdge
            push!(expanded, e)
        else
            push!(expanded, EdgeMap(e.name * "_in", e, :in, m), e, EdgeMap(e.name * "_out", e, :out, m))
        end
    end
    return expanded
end

# ----------------------------------------------------------------------------- hard edges
function _hard(e::SBend, side, m::HardEdge, ref)
    x, px, y, py, l, d = Z
    h = curvature(e)
    ang = side === :in ? e.e1 : e.e2
    g = m.gap === nothing ? e.gap : m.gap
    fi = m.fint === nothing ? e.fint : m.fint
    psi = ang - h * g * fi * (1 + sin(ang)^2) / cos(ang)
    te, tp = tan(ang), tan(psi)
    K = -h * te * x^2 / 2
    if m.order <= 1
        K += h * tp * y^2 / 2
    else
        ps = sqrt((1 + d)^2 - px^2 - py^2)
        t = px / ps                                 # tan of the slope angle x'
        T = side === :in ? (tp + t) / (1 - tp * t) : (tp - t) / (1 + tp * t)
        K += h * y^2 / 2 * T
        m.order == 2 && (K = truncate_to(K, 3))
    end
    K -= constant_term(K)
    return lie_map(K, 1.0)
end

function _hard(e::Quadrupole, side, m::HardEdge, ref)
    m.order <= 1 && return collect(Z)               # no linear hard-edge effect
    x, px, y, py, l, d = Z
    s = side === :in ? 1.0 : -1.0
    K = s * e.k1 / (12 * (1 + d)) * (x^3 * px + 3x * y^2 * px - 3x^2 * y * py - y^3 * py)
    return lie_map(K, 1.0)
end

_hard(e::AbstractElement, side, m::HardEdge, ref) =
    error("HardEdge is implemented for SBend and Quadrupole, not for $(kind_string(e)) $(e.name)")

# ----------------------------------------------------------------------------- soft edges
"""
Soft-edge Hamiltonian at field level F, slope F' (per metre), curvature h0.
b1, b2: gradient and sextupole-like terms of the midplane field
h0 F(s + x theta) of a straight EFB seen from the curved frame.
"""
function fringe_hamiltonian(h0::Real, F::Real, Fp::Real, ref::RefParticle; b1::Real = 0.0, b2::Real = 0.0)
    x, px, y, py, l, d = Z
    h = h0 * F
    ax = -y^2 * (h0 * Fp) / (2 * (1 + h * x))
    ps = sqrt((1 + d)^2 - (px - ax)^2 - py^2)
    # field terms of the curved frame: the force carries the metric factor (1 + h x), so the gradient
    # b1 contributes h b1 x^2 (horizontal) and h b1 x y (vertical) as well; h b1 (x^3/3 - x y^2/2) is
    # the lowest-order curl-free completion (A_s = -int (1 + h x) B_y dx).  Without it the soft edge
    # had a spurious sextupole-like term: int h b1 ds = -int b2/2 ds exactly (theta' = h0 F).
    H = energy_term(d, ref) - (1 + h * x) * ps + h * x + h^2 * x^2 / 2 +
        b1 * (x^2 - y^2) / 2 + h * b1 * (x^3 / 3 - x * y^2 / 2) + b2 * (x^3 - 3x * y^2) / 6
    return H - constant_term(H)
end

_along(F, side) = side === :in ? F : (s -> F(-s))

# Field seen by the curved frame: F(zeta), zeta = s + x theta(s).  The EFB is a
# straight line perpendicular to the straight (field-free) orbit outside the
# magnet; theta(s) is the angle of the local orbit to that straight orbit:
#   entrance  theta(s) =  h0 int_{-W}^{s} F,   exit  theta(s) = -h0 int_{s}^{W} F.
# The integrated gradient int h0 F' theta = h0^2 int F (1 - F) then reproduces the
# vertical FINT focusing of the hard-edge limit and cancels the horizontal deficit
# int (F^2 - step) of the weak focusing.
struct _Edge
    h0::Float64
    F::Function
    W::Float64
    side::Symbol
    ep::Float64
end
_Edge(h0, F, W, side) = _Edge(h0, F, W, side, 1e-5 * W)
# Derivatives as averages over a window [s - w, s + w]: exact for smooth F as
# w -> 0, and integrable for profiles with kinks (LinearRamp, Tabulated).
_dF(E::_Edge, s, w = E.ep) = (E.F(s + w) - E.F(s - w)) / (2w)
_d2F(E::_Edge, s, w = E.ep) = (E.F(s + 2w) - 2E.F(s) + E.F(s - 2w)) / (4w^2)
_theta(E::_Edge, s) = E.side === :in ? E.h0 * _simpson(E.F, -E.W, s; n = 400) :
                                       -E.h0 * _simpson(E.F, s, E.W; n = 400)
"H at a node s (F', F'' averaged over a window of half-width w)."
_Hnode(E::_Edge, s, ref, w = E.ep) = (th = _theta(E, s); Fp = _dF(E, s, w);
    fringe_hamiltonian(E.h0, E.F(s), Fp, ref; b1 = E.h0 * Fp * th, b2 = E.h0 * _d2F(E, s, w) * th^2))
"H of a slice [a,b] (slice averages of F', F'')."
function _Hslice(E::_Edge, a, b, ref)
    m = (a + b) / 2; th = _theta(E, m)
    Fp = (E.F(b) - E.F(a)) / (b - a)
    Fpp = (_dF(E, b) - _dF(E, a)) / (b - a)
    return fringe_hamiltonian(E.h0, E.F(m), Fp, ref; b1 = E.h0 * Fp * th, b2 = E.h0 * Fpp * th^2)
end

function _soft(E::_Edge, m::Splitting, ref)
    W = E.W
    ws = splitting_weights(m.order >= 4 ? Yoshida(4) : Leapfrog())
    hs = 2W / m.n
    M = collect(Z)
    for k in 0:m.n-1
        a = -W + k * hs
        for w in ws
            b = a + w * hs
            M = compose_maps(lie_map(_Hslice(E, a, b, ref), b - a), M)
            a = b
        end
    end
    return M
end

function _soft(E::_Edge, m::Magnus, ref)
    W = E.W
    hs = 2W / m.n
    Hs(s) = _Hnode(E, s, ref, max(E.ep, hs * sqrt(3) / 6))
    M = collect(Z)
    for k in 0:m.n-1
        mid = -W + (k + 0.5) * hs
        if m.order <= 2
            Om = hs * Hs(mid)
        else
            H1 = Hs(mid - hs * sqrt(3) / 6); H2 = Hs(mid + hs * sqrt(3) / 6)
            Om = hs / 2 * (H1 + H2) - sqrt(3) / 12 * hs^2 * poisson(H1, H2)
        end
        M = compose_maps(lie_map(Om, 1.0), M)
    end
    return M
end

function _rk4(E::_Edge, n, ref)
    W = E.W
    hs = 2W / n
    function f(z, s)
        H = _Hnode(E, s, ref, max(E.ep, hs / 4))
        V = hamiltonian_field(H)
        return [truncate_to(evaluate(V[i], z), MAP_ORDER) for i in 1:6]
    end
    r = integrate_fixed(f, collect(Z), -W, W, n; method = RK4())
    return [truncate_to(zi, MAP_ORDER) for zi in r.z]
end

const _RK_STEPS = Dict{Any,Int}()

function _soft(E::_Edge, m::TaylorIntegrate, ref)
    n = 16
    M = _rk4(E, n, ref)
    while 2n <= m.nmax
        M2 = _rk4(E, 2n, ref)
        diff = maximum(maxabs(M2[i] - M[i]) for i in 1:6)
        M, n = M2, 2n
        diff < m.tol && break
    end
    _RK_STEPS[(E.h0, E.W, m)] = n
    return M
end

function _soft_edge(e::SBend, side, m::SoftEdge, ref)
    F0, W, _ = profile_function(m.profile)
    2W < e.L || error("fringe region 2W = $(2W) m is longer than $(e.name) (L = $(e.L) m)")
    F = _along(F0, side)
    h0 = curvature(e)
    S = _soft(_Edge(h0, F, W, side), m.method, ref)
    ang = side === :in ? e.e1 : e.e2
    if ang != 0                                  # pole-face rotation: linear wedge at the EFB
        x, px, y, py, l, d = Z
        Wg = lie_map(-h0 * tan(ang) * (x^2 - y^2) / 2, 1.0)
        S = side === :in ? compose_maps(S, Wg) : compose_maps(Wg, S)
    end
    Bm = lie_map(hamiltonian(e, ref), -W)
    Dm = lie_map(hamiltonian(Drift("_fringe_drift", 1.0), ref), -W)
    return side === :in ? compose_maps(Bm, compose_maps(S, Dm)) : compose_maps(Dm, compose_maps(S, Bm))
end
_soft_edge(e::AbstractElement, side, m::SoftEdge, ref) =
    error("SoftEdge is implemented for SBend only, not for $(kind_string(e)) $(e.name)")

const _ECACHE = Dict{Any,Any}()

# Cache key of an edge model: repr() for analytic models, plus the table contents for a
# Tabulated profile (its repr is only "Tabulated(N points)", so two different tables of the
# same length used to share one cache entry).
_model_key(m::EdgeModel) = m isa SoftEdge && m.profile isa Tabulated ?
    (repr(m), Tuple(m.profile.z), Tuple(m.profile.B)) : repr(m)

"""
    edge_map(E::EdgeMap, ref) -> Vector{TaylorN}

Polynomial map (order MAP_ORDER) of a thin edge item: identity for `NoEdge`, the
hard-edge Lie map for `HardEdge`, or the composed soft-edge map for `SoftEdge`.
Cached per (parent magnet, side, model, p0).
"""
edge_map(E::EdgeMap, ref::RefParticle) =
    get!(_ECACHE, (E.parent, E.side, _model_key(E.model), ref.p0_MeV)) do
        E.model isa NoEdge && return collect(Z)
        E.model isa HardEdge && return _hard(E.parent, E.side, E.model, ref)
        return _soft_edge(E.parent, E.side, E.model, ref)
    end

# ----------------------------------------------------------------------------- ring
"One factor per item of fringe_sequence (edges are separate factors)."
function build_factors(lattice::Lattice, model::FringeModel)
    cache = Dict{Any,LieFactor}()
    factors = Any[]
    for it in fringe_sequence(lattice, model)
        if it isa RFCavity
            push!(factors, it)
        elseif it isa EdgeMap
            push!(factors, get!(() -> dragt_finn_factor(edge_map(it, lattice.ref); name = it.name), cache,
                            (it.parent, it.side, _model_key(it.model))))
        else
            push!(factors, get!(() -> dragt_finn_factor(element_lie_map(it, lattice.ref); name = it.name), cache, it))
        end
    end
    return factors
end

"One-turn Lie factor of the ring with the given fringe model."
ring_lie_factor(lattice::Lattice, model::FringeModel; name = lattice.name) =
    product_factors(LieFactor[F for F in build_factors(lattice, model) if F isa LieFactor]; name)

"""
    fringe_observables(lattice, model) -> NamedTuple

Fractional tunes, first/second-order chromaticity of the one-turn map and the
symplecticity of the edge maps (linear part, and Jacobian at a test point z0).
"""
function fringe_observables(lattice::Lattice, model::FringeModel = FringeModel(); dd = 1e-4,
                            z0 = [1e-3, 1e-4, 1e-3, 1e-4, 0.0, 1e-3])
    t0 = time()
    F = ring_lie_factor(lattice, model)
    M = factor_to_map(F)
    ch = chromaticity_orders(M; dd)
    se_lin, se_pt = 0.0, 0.0
    for it in fringe_sequence(lattice, model)
        it isa EdgeMap || continue
        E = edge_map(it, lattice.ref)
        se_lin = max(se_lin, symplectic_error(map_linear_part(E)))
        se_pt = max(se_pt, symplectic_error(jacobian_of_map(E, z0)))
    end
    return (Qx = ch.x[1], Qy = ch.y[1], xix = ch.x[2], xiy = ch.y[2], xi2x = ch.x[3], xi2y = ch.y[3],
            symp_edge_lin = se_lin, symp_edge_pt = se_pt, symp_ring = symplectic_error(F.R),
            seconds = time() - t0, model = model)
end

end # module
