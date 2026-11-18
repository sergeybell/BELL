# =============================================================================
# fringe.jl - fringe-field model of a magnet edge
# =============================================================================
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# The lattice stays physical: an SBend carries L, angle, e1, e2 and optional
# gap / FINT. No element carries a fringe flag. The model is chosen when the
# ring is built. FringeModel() is NoEdge everywhere and leaves the ring as is.
#
# This module is the model only: edge types, longitudinal profiles, the
# s-dependent fringe Hamiltonian, and the sequence of edge items. It does not
# build a map. Hard-edge generators, the soft-edge flow, caches, and
# fringe_observables live in maps/lie/fringe_maps.jl.
#
# HardEdge is a thin lens at each end. The body Hamiltonian is not changed.
# SoftEdge replaces a slice of body and a slice of drift by the s-dependent
# Hamiltonian. The body gradient k1 is carried by that Hamiltonian as k1*F(s),
# so F = 1 reproduces the body and F = 0 is a drift. Half-width W is read from
# the profile. The magnet must satisfy 2W < L.
#
# enge_default is an unverified default coefficient set stored in this
# repository. There is no citation for it; it is not a standard.
# =============================================================================
module Fringe

using ..Coordinates, ..Elements, ..Structure

export EdgeModel, NoEdge, HardEdge, SoftEdge,
       FringeProfile, Enge, Tanh, LinearRamp, Tabulated, enge_default,
       FringeMethod, Splitting, Magnus, TaylorIntegrate,
       FringeModel, EdgeMap, edge_model, fringe_sequence,
       profile_function, fringe_integral, fringe_hamiltonian

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

Thin hard-edge model at the magnet ends. The map is built in `FringeMaps`.

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

Soft-edge integration of ``dz/ds = S\\nabla H(s)`` by [`element_flow`](@ref).
Starts with 16 steps and doubles the step count until the map changes by less
than `tol` (max coefficient) or `nmax` steps are reached. Not symplectic by
construction; it is the reference for the other methods.
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
`method::FringeMethod`. The body gradient enters as `k1 * F(s)`.

# Layout
The soft region replaces "drift W, then body W" of the hard-edge layout, so the
lattice length is unchanged:
- EdgeIn = Body(−W) ∘ Soft ∘ Drift(−W);
- EdgeOut = Drift(−W) ∘ Soft_out ∘ Body(−W).

A pole-face angle is a linear wedge at the EFB, local s = 0, between the
upstream and downstream halves of the fringe. The magnet must satisfy 2W < L.
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
    return get(m.override, e.name, _family_edge(m, e))
end
_family_edge(m::FringeModel, ::SBend) = m.dipole
_family_edge(m::FringeModel, ::Quadrupole) = m.quad
_family_edge(m::FringeModel, ::Sextupole) = m.sext
_family_edge(::FringeModel, ::AbstractElement) = NoEdge()

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
_no_edge(::NoEdge) = true
_no_edge(::EdgeModel) = false
function _announce_fringe(model::FringeModel)
    _no_edge(model.dipole) && _no_edge(model.quad) && _no_edge(model.sext) &&
        isempty(model.override) && return
    key = repr((model.dipole, model.quad, model.sext, sort(collect(keys(model.override)))))
    _LAST_FRINGE[] == key && return
    _LAST_FRINGE[] = key
    @info "fringe model enabled" dipole = model.dipole quad = model.quad sext = model.sext override = model.override
end

_expand_edge!(expanded, e::AbstractElement, ::NoEdge) = push!(expanded, e)
function _expand_edge!(expanded, e::AbstractElement, m::EdgeModel)
    push!(expanded, EdgeMap(e.name * "_in", e, :in, m), e, EdgeMap(e.name * "_out", e, :out, m))
    return expanded
end

"Element list with EdgeMap items around every magnet whose model is not NoEdge."
function fringe_sequence(lattice::Lattice, model::FringeModel)
    _announce_fringe(model)
    expanded = Any[]
    for e in lattice.elements
        _expand_edge!(expanded, e, edge_model(model, e))
    end
    return expanded
end

# ----------------------------------------------------------------------------- Hamiltonian
@doc raw"""
    fringe_hamiltonian(h0, F, Fp, ref, z; b1 = 0, b2 = 0, k1 = 0)

Soft-edge s-Hamiltonian at field level `F` and slope `Fp` [1/m], evaluated on
the canonical state `z` (numbers, `TaylorN`, or DA).

`h0` is the body curvature. `h = h0 * F`. `k1` is the local normalised
gradient. On a soft edge the caller passes the body gradient times `F`, so
`F = 1` and `Fp = 0` reproduce the body Hamiltonian of an [`SBend`](@ref)
(including its `k1`) and `F = 0` is a drift.

`b1` and `b2` are the gradient and sextupole-like terms of a straight effective
field boundary seen from the curved frame. They are not the design `k1`.
The metric factor `(1 + h x)` puts `h b1 (x^3/3 - x y^2/2)` into `H`; without
it the straight-boundary gradient leaves a spurious sextupole-like term.
The Lie path evaluates this function on [`Z`](@ref).
"""
function fringe_hamiltonian(h0::Real, F::Real, Fp::Real, ref::RefParticle, z::AbstractVector;
                            b1::Real = 0.0, b2::Real = 0.0, k1::Real = 0.0)
    x, px, y, py, _, d = z
    h = h0 * F
    ax = -y^2 * (h0 * Fp) / (2 * (1 + h * x))
    ps = sqrt((1 + d)^2 - (px - ax)^2 - py^2)
    # k1 is only the quadrupole term. Folding it into b1 would add the curved-frame
    # cubic that the body Hamiltonian does not carry, and Body(-W) would no longer cancel.
    return energy_term(d, ref) - (1 + h * x) * ps + h * x + h^2 * x^2 / 2 +
           (k1 + b1) * (x^2 - y^2) / 2 + h * b1 * (x^3 / 3 - x * y^2 / 2) +
           b2 * (x^3 - 3x * y^2) / 6
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
_Edge(h0, F, W, side) = _Edge(Float64(h0), F, Float64(W), side, 1e-5 * Float64(W))

"Fringe field of one `SBend` end. `F` is already oriented along the integration coordinate."
function edge_field(e::SBend, profile::FringeProfile, side::Symbol)
    F0, W, _ = profile_function(profile)
    2W < e.L || error("fringe region 2W = $(2W) m is longer than $(e.name) (L = $(e.L) m)")
    return _Edge(curvature(e), _along(F0, side), W, side)
end

# Derivatives as averages over a window [s - w, s + w]: exact for smooth F as
# w -> 0, and integrable for profiles with kinks (LinearRamp, Tabulated).
_dF(E::_Edge, s, w = E.ep) = (E.F(s + w) - E.F(s - w)) / (2w)
_d2F(E::_Edge, s, w = E.ep) = (E.F(s + 2w) - 2E.F(s) + E.F(s - 2w)) / (4w^2)
_theta(E::_Edge, s) = E.side === :in ? E.h0 * _simpson(E.F, -E.W, s; n = 400) :
                                       -E.h0 * _simpson(E.F, s, E.W; n = 400)

"""
    fringe_hamiltonian(E, ref, s, z; k1 = 0, w = E.ep)

Hamiltonian of an edge field at station `s`. `k1` is the body gradient; the
local gradient is `k1 * F(s)`.
"""
function fringe_hamiltonian(E::_Edge, ref::RefParticle, s::Real, z::AbstractVector;
                            k1::Real = 0.0, w::Real = E.ep)
    th = _theta(E, s)
    Fp = _dF(E, s, w)
    F = E.F(s)
    return fringe_hamiltonian(E.h0, F, Fp, ref, z; b1 = E.h0 * Fp * th,
                              b2 = E.h0 * _d2F(E, s, w) * th^2, k1 = k1 * F)
end

"Hamiltonian of a slice [a, b] (slice averages of F' and F''). `k1` is the body gradient."
function _slice_hamiltonian(E::_Edge, a::Real, b::Real, ref::RefParticle, z::AbstractVector;
                            k1::Real = 0.0)
    m = (a + b) / 2
    th = _theta(E, m)
    Fp = (E.F(b) - E.F(a)) / (b - a)
    Fpp = (_dF(E, b) - _dF(E, a)) / (b - a)
    F = E.F(m)
    return fringe_hamiltonian(E.h0, F, Fp, ref, z; b1 = E.h0 * Fp * th,
                              b2 = E.h0 * Fpp * th^2, k1 = k1 * F)
end

end # module
