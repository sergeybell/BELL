# =============================================================================
# Hamiltonians.jl - the s-Hamiltonian of every element, expanded to order 3
# =============================================================================
#
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# WHAT
#   hamiltonian(element, ref) returns H(z) as a TaylorN polynomial in the
#   canonical variables z = (x, px, y, py, l, delta), expanded to total order
#   3 (constant dropped).  It is computed by TaylorSeries arithmetic directly
#   from the EXACT expressions below, so no hand expansion can go wrong.
#
#   With ps = sqrt((1+delta)^2 - px^2 - py^2) and the energy term
#   e(delta) = sqrt((1+delta)^2 + 1/(beta0 gamma0)^2)/beta0 - 1/beta0^2 :
#
#     Drift       H = e(delta) - ps
#     Quadrupole  H = e(delta) - ps + k1/2 (x^2 - y^2)
#     Sextupole   H = e(delta) - ps + k2/6 (x^3 - 3 x y^2)
#     SBend       H = e(delta) - (1 + h x) ps + h x + h^2 x^2/2 + k1/2 (x^2 - y^2)
#
#   Origin of the bend terms: in the curved frame (radius 1/h) the
#   Hamiltonian is  e(delta) - (1+hx) ps - (1+hx) a_s, with a_s = q A_s/p0.
#   For the uniform design field B0 = (B rho) h we need
#   b_y = (1/(1+hx)) d/dx[ -(1+hx) a_s ] = h, i.e. -(1+hx) a_s = h x + h^2 x^2/2.
#   Check: dH/dx at z = 0 equals -h + h = 0, so the reference arc is a
#   solution.  The k1 term in a bend uses the straight-magnet potential
#   (curvature corrections to a gradient in a bend are O(h k1 x^3) and are
#   not needed here because all OptiM bends have G = 0).
#
#   Hamilton's equations:  dz/ds = [z, H] = S grad H   with
#   [f, g] = sum_k (df/dq_k dg/dp_k - df/dp_k dg/dq_k),  q = (x, y, l),
#   p = (px, py, delta).
#
# WHY
#   Everything else (linear matrices = exp(L S Hess H), Lie generators,
#   Taylor maps) is derived from these H; this file is the single place where
#   element physics is defined.
#
# WHAT YOU SEE
#   Nothing printed; scripts/02 prints H of every element.
# =============================================================================
module Hamiltonians

using TaylorSeries
using ..Coordinates, ..Elements, ..LieAlgebra

export hamiltonian, hamiltonian_formula

function _clean(H::TaylorN, e)
    H = H - constant_term(H)
    lin = maxabs(TaylorN(H[1], JET_ORDER))
    lin > 1e-12 && @warn "Hamiltonian of $(e.name) has linear terms (max coeff $lin): reference orbit is not a solution"
    return H
end

function _parts(ref)
    x, px, y, py, l, d = Z
    ps = sqrt((1 + d)^2 - px^2 - py^2)
    return x, y, d, ps, energy_term(d, ref)
end

"""
    hamiltonian(e::AbstractElement, ref::RefParticle) -> TaylorN

s-Hamiltonian of element `e` in canonical variables z = (x, px, y, py, l, δ), as a
polynomial of degree [`JET_ORDER`](@ref) (constant term dropped).

With ps = √((1+δ)² − px² − py²) and e(δ) = [`energy_term`](@ref)(δ, ref):

| element | H |
|:--|:--|
| `Drift`, `Cavity` (RF off) | e(δ) − ps |
| `Quadrupole` | e(δ) − ps + k1(x² − y²)/2 |
| `Sextupole` | e(δ) − ps + k2(x³ − 3xy²)/6 |
| `SBend` | e(δ) − (1+hx) ps + hx + h²x²/2 + k1(x² − y²)/2 |
| `Marker`, thin `RFCavity` | 0 |

# Algorithm
The exact expressions are evaluated in `TaylorN` arithmetic, so the expansion
(square root included) is done by the polynomial engine, not by hand. Then the
constant is removed, and the result is checked to have no linear terms. Linear
terms would mean the reference orbit is not a solution, and a warning is issued.

The bend terms come from the curved-frame vector potential
−(1+hx)a_s = hx + h²x²/2 (uniform field h·Bρ). ∂H/∂x = 0 at z = 0.
This is the only place in Bell where element physics is defined: linear
matrices, Lie maps, DA maps and the exact integrator all derive from these H.
"""
hamiltonian(e::Drift, ref::RefParticle) = (p = _parts(ref); _clean(p[5] - p[4], e))

function hamiltonian(e::Quadrupole, ref::RefParticle)
    x, y, d, ps, en = _parts(ref)
    _clean(en - ps + e.k1 / 2 * (x^2 - y^2), e)
end

function hamiltonian(e::Sextupole, ref::RefParticle)
    x, y, d, ps, en = _parts(ref)
    _clean(en - ps + e.k2 / 6 * (x^3 - 3x * y^2), e)
end

function hamiltonian(e::SBend, ref::RefParticle)
    x, y, d, ps, en = _parts(ref)
    h = curvature(e)
    _clean(en - (1 + h * x) * ps + h * x + h^2 * x^2 / 2 + e.k1 / 2 * (x^2 - y^2), e)
end

hamiltonian(e::Marker, ref::RefParticle) = zero(Z[1])
# The RF part of a finite cavity is s- and time-dependent -> ExactIntegrator.jl.
# For maps / linear optics the cavity is taken with RF OFF, i.e. as a drift.
hamiltonian(e::Cavity, ref::RefParticle) = hamiltonian(Drift(e.name, e.L), ref)
hamiltonian(e::RFCavity, ref::RefParticle) =
    error("RFCavity is a thin kick (see Tracking.rf_kick); it has no s-Hamiltonian in phase 1")

"""
    hamiltonian_formula(e) -> String

The formula of [`hamiltonian`](@ref) for the type of `e` as plain text (for printing).
"""
hamiltonian_formula(::Drift) = "H = e(delta) - ps"
hamiltonian_formula(::Quadrupole) = "H = e(delta) - ps + k1/2 (x^2 - y^2)"
hamiltonian_formula(::Sextupole) = "H = e(delta) - ps + k2/6 (x^3 - 3 x y^2)"
hamiltonian_formula(::SBend) = "H = e(delta) - (1+hx) ps + h x + h^2 x^2/2 + k1/2 (x^2-y^2)"
hamiltonian_formula(::AbstractElement) = "H = 0"

end # module
