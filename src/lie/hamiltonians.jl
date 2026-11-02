# =============================================================================
# hamiltonians.jl - s-Hamiltonian of every element (Lie path only)
# =============================================================================
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# Polynomial H truncated at JET_ORDER for the Lie / Dragt–Finn path.
# Field ODEs for Float64 / DA tracking live in da/fields.jl and do not
# call this module.
# =============================================================================
module Hamiltonians

using TaylorSeries
using ..Coordinates, ..Elements, ..LieAlgebra

export hamiltonian, hamiltonian_formula, hamiltonian_field, rf_potential, rf_wavenumber

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

@doc raw"""
    hamiltonian(e::AbstractElement, ref::RefParticle) -> TaylorN

s-Hamiltonian of element `e` in canonical variables ``z = (x, p_x, y, p_y, l, \delta)``, as a
polynomial of degree [`JET_ORDER`](@ref) = 4 (constant term dropped).

With ``p_s = \sqrt{(1+\delta)^2 - p_x^2 - p_y^2}`` and [`energy_term`](@ref):

| element | ``H`` |
|:--|:--|
| `Drift` | ``\mathrm{energy\_term}(\delta) - p_s`` |
| `Quadrupole` | ``\mathrm{energy\_term}(\delta) - p_s + \tfrac{k_1}{2}(x^2 - y^2)`` |
| `Sextupole` | ``\mathrm{energy\_term}(\delta) - p_s + \tfrac{k_2}{6}(x^3 - 3xy^2)`` |
| `SBend` | ``\mathrm{energy\_term}(\delta) - (1+hx)p_s + hx + \tfrac{h^2x^2}{2} + \tfrac{k_1}{2}(x^2 - y^2)`` |
| `Marker` | ``0`` (thin `RFCavity`: error, it is a kick) |

The exact expressions are evaluated in `TaylorN` arithmetic. The constant is
removed, and linear terms warn that the reference orbit is not a solution.
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
hamiltonian(e::RFCavity, ref::RefParticle) =
    error("hamiltonian(::RFCavity, ref, path_length): the RF frequency needs the path length")
hamiltonian(e::AbstractElement, ref::RefParticle, path_length::Real) = hamiltonian(e, ref)

@doc raw"""
    rf_wavenumber(cav, ref, path_length) -> k [1/m]

``k = h\omega_0/c = 2\pi h\beta_0/C`` (``C`` = `path_length`).
"""
rf_wavenumber(cav::RFCavity, ref::RefParticle, path_length::Real) = 2pi * cav.harmonic * ref.beta / path_length

@doc raw"""
    rf_potential(cav, ref, path_length; z = Z) -> G

Integrated RF generator ``G = \int H_{rf}\,ds`` of a cavity (exact expression,
numbers / `TaylorN` / DA):

```math
G(l,\delta) = \frac{ZV}{p_0 c\,k}\cos\!\Big(\varphi_s + \frac{k\,l}{\beta_0\,v(\delta)/v_0}\Big),
\qquad k = \frac{2\pi h\beta_0}{C}.
```

Derivation: in the pair ``(q, P) = (-c\Delta t,\ E/(p_0c))`` the force of a
uniform gap is ``dP/ds = -\partial H/\partial q = ZV\sin(\varphi_s + kq)/(p_0cL)``.
Bell's ``(l, \delta)`` is the point transformation ``\delta = \delta(P)``,
``l = q\,dP/d\delta = q\beta = -v\Delta t``, so ``H`` is unchanged with
``q = l/(\beta_0 v/v_0)``; ``kq = 2\pi h\,l/(C\,v/v_0)`` is the RF phase of
[`rf_kick`](@ref) and [`canonical_to_laboratory`](@ref). Thick cavity:
``H_{rf} = G/L``. Thin cavity: the map is ``\exp(-{:}G{:})``, which is exactly
[`rf_kick`](@ref) (``\Delta t`` fixed, ``E \mathrel{+}= ZV\sin(\varphi_s+\varphi)``).
"""
function rf_potential(cav::RFCavity, ref::RefParticle, path_length::Real; z = Z)
    l, d = z[5], z[6]
    k = rf_wavenumber(cav, ref, path_length)
    q = l / (ref.beta * velocity_ratio(d, ref))
    return ref.particle.charge * cav.voltage_MV / (ref.p0_MeV * k) * cos(cav.phi_s + k * q)
end

@doc raw"""
    hamiltonian(cav::RFCavity, ref, path_length) -> TaylorN

Thick cavity (``L > 0``): ``H = \mathrm{energy\_term}(\delta) - p_s + G/L`` with
``G`` = [`rf_potential`](@ref). For ``\sin\varphi_s \ne 0`` ``H`` has a linear
term in ``l`` (net acceleration of the reference); see `dev/rf_acceleration`.
"""
function hamiltonian(e::RFCavity, ref::RefParticle, path_length::Real)
    e.L > 0 || error("thin RFCavity: use rf_potential (map exp(-:G:)), not hamiltonian")
    x, y, d, ps, en = _parts(ref)
    _clean(en - ps + rf_potential(e, ref, path_length) / e.L, e)
end



"""
    hamiltonian_formula(e) -> String

The formula of [`hamiltonian`](@ref) for the type of `e` as plain text (for printing).
"""
hamiltonian_formula(::Drift) = "H = energy_term(delta) - ps"
hamiltonian_formula(::Quadrupole) = "H = energy_term(delta) - ps + k1/2 (x^2 - y^2)"
hamiltonian_formula(::Sextupole) = "H = energy_term(delta) - ps + k2/6 (x^3 - 3 x y^2)"
hamiltonian_formula(::SBend) = "H = energy_term(delta) - (1+hx) ps + h x + h^2 x^2/2 + k1/2 (x^2-y^2)"
hamiltonian_formula(::RFCavity) = "H = energy_term(delta) - ps + ZV/(p0 c k L) cos(phi_s + k l/(beta0 v/v0)),  k = 2 pi h beta0/C"
hamiltonian_formula(::AbstractElement) = "H = 0"

"""
    hamiltonian_field(H) -> Vector

Polynomial vector field ``dz/ds = S\\nabla H`` of a TaylorN Hamiltonian.
"""
function hamiltonian_field(H)
    g = [differentiate(H, j) for j in 1:6]
    return [sum(SYMPLECTIC_S[i, j] * g[j] for j in 1:6) for i in 1:6]
end

end # module
