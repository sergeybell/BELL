# =============================================================================
# hamiltonians.jl - s-Hamiltonian of every element
# =============================================================================
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# One formula per element. hamiltonian(e, ref) evaluates that formula on Z
# and returns a TaylorN polynomial (JET_ORDER, constant removed) for the Lie path.
# hamiltonian(e, ref, z) evaluates the same expression on a 6-vector z
# (numbers, TaylorN, or DA). Magnet maps and a thick RF cavity integrate that
# value. A thin cavity has no continuous Hamiltonian. The field equations in
# physics/fields.jl are the comparison formulation and do not call this module.
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

"Drop the constant and warn on linear terms only for the Lie-path polynomial."
_finish(H::TaylorN, e::AbstractElement) = _clean(H, e)
_finish(H, ::AbstractElement) = H

function _parts(z, ref::RefParticle)
    x, px, y, py, _, d = z
    ps = sqrt((1 + d)^2 - px^2 - py^2)
    return x, y, ps, energy_term(d, ref)
end

@doc raw"""
    hamiltonian(e::AbstractElement, ref::RefParticle) -> TaylorN
    hamiltonian(e::AbstractElement, ref::RefParticle, z::AbstractVector)

s-Hamiltonian of element `e` in canonical variables ``z = (x, p_x, y, p_y, l, \delta)``.

With ``p_s = \sqrt{(1+\delta)^2 - p_x^2 - p_y^2}`` and [`energy_term`](@ref):

| element | ``H`` |
|:--|:--|
| `Drift` | ``\mathrm{energy\_term}(\delta) - p_s`` |
| `Quadrupole` | ``\mathrm{energy\_term}(\delta) - p_s + \tfrac{k_1}{2}(x^2 - y^2)`` |
| `Sextupole` | ``\mathrm{energy\_term}(\delta) - p_s + \tfrac{k_2}{6}(x^3 - 3xy^2)`` |
| `SBend` | ``\mathrm{energy\_term}(\delta) - (1+hx)p_s + hx + \tfrac{h^2x^2}{2} + \tfrac{k_1}{2}(x^2 - y^2)`` |
| `Marker` | ``0`` |

`k_1 > 0` focuses in ``x``. `z` is a 6-vector of numbers, `TaylorN`, or DA;
`H` has the same scalar type. Calling this method on an `RFCavity` without
`path_length` throws: the RF frequency needs the ring length. A thin cavity
has no continuous Hamiltonian; its map is [`rf_kick`](@ref).
Without `z`, the Lie path evaluates the formula on [`Z`](@ref) and returns a
polynomial of degree [`JET_ORDER`](@ref) = 4. Only that polynomial has its
constant removed; linear terms warn that the reference orbit is not a solution.
"""
function hamiltonian(e::Drift, ref::RefParticle, z::AbstractVector)
    _, _, ps, en = _parts(z, ref)
    return _finish(en - ps, e)
end

function hamiltonian(e::Quadrupole, ref::RefParticle, z::AbstractVector)
    x, y, ps, en = _parts(z, ref)
    return _finish(en - ps + e.k1 / 2 * (x^2 - y^2), e)
end

function hamiltonian(e::Sextupole, ref::RefParticle, z::AbstractVector)
    x, y, ps, en = _parts(z, ref)
    return _finish(en - ps + e.k2 / 6 * (x^3 - 3x * y^2), e)
end

function hamiltonian(e::SBend, ref::RefParticle, z::AbstractVector)
    x, y, ps, en = _parts(z, ref)
    h = curvature(e)
    return _finish(en - (1 + h * x) * ps + h * x + h^2 * x^2 / 2 + e.k1 / 2 * (x^2 - y^2), e)
end

hamiltonian(e::Marker, ref::RefParticle, z::AbstractVector) = _finish(zero(z[1]), e)
hamiltonian(e::RFCavity, ref::RefParticle, ::AbstractVector) =
    error("hamiltonian(::RFCavity, ref, path_length): the RF frequency needs the path length")

hamiltonian(e::Drift, ref::RefParticle) = hamiltonian(e, ref, Z)
hamiltonian(e::Quadrupole, ref::RefParticle) = hamiltonian(e, ref, Z)
hamiltonian(e::Sextupole, ref::RefParticle) = hamiltonian(e, ref, Z)
hamiltonian(e::SBend, ref::RefParticle) = hamiltonian(e, ref, Z)
hamiltonian(e::Marker, ref::RefParticle) = hamiltonian(e, ref, Z)
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
    hamiltonian(cav::RFCavity, ref, path_length, z) -> scalar

Thick cavity (``L > 0``): ``H = \mathrm{energy\_term}(\delta) - p_s + G/L`` with
``G`` = [`rf_potential`](@ref) at `z` (default [`Z`](@ref)). For ``\sin\varphi_s \ne 0``
the TaylorN polynomial has a linear term in ``l`` (net acceleration of the
reference). A thin cavity (``L = 0``) has no continuous Hamiltonian.
"""
function hamiltonian(e::RFCavity, ref::RefParticle, path_length::Real, z::AbstractVector)
    e.L > 0 || error("thin RFCavity: use rf_potential (map exp(-:G:)), not hamiltonian")
    _, _, ps, en = _parts(z, ref)
    return _finish(en - ps + rf_potential(e, ref, path_length; z = z) / e.L, e)
end

hamiltonian(e::RFCavity, ref::RefParticle, path_length::Real) = hamiltonian(e, ref, path_length, Z)



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
