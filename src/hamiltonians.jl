# =============================================================================
# hamiltonians.jl - s-Hamiltonian of every element
# =============================================================================
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# Two formulations of the same equations live here:
#   hamiltonian(e, ref)          polynomial (TaylorN), truncated at JET_ORDER
#   magnetic_field / cavity_field  exact (unexpanded) vector field dz/ds = S ∇H
# The tracker does not rewrite these expressions.
# =============================================================================
module Hamiltonians

using TaylorSeries
using ..Coordinates, ..Elements, ..LieAlgebra

export hamiltonian, hamiltonian_formula, hamiltonian_field, magnetic_field, cavity_field,
       magnetic_parameters, velocity_derivative

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
| `Drift`, `LongCavity` (RF off) | ``\mathrm{energy\_term}(\delta) - p_s`` |
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
hamiltonian(e::LongCavity, ref::RefParticle) = hamiltonian(Drift(e.name, e.L), ref)
hamiltonian(e::RFCavity, ref::RefParticle) =
    error("RFCavity is a thin kick (see rf_kick); it has no s-Hamiltonian")

"""
    hamiltonian_formula(e) -> String

The formula of [`hamiltonian`](@ref) for the type of `e` as plain text (for printing).
"""
hamiltonian_formula(::Drift) = "H = energy_term(delta) - ps"
hamiltonian_formula(::Quadrupole) = "H = energy_term(delta) - ps + k1/2 (x^2 - y^2)"
hamiltonian_formula(::Sextupole) = "H = energy_term(delta) - ps + k2/6 (x^3 - 3 x y^2)"
hamiltonian_formula(::SBend) = "H = energy_term(delta) - (1+hx) ps + h x + h^2 x^2/2 + k1/2 (x^2-y^2)"
hamiltonian_formula(::AbstractElement) = "H = 0"

"""
    hamiltonian_field(H) -> Vector

Polynomial vector field ``dz/ds = S\\nabla H`` of a TaylorN Hamiltonian.
"""
function hamiltonian_field(H)
    g = [differentiate(H, j) for j in 1:6]
    return [sum(SYMPLECTIC_S[i, j] * g[j] for j in 1:6) for i in 1:6]
end

"Curvature and normalised strengths (h, k1, k2) of a magnetic element."
magnetic_parameters(e::Drift) = (0.0, 0.0, 0.0)
magnetic_parameters(e::Quadrupole) = (0.0, e.k1, 0.0)
magnetic_parameters(e::Sextupole) = (0.0, 0.0, e.k2)
magnetic_parameters(e::SBend) = (curvature(e), e.k1, 0.0)
magnetic_parameters(e::AbstractElement) = (0.0, 0.0, 0.0)

"d(v/v0)/dδ of the reference particle."
velocity_derivative(d, ref::RefParticle) = (a2 = 1 / ref.betagamma^2; a2 / ((1 + d)^2 + a2)^1.5 / ref.beta)

@doc raw"""
    magnetic_field(z, ref, h, k1, k2)

Exact vector field ``dz/ds = S\nabla H`` of the magnetic Hamiltonian
``H = \mathrm{energy\_term}(\delta) - (1+hx)p_s + hx + h^2x^2/2 + k_1(x^2-y^2)/2 + k_2(x^3-3xy^2)/6``
(no expansion in ``p_x, p_y, \delta``).
"""
function magnetic_field(z, ref, h, k1, k2)
    x, px, y, py, l, d = z
    ps = sqrt((1 + d)^2 - px^2 - py^2)
    hx1 = 1 + h * x
    Hx = -h * ps + h + h^2 * x + k1 * x + k2 / 2 * (x^2 - y^2)
    Hy = -k1 * y - k2 * x * y
    return [hx1 * px / ps, -Hx, hx1 * py / ps, -Hy, velocity_ratio(d, ref) - hx1 * (1 + d) / ps, zero(x)]
end

@doc raw"""
    cavity_field(z, s, ref, L, ω, A, v0, φs)

Exact s-dependent vector field of the long-cavity Hamiltonian
``H = \mathrm{energy\_term}(\delta) - p_s + A\cos\Phi``,
``\Phi = \omega(s - L/2)/v_0 - \omega l / v(\delta) - \varphi_s``.
"""
function cavity_field(z, s, ref, L, omega, A, v0, phis)
    x, px, y, py, l, d = z
    ps = sqrt((1 + d)^2 - px^2 - py^2)
    vr = velocity_ratio(d, ref); v = vr * v0
    Phi = omega * (s - L / 2) / v0 - omega * l / v - phis
    sP = sin(Phi)
    dPhi_dd = omega * l * velocity_derivative(d, ref) * v0 / v^2
    Hl = A * sP * omega / v
    Hd = vr - (1 + d) / ps - A * sP * dPhi_dd
    return [px / ps, zero(x), py / ps, zero(x), Hd, -Hl]
end

end # module
