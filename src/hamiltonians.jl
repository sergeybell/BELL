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

export hamiltonian, hamiltonian_formula, hamiltonian_field

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

end # module
