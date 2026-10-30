# =============================================================================
# longitudinal.jl - averaged turn-by-turn motion in (phi, delta)
# =============================================================================
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# WHAT
#   Smooth longitudinal Hamiltonian from already-read slip coefficients:
#     dphi/dn   = -2 pi h (eta0 delta + eta1 delta^2)
#     ddelta/dn = (Z V T /(beta^2 E)) (sin(phi_s + phi) - sin phi_s)
#   Synchrotron tune, fixed points, separatrix, bucket height.
#   This file does not track the ring in 6D.  eta comes from the higher-order
#   reading of the one-turn Lie factor.
#
# =============================================================================
module Longitudinal

using ..Coordinates

export longitudinal_hamiltonian_text, synchrotron_tune,
       LongitudinalParams, long_H, fixed_points, separatrix, bucket_height

@doc raw"""
    synchrotron_tune(ref, eta0, V_MV, h, phi_s; T = 1.0) -> Q_s

Small-amplitude synchrotron tune of the smooth model (analysis only):
``Q_s = \sqrt{\dfrac{h\,\eta_0\,ZVT\cos\varphi_s}{2\pi\beta_0^2E}}``, ``E = \gamma_0mc^2`` [MeV], V [MV],
``T`` the transit factor. Returns NaN if ``\eta_0\cos\varphi_s \le 0`` (unstable).
"""
function synchrotron_tune(ref::RefParticle, eta0, V_MV, h, phi_s; T = 1.0)
    E = ref.gamma * ref.mass_MeV
    a = h * eta0 * ref.particle.charge * V_MV * T * cos(phi_s) / (2pi * ref.beta^2 * E)
    return a > 0 ? sqrt(a) : NaN
end

"""
    longitudinal_hamiltonian_text() -> String

Plain-text statement of the smooth turn-by-turn longitudinal model (equations,
H(φ, δ), Q_s, stability condition) for printing in examples.
"""
longitudinal_hamiltonian_text() = """
    turn-by-turn longitudinal motion (n = turn number, phi = -h w0 dt, delta = dp/p):
      dphi/dn   = -2 pi h (eta0 delta + eta1 delta^2)
      ddelta/dn = (Z V T/(beta^2 E)) (sin(phi_s + phi) - sin(phi_s))
      H(phi,delta) = -2 pi h (eta0 delta^2/2 + eta1 delta^3/3)
                     + (Z V T/(beta^2 E)) (cos(phi_s + phi) + phi sin(phi_s))
      small amplitude:  Q_s = sqrt( h eta0 Z V T cos(phi_s) / (2 pi beta^2 E) ),
      stable iff eta0*cos(phi_s) > 0.  eta1 adds a second fixed point at
      delta ~ -eta0/eta1 (alpha-bucket) and makes the bucket asymmetric."""

@doc raw"""
    LongitudinalParams(h, eta0, eta1, A, phi_s)
    LongitudinalParams(ref, alpha0, alpha1, V_MV, h, phi_s; T = 1)

Parameters of the averaged (smooth) longitudinal model in ``(\varphi, \delta)``, ``n`` = turn.
This model is an analysis object, truncated at ``\delta^2``; no tracking routine uses it.

```math
\frac{d\varphi}{dn} = -2\pi h(\eta_0\delta + \eta_1\delta^2),\qquad
\frac{d\delta}{dn} = A\big(\sin(\varphi_s+\varphi) - \sin\varphi_s\big).
```

# Fields
- `h::Int` – harmonic number.
- `eta0`, `eta1` – ``\eta_0 = \alpha_0 - 1/\gamma^2``,
  ``\eta_1 = \alpha_1 - \alpha_0/\gamma^2 + 3\beta^2/(2\gamma^2) + 1/\gamma^4``.
- `A` – per-turn ``\delta``-kick amplitude ``ZVT/(\beta^2E)`` [1].
- `phi_s` – synchronous phase [rad].

The second constructor computes ``\eta_0, \eta_1, A`` from the momentum compaction
``\alpha_0, \alpha_1`` at the energy of `ref`.
"""
struct LongitudinalParams
    h::Int; eta0::Float64; eta1::Float64; A::Float64; phi_s::Float64
end
function LongitudinalParams(ref::RefParticle, alpha0, alpha1, V_MV, h, phi_s; T = 1.0)
    g2 = ref.gamma^2; b2 = ref.beta^2
    eta0 = alpha0 - 1 / g2
    eta1 = alpha1 - alpha0 / g2 + 1.5 * b2 / g2 + 1 / g2^2
    A = ref.particle.charge * V_MV * T / (b2 * ref.gamma * ref.mass_MeV)
    return LongitudinalParams(h, eta0, eta1, A, phi_s)
end

@doc raw"""
    long_H(p, φ, δ) -> Float64

Longitudinal Hamiltonian (invariant of the smooth motion):
``H = -2\pi h\big(\eta_0\delta^2/2 + \eta_1\delta^3/3\big) + A\big(\cos(\varphi_s+\varphi) + \varphi\sin\varphi_s\big)``.
"""
long_H(p::LongitudinalParams, phi, d) =
    -2pi * p.h * (p.eta0 * d^2 / 2 + p.eta1 * d^3 / 3) + p.A * (cos(p.phi_s + phi) + phi * sin(p.phi_s))

"""
    fixed_points(p) -> (stable, unstable, delta2)

Fixed points of [`long_H`](@ref): stable at φ = 0, unstable at φ_u = π − 2φs
(both at δ = 0), and the second (α-bucket) momentum δ₂ = −η0/η1.
"""
function fixed_points(p::LongitudinalParams)
    phu = pi - 2p.phi_s
    return (stable = 0.0, unstable = phu, delta2 = -p.eta0 / p.eta1)
end

"""
    separatrix(p) -> (H_sep, inside)

Separatrix level H_sep = H(φ_u, 0) and a predicate `inside(φ, δ)` that is true
when the point lies on the same side of H_sep as the stable fixed point.
"""
function separatrix(p::LongitudinalParams)
    Hs = long_H(p, fixed_points(p).unstable, 0.0)
    H0 = long_H(p, 0.0, 0.0)
    inside(phi, d) = sign(long_H(p, phi, d) - Hs) == sign(H0 - Hs)
    return Hs, inside
end

"""
    bucket_height(p) -> (δ_plus, δ_minus)

Bucket half-heights at φ = 0: the roots of H(0, δ) = H_sep for δ > 0 and δ < 0,
found by bisection (200 steps) in |δ| < |η0/η1|/2, so the α-bucket is excluded.
η1 makes the bucket asymmetric.
"""
function bucket_height(p::LongitudinalParams)
    Hs, _ = separatrix(p)
    f(d) = long_H(p, 0.0, d) - Hs
    function root(a, b)
        for _ in 1:200
            m = (a + b) / 2
            sign(f(m)) == sign(f(a)) ? (a = m) : (b = m)
        end
        (a + b) / 2
    end
    dmax = 0.5 * abs(p.eta0 / p.eta1)
    return (root(0.0, dmax), root(0.0, -dmax))
end

end # module
