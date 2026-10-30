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

"""
    synchrotron_tune(ref, eta0, V_MV, h, phi_s; T = 1.0) -> Q_s

Small-amplitude synchrotron tune (oscillations per turn):
Q_s = √(h η0 Z V T cos φs / (2π β0² E)), where E = γ0 m c² [MeV], V [MV],
T is the transit factor. Returns NaN if η0 cos φs ≤ 0 (unstable).
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

"""
    LongitudinalParams(h, eta0, eta1, A, phi_s)
    LongitudinalParams(ref, alpha0, alpha1, V_MV, h, phi_s; T = 1)

Parameters of the averaged longitudinal motion in (φ, δ), with n = turn:

dφ/dn = −2π h (η0 δ + η1 δ²),  dδ/dn = A (sin(φs + φ) − sin φs).

# Fields
- `h::Int` – harmonic number.
- `eta0`, `eta1` – slip factors η0 = α0 − 1/γ², η1 = α1 − α0/γ² + 3β²/(2γ²) + 1/γ⁴.
- `A` – per-turn δ-kick amplitude Z V T/(β² E) [1].
- `phi_s` – synchronous phase [rad].

The second constructor computes η0, η1, A from momentum compaction α0, α1 at the
energy of `ref`.
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

"""
    long_H(p, φ, δ) -> Float64

Longitudinal Hamiltonian (invariant of the smooth motion):
H = −2π h (η0 δ²/2 + η1 δ³/3) + A (cos(φs + φ) + φ sin φs).
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
