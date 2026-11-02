# =============================================================================
# smooth.jl - averaged turn-by-turn motion in (phi, delta)
#
# NOTE (v0.17): approximate / expanded longitudinal model, kept for comparison.
# The reference 6D RF model is the RFCavity s-Hamiltonian (src/lie/hamiltonians.jl)
# and its exact thin limit rf_kick (src/da/field_track.jl).
# =============================================================================
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# Smooth is a separate averaged equation in (φ, δ). It is not 6D RF tracking
# and not a replacement for thin rf_kick. Acceleration uses FieldTrack +
# RFCavity. Higher order of this model is further η from map readings.
# =============================================================================
module Smooth

using ..Coordinates

export longitudinal_hamiltonian_text, synchrotron_tune,
       LongitudinalParams, long_H, fixed_points, separatrix, bucket_height

@doc raw"""
    synchrotron_tune(ref, eta0, V_MV, h, phi_s; T = 1.0) -> Q_s

Small-amplitude synchrotron tune of the smooth model (analysis only):
``Q_s = \sqrt{\dfrac{h\,\eta_0\,ZVT\cos\varphi_s}{2\pi\beta_0^2E}}``.
"""
function synchrotron_tune(ref::RefParticle, eta0, V_MV, h, phi_s; T = 1.0)
    E = ref.gamma * ref.mass_MeV
    a = h * eta0 * ref.particle.charge * V_MV * T * cos(phi_s) / (2pi * ref.beta^2 * E)
    return a > 0 ? sqrt(a) : NaN
end

"""
    longitudinal_hamiltonian_text() -> String

Plain-text statement of the smooth turn-by-turn model.
"""
longitudinal_hamiltonian_text() = """
    smooth turn-by-turn motion (n = turn number, phi = -h w0 dt, delta = dp/p):
      dphi/dn   = -2 pi h (eta0 delta + eta1 delta^2)
      ddelta/dn = (Z V T/(beta^2 E)) (sin(phi_s + phi) - sin(phi_s))
      H(phi,delta) = -2 pi h (eta0 delta^2/2 + eta1 delta^3/3)
                     + (Z V T/(beta^2 E)) (cos(phi_s + phi) + phi sin(phi_s))
      small amplitude:  Q_s = sqrt( h eta0 Z V T cos(phi_s) / (2 pi beta^2 E) ),
      stable iff eta0*cos(phi_s) > 0.  eta1 adds a second fixed point at
      delta ~ -eta0/eta1 (alpha-bucket) and makes the bucket asymmetric.
      This is not an element Hamiltonian and does not track the ring in 6D."""

@doc raw"""
    LongitudinalParams(h, eta0, eta1, A, phi_s)
    LongitudinalParams(ref, alpha0, alpha1, V_MV, h, phi_s; T = 1)

Parameters of the averaged (smooth) model in ``(\varphi, \delta)``, ``n`` = turn.
Truncated at ``\eta_1``; further η would be the same model with more terms.
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

long_H(p::LongitudinalParams, phi, d) =
    -2pi * p.h * (p.eta0 * d^2 / 2 + p.eta1 * d^3 / 3) + p.A * (cos(p.phi_s + phi) + phi * sin(p.phi_s))

function fixed_points(p::LongitudinalParams)
    phu = pi - 2p.phi_s
    return (stable = 0.0, unstable = phu, delta2 = -p.eta0 / p.eta1)
end

function separatrix(p::LongitudinalParams)
    Hs = long_H(p, fixed_points(p).unstable, 0.0)
    H0 = long_H(p, 0.0, 0.0)
    inside(phi, d) = sign(long_H(p, phi, d) - Hs) == sign(H0 - Hs)
    return Hs, inside
end

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
