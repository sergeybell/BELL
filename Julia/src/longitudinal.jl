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

"Linear synchrotron tune; returns NaN if unstable (eta0 cos(phi_s) <= 0)."
function synchrotron_tune(ref::RefParticle, eta0, V_MV, h, phi_s; T = 1.0)
    E = ref.gamma * ref.mass_MeV
    a = h * eta0 * ref.particle.charge * V_MV * T * cos(phi_s) / (2pi * ref.beta^2 * E)
    return a > 0 ? sqrt(a) : NaN
end

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
LongitudinalParams(ref, alpha0, alpha1, V_MV, h, phi_s; T = 1)
eta0 = alpha0 - 1/g^2, eta1 = alpha1 - alpha0/g^2 + 1.5 b^2/g^2 + 1/g^4 at the
energy of `ref`; A = Z V T/(beta^2 E) is the per-turn delta-kick amplitude.
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

"H(phi, delta) = -2 pi h (eta0 d^2/2 + eta1 d^3/3) + A (cos(phi_s + phi) + phi sin(phi_s))."
long_H(p::LongitudinalParams, phi, d) =
    -2pi * p.h * (p.eta0 * d^2 / 2 + p.eta1 * d^3 / 3) + p.A * (cos(p.phi_s + phi) + phi * sin(p.phi_s))

function fixed_points(p::LongitudinalParams)
    phu = pi - 2p.phi_s
    return (stable = 0.0, unstable = phu, delta2 = -p.eta0 / p.eta1)
end

"Separatrix level H_sep = H(phi_u, 0) and a function inside(phi, d) for the bucket."
function separatrix(p::LongitudinalParams)
    Hs = long_H(p, fixed_points(p).unstable, 0.0)
    H0 = long_H(p, 0.0, 0.0)
    inside(phi, d) = sign(long_H(p, phi, d) - Hs) == sign(H0 - Hs)
    return Hs, inside
end

"Bucket half-height in delta at phi = 0 (root of H(0, d) = H_sep, d > 0 and d < 0)."
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
