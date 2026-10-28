# =============================================================================
# Constants.jl - physical constants, reference particle, phase-space coordinates
# =============================================================================
#
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# WHAT
#   * Physical constants (speed of light) and the Particle type (deuteron).
#   * RefParticle: everything derived from the reference momentum
#     (gamma0, beta0, p0, magnetic rigidity B*rho).
#   * The definition of the 6D coordinates used everywhere in the program and
#     the exact transformation between the CANONICAL variables (used inside
#     all Hamiltonians and Lie maps) and the USER variables (printed / input).
#
# WHY
#   Every other module must agree on one set of coordinates and sign
#   conventions; putting them in one file is the only way to avoid a silent
#   sign/convention mix-up between elements.
#
# COORDINATES (independent variable: arc length s along the reference orbit)
#
#   USER coordinates  u = (x, x', y, y', phi, delta)
#     x, y   [m]   transverse offsets (x radially OUTWARD in bends, y vertical)
#     x', y' [rad] slopes relative to the local reference direction,
#                  defined exactly as x' = px/ps, y' = py/ps (see below);
#                  outside bends this is exactly dx/ds, dy/ds.
#     phi    [rad] RF phase of arrival relative to the reference particle:
#                      phi = -h * omega0 * dt ,   dt = t - t_ref(s)
#                  omega0 = 2*pi*beta0*c/C (revolution angular frequency),
#                  h = harmonic number (h = 1 if no RF is installed).
#                  SIGN: a particle arriving LATER than the reference
#                  (dt > 0) has phi < 0;  an EARLY particle has phi > 0.
#     delta  [1]   relative MOMENTUM deviation  delta = (p - p0)/p0
#                  (momentum, not energy).
#
#   CANONICAL coordinates  z = (x, px, y, py, l, delta)
#     px, py = Px/p0, Py/p0           (normalised canonical momenta)
#     l      = -v(delta) * dt   [m]   conjugate to delta
#     Derivation (done in PLAN.md): start from (t, -E) as the longitudinal
#     pair of the s-Hamiltonian, use the type-2 generating function
#        F2(t, delta; s) = -(t - s/v0) * E(delta) / p0 .
#     Then  l = dF2/d(delta) = -(t - s/v0) * v(delta) and the new Hamiltonian
#     gains  dF2/ds = E(delta)/(p0 v0), i.e. the "energy term"
#        e(delta) = sqrt((1+delta)^2 + 1/(beta0 gamma0)^2)/beta0 - 1/beta0^2
#     (constant dropped).  Because l is conjugate to delta EXACTLY, delta is
#     the true dp/p, as requested.
#
#   Relation  phi <-> l :   phi = h*omega0*l/v = 2*pi*h*l / (C * v/v0)
#   Relation  x' <-> px :   x' = px/ps,  ps = sqrt((1+delta)^2 - px^2 - py^2)
#   Both are exact; at linear order x' = px and phi = 2*pi*h*l/C.
#
# WHAT YOU SEE
#   Nothing is printed by this file; scripts print the reference particle.
# =============================================================================
module Coordinates

export C_LIGHT, Particle, DEUTERON, RefParticle, reference_from_kinetic,
       IX, IPX, IY, IPY, IL, IDELTA, VAR_NAMES, USER_NAMES,
       energy_term, velocity_ratio, canonical_to_user, user_to_canonical,
       revolution_frequency, describe_reference

"Speed of light in vacuum [m/s] (exact)."
const C_LIGHT = 299_792_458.0

# indices of the canonical (and user) coordinates in a 6-vector
const IX, IPX, IY, IPY, IL, IDELTA = 1, 2, 3, 4, 5, 6
const VAR_NAMES  = ["x", "px", "y", "py", "l", "delta"]
const USER_NAMES = ["x", "x'", "y", "y'", "phi", "delta"]

"A particle species: rest energy [MeV], charge number Z, anomalous magnetic moment G."
struct Particle
    name::String
    mass_MeV::Float64
    charge::Float64
    G::Float64
end

"Deuteron (CODATA 2018 rest energy; G = -0.142987 for later spin work)."
const DEUTERON = Particle("deuteron", 1875.61294257, 1.0, -0.142987272)

"Reference particle: all quantities derived from the design momentum p0."
struct RefParticle
    particle::Particle
    mass_MeV::Float64     # rest energy actually used (may differ from particle.mass_MeV, e.g. OptiM 2*938.272)
    Ekin_MeV::Float64     # kinetic energy of the whole ion
    gamma::Float64
    beta::Float64
    betagamma::Float64
    p0_MeV::Float64       # p0*c [MeV]
    brho_Tm::Float64      # magnetic rigidity [T m]
end

"Build the reference particle from its kinetic energy [MeV] (optionally overriding the mass)."
function reference_from_kinetic(p::Particle, Ekin_MeV::Real; mass_MeV::Real = p.mass_MeV)
    gamma = 1 + Ekin_MeV / mass_MeV
    bg    = sqrt(gamma^2 - 1)
    beta  = bg / gamma
    p0    = mass_MeV * bg
    brho  = p0 * 1e6 / C_LIGHT / p.charge
    return RefParticle(p, mass_MeV, Ekin_MeV, gamma, beta, bg, p0, brho)
end

"""
    energy_term(delta, ref)

e(delta) = E(delta)/(p0 v0) - E0/(p0 v0)
         = sqrt((1+delta)^2 + 1/(beta0 gamma0)^2)/beta0 - 1/beta0^2.
Its expansion is  delta + delta^2/(2 gamma0^2) - delta^3/(2 beta0^2 gamma0^2) + ...
Works for numbers and TaylorN polynomials.
"""
energy_term(d, ref::RefParticle) = sqrt((1 + d)^2 + 1 / ref.betagamma^2) / ref.beta - 1 / ref.beta^2

"v(delta)/v0 for a particle of relative momentum deviation delta."
velocity_ratio(d, ref::RefParticle) = (1 + d) / sqrt((1 + d)^2 + 1 / ref.betagamma^2) / ref.beta

"Revolution frequency f0 = beta0 c / C [Hz]."
revolution_frequency(ref::RefParticle, C::Real) = ref.beta * C_LIGHT / C

"""
    canonical_to_user(z, ref, C, h) -> (x, x', y, y', phi, delta)

Exact transformation (works on numbers and on TaylorN polynomials).
"""
function canonical_to_user(z, ref::RefParticle, C::Real, h::Real)
    x, px, y, py, l, d = z
    ps  = sqrt((1 + d)^2 - px^2 - py^2)
    phi = 2pi * h * l / (C * velocity_ratio(d, ref))
    return [x, px / ps, y, py / ps, phi, d]
end

"""
    user_to_canonical(u, ref, C, h) -> (x, px, y, py, l, delta)

Inverse of canonical_to_user:  px = x'(1+delta)/sqrt(1+x'^2+y'^2),
l = phi * C * (v/v0) / (2 pi h).
"""
function user_to_canonical(u, ref::RefParticle, C::Real, h::Real)
    x, xp, y, yp, phi, d = u
    n = sqrt(1 + xp^2 + yp^2)
    return [x, xp * (1 + d) / n, y, yp * (1 + d) / n, phi * C * velocity_ratio(d, ref) / (2pi * h), d]
end

"Multi-line human readable description of the reference particle."
function describe_reference(ref::RefParticle)
    return """
      particle            : $(ref.particle.name) (Z = $(ref.particle.charge))
      rest energy used    : $(ref.mass_MeV) MeV
      kinetic energy      : $(ref.Ekin_MeV) MeV
      gamma0, beta0       : $(round(ref.gamma, digits=8)), $(round(ref.beta, digits=8))
      p0*c                : $(round(ref.p0_MeV, digits=5)) MeV
      B*rho               : $(round(ref.brho_Tm, digits=7)) T*m
      1/gamma0^2          : $(round(1/ref.gamma^2, digits=8))"""
end

end # module
