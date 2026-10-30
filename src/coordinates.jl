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

"""
    C_LIGHT :: Float64

Speed of light in vacuum, 299 792 458 m/s (exact SI value).
Used for v0 = β0·c, revolution frequency and magnetic rigidity Bρ = p0/(Z e c).
"""
const C_LIGHT = 299_792_458.0

# indices of the canonical (and user) coordinates in a 6-vector
const IX, IPX, IY, IPY, IL, IDELTA = 1, 2, 3, 4, 5, 6
const VAR_NAMES  = ["x", "px", "y", "py", "l", "delta"]
const USER_NAMES = ["x", "x'", "y", "y'", "phi", "delta"]
@doc """
    IX, IPX, IY, IPY, IL, IDELTA :: Int  (= 1, 2, 3, 4, 5, 6)

Indices of the coordinates in every 6-vector of Bell (canonical `z` and user `u`):
`z[IX] = x`, `z[IPX] = px`, `z[IY] = y`, `z[IPY] = py`, `z[IL] = l`, `z[IDELTA] = δ`.
The same order is used for the rows/columns of all 6×6 matrices and for
the variables `Z[1..6]` of the polynomial space.
""" IX
@doc "Index of px in a 6-vector (= 2). See [`IX`](@ref)." IPX
@doc "Index of y in a 6-vector (= 3). See [`IX`](@ref)." IY
@doc "Index of py in a 6-vector (= 4). See [`IX`](@ref)." IPY
@doc "Index of l (path-length / time-of-flight variable) in a 6-vector (= 5). See [`IX`](@ref)." IL
@doc "Index of δ = Δp/p0 in a 6-vector (= 6). See [`IX`](@ref)." IDELTA
@doc """
    VAR_NAMES :: Vector{String} = ["x", "px", "y", "py", "l", "delta"]

Printable names of the canonical coordinates (used by `poly_string`, `map_string`).
""" VAR_NAMES
@doc """
    USER_NAMES :: Vector{String} = ["x", "x'", "y", "y'", "phi", "delta"]

Printable names of the user coordinates returned by [`canonical_to_user`](@ref).
""" USER_NAMES

"""
    Particle(name, mass_MeV, charge, G)

Particle species (immutable data).

# Fields
- `name::String` – species name, e.g. `"deuteron"` (used in structure files).
- `mass_MeV::Float64` – rest energy m c² [MeV].
- `charge::Float64` – charge number Z (charge in units of e).
- `G::Float64` – anomalous magnetic moment (reserved for spin work).

Predefined: [`DEUTERON`](@ref); `proton` and `electron` exist in the structure-file parser.
"""
struct Particle
    name::String
    mass_MeV::Float64
    charge::Float64
    G::Float64
end

"""
    DEUTERON :: Particle

Deuteron: rest energy 1875.61294257 MeV (CODATA 2018), Z = 1, G = −0.142987272.
Default particle of `convert_optim` and of structure files with `particle deuteron`.
"""
const DEUTERON = Particle("deuteron", 1875.61294257, 1.0, -0.142987272)

"""
    RefParticle

Reference (design) particle. All normalised quantities of Bell (k1, k2, h,
px = Px/p0, δ = Δp/p0) refer to its momentum p0. Build it with
[`reference_from_kinetic`](@ref); do not fill the fields by hand.

# Fields
- `particle::Particle` – species.
- `mass_MeV::Float64` – rest energy actually used [MeV] (normally `particle.mass_MeV`).
- `Ekin_MeV::Float64` – kinetic energy of the whole ion [MeV].
- `gamma::Float64`, `beta::Float64`, `betagamma::Float64` – γ0, β0, β0γ0.
- `p0_MeV::Float64` – p0·c [MeV].
- `brho_Tm::Float64` – magnetic rigidity Bρ [T·m] = p0/(Z e).

`RefParticle` values are used as cache keys (via `p0_MeV`) by the map builders,
so a lattice at another energy gets new maps automatically.
"""
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

"""
    reference_from_kinetic(p::Particle, Ekin_MeV; mass_MeV = p.mass_MeV) -> RefParticle

Reference particle of species `p` with kinetic energy `Ekin_MeV` [MeV].

# Algorithm
γ0 = 1 + Ekin/m, β0γ0 = √(γ0² − 1), β0 = β0γ0/γ0, p0c = m·β0γ0 [MeV],
Bρ = p0c·10⁶/c/Z [T·m].

# Example
```julia
ref = reference_from_kinetic(DEUTERON, 270.0)
ref.brho_Tm      # ≈ 2.4 T·m
```
"""
function reference_from_kinetic(p::Particle, Ekin_MeV::Real; mass_MeV::Real = p.mass_MeV)
    gamma = 1 + Ekin_MeV / mass_MeV
    bg    = sqrt(gamma^2 - 1)
    beta  = bg / gamma
    p0    = mass_MeV * bg
    brho  = p0 * 1e6 / C_LIGHT / p.charge
    return RefParticle(p, mass_MeV, Ekin_MeV, gamma, beta, bg, p0, brho)
end

@doc raw"""
    energy_term(δ, ref) -> same type as δ

Energy term of the s-Hamiltonian in the canonical pair ``(l, \delta)``:

```math
e(\delta) = \frac{1}{\beta_0}\sqrt{(1+\delta)^2 + \frac{1}{\beta_0^2\gamma_0^2}} - \frac{1}{\beta_0^2}
\qquad (e(0) = 0).
```

Expansion:

```math
e(\delta) = \delta + \frac{\delta^2}{2\gamma_0^2} - \frac{\beta_0^2\delta^3}{2\gamma_0^2}
+ \Big(-\frac{\beta_0^2}{8} + \frac{3\beta_0^4}{4} - \frac{5\beta_0^6}{8}\Big)\delta^4 + \dots
```

and ``\partial e/\partial\delta = v(\delta)/v_0``. Generic: `δ` may be a number or a
`TaylorN` polynomial (this is how every element Hamiltonian gets it).
"""
energy_term(d, ref::RefParticle) = sqrt((1 + d)^2 + 1 / ref.betagamma^2) / ref.beta - 1 / ref.beta^2

@doc raw"""
    velocity_ratio(δ, ref) -> v(δ)/v0

Velocity relative to the reference velocity,

```math
\frac{v}{v_0} = \frac{1+\delta}{\beta_0\sqrt{(1+\delta)^2 + 1/(\beta_0\gamma_0)^2}}
= 1 + \frac{\delta}{\gamma_0^2} - \frac{3\beta_0^2\delta^2}{2\gamma_0^2} + \dots
```

Numbers or `TaylorN`. Used to convert ``l \leftrightarrow`` time of flight
(``l = -v\,\Delta t``) and RF phase.
"""
velocity_ratio(d, ref::RefParticle) = (1 + d) / sqrt((1 + d)^2 + 1 / ref.betagamma^2) / ref.beta

"""
    revolution_frequency(ref, C) -> f0 [Hz]

Revolution frequency f0 = β0 c / C of the reference particle on a ring of
circumference `C` [m].
"""
revolution_frequency(ref::RefParticle, C::Real) = ref.beta * C_LIGHT / C

"""
    canonical_to_user(z, ref, C, h) -> u = [x, x', y, y', φ, δ]

Exact transformation from canonical z = (x, px, y, py, l, δ) to user
coordinates.

# Arguments
- `z` – 6-vector (numbers or `TaylorN`).
- `ref::RefParticle`, `C` – circumference [m], `h` – RF harmonic number.

# Formulas
ps = √((1+δ)² − px² − py²); x' = px/ps, y' = py/ps;
φ = 2π h l / (C · v(δ)/v0) [rad] (φ > 0: particle early).
"""
function canonical_to_user(z, ref::RefParticle, C::Real, h::Real)
    x, px, y, py, l, d = z
    ps  = sqrt((1 + d)^2 - px^2 - py^2)
    phi = 2pi * h * l / (C * velocity_ratio(d, ref))
    return [x, px / ps, y, py / ps, phi, d]
end

"""
    user_to_canonical(u, ref, C, h) -> z = [x, px, y, py, l, δ]

Inverse of [`canonical_to_user`](@ref):
px = x'(1+δ)/√(1+x'²+y'²), py likewise, l = φ C (v/v0)/(2π h).
Used to start tracking from user-specified particles.
"""
function user_to_canonical(u, ref::RefParticle, C::Real, h::Real)
    x, xp, y, yp, phi, d = u
    n = sqrt(1 + xp^2 + yp^2)
    return [x, xp * (1 + d) / n, y, yp * (1 + d) / n, phi * C * velocity_ratio(d, ref) / (2pi * h), d]
end

"""
    describe_reference(ref) -> String

Multi-line human-readable summary of a [`RefParticle`](@ref)
(species, mass, E_kin, γ0, β0, p0c, Bρ, 1/γ0²) for printing in examples.
"""
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
