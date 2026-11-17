# =============================================================================
# coordinates.jl - physical constants, particle species, reference particle,
#                  canonical and laboratory coordinates
# =============================================================================
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# WHAT
#   Particle species (name, rest energy, charge number, anomalous moment).
#   RefParticle: quantities derived from the kinetic energy of one such
#   species (gamma, beta, p0, magnetic rigidity).
#   Canonical coordinates z = (x, px, y, py, l, delta) used in every
#   Hamiltonian, and laboratory coordinates (x, x', y, y', phi, delta)
#   used as input and output.
#
# =============================================================================
module Coordinates

export c_light, Particle, deuteron, proton, electron, RefParticle, reference_from_kinetic,
       index_x, index_px, index_y, index_py, index_l, index_delta, VAR_NAMES, laboratory_names,
       energy_term, velocity_ratio, canonical_to_laboratory, laboratory_to_canonical,
       revolution_frequency, describe_reference

"""
    c_light :: Float64

Speed of light in vacuum, 299 792 458 m/s (exact SI value).
Used for v0 = β0·c, revolution frequency and magnetic rigidity Bρ = p0/(Z e c).
"""
const c_light = 299_792_458.0

const index_x, index_px, index_y, index_py, index_l, index_delta = 1, 2, 3, 4, 5, 6
const VAR_NAMES = ["x", "px", "y", "py", "l", "delta"]
const laboratory_names = ["x", "x'", "y", "y'", "phi", "delta"]
@doc """
    index_x, index_px, index_y, index_py, index_l, index_delta :: Int  (= 1, 2, 3, 4, 5, 6)

Indices of the coordinates in every 6-vector of Bell (canonical `z` and laboratory `u`):
`z[index_x] = x`, `z[index_px] = px`, `z[index_y] = y`, `z[index_py] = py`,
`z[index_l] = l`, `z[index_delta] = δ`.
The same order is used for the rows/columns of all 6×6 matrices and for
the variables `Z[1..6]` of the polynomial space.
""" index_x
@doc "Index of px in a 6-vector (= 2). See [`index_x`](@ref)." index_px
@doc "Index of y in a 6-vector (= 3). See [`index_x`](@ref)." index_y
@doc "Index of py in a 6-vector (= 4). See [`index_x`](@ref)." index_py
@doc "Index of l (path-length / time-of-flight variable) in a 6-vector (= 5). See [`index_x`](@ref)." index_l
@doc "Index of δ = Δp/p0 in a 6-vector (= 6). See [`index_x`](@ref)." index_delta
@doc """
    VAR_NAMES :: Vector{String} = ["x", "px", "y", "py", "l", "delta"]

Printable names of the canonical coordinates (used by `poly_string`, `map_string`).
""" VAR_NAMES
@doc """
    laboratory_names :: Vector{String} = ["x", "x'", "y", "y'", "phi", "delta"]

Printable names of the laboratory coordinates returned by [`canonical_to_laboratory`](@ref).
""" laboratory_names

"""
    Particle(name, mass_MeV, charge, G)

Particle species (immutable data).

# Fields
- `name::String` – species name, e.g. `"deuteron"` (used in structure files).
- `mass_MeV::Float64` – rest energy m c² [MeV].
- `charge::Float64` – charge number Z (charge in units of e).
- `G::Float64` – anomalous magnetic moment (reserved for spin work).

Predefined: [`deuteron`](@ref), [`proton`](@ref), [`electron`](@ref).
"""
struct Particle
    name::String
    mass_MeV::Float64
    charge::Float64
    G::Float64
end

"""
    deuteron :: Particle

Deuteron: rest energy 1875.61294257 MeV (CODATA 2018), Z = 1, G = −0.142987272.
Default particle of `convert_optim` and of structure files with `particle deuteron`.
"""
const deuteron = Particle("deuteron", 1875.61294257, 1.0, -0.142987272)

"""
    proton :: Particle

Proton: rest energy 938.27208816 MeV (CODATA 2018), Z = 1, G = 1.7928473443.
"""
const proton = Particle("proton", 938.27208816, 1.0, 1.7928473443)

"""
    electron :: Particle

Electron: rest energy 0.51099895000 MeV, Z = −1, G = 0.001159652181.
"""
const electron = Particle("electron", 0.51099895000, -1.0, 0.001159652181)

"""
    RefParticle

Reference (design) particle. All normalised quantities of Bell (k1, k2, h,
px = Px/p0, δ = Δp/p0) refer to its momentum p0. The species is a [`Particle`](@ref);
`RefParticle` stores that species, the kinetic energy, and the quantities
computed from it. Build it with [`reference_from_kinetic`](@ref); do not fill
the fields by hand.

# Fields
- `particle::Particle` – species.
- `mass_MeV::Float64` – rest energy actually used [MeV] (normally `particle.mass_MeV`).
- `Ekin_MeV::Float64` – kinetic energy of the whole ion [MeV].
- `gamma::Float64`, `beta::Float64`, `betagamma::Float64` – γ0, β0, β0γ0.
- `p0_MeV::Float64` – p0·c [MeV].
- `brho_Tm::Float64` – magnetic rigidity Bρ [T·m] = p0/(Z e).
"""
struct RefParticle
    particle::Particle
    mass_MeV::Float64
    Ekin_MeV::Float64
    gamma::Float64
    beta::Float64
    betagamma::Float64
    p0_MeV::Float64
    brho_Tm::Float64
end

"""
    reference_from_kinetic(p::Particle, Ekin_MeV; mass_MeV = p.mass_MeV) -> RefParticle

Reference particle of species `p` with kinetic energy `Ekin_MeV` [MeV].
We study the reference particle: γ0 = 1 + Ekin/m, β0γ0 = √(γ0² − 1),
β0 = β0γ0/γ0, p0c = m·β0γ0 [MeV], Bρ = p0c·10⁶/c/Z [T·m].

# Example
```julia
ref = reference_from_kinetic(deuteron, 270.0)
ref.brho_Tm      # ≈ 2.4 T·m
```
"""
function reference_from_kinetic(p::Particle, Ekin_MeV::Real; mass_MeV::Real = p.mass_MeV)
    gamma = 1 + Ekin_MeV / mass_MeV
    bg    = sqrt(gamma^2 - 1)
    beta  = bg / gamma
    p0    = mass_MeV * bg
    brho  = p0 * 1e6 / c_light / p.charge
    return RefParticle(p, mass_MeV, Ekin_MeV, gamma, beta, bg, p0, brho)
end

@doc raw"""
    energy_term(δ, ref) -> same type as δ

Energy term of the s-Hamiltonian in the canonical pair ``(l, \delta)``.
This is not the kinetic energy (that is `Ekin_MeV` of the reference). On the
reference orbit ``\mathrm{energy\_term}(0) = 0``.

```math
\mathrm{energy\_term}(\delta) =
\frac{1}{\beta_0}\sqrt{(1+\delta)^2 + \frac{1}{\beta_0^2\gamma_0^2}} - \frac{1}{\beta_0^2}.
```

and ``\partial(\mathrm{energy\_term})/\partial\delta = v(\delta)/v_0``. Generic: `δ` may be a number or a
`TaylorN` polynomial (this is how every element Hamiltonian gets it).
"""
energy_term(d, ref::RefParticle) = sqrt((1 + d)^2 + 1 / ref.betagamma^2) / ref.beta - 1 / ref.beta^2

@doc raw"""
    velocity_ratio(δ, ref) -> v(δ)/v0

Velocity relative to the reference velocity,

```math
\frac{v}{v_0} = \frac{1+\delta}{\beta_0\sqrt{(1+\delta)^2 + 1/(\beta_0\gamma_0)^2}}.
```

Numbers or `TaylorN`. Used to convert ``l \leftrightarrow`` time of flight
(``l = -v\,\Delta t``) and RF phase.
"""
velocity_ratio(d, ref::RefParticle) = (1 + d) / sqrt((1 + d)^2 + 1 / ref.betagamma^2) / ref.beta

"""
    revolution_frequency(ref, path_length) -> f0 [Hz]

Revolution frequency f0 = β0 c / path_length of the reference particle on a
line of length `path_length` [m] (the sum of element lengths; not necessarily a closed ring).
"""
revolution_frequency(ref::RefParticle, path_length::Real) = ref.beta * c_light / path_length

"""
    canonical_to_laboratory(z, ref, path_length, h) -> u = [x, x', y, y', φ, δ]

Exact transformation from canonical z = (x, px, y, py, l, δ) to laboratory
coordinates (x, x', y, y', φ, δ).

# Arguments
- `z` – 6-vector (numbers or `TaylorN`).
- `ref::RefParticle`, `path_length` – length of the reference orbit [m], `h` – RF harmonic number.

# Formulas
ps = √((1+δ)² − px² − py²); x' = px/ps, y' = py/ps;
φ = 2π h l / (path_length · v(δ)/v0) [rad] (φ > 0: particle early).
"""
function canonical_to_laboratory(z, ref::RefParticle, path_length::Real, h::Real)
    x, px, y, py, l, d = z
    ps  = sqrt((1 + d)^2 - px^2 - py^2)
    phi = 2pi * h * l / (path_length * velocity_ratio(d, ref))
    return [x, px / ps, y, py / ps, phi, d]
end

"""
    laboratory_to_canonical(u, ref, path_length, h) -> z = [x, px, y, py, l, δ]

Inverse of [`canonical_to_laboratory`](@ref):
px = x'(1+δ)/√(1+x'²+y'²), py likewise, l = φ path_length (v/v0)/(2π h).
Used to start tracking from laboratory-specified particles.
"""
function laboratory_to_canonical(u, ref::RefParticle, path_length::Real, h::Real)
    x, xp, y, yp, phi, d = u
    n = sqrt(1 + xp^2 + yp^2)
    return [x, xp * (1 + d) / n, y, yp * (1 + d) / n, phi * path_length * velocity_ratio(d, ref) / (2pi * h), d]
end

"""
    describe_reference(ref) -> String

We study the reference particle. Multi-line summary of a [`RefParticle`](@ref)
(species, mass, E_kin, γ0, β0, p0c, Bρ, 1/γ0²).
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
