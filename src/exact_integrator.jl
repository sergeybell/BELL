# =============================================================================
# ExactIntegrator.jl - symplectic integration of the EXACT (unexpanded) H
# =============================================================================
#
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# WHAT
#   Element-by-element integration of Hamilton's equations dz/ds = S grad H
#   with the full, non-truncated Hamiltonians of Hamiltonians.jl:
#     H_mag = e(delta) - (1+hx) ps + h x + h^2 x^2/2 + k1/2 (x^2-y^2)
#             + k2/6 (x^3 - 3 x y^2)
#   (drift, quad, sextupole, sector bend are special cases) and, for the
#   finite-length cavity, H_cav = e(delta) - ps + A cos(Phi(l, delta, s)):
#       Phi = omega (s - L/2)/v0 - omega l / v(delta) - phi_s ,
#       A   = (Z V/(p0 c)) * c/(omega L),   omega = h * omega0 = h 2 pi v0 / C.
#   Derivation of H_cav: on-axis field E_s = -dA_s/dt with
#   q A_s/p0 = -A cos(omega t + psi); with t = s/v0 - l/v (definition of l)
#   the term -q A_s/p0 of the s-Hamiltonian becomes +A cos(Phi).  Then
#       d delta/ds = -dH/dl = -(Z V/(p0 c L)) (c/v) sin(Phi)
#   which is exactly dp/ds = dE/(v ds) with dE/ds = -(Z V/L) sin(Phi); the
#   energy gain over the cavity is Z V T sin(phi_s + phi), T = sin(th/2)/(th/2),
#   th = omega L / v0 (transit-time factor).  Transverse RF fields are
#   neglected (pure longitudinal model - not Maxwell-complete off axis).
#
#   Integrator: implicit midpoint rule (symplectic for any, also
#   s-dependent, H) composed as Yoshida's 4th-order triple jump; the implicit
#   equation is solved by fixed-point iteration (to 1e-16).  Generic code:
#   runs on Float64 and on TaylorN (for exact Jacobians / symplecticity).
#
#   exact_track_element(z, e, ref, C), exact_track_turn(z, lat, ref, C)
#   renormalize(zs, ref, isync): after a cavity, make the synchronous particle
#     the new reference: p0' = p0 (1 + delta_s),  delta' = (1+delta)/(1+delta_s)-1,
#     px' = px/(1+delta_s), py' likewise (this IS adiabatic damping),
#     dt' = dt - dt_s,  l' = -v dt'.  Fields are assumed to ramp with p0, so
#     normalised strengths (k1, k2, h) stay fixed.
#
# WHY
#   (1) independent check of the truncated Lie / Taylor maps, (2) tracking
#   at ANY energy without rebuilding maps (needed for acceleration),
#   (3) the RF cavity is s- and time-dependent, which the s-independent Lie
#   maps do not cover.
#
# WHAT YOU SEE
#   Nothing printed here.
# =============================================================================
module ExactIntegrator

using ..Coordinates, ..Elements, ..Structure

export exact_track_element, exact_track_turn, renormalize, cavity_parameters, transit_factor

"""
Yoshida triple-jump weights: W1 = 1/(2 − 2^(1/3)), W0 = 1 − 2 W1. Composing
the 2nd-order implicit midpoint step as (W1 h, W0 h, W1 h) gives a 4th-order
symplectic method.
"""
const W1 = 1 / (2 - 2^(1 / 3))
const W0 = 1 - 2W1

_vel(d, ref) = velocity_ratio(d, ref)                       # v/v0
_dvel(d, ref) = (a2 = 1 / ref.betagamma^2; a2 / ((1 + d)^2 + a2)^1.5 / ref.beta)   # d(v/v0)/ddelta

@doc raw"""
    cavity_parameters(cav::Cavity, ref, C) -> (ω, A, v0)

Constants of the finite-cavity Hamiltonian, recomputed from the current `ref`
at every call (so they follow the energy ramp):

- ``\omega = h\,2\pi v_0/C`` [rad/s];
- ``A = \dfrac{ZV}{p_0c}\,\dfrac{c}{\omega L}`` [1];
- ``v_0 = \beta_0 c`` [m/s].
"""
function cavity_parameters(cav::Cavity, ref::RefParticle, C::Real)
    v0 = ref.beta * C_LIGHT
    omega = cav.harmonic * 2pi * v0 / C
    A = ref.particle.charge * cav.voltage_MV / ref.p0_MeV * C_LIGHT / (omega * cav.L)
    return omega, A, v0
end

@doc raw"""
    transit_factor(cav, ref, C) -> T

Transit-time factor of a finite cavity: ``T = \dfrac{\sin(\theta/2)}{\theta/2}``,
``\theta = 2\pi hL/C = \omega L/v_0``. The energy gain per passage is ``ZVT\sin(\varphi_s + \varphi)``.
"""
transit_factor(cav::Cavity, ref::RefParticle, C::Real) =
    (th = 2pi * cav.harmonic * cav.L / C; sin(th / 2) / (th / 2))

# vector field dz/ds = S grad H for magnetic elements
@doc raw"""
    _field_mag(z, ref, h, k1, k2)

Exact vector field ``dz/ds = S\nabla H`` of the magnetic Hamiltonian
``H = e(\delta) - (1+hx)p_s + hx + h^2x^2/2 + k_1(x^2 - y^2)/2 + k_2(x^3 - 3xy^2)/6``
(no expansion in ``p_x, p_y, \delta``):

``x' = (1+hx)p_x/p_s``, ``p_x' = hp_s - h - h^2x - k_1x - \tfrac{k_2}{2}(x^2-y^2)``,
``y' = (1+hx)p_y/p_s``, ``p_y' = k_1y + k_2xy``,
``l' = v/v_0 - (1+hx)(1+\delta)/p_s``, ``\delta' = 0``.
"""
function _field_mag(z, ref, h, k1, k2)
    x, px, y, py, l, d = z
    ps = sqrt((1 + d)^2 - px^2 - py^2)
    hx1 = 1 + h * x
    Hx = -h * ps + h + h^2 * x + k1 * x + k2 / 2 * (x^2 - y^2)
    Hy = -k1 * y - k2 * x * y
    return [hx1 * px / ps, -Hx, hx1 * py / ps, -Hy, _vel(d, ref) - hx1 * (1 + d) / ps, zero(x)]
end

@doc raw"""
    _field_cav(z, s, ref, L, ω, A, v0, φs)

Exact s-dependent vector field of the cavity Hamiltonian
``H = e(\delta) - p_s + A\cos\Phi``, ``\Phi = \omega(s - L/2)/v_0 - \omega l/v(\delta) - \varphi_s``:

``x' = p_x/p_s``, ``y' = p_y/p_s``, ``p_x' = p_y' = 0``,
``l' = v/v_0 - (1+\delta)/p_s - A\sin\Phi\,\partial_\delta\Phi``,
``\delta' = -A\,(\omega/v)\sin\Phi``.
"""
function _field_cav(z, s, ref, L, omega, A, v0, phis)
    x, px, y, py, l, d = z
    ps = sqrt((1 + d)^2 - px^2 - py^2)
    vr = _vel(d, ref); v = vr * v0
    Phi = omega * (s - L / 2) / v0 - omega * l / v - phis
    sP = sin(Phi)
    dPhi_dd = omega * l * _dvel(d, ref) * v0 / v^2
    Hl = A * sP * omega / v
    Hd = vr - (1 + d) / ps - A * sP * dPhi_dd
    return [px / ps, zero(x), py / ps, zero(x), Hd, -Hl]
end

"""
    _midpoint(F, z, s, hs)

One implicit midpoint step z' = z + hs F((z + z')/2, s + hs/2). Fixed-point
iterations (≤ 40) stop at |Δ| ≤ 1e-16 for numbers. Symplectic for any, even
s-dependent, Hamiltonian; order 2.
"""
function _midpoint(F, z, s, hs)
    zn = z .+ hs .* F(z, s + hs / 2)
    for _ in 1:40
        znew = z .+ hs .* F((z .+ zn) ./ 2, s + hs / 2)
        if eltype(znew) <: Real
            done = maximum(abs.(znew .- zn)) <= 1e-16
            zn = znew
            done && break
        else
            zn = znew
        end
    end
    return zn
end

"""
    _integrate(F, z, L, n)

n steps of the 4th-order Yoshida composition of [`_midpoint`](@ref Bell.ExactIntegrator._midpoint)
over length L (substeps W1 h, W0 h, W1 h; h = L/n), with the s value advanced
correctly for s-dependent fields.
"""
function _integrate(F, z, L, n)
    hs = L / n; s = 0.0; w = collect(z)
    for _ in 1:n
        w = _midpoint(F, w, s, W1 * hs);  s += W1 * hs
        w = _midpoint(F, w, s, W0 * hs);  s += W0 * hs
        w = _midpoint(F, w, s, W1 * hs);  s += W1 * hs
    end
    return w
end

_params(e::Drift) = (0.0, 0.0, 0.0)
_params(e::Quadrupole) = (0.0, e.k1, 0.0)
_params(e::Sextupole) = (0.0, 0.0, e.k2)
_params(e::SBend) = (curvature(e), e.k1, 0.0)

@doc raw"""
    exact_track_element(z, e, ref, C; ds = 0.1) -> Vector

Track a canonical point through one element with the exact (non-expanded)
Hamiltonian; nothing is truncated, longitudinal motion included.

# Algorithm
- `Marker`: identity.
- `Cavity`: the field [`_field_cav`](@ref Bell.ExactIntegrator._field_cav) of the RF
  Hamiltonian, ``n = \max(4, \lceil L/0.02\,\text{m}\rceil)`` Yoshida-4 steps.
- Magnets: [`_field_mag`](@ref Bell.ExactIntegrator._field_mag) with ``(h, k_1, k_2)``
  from the element type, ``n = \lceil L/ds\rceil`` steps (one for a drift).
- A thin `RFCavity` raises an error: use `rf_kick` for it.

# Properties
4th order in the step and symplectic. Generic: works on Float64 and `TaylorN`,
so `jacobian_of_function` gives exact Jacobians for symplecticity checks.
"""
function exact_track_element(z, e::AbstractElement, ref::RefParticle, C::Real; ds = 0.1)
    e isa Marker && return collect(z)
    if e isa Cavity
        omega, A, v0 = cavity_parameters(e, ref, C)
        n = max(4, ceil(Int, e.L / 0.02))
        return _integrate((w, s) -> _field_cav(w, s, ref, e.L, omega, A, v0, e.phi_s), z, e.L, n)
    end
    e isa RFCavity && error("thin RFCavity: use Tracking.rf_kick")
    h, k1, k2 = _params(e)
    n = e isa Drift ? 1 : max(1, ceil(Int, e.L / ds))
    return _integrate((w, s) -> _field_mag(w, ref, h, k1, k2), z, e.L, n)
end

"""
    exact_track_turn(z, lat, ref = lat.ref, C = circumference(lat); ds = 0.1) -> Vector

One turn of [`exact_track_element`](@ref) through all elements. Being exact in
the Hamiltonian, it is the reference against which the truncated maps are checked,
and it can track at any energy without rebuilding maps.
"""
function exact_track_turn(z, lat::Lattice, ref::RefParticle = lat.ref, C::Real = circumference(lat); ds = 0.1)
    w = collect(z)
    for e in lat.elements
        w = exact_track_element(w, e, ref, C; ds = ds)
    end
    return w
end

@doc raw"""
    renormalize(zs, ref, isync) -> (zs', ref')

Make particle `isync` the new reference, i.e. ramp the reference energy after
a cavity.

# Algorithm
With ``(x_s, p_{xs}, y_s, p_{ys}, l_s, \delta_s)`` the coordinates of the synchronous particle:
- ``p_0' = p_0(1 + \delta_s)``, and `ref'` is rebuilt from the new kinetic energy
  ``\sqrt{p_0'^2 + m^2} - m`` (new ``\gamma_0, \beta_0, B\rho``);
- ``\delta' = (1+\delta)/(1+\delta_s) - 1``;
- ``p_x' = (p_x - p_{xs})/(1+\delta_s)``, ``p_y'`` likewise: this is the adiabatic damping;
- ``x' = x - x_s``, ``y' = y - y_s``;
- ``\Delta t' = \Delta t - \Delta t_s`` with ``\Delta t = -l/v(\delta)``, and ``l' = -v(\delta)\,\Delta t'``.

The magnet fields are assumed to ramp with ``p_0``, so normalised strengths stay fixed.
"""
function renormalize(zs::Vector{<:AbstractVector}, ref::RefParticle, isync::Int)
    zsy = zs[isync]
    ds_ = zsy[6]
    v0 = ref.beta * C_LIGHT
    dts = -zsy[5] / (_vel(ds_, ref) * v0)
    p0n = ref.p0_MeV * (1 + ds_)
    Ekn = sqrt(p0n^2 + ref.mass_MeV^2) - ref.mass_MeV
    refn = reference_from_kinetic(ref.particle, Ekn; mass_MeV = ref.mass_MeV)
    out = similar(zs)
    for (k, z) in enumerate(zs)
        x, px, y, py, l, d = z
        v = _vel(d, ref) * v0
        dt = -l / v - dts
        dn = (1 + d) / (1 + ds_) - 1
        out[k] = [x - zsy[1], px / (1 + ds_) - zsy[2] / (1 + ds_), y - zsy[3], py / (1 + ds_) - zsy[4], -v * dt, dn]
    end
    return out, refn
end

end # module
