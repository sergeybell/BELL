# =============================================================================
# exact_integrator.jl - tracking with the unexpanded Hamiltonian
# =============================================================================
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# The integrator is chosen by `method` (default ImplicitYoshida(4)).
# The vector field comes from hamiltonians.jl. Thin RF is rf_kick.
# renormalize updates the reference after a long cavity.
# =============================================================================
module ExactIntegrator

using ..Coordinates, ..Elements, ..Structure, ..Hamiltonians, ..Integrators

export exact_track_element, exact_track_turn, renormalize, cavity_parameters, transit_factor, rf_kick

@doc raw"""
    cavity_parameters(cav::LongCavity, ref, path_length) -> (ω, A, v0)

Constants of the finite-cavity Hamiltonian, recomputed from the current `ref`
at every call (so they follow the energy ramp):

- ``\omega = h\,2\pi v_0/\mathrm{path\_length}`` [rad/s];
- ``A = \dfrac{ZV}{p_0c}\,\dfrac{c}{\omega L}`` [1];
- ``v_0 = \beta_0 c`` [m/s].
"""
function cavity_parameters(cav::LongCavity, ref::RefParticle, path_length::Real)
    v0 = ref.beta * c_light
    omega = cav.harmonic * 2pi * v0 / path_length
    A = ref.particle.charge * cav.voltage_MV / ref.p0_MeV * c_light / (omega * cav.L)
    return omega, A, v0
end

@doc raw"""
    transit_factor(cav, ref, path_length) -> T

Transit-time factor of a finite cavity: ``T = \dfrac{\sin(\theta/2)}{\theta/2}``,
``\theta = 2\pi hL/\mathrm{path\_length} = \omega L/v_0``.
"""
transit_factor(cav::LongCavity, ref::RefParticle, path_length::Real) =
    (th = 2pi * cav.harmonic * cav.L / path_length; sin(th / 2) / (th / 2))

@doc raw"""
    rf_kick(z, cav::RFCavity, ref, path_length) -> Vector

Thin RF cavity acting on a canonical point (exact formulas, no expansion).

1. Arrival-time offset ``\Delta t = -l/v(\delta)``, RF phase ``\varphi = -h\omega_0\Delta t``.
2. Energy kick ``E' = E + ZV\sin(\varphi_s + \varphi)`` [MeV].
3. New momentum ``\delta' = \sqrt{E'^2 - m^2}/p_0 - 1``.
4. Same arrival time: ``l' = -v(\delta')\,\Delta t``.

``x, p_x, y, p_y`` are unchanged. The reference momentum is not changed.
"""
function rf_kick(z, cav::RFCavity, ref::RefParticle, path_length::Real)
    x, px, y, py, l, d = z
    v0 = ref.beta * c_light
    v = velocity_ratio(d, ref) * v0
    dt = -l / v
    phi = -cav.harmonic * 2pi * v0 / path_length * dt
    E = sqrt((ref.p0_MeV * (1 + d))^2 + ref.mass_MeV^2)
    E2 = E + ref.particle.charge * cav.voltage_MV * sin(cav.phi_s + phi)
    d2 = sqrt(E2^2 - ref.mass_MeV^2) / ref.p0_MeV - 1
    l2 = -velocity_ratio(d2, ref) * v0 * dt
    return [x, px, y, py, l2, d2]
end

@doc raw"""
    exact_track_element(z, e, ref, path_length; ds = 0.1, method = ImplicitYoshida(4)) -> Vector

Track a canonical point through one element with the exact (non-expanded)
Hamiltonian. The integrator is `method`; the field is [`magnetic_field`](@ref)
or [`cavity_field`](@ref). Nothing is truncated; the cut is the step in s.
"""
function exact_track_element(z, e::AbstractElement, ref::RefParticle, path_length::Real;
                             ds = 0.1, method::ODEMethod = ImplicitYoshida(4))
    e isa Marker && return collect(z)
    if e isa LongCavity
        omega, A, v0 = cavity_parameters(e, ref, path_length)
        n = max(4, ceil(Int, e.L / 0.02))
        f = (w, s) -> cavity_field(w, s, ref, e.L, omega, A, v0, e.phi_s)
        return integrate_fixed(f, collect(z), 0.0, e.L, n; method).z
    end
    e isa RFCavity && error("thin RFCavity: use rf_kick")
    h, k1, k2 = magnetic_parameters(e)
    n = e isa Drift ? 1 : max(1, ceil(Int, e.L / ds))
    f = (w, s) -> magnetic_field(w, ref, h, k1, k2)
    return integrate_fixed(f, collect(z), 0.0, e.L, n; method).z
end

"""
    exact_track_turn(z, lattice; ds = 0.1, method = ImplicitYoshida(4)) -> Vector

One turn of [`exact_track_element`](@ref) through all elements.
"""
function exact_track_turn(z, lattice::Lattice, ref::RefParticle = lattice.ref,
                          Lpath::Real = path_length(lattice); ds = 0.1,
                          method::ODEMethod = ImplicitYoshida(4))
    w = collect(z)
    for e in lattice.elements
        w = exact_track_element(w, e, ref, Lpath; ds = ds, method = method)
    end
    return w
end

@doc raw"""
    renormalize(particles, ref, isync) -> (particles', ref')

Make particle `isync` the new reference, i.e. ramp the reference energy after
a cavity.

With ``(x_s, p_{xs}, y_s, p_{ys}, l_s, \delta_s)`` the coordinates of the synchronous particle:
- ``p_0' = p_0(1 + \delta_s)``, and `ref'` is rebuilt from the new kinetic energy;
- ``\delta' = (1+\delta)/(1+\delta_s) - 1``;
- ``p_x' = (p_x - p_{xs})/(1+\delta_s)``, ``p_y'`` likewise (adiabatic damping);
- ``x' = x - x_s``, ``y' = y - y_s``;
- ``\Delta t' = \Delta t - \Delta t_s`` with ``\Delta t = -l/v(\delta)``, and ``l' = -v(\delta)\,\Delta t'``.

Magnet fields are assumed to ramp with ``p_0``, so normalised strengths stay fixed.
"""
function renormalize(particles::Vector{<:AbstractVector}, ref::RefParticle, isync::Int)
    zsy = particles[isync]
    ds_ = zsy[6]
    v0 = ref.beta * c_light
    dts = -zsy[5] / (velocity_ratio(ds_, ref) * v0)
    p0n = ref.p0_MeV * (1 + ds_)
    Ekn = sqrt(p0n^2 + ref.mass_MeV^2) - ref.mass_MeV
    refn = reference_from_kinetic(ref.particle, Ekn; mass_MeV = ref.mass_MeV)
    shifted = similar(particles)
    for (i, z) in enumerate(particles)
        x, px, y, py, l, d = z
        v = velocity_ratio(d, ref) * v0
        dt = -l / v - dts
        dn = (1 + d) / (1 + ds_) - 1
        shifted[i] = [x - zsy[1], px / (1 + ds_) - zsy[2] / (1 + ds_), y - zsy[3], py / (1 + ds_) - zsy[4], -v * dt, dn]
    end
    return shifted, refn
end

end # module
