# =============================================================================
# field_track.jl - 6D tracking with field ODEs and thin RF kick
# =============================================================================
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# Magnets: integrate Fields.element_rhs with a chosen ODE method (no polynomial
# cut; the cut is the step in s). Thin RFCavity: rf_kick (closed-form formulas
# in canonical z). renormalize updates the reference after a kick on a ramp.
# Independent of Hamiltonians / Lie algebra.
# =============================================================================
module FieldTrack

using ..Coordinates, ..Elements, ..Structure, ..Fields, ..Integrators

export exact_track_element, exact_track_turn, renormalize, rf_kick

@doc raw"""
    rf_kick(z, cav::RFCavity, ref, path_length) -> Vector

Thin RF cavity acting on a canonical point (exact formulas, no series cut).
The same expressions accept numbers or DA.

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

Track a canonical point through one element with the field equations of motion
([`element_rhs`](@ref Bell.Fields.element_rhs)). Thin [`RFCavity`](@ref) uses
[`rf_kick`](@ref). Nothing is truncated as a polynomial; the cut is the step
in ``s``.
"""
function exact_track_element(z, e::AbstractElement, ref::RefParticle, path_length::Real;
                             ds = 0.1, method::ODEMethod = ImplicitYoshida(4))
    e isa Marker && return collect(z)
    if e isa RFCavity
        e.L == 0 && return rf_kick(z, e, ref, path_length)
        # thick: Strang splitting of H = H_drift + G/L (2nd order in ds)
        n = max(1, ceil(Int, e.L / ds)); hs = e.L / n
        kick = RFCavity(e.name, e.voltage_MV / n, e.harmonic, e.phi_s, 0.0)
        fd = element_rhs(Drift(e.name, hs), ref)
        w = collect(z)
        for _ in 1:n
            w = integrate_fixed(fd, w, 0.0, hs / 2, 1; method).z
            w = rf_kick(w, kick, ref, path_length)
            w = integrate_fixed(fd, w, 0.0, hs / 2, 1; method).z
        end
        return w
    end
    L = elem_length(e)
    L == 0 && return collect(z)
    n = e isa Drift ? 1 : max(1, ceil(Int, L / ds))
    f = element_rhs(e, ref)
    return integrate_fixed(f, collect(z), 0.0, L, n; method).z
end

"""
    exact_track_turn(z, lattice; ds = 0.1, method = ImplicitYoshida(4)) -> Vector

One turn of [`exact_track_element`](@ref) through all elements (magnets and
thin RF kicks).
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

Make particle `isync` the new reference after a thin RF kick on a ramp.

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
        shifted[i] = [x - zsy[1], px / (1 + ds_) - zsy[2] / (1 + ds_), y - zsy[3],
                      py / (1 + ds_) - zsy[4] / (1 + ds_), -v * dt, dn]
    end
    return shifted, refn
end

end # module
