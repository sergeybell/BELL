# =============================================================================
# field_track.jl - 6D tracking of one element
# =============================================================================
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# The model of an element is its Hamiltonian H, and dz/ds = S ∇H.
# A magnet and a thick RF cavity go through ElementFlow, which integrates H.
# There is no polynomial truncation on a point: the cut is the step in s.
# A Marker is the identity. A thin RF cavity applies closed-form formulas.
# The field equations remain for comparison and for tests. After a kick on a
# ramp the reference is updated.
# =============================================================================
module FieldTrack

using ..Coordinates, ..Elements, ..Structure, ..Integrators, ..ElementFlow

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

Track a canonical point through one element. A magnet and a thick
[`RFCavity`](@ref) are [`element_flow`](@ref): ``dz/ds = S\nabla H``.
A [`Marker`](@ref) is the identity. A thin cavity uses [`rf_kick`](@ref).
A [`Drift`](@ref) is one step of `method`. Nothing is truncated as a
polynomial; the cut is the step in ``s``.
"""
function exact_track_element(z, ::Marker, ::RefParticle, ::Real;
                             ds = 0.1, method::ODEMethod = ImplicitYoshida(4))
    return collect(z)
end

function exact_track_element(z, e::RFCavity, ref::RefParticle, path_length::Real;
                             ds = 0.1, method::ODEMethod = ImplicitYoshida(4))
    e.L == 0 && return rf_kick(z, e, ref, path_length)
    n = max(1, ceil(Int, e.L / ds))
    return element_flow(e, ref, z; method = method, nsteps = n, path_length = path_length)
end

function exact_track_element(z, e::Drift, ref::RefParticle, path_length::Real;
                             ds = 0.1, method::ODEMethod = ImplicitYoshida(4))
    return element_flow(e, ref, z; method = method, nsteps = 1, path_length = path_length)
end

function exact_track_element(z, e::AbstractElement, ref::RefParticle, path_length::Real;
                             ds = 0.1, method::ODEMethod = ImplicitYoshida(4))
    L = elem_length(e)
    L == 0 && return collect(z)
    n = max(1, ceil(Int, L / ds))
    return element_flow(e, ref, z; method = method, nsteps = n, path_length = path_length)
end

"""
    exact_track_turn(z, lattice; ds = 0.1, method = ImplicitYoshida(4)) -> Vector

One turn of [`exact_track_element`](@ref) through every element.
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
