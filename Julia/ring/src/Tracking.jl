# =============================================================================
# Tracking.jl - symplectic tracking with factorised Lie maps
# =============================================================================
#
# WHAT
#   build_factors(lat)           one LieFactor (R, f3) per element: the order-2
#                                Lie map exp(-L:H:) of the element is computed
#                                by LieOperators.lie_map and factorised.
#   apply_factor(F, z)           z -> R * N3(P4(z)) (f4 flow, then f3 flow), where each is
#                                applied as the time-1 flow of the Hamiltonian
#                                -f3, integrated with the IMPLICIT MIDPOINT
#                                rule (symplectic for any Hamiltonian), so each
#                                factor is a canonical transformation to
#                                machine precision at any amplitude.
#   track(factors, z; turns)     factor-by-factor tracking, z in CANONICAL
#                                coordinates; RF cavities are applied as thin
#                                kicks (rf_kick).
#   rf_kick(z, cav, ref, C)      thin RF cavity, own derivation:
#                                 dt = -l/v,  phi = -h omega0 dt,
#                                 E -> E + Z V sin(phi_s + phi),
#                                 delta from the new E, l -> -v(delta_new) dt.
#                                (Energy kick at fixed arrival time = canonical
#                                in (t,-E); the change of variables is canonical
#                                too, so the kick is symplectic.)  Phase 1: the
#                                kick is implemented but NOT used/validated.
#
# WHY
#   This is the "Lie" tracking engine of the program: generic code (works on
#   Float64 and on TaylorN, so Symplectic.jacobian_of_function can check it).
#
# WHAT YOU SEE
#   Nothing printed here; scripts/02 prints the tracking checks.
# =============================================================================
module Tracking

using TaylorSeries, LinearAlgebra
using ..Constants, ..Elements, ..Hamiltonians, ..Jets, ..LieOperators

export build_factors, apply_factor, track, rf_kick, element_lie_map

const _MCACHE = Dict{Any,Any}()

"Order-2 Lie map exp(-L:H:) z of an element (cached by element value)."
element_lie_map(e::AbstractElement, ref::RefParticle) =
    get!(_MCACHE, (e, ref.p0_MeV)) do
        (e isa Marker) ? collect(Z) : lie_map(hamiltonian(e, ref), elem_length(e))
    end

"One LieFactor per element (RFCavity entries are kept as the cavity itself)."
function build_factors(lat::Lattice)
    cache = Dict{Any,LieFactor}()
    out = Any[]
    for e in lat.elements
        if e isa RFCavity
            push!(out, e)
        else
            push!(out, get!(() -> dragt_finn_factor(element_lie_map(e, lat.ref); name = e.name), cache, e))
        end
    end
    return out
end

"exp(:f:) z by implicit midpoint on dz/dt = field(z), t in [0,1]."
function _flow(fieldvec, z::AbstractVector; nsteps = 2, iters = 12)
    h = 1 / nsteps
    w = collect(z)
    ev(u) = [evaluate(fieldvec[i], u) for i in 1:6]
    for _ in 1:nsteps
        wn = w .+ h .* ev(w)
        for _ in 1:iters
            wnew = w .+ h .* ev((w .+ wn) ./ 2)
            if eltype(wnew) <: Real
                done = maximum(abs.(wnew .- wn)) <= 1e-17
                wn = wnew
                done && break
            else
                wn = wnew
            end
        end
        w = wn
    end
    return w
end

"N3 o P4 applied to z: first the f4 flow, then the f3 flow (each symplectic)."
apply_nonlinear(F::LieFactor, z::AbstractVector) =
    _flow(F.field, iszero(F.f4) ? collect(z) : _flow(F.field4, z))
apply_factor(F::LieFactor, z) = (w = apply_nonlinear(F, z); [sum(F.R[i, j] * w[j] for j in 1:6) for i in 1:6])

function rf_kick(z, cav::RFCavity, ref::RefParticle, C::Real)
    x, px, y, py, l, d = z
    v0 = ref.beta * C_LIGHT
    v = velocity_ratio(d, ref) * v0
    dt = -l / v
    phi = -cav.harmonic * 2pi * v0 / C * dt
    E = sqrt((ref.p0_MeV * (1 + d))^2 + ref.mass_MeV^2)
    E2 = E + ref.particle.charge * cav.voltage_MV * sin(cav.phi_s + phi)
    d2 = sqrt(E2^2 - ref.mass_MeV^2) / ref.p0_MeV - 1
    l2 = -velocity_ratio(d2, ref) * v0 * dt
    return [x, px, y, py, l2, d2]
end

"Track z (canonical) through the factor list for `turns` turns; optional RF."
function track(factors, z; turns::Int = 1, ref = nothing, C = 0.0)
    w = collect(z)
    for _ in 1:turns, F in factors
        w = F isa LieFactor ? apply_factor(F, w) : rf_kick(w, F, ref, C)
    end
    return w
end

end # module
