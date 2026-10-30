# =============================================================================
# Tracking.jl - symplectic tracking with factorised Lie maps
# =============================================================================
#
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# WHAT
#   build_factors(lat)           one Dragt-Finn factor per element
#   ring_lie_factor(lat)         product of those factors in the Lie algebra
#   apply_factor(F, z)           flow of f4, then f3, then R (implicit midpoint)
#   track(factors, z; turns)     apply each factor in order; RF as a thin kick
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
using ..Coordinates, ..Elements, ..Structure, ..Hamiltonians, ..LieAlgebra

export build_factors, apply_factor, track, rf_kick, element_lie_map, ring_lie_factor

const _MCACHE = Dict{Any,Any}()

"""
    element_lie_map(e, ref) -> Vector{TaylorN}

Map exp(−L :H:) z of one element ([`lie_map`](@ref) of its [`hamiltonian`](@ref)),
truncated to MAP_ORDER; identity for markers. Cached per (element, p0).
"""
element_lie_map(e::AbstractElement, ref::RefParticle) =
    get!(_MCACHE, (e, ref.p0_MeV)) do
        (e isa Marker) ? collect(Z) : lie_map(hamiltonian(e, ref), elem_length(e))
    end

"""
    ring_lie_factor(lat; name = lat.name) -> LieFactor

One-turn Dragt-Finn factor of the magnetic ring: [`product_factors`](@ref) of
[`build_factors`](@ref)(lat), with thin RF cavities skipped. The RF factor is
concatenated separately (`one_turn_factor_rf`). This is the primary one-turn
object of Bell; `observables` reads all ring quantities from it.
"""
function ring_lie_factor(lat::Lattice; name = lat.name)
    fs = LieFactor[F for F in build_factors(lat) if F isa LieFactor]
    return product_factors(fs; name)
end

"""
    build_factors(lat) -> Vector{Any}

One [`LieFactor`](@ref) per element: `dragt_finn_factor(element_lie_map(e, ref))`,
computed once per distinct element value. A thin [`RFCavity`](@ref) is kept
as the cavity object itself, so that [`track`](@ref) applies [`rf_kick`](@ref) there.
"""
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

"""
    _flow(field, z; nsteps = 2, iters = 12)

Time-1 flow of the polynomial vector field `field` (that of a generator f, i.e.
exp(:f:) z) by the implicit midpoint rule.

# Algorithm
`nsteps` steps of h = 1/nsteps:
w_{n+1} = w_n + h F((w_n + w_{n+1})/2).
The implicit equation is solved by fixed-point iteration from an explicit Euler
guess. With numbers it stops at |Δ| ≤ 1e-17 or after `iters` iterations; with
`TaylorN` all `iters` iterations are done.

# Properties
- The implicit midpoint rule is symplectic for any Hamiltonian vector field and
  of order 2.
- Because f3 and f4 are small (non-linear) generators, two steps reproduce the
  Lie series to well below the truncation error, while staying exactly
  symplectic.
- Generic: works on Float64 and on TaylorN (exact Jacobians).
"""
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

"""
    apply_nonlinear(F, z) -> Vector

Apply the non-linear part of a factor, N3 ∘ P4: first the f4 flow (skipped if
f4 = 0), then the f3 flow, each by [`_flow`](@ref Bell.Tracking._flow).
"""
apply_nonlinear(F::LieFactor, z::AbstractVector) =
    _flow(F.field, iszero(F.f4) ? collect(z) : _flow(F.field4, z))
"""
    apply_factor(F::LieFactor, z) -> Vector

Track one point (numbers or `TaylorN`) through a factor: w = N3(P4(z)), then R·w.
Each step is symplectic, so the composition is symplectic to round-off. It is not
the truncated polynomial `factor_to_map(F)`.
"""
apply_factor(F::LieFactor, z) = (w = apply_nonlinear(F, z); [sum(F.R[i, j] * w[j] for j in 1:6) for i in 1:6])

"""
    rf_kick(z, cav::RFCavity, ref, C) -> Vector

Thin RF cavity acting on a canonical point.

# Algorithm
1. Arrival-time offset dt = −l/v(δ), RF phase φ = −h ω0 dt with ω0 = 2π v0/C.
2. Energy kick E → E + Z V sin(φs + φ) [MeV].
3. New δ from the new energy: δ' = √(E'² − m²)/p0 − 1.
4. Keep the arrival time: l' = −v(δ') dt.

An energy kick at fixed arrival time is canonical in (t, −E), and the change of
variables to (l, δ) is canonical too, so the kick is symplectic. x, px, y, py are
unchanged (no transverse RF field). Works on numbers and `TaylorN`.
"""
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

"""
    track(factors, z; turns = 1, ref = nothing, C = 0.0) -> Vector

Track a canonical point `z` for `turns` turns through the list from
[`build_factors`](@ref). Entries that are `LieFactor`s use [`apply_factor`](@ref);
`RFCavity` entries use [`rf_kick`](@ref) and then need `ref` and `C`.
"""
function track(factors, z; turns::Int = 1, ref = nothing, C = 0.0)
    w = collect(z)
    for _ in 1:turns, F in factors
        w = F isa LieFactor ? apply_factor(F, w) : rf_kick(w, F, ref, C)
    end
    return w
end

end # module
