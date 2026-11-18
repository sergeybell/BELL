# =============================================================================
# tracking/trackers.jl - dispatch hooks over the existing tracking paths
# =============================================================================
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# FieldTracker calls exact_track_turn. LieTracker calls track on factors it
# already holds. MapTracker calls a one-turn map. None of them reimplements
# the physics, and none of them transports spin.
# =============================================================================
"""
    Trackers

[`track!`](@ref) for a [`Beam`](@ref) through a lattice with one of
[`FieldTracker`](@ref), [`LieTracker`](@ref) or [`MapTracker`](@ref).

The field and Lie implementations stay in `FieldTrack` and `Tracking`.
A tracker only chooses which of those functions to call. Spin on a
[`ParticleState`](@ref) is left as stored.
"""
module Trackers

using ..BeamState: Beam
using ..Structure: Lattice
using ..Integrators: ODEMethod, ImplicitYoshida
using ..FieldTrack: exact_track_turn
using ..Tracking: build_factors, track
using ..LieAlgebra: LieFactor

export AbstractTracker, FieldTracker, LieTracker, MapTracker, track!

"""
    AbstractTracker

Supertype of objects [`track!`](@ref) can run. The lattice is always passed
to `track!`; each concrete tracker documents whether it reads that lattice
or uses data it already holds.
"""
abstract type AbstractTracker end

"""
    FieldTracker(; ds = 0.1, method = ImplicitYoshida(4))
    FieldTracker(ds, method = ImplicitYoshida(4))

[`track!`](@ref) calls [`exact_track_turn`](@ref) with step `ds` and ODE
`method`. That turn integrates the element s-Hamiltonian. The default method
is the default of `exact_track_turn`. `ds` must be positive.
"""
struct FieldTracker{M<:ODEMethod} <: AbstractTracker
    ds::Float64
    method::M
end

function FieldTracker(; ds::Real = 0.1, method::ODEMethod = ImplicitYoshida(4))
    ds > 0 || throw(ArgumentError("FieldTracker: ds must be positive, got $ds"))
    return FieldTracker{typeof(method)}(Float64(ds), method)
end

FieldTracker(ds::Real, method::ODEMethod = ImplicitYoshida(4)) = FieldTracker(; ds, method)

function Base.show(io::IO, tracker::FieldTracker)
    print(io, "FieldTracker(ds=", tracker.ds, ", method=", tracker.method, ")")
    return nothing
end

"""
    LieTracker(factors)
    LieTracker(lattice) -> LieTracker

Lie tracker holding one [`LieFactor`](@ref) per element, as returned by
[`build_factors`](@ref). The lattice constructor stores `build_factors(lattice)`.

[`track!`](@ref) calls [`track`](@ref) on these factors. The `lattice` argument
of `track!` is not read: rebuild a tracker to follow a different lattice.
"""
struct LieTracker <: AbstractTracker
    factors::Vector{LieFactor}
end

LieTracker(lattice::Lattice) = LieTracker(build_factors(lattice))

Base.show(io::IO, tracker::LieTracker) =
    print(io, "LieTracker(", length(tracker.factors), " factors)")

"""
    MapTracker(one_turn)

Tracker that applies a one-turn callable. `one_turn(z)` must return a 6-vector.
A [`LinearMap`](@ref) or [`DAMap`](@ref) is such a callable.

[`track!`](@ref) does not read the lattice; the callable is already the turn.
"""
struct MapTracker{F} <: AbstractTracker
    one_turn::F
end

Base.show(io::IO, tracker::MapTracker) = print(io, "MapTracker(", tracker.one_turn, ")")

function _write_turns!(beam::Beam, turns::Int, oneturn)
    for particle in beam.particles
        z = particle.z
        for _ in 1:turns
            z = oneturn(z)
        end
        particle.z .= z
    end
    return beam
end

"""
    track!(beam::Beam, lattice, tracker; turns = 1) -> beam
    track!(z::AbstractVector, lattice, tracker; turns = 1) -> z

Track `turns` turns and write the new canonical coordinates in place.

- [`FieldTracker`](@ref) uses [`exact_track_turn`](@ref) on `lattice`.
- [`LieTracker`](@ref) uses [`track`](@ref) on the factors it holds.
- [`MapTracker`](@ref) calls its one-turn map once per turn.

Spin is not modified. `z` must be a mutable 6-vector; the result is written
back into it. `turns <= 0` leaves the coordinates unchanged, as `track` does.
"""
function track!(beam::Beam, lattice::Lattice, tracker::FieldTracker; turns::Int = 1)
    ds, method = tracker.ds, tracker.method
    return _write_turns!(beam, turns, z -> exact_track_turn(z, lattice; ds, method))
end

function track!(beam::Beam, ::Lattice, tracker::LieTracker; turns::Int = 1)
    for particle in beam.particles
        particle.z .= track(tracker.factors, particle.z; turns = turns)
    end
    return beam
end

function track!(beam::Beam, ::Lattice, tracker::MapTracker; turns::Int = 1)
    return _write_turns!(beam, turns, tracker.one_turn)
end

function track!(z::AbstractVector{<:Real}, lattice::Lattice, tracker::AbstractTracker;
                turns::Int = 1)
    length(z) == 6 ||
        throw(ArgumentError("track!: expected 6 canonical coordinates, got $(length(z))"))
    beam = Beam(z)
    track!(beam, lattice, tracker; turns = turns)
    z .= beam.particles[1].z
    return z
end

end # module
