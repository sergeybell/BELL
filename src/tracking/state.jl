# =============================================================================
# tracking/state.jl - canonical particle state and a beam
# =============================================================================
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# ParticleState is a phase-space point. Coordinates.Particle is the species.
# Spin is stored only; nothing transports it.
# =============================================================================
"""
    BeamState

Canonical phase-space state of one particle and a beam of such particles.

[`Particle`](@ref) in `Coordinates` is the species (mass, charge, anomalous
moment). [`ParticleState`](@ref) is a point in the 6D canonical coordinates,
plus an optional spin slot. No Thomas–Bargmann–Michel–Telegdi equation is
applied here.
"""
module BeamState

export ParticleState, Beam

"""
    ParticleState(z, spin = nothing)

One particle in canonical coordinates `z = (x, px, y, py, l, δ)`.

# Fields
- `z::Vector{Float64}` – length-6 canonical vector. The constructor copies it.
- `spin::Union{Nothing,NTuple{3,Float64}}` – optional spin, default `nothing`.
  Stored only. Tracking does not change it.

`z` is mutated in place by [`track!`](@ref).
"""
struct ParticleState
    z::Vector{Float64}
    spin::Union{Nothing,NTuple{3,Float64}}
    function ParticleState(z::AbstractVector{<:Real},
                           spin::Union{Nothing,NTuple{3,Real}} = nothing)
        length(z) == 6 ||
            throw(ArgumentError("ParticleState: expected 6 canonical coordinates, got $(length(z))"))
        return new(Vector{Float64}(z), _spin64(spin))
    end
end

_spin64(::Nothing) = nothing
_spin64(s::NTuple{3,Real}) = (Float64(s[1]), Float64(s[2]), Float64(s[3]))

Base.show(io::IO, p::ParticleState) = print(io, "ParticleState(z=", p.z, ", spin=", p.spin, ")")

"""
    Beam(particles::Vector{ParticleState})
    Beam(p::ParticleState)
    Beam(z::AbstractVector{<:Real})
    Beam(zs::AbstractVector{<:AbstractVector{<:Real}})

A beam as a vector of [`ParticleState`](@ref).

A single 6-vector or a single [`ParticleState`](@ref) builds a one-particle beam.
A vector of 6-vectors builds one [`ParticleState`](@ref) per entry (each `z` is
copied). Passing an existing `Vector{ParticleState}` keeps those objects:
[`track!`](@ref) then updates their `z` vectors.
"""
struct Beam
    particles::Vector{ParticleState}
end

Beam(p::ParticleState) = Beam([p])
Beam(z::AbstractVector{<:Real}) = Beam([ParticleState(z)])

function Beam(zs::AbstractVector{<:AbstractVector{<:Real}})
    return Beam(ParticleState[ParticleState(z) for z in zs])
end

Base.show(io::IO, beam::Beam) = print(io, "Beam(", length(beam.particles), " particles)")

end # module
