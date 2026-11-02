# =============================================================================
# acceleration.jl - ramping the reference energy with a finite RF cavity
# =============================================================================
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# accelerate tracks several particles in 6D with exact_track_element.
# There is no polynomial-order cut; the cut is the step ds.
# Particle 1 is the synchronous (reference) particle: after each LongCavity,
# renormalize makes it the new reference. Smooth (φ, δ) is not used here.
#
# turn_snapshots records one turn element by element. Cavity parameters
# (beta, p0, A) are recomputed from the current reference after renormalize,
# which is required on a ramp. They are not needed for a many-turn ramp
# summary, so the snapshot path is a separate function.
# =============================================================================
module Acceleration

using ..Coordinates, ..Elements, ..Structure, ..ExactIntegrator

export accelerate, measure_tune, turn_snapshots

"""
    accelerate(lattice, injection_reference, laboratory_particles;
               Ekin_target = 270.0, max_turns = 10_000,
               flat_turns = 1000, phi_flat = π, ds = 0.2) -> (record, flat_start)

Turn-by-turn acceleration with a finite [`LongCavity`](@ref) in `lattice`.

We study the reference particle. It is the first entry of `laboratory_particles`
and starts on the design orbit. After every long cavity, [`renormalize`](@ref)
makes that particle the new reference, which is how the energy ramps.

# Arguments
- `lattice` – lattice containing a `LongCavity` (see [`insert_cavity`](@ref)).
- `injection_reference` – reference particle at injection.
- `laboratory_particles` – initial particles in laboratory coordinates
  (x, x', y, y', φ, δ). Particle 1 must be the synchronous particle.
- `Ekin_target` [MeV] – when the reference reaches it, the cavity phase is set to
  `phi_flat` (default π: stationary bucket) for `flat_turns` more turns.
- `ds` – integrator step [m].

# Returns
- `record` – NamedTuple of per-turn vectors `turn, Ekin, gamma, beta, p0` and
  `u[k]` = laboratory coordinates of particle k on turn n.
- `flat_start` – turn at which the flat top started (0 if never).
"""
function accelerate(lattice::Lattice, injection_reference::RefParticle, laboratory_particles::Vector{<:AbstractVector};
                    Ekin_target = 270.0, max_turns = 10_000, flat_turns = 1000, phi_flat = pi, ds = 0.2)
    Lpath = path_length(lattice)
    h = lattice.harmonic
    ref = injection_reference
    zs = [laboratory_to_canonical(u, ref, Lpath, h) for u in laboratory_particles]
    record = (turn = Int[], Ekin = Float64[], gamma = Float64[], beta = Float64[], p0 = Float64[],
              u = [Vector{Vector{Float64}}() for _ in laboratory_particles])
    function record!(n)
        push!(record.turn, n); push!(record.Ekin, ref.Ekin_MeV); push!(record.gamma, ref.gamma); push!(record.beta, ref.beta)
        push!(record.p0, ref.p0_MeV)
        for k in eachindex(zs)
            push!(record.u[k], canonical_to_laboratory(zs[k], ref, Lpath, h))
        end
    end
    record!(0)
    curlat = lattice
    flat_start = 0
    for n in 1:max_turns
        for e in curlat.elements
            zs = [exact_track_element(z, e, ref, Lpath; ds = ds) for z in zs]
            if e isa LongCavity
                zs, ref = renormalize(zs, ref, 1)
            end
        end
        record!(n)
        if flat_start == 0 && ref.Ekin_MeV >= Ekin_target
            flat_start = n
            curlat = set_cavity_phase(curlat, phi_flat)
        end
        flat_start > 0 && n >= flat_start + flat_turns && break
    end
    return record, flat_start
end

"""
    measure_tune(x) -> Q ∈ [0, 0.5]

Tune of turn-by-turn data of linear motion about an unknown fixed point.
"""
function measure_tune(x::AbstractVector)
    n = length(x)
    A = hcat(2 .* x[2:n-1], ones(n - 2))
    y = x[3:n] .+ x[1:n-2]
    c = (A \ y)[1]
    return acos(clamp(c, -1, 1)) / 2pi
end

"""
    turn_snapshots(lattice, ref, zs; ds = 0.2) -> (snapshots, names, s, zs_end, ref_end)

One turn element by element. `snapshots[k]` holds the laboratory coordinates of
all particles after element k−1 (`snapshots[1]` = start). After a long cavity
the same [`renormalize`](@ref) as in [`accelerate`](@ref) is applied, so
`cavity_parameters` (β, p0, A) are taken from the current reference. That
recomputation is required on a ramp; a many-turn `accelerate` call does not
need the per-element snapshots, which is why this function is separate.
"""
function turn_snapshots(lattice::Lattice, ref::RefParticle, zs; ds = 0.2)
    Lpath = path_length(lattice); h = lattice.harmonic
    snapshots = [[canonical_to_laboratory(z, ref, Lpath, h) for z in zs]]
    names = ["start"]; svals = [0.0]; s = 0.0
    for e in lattice.elements
        zs = [exact_track_element(z, e, ref, Lpath; ds = ds) for z in zs]
        e isa LongCavity && ((zs, ref) = renormalize(zs, ref, 1))
        s += elem_length(e)
        push!(snapshots, [canonical_to_laboratory(z, ref, Lpath, h) for z in zs]); push!(names, e.name); push!(svals, s)
    end
    return snapshots, names, svals, zs, ref
end

end # module
