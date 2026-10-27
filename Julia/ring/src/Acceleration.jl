# =============================================================================
# Acceleration.jl - ramping the reference energy with a finite RF cavity
# =============================================================================
#
# WHAT
#   insert_cavity(lat, cav; in_element) splits one long drift (default: the
#   first "OBM" straight-section drift) as  drift - cavity - drift.
#   accelerate(lat, ref0, users; ...) tracks several particles turn by turn
#   with ExactIntegrator; directly after the cavity the synchronous particle
#   (particle 1, starts on the reference orbit) becomes the new reference
#   (renormalize): this is how the reference energy ramps, and it is where
#   adiabatic damping of x', y', delta appears (momenta are normalised to p0).
#   When Ekin >= target the cavity phase is switched to phi_flat (default
#   pi = stationary bucket, no net acceleration) for `flat_turns` turns.
#   Records per turn: kinetic energy of the reference, and all particles in
#   user coordinates (x, x', y, y', phi, delta).
#   measure_tune(series): linear-motion tune from x_{n+1} + x_{n-1} = 2 cos(2 pi Q) x_n + b
#   (least squares; exact for linear motion about a KNOWN fixed point), Q in [0, 0.5].
#
# ASSUMPTIONS (stated in the scripts as well)
#   * magnet fields ramp in proportion to p0 (fixed normalised strengths);
#   * RF frequency follows h * f0(p0) (phase-locked to the synchronous particle);
#   * cavity field purely longitudinal on axis (no transverse RF fields).
#
# WHAT YOU SEE
#   Nothing printed here.
# =============================================================================
module Acceleration

using ..Constants, ..Elements, ..ExactIntegrator

export accelerate, measure_tune, turn_snapshots

function accelerate(lat::Lattice, ref0::RefParticle, users::Vector{<:AbstractVector};
                    Ekin_target = 270.0, max_turns = 10_000, flat_turns = 1000, phi_flat = pi, ds = 0.2)
    C = circumference(lat)
    h = lat.harmonic
    ref = ref0
    zs = [user_to_canonical(u, ref, C, h) for u in users]
    rec = (turn = Int[], Ekin = Float64[], gamma = Float64[], beta = Float64[], p0 = Float64[],
           u = [Vector{Vector{Float64}}() for _ in users])
    function record!(n)
        push!(rec.turn, n); push!(rec.Ekin, ref.Ekin_MeV); push!(rec.gamma, ref.gamma); push!(rec.beta, ref.beta)
        push!(rec.p0, ref.p0_MeV)
        for k in eachindex(zs)
            push!(rec.u[k], canonical_to_user(zs[k], ref, C, h))
        end
    end
    record!(0)
    curlat = lat
    flat_start = 0
    for n in 1:max_turns
        for e in curlat.elements
            zs = [exact_track_element(z, e, ref, C; ds = ds) for z in zs]
            if e isa Cavity
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
    return rec, flat_start
end

"""
Tune of turn-by-turn data x_n of linear motion about a fixed point e:
x_{n+1} + x_{n-1} = 2 cos(2 pi Q) x_n + b,  b = (2 - 2cos) e  (least squares in
(cos, b), so an imperfectly known fixed point does not bias Q). Q in [0, 0.5].
"""
function measure_tune(x::AbstractVector)
    n = length(x)
    A = hcat(2 .* x[2:n-1], ones(n - 2))
    y = x[3:n] .+ x[1:n-2]
    c = (A \ y)[1]
    return acos(clamp(c, -1, 1)) / 2pi
end
"""
turn_snapshots(lat, ref, zs) -> (snaps, names, s, zs_end, ref_end)
One turn, element by element: user coordinates of all particles after
every element (renormalisation right after the cavity, as in accelerate).
"""
function turn_snapshots(lat::Lattice, ref::RefParticle, zs; ds = 0.2)
    C = circumference(lat); h = lat.harmonic
    snaps = [[canonical_to_user(z, ref, C, h) for z in zs]]
    names = ["start"]; svals = [0.0]; s = 0.0
    for e in lat.elements
        zs = [exact_track_element(z, e, ref, C; ds = ds) for z in zs]
        e isa Cavity && ((zs, ref) = renormalize(zs, ref, 1))
        s += elem_length(e)
        push!(snaps, [canonical_to_user(z, ref, C, h) for z in zs]); push!(names, e.name); push!(svals, s)
    end
    return snaps, names, svals, zs, ref
end

end # module