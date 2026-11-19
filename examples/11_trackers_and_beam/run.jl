# 11_trackers_and_beam/run.jl
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# ParticleState, Beam, and the three trackers on lattices/fodo.bell.
# FieldTracker integrates the element s-Hamiltonian. LieTracker applies the
# stored Lie factors. MapTracker repeats the linear one-turn matrix.
#
#   julia --project=examples examples/11_trackers_and_beam/run.jl

using Bell, Printf

include(joinpath(@__DIR__, "..", "common", "print_util.jl"))

const DS = 0.5
const N_TURNS = 3
const STARTED = time()

function beam_from(rows)
    particles = Vector{ParticleState}(undef, length(rows))
    for i in eachindex(rows)
        spin = i == 1 ? (0.0, 0.0, 1.0) : nothing
        particles[i] = ParticleState(rows[i], spin)
    end
    return Beam(particles)
end

function snapshot(beam)
    return [copy(p.z) for p in beam.particles]
end

function print_rows(tag, rows)
    @printf("  %-22s %12s %12s %12s %12s %12s %12s\n",
            tag, VAR_NAMES[1], VAR_NAMES[2], VAR_NAMES[3], VAR_NAMES[4], VAR_NAMES[5], VAR_NAMES[6])
    for (i, z) in enumerate(rows)
        @printf("  %-22s %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e\n",
                "particle $i", z[1], z[2], z[3], z[4], z[5], z[6])
    end
    return nothing
end

function max_sep(a, b)
    return maximum(maximum(abs, u - v) for (u, v) in zip(a, b))
end


banner("Trackers and a beam on the FODO ring")

lat = load_structure(joinpath(@__DIR__, "..", "lattices", "fodo.bell"))
ref = lat.ref
C = path_length(lat)
starts = [
    [1.0e-3, 0.0, 1.2e-3, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.0, 1.0e-3],
    [-5.0e-4, 1.0e-4, 4.0e-4, -8.0e-5, 2.0e-4, -5.0e-4],
]


section("1. Lattice and initial beam",
        "Three particles. Particle 1 carries a stored spin; tracking does not change it.",
        "Coordinates are canonical. l = -v dt, delta = dp/p0.")

println(describe_reference(ref))
println()
@printf("  elements = %d    C = %.6f m    harmonic = %d\n", length(lat.elements), C, lat.harmonic)
println()
print_rows("initial", starts)
probe = ParticleState(starts[1], (0.0, 0.0, 1.0))
@printf("\n  ParticleState copies z and stores spin = %s\n", probe.spin)
@printf("  Beam of the three starts: %s\n", beam_from(starts))


section("2. Three turns",
        "FieldTracker(ds = $(DS) m) calls exact_track_turn.",
        "LieTracker holds one Lie factor per element from this lattice.",
        "MapTracker repeats LinearMap(ring_matrix). That map is linear, so it",
        "drops the energy-dependent and chromatic terms the other two keep.")

field = FieldTracker(ds = DS)
t_lie_build = @elapsed lie = LieTracker(lat)
R = ring_matrix(lat)
mapper = MapTracker(LinearMap(R))
@printf("  %s\n", field)
@printf("  %s    built in %.3f s\n", lie, t_lie_build)
@printf("  %s\n", mapper)
@printf("  symplectic error of the one-turn matrix = %.3e\n", symplectic_error(R))

beam_field = beam_from(starts)
beam_lie = beam_from(starts)
beam_map = beam_from(starts)
before = snapshot(beam_field)

t_field = @elapsed track!(beam_field, lat, field; turns = N_TURNS)
t_lie = @elapsed track!(beam_lie, lat, lie; turns = N_TURNS)
t_map = @elapsed track!(beam_map, lat, mapper; turns = N_TURNS)
after_field = snapshot(beam_field)
after_lie = snapshot(beam_lie)
after_map = snapshot(beam_map)

@printf("\n  %d turns:  Field %.3f s    Lie %.3f s    linear map %.3f s\n",
        N_TURNS, t_field, t_lie, t_map)
println()
print_rows("FieldTracker", after_field)
println()
print_rows("LieTracker", after_lie)
println()
print_rows("MapTracker", after_map)
println()
@printf("  |Field - initial|  = %.3e\n", max_sep(after_field, before))
@printf("  |Lie - initial|    = %.3e\n", max_sep(after_lie, before))
@printf("  |linear - initial| = %.3e\n", max_sep(after_map, before))
@printf("  |Field - Lie|      = %.3e\n", max_sep(after_field, after_lie))
@printf("  |Field - linear|   = %.3e\n", max_sep(after_field, after_map))
@printf("  spin after Field / Lie / linear = %s / %s / %s\n",
        beam_field.particles[1].spin, beam_lie.particles[1].spin, beam_map.particles[1].spin)


section("3. One turn on a mutable 6-vector",
        "track!(z, lattice, tracker) writes the new coordinates into z.",
        "The field result is exact_track_turn with the same ds and method.",
        "The Lie result is track on the factors the tracker already holds.")

z_field = copy(starts[1])
z_exact = exact_track_turn(starts[1], lat; ds = field.ds, method = field.method)
track!(z_field, lat, field; turns = 1)
z_lie = copy(starts[1])
track!(z_lie, lat, lie; turns = 1)
z_factors = track(lie.factors, starts[1]; turns = 1)
@printf("  |track!(z, FieldTracker) - exact_track_turn| = %.3e\n", maximum(abs, z_field - z_exact))
@printf("  |track!(z, LieTracker) - track(factors)|     = %.3e\n", maximum(abs, z_lie - z_factors))
@printf("  z in  = %s\n", join([@sprintf("% .4e", v) for v in starts[1]], " "))
@printf("  z out = %s\n", join([@sprintf("% .4e", v) for v in z_field], " "))


section("4. Symplecticity of one step",
        "symplectic_error(J) is max |J' S J - S|.",
        "The one-turn matrix was printed above. Here the same check is the",
        "quadrupole field step at the first particle, and one Lie factor of that quadrupole.",
        "A full-turn Jacobian in series arithmetic is the same call on exact_track_turn;",
        "it is not run on this 16-cell ring.")

qf = first(e for e in lat.elements if e isa Quadrupole && e.k1 > 0)
t_sym = @elapsed err_step = symplectic_error(
    z -> exact_track_element(z, qf, ref, C; ds = DS), starts[1])
Fq = first(F for F in lie.factors if F.name == qf.name)
t_sym_lie = @elapsed err_factor = symplectic_error(z -> apply_factor(Fq, z), starts[1])
@printf("  %s    L = %.3f m    k1 = %.4f 1/m^2\n", qf, qf.L, qf.k1)
@printf("  symplectic error of one field step     = %.3e    (%.3f s)\n", err_step, t_sym)
@printf("  symplectic error of Lie factor %s     = %.3e    (%.3f s)\n", Fq.name, err_factor, t_sym_lie)
@printf("  symplectic error of that factor's R    = %.3e\n", symplectic_error(Fq.R))

@printf("\n  wall time %.2f s\n", time() - STARTED)
println("Done.")
