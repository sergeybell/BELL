# da_basics.jl
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# Differential-algebra arithmetic, then element and one-turn maps of
# lattices/fodo.bell. The one-turn map uses a fixed RK4 step count.
# chromaticity_series reads that map. The lattice method of chromaticity_series
# builds one_turn_map_da with the default adaptive stepper, so it is not called.
#
#   julia --project=examples examples/da_basics/run.jl

using Bell, Bell.DiffAlgebra, Bell.Integrators, Printf

include(joinpath(@__DIR__, "..", "common", "print_util.jl"))

const NSTEPS = 4
const METHOD = RK4()
const Z0 = [1.0e-3, -2.0e-4, 1.2e-3, 1.0e-4, 5.0e-4, 1.0e-3]
const STARTED = time()

function show_mono(label, p, powers)
    @printf("  %-28s = % .6e\n", label, coefficient(p, powers))
    return nothing
end


banner("Differential algebra, from arithmetic to a FODO one-turn map")

section("1. Arithmetic in a 2-variable space",
        "DASpace(nv, no) is the algebra of polynomials in nv variables, truncated at total degree no.",
        "da_vars are the coordinate functions. deriv differentiates one coordinate. compose substitutes a map.")

sp2 = DASpace(2, 3)
x, y = da_vars(sp2)
a = x^2 + 3 * x * y - y^3 / 6
dx = deriv(a, 1)
s = sin(x)
ac = compose([a], [x + y, y])[1]

@printf("  DASpace(%d, %d) has %d monomials\n", sp2.nv, sp2.no, nmonomials(sp2))
show_mono("a: x^2", a, [2, 0])
show_mono("a: x*y", a, [1, 1])
show_mono("a: y^3", a, [0, 3])
show_mono("da/dx: x", dx, [1, 0])
show_mono("da/dx: y", dx, [0, 1])
show_mono("sin(x): x^3", s, [3, 0])
@printf("  -1/6                         = % .6e\n", -1 / 6)
show_mono("a(x+y, y): x^2", ac, [2, 0])
show_mono("a(x+y, y): x*y", ac, [1, 1])
show_mono("a(x+y, y): y^2", ac, [0, 2])
show_mono("a(x+y, y): y^3", ac, [0, 3])
note("a(x+y, y) = x^2 + 5 x y + 4 y^2 - y^3/6.")


section("2. Shared 6-variable space",
        "map_space(N) is the space used by element_da_map: z = (x, px, y, py, l, delta).")

for N in (2, 3)
    sp = map_space(N)
    z = da_vars(sp)
    id = compose(z, z)
    gap = maximum(coeffnorm(id[i] - z[i]) for i in eachindex(z))
    @printf("  map_space(%d): %d variables, order %d, %d monomials, |compose(id, id) - id| = %.3e\n",
            N, sp.nv, sp.no, nmonomials(sp), gap)
end


section("3. Quadrupole map",
        "element_da_map integrates dz/ds = S grad H on the DA identity, RK4, $(NSTEPS) steps.")

lat = load_structure(joinpath(@__DIR__, "..", "lattices", "fodo.bell"))
ref = lat.ref
C = path_length(lat)
qf = first(e for e in lat.elements if e isa Quadrupole && e.k1 > 0)
R = linear_matrix(qf, ref)
println(describe_reference(ref))
@printf("  %s  L = %.3f m  k1 = %.4f 1/m^2\n", qf.name, qf.L, qf.k1)

for N in (2, 3)
    t = @elapsed M = element_da_map(qf, ref, N; method = METHOD, nsteps = NSTEPS)
    lin = linear_part(M)
    zf = [Float64(M[i](Z0)) for i in eachindex(M)]
    zn = element_flow(qf, ref, Z0; method = METHOD, nsteps = NSTEPS)
    @printf("  order %d in %.4f s    |linear_part - linear_matrix| = %.3e    symplectic error = %.3e\n",
            N, t, maximum(abs, lin - R), symplectic_error(M))
    @printf("           |DA(z) - element_flow| = %.3e\n", maximum(abs, zf - zn))
end


section("4. FODO one-turn map",
        "one_turn_map_da composes element maps. $(length(lat.elements)) elements, RK4, $(NSTEPS) steps.",
        "Orders 2 and 3. The heavy one-turn stays on this FODO file.")

@printf("  elements %d    C = %.7f m\n", length(lat.elements), C)
Rturn = ring_matrix(lat)
maps = Dict{Int,Any}()
for N in (2, 3)
    t = @elapsed M = one_turn_map_da(lat, N; method = METHOD, nsteps = NSTEPS)
    maps[N] = M
    lin = linear_part(M)
    @printf("  order %d in %.3f s    |linear_part - ring_matrix| = %.3e    symplectic error = %.3e\n",
            N, t, maximum(abs, lin - Rturn), symplectic_error(M))
    @printf("           constant part max = %.3e\n", maximum(abs, constant_part.(M)))
end


section("5. Chromaticity series",
        "chromaticity_series(M) reads the order-3 one-turn DA map already built.",
        "xi_x and xi_y are the derivatives d^n Q / d delta^n. The Lie column is chromaticity of the factor map.")

M3 = maps[3]
t_ch = @elapsed ch = chromaticity_series(M3)
@printf("  chromaticity_series in %.3f s    nmax = %d\n", t_ch, ch.nmax)
@printf("  Qx series (powers of delta): %s\n", join((@sprintf("%.6e", q) for q in ch.x), "  "))
@printf("  Qy series (powers of delta): %s\n", join((@sprintf("%.6e", q) for q in ch.y), "  "))
@printf("  xi_x: %s\n", join((@sprintf("%.6e", q) for q in ch.xi_x), "  "))
@printf("  xi_y: %s\n", join((@sprintf("%.6e", q) for q in ch.xi_y), "  "))

t_lie = @elapsed Mlie = factor_to_map(ring_lie_factor(lat))
xi = chromaticity(Mlie)
@printf("  Lie factor map in %.3f s    chromaticity xi_x = %.6f    xi_y = %.6f\n", t_lie, xi[1], xi[2])
@printf("  |xi_x DA - xi_x Lie| = %.3e    |xi_y DA - xi_y Lie| = %.3e\n",
        abs(ch.xi_x[1] - xi[1]), abs(ch.xi_y[1] - xi[2]))

@printf("\n  wall time %.2f s\n", time() - STARTED)
println("Done.")
