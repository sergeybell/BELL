# element_hamiltonians.jl
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# s-Hamiltonian of a drift, a quadrupole, a sector bend, a sextupole,
# a thin RF cavity, and a thick RF cavity. The FODO file supplies the
# drift, the quadrupoles, and the bends. The sextupole and the cavities
# are constructed here.
#
#   julia --project=examples examples/element_hamiltonians/run.jl

using Bell, Printf

include(joinpath(@__DIR__, "..", "common", "print_util.jl"))

banner("Element s-Hamiltonians (FODO, plus a sextupole and RF)")

lat = load_structure(joinpath(@__DIR__, "..", "lattices", "fodo.bell"))
ref = lat.ref
C = path_length(lat)
z0 = zeros(6)

function mono(powers...)
    e = zeros(Int, 6)
    for (i, p) in pairs(powers)
        e[i] = p
    end
    return e
end

function print_polynomial(e, ref)
    println("  ", e)
    println("  ", hamiltonian_formula(e))
    if e isa RFCavity
        e.L > 0 || error("print_polynomial is for elements that have a continuous Hamiltonian")
        H0 = hamiltonian(e, ref, C, z0)
        Hp = hamiltonian(e, ref, C)
    else
        H0 = hamiltonian(e, ref, z0)
        Hp = hamiltonian(e, ref)
    end
    @printf("  H(z = 0) = %.8f     Lie polynomial (constant removed):\n", H0)
    println("  ", poly_string(Hp))
    println()
    return Hp
end


section("1. Reference orbit",
        "energy_term is the longitudinal term in the pair (l, delta). It is 0 at delta = 0.",
        "At z = 0, ps = 1, so a pure drift evaluates to energy_term(0) - 1 = -1.",
        "hamiltonian(e, ref) returns the Lie-path polynomial and drops that constant.")

@printf("  species %s, Ekin = %.1f MeV, beta = %.8f, gamma = %.8f\n",
        ref.particle.name, ref.Ekin_MeV, ref.beta, ref.gamma)
@printf("  energy_term(0) = %.3e\n", energy_term(0.0, ref))
@printf("  velocity_ratio(0) = %.16f\n", velocity_ratio(0.0, ref))
@printf("  1/gamma^2 = %.8f    (quadratic coefficient of energy_term is half of this)\n",
        1 / ref.gamma^2)


section("2. Drift, quadrupoles, and sector bends from the FODO file",
        "H does not contain the element length. The length is the interval over which the map integrates H.",
        "k1 > 0 focuses in x: the quadrupole term is (k1/2) (x^2 - y^2).")

by_kind = Dict{DataType,Any}()
for e in unique_elements(lat)
    get!(by_kind, typeof(e), e)
end

Hd = print_polynomial(by_kind[Drift], ref)
Hq = print_polynomial(by_kind[Quadrupole], ref)
Hb = print_polynomial(by_kind[SBend], ref)

qd = lat.elements[findfirst(e -> e isa Quadrupole && e.k1 < 0, lat.elements)]
bend = by_kind[SBend]
@printf("  drift delta^2 coefficient = %.8f    1/(2 gamma^2) = %.8f\n",
        coeff(Hd, mono(0, 0, 0, 0, 0, 2)), 1 / (2 * ref.gamma^2))
@printf("  %s k1 = %.4f 1/m^2,  x^2 coefficient = %.4f,  y^2 coefficient = %.4f\n",
        by_kind[Quadrupole].name, by_kind[Quadrupole].k1,
        coeff(Hq, mono(2, 0, 0, 0, 0, 0)), coeff(Hq, mono(0, 0, 2, 0, 0, 0)))
@printf("  %s k1 = %.4f 1/m^2  (defocusing in x; the x^2 sign flips with k1)\n", qd.name, qd.k1)
@printf("  %s curvature h = angle/L = %.6f 1/m,  x*delta coefficient = %.6f\n",
        bend.name, curvature(bend), coeff(Hb, mono(1, 0, 0, 0, 0, 1)))


section("3. Sextupole",
        "The FODO file has no sextupole. This element is constructed for the formula.",
        "The term is (k2/6) (x^3 - 3 x y^2).")

sx = Sextupole("SF", 0.20, 1.2)
Hs = print_polynomial(sx, ref)
@printf("  k2 = %.4f 1/m^3,  x^3 coefficient = %.6f,  k2/6 = %.6f\n",
        sx.k2, coeff(Hs, mono(3, 0, 0, 0, 0, 0)), sx.k2 / 6)
@printf("  x*y^2 coefficient = %.6f,  -k2/2 = %.6f\n",
        coeff(Hs, mono(1, 0, 2, 0, 0, 0)), -sx.k2 / 2)


section("4. Thin RF cavity",
        "A thin cavity has no continuous Hamiltonian. Its generator is rf_potential,",
        "and the map is the kick. hamiltonian on a thin cavity throws.")

thin = RFCavity("RFthin", 0.05, lat.harmonic, 0.0)
println("  ", thin)
println("  ", hamiltonian_formula(thin))
k = rf_wavenumber(thin, ref, C)
@printf("  rf_wavenumber k = 2 pi h beta / C = %.8e 1/m\n", k)
@printf("  rf_potential at z = 0 : %.8e\n", rf_potential(thin, ref, C; z = z0))
z_l = [0.0, 0.0, 0.0, 0.0, 0.25, 0.0]
phase = k * z_l[index_l] / (ref.beta * velocity_ratio(z_l[index_delta], ref))
@printf("  at l = %.2f m, delta = 0: phase advance k l / (beta v/v0) = %.6f rad\n",
        z_l[index_l], phase)
@printf("  rf_potential there         : %.8e\n", rf_potential(thin, ref, C; z = z_l))
print("  hamiltonian(thin, ref, C, z) : ")
try
    hamiltonian(thin, ref, C, z0)
    println("returned a value")
catch err
    println(sprint(showerror, err))
end


section("5. Thick RF cavity",
        "L > 0: H = energy_term(delta) - ps + rf_potential/L.",
        "phi_s = 0, so the reference particle is not accelerated and the polynomial has no linear term.",
        "The constant cos(phi_s) shifts H(0) away from -1. The Lie polynomial drops that constant.")

thick = RFCavity("RFthick", 0.05, lat.harmonic, 0.0, 0.40)
print_polynomial(thick, ref)
@printf("  rf_potential(0)/L = %.8e 1/m\n", rf_potential(thick, ref, C; z = z0) / thick.L)
@printf("  H(0) + 1          = %.8e   (the constant RF piece)\n",
        hamiltonian(thick, ref, C, z0) + 1)

println()
println("Done.")
