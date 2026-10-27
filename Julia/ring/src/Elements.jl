# =============================================================================
# Elements.jl - element types and the Lattice container
# =============================================================================
#
# WHAT
#   Plain data types describing beam-line elements, in SI units and in the
#   normalised strengths that enter the Hamiltonians:
#     Drift(name, L)
#     SBend(name, L, angle, k1)      sector dipole, no fringe fields,
#                                    h = angle/L = 1/rho, optional gradient k1
#     Quadrupole(name, L, k1)        k1 = (dBy/dx)/(B rho)   [1/m^2], k1>0 focuses x
#     Sextupole(name, L, k2)         k2 = (d2By/dx2)/(B rho) [1/m^3]
#     RFCavity(name, V, h, phi_s)    THIN cavity, voltage V [MV] (per unit charge),
#                                    harmonic h, synchronous phase phi_s [rad]
#     Cavity(name, L, V, h, phi_s)  finite-length RF cavity (phase 2)
#     Marker(name)
#   and Lattice = ordered vector of elements + reference particle.
#
# WHY
#   The element types carry NO physics; the physics (Hamiltonian) of each type
#   lives in Hamiltonians.jl.  Linear optics, Lie maps, Taylor maps and
#   tracking all dispatch on these types, so adding an element = one new
#   struct here + one Hamiltonian method.
#
# WHAT YOU SEE
#   Nothing printed; Base.show methods give one-line summaries.
# =============================================================================
module Elements

using ..Constants

export AbstractElement, Drift, SBend, Quadrupole, Sextupole, RFCavity, Cavity, Marker,
       Lattice, curvature, elem_length, circumference, total_bend_angle,
       kind_string, unique_elements, SEXT_FAMILIES, find_optics_maxima, place_sextupole_families,
       set_family_strengths, insert_cavity, set_cavity_phase, with_reference

abstract type AbstractElement end

struct Drift <: AbstractElement
    name::String
    L::Float64
end

"Sector bend: reference orbit is an arc of radius rho = L/angle; no fringe fields."
struct SBend <: AbstractElement
    name::String
    L::Float64
    angle::Float64
    k1::Float64
end
SBend(name, L, angle) = SBend(name, L, angle, 0.0)

struct Quadrupole <: AbstractElement
    name::String
    L::Float64
    k1::Float64
end

struct Sextupole <: AbstractElement
    name::String
    L::Float64
    k2::Float64
end

"""
Thin RF cavity (phase 1: type designed, kick implemented in Tracking.jl,
not used in the phase-1 optics).  Energy gain per passage
    dE = Z * e * V * sin(phi_s + phi)    [MeV for V in MV]
which for small kicks is  d(delta) ~ Z e V/(beta0^2 E0) sin(phi_s + phi).
"""
struct RFCavity <: AbstractElement
    name::String
    voltage_MV::Float64
    harmonic::Int
    phi_s::Float64
end

"""
Finite-length RF cavity (phase 2), longitudinal on-axis field only:
    E_s(s,t) = -(V/L) sin(Phi),  Phi = omega (s - L/2)/v0 - omega l/v - phi_s,
omega = h*omega0.  s-Hamiltonian (derivation in PLAN.md, code in ExactIntegrator.jl):
    H = e(delta) - ps + A cos(Phi),   A = (Z V/(p0 c)) * c/(omega L).
In the thin limit the energy gain is Z V sin(phi_s + phi) (same as RFCavity).
For linear optics / Lie maps (RF off) it behaves as a drift of length L.
"""
struct Cavity <: AbstractElement
    name::String
    L::Float64
    voltage_MV::Float64
    harmonic::Int
    phi_s::Float64
end
struct Marker <: AbstractElement
    name::String
end

elem_length(e::AbstractElement) = hasproperty(e, :L) ? e.L : 0.0
curvature(b::SBend) = b.angle / b.L
curvature(::AbstractElement) = 0.0

kind_string(::Drift) = "Drift"
kind_string(::SBend) = "SBend"
kind_string(::Quadrupole) = "Quadrupole"
kind_string(::Sextupole) = "Sextupole"
kind_string(::RFCavity) = "RFCavity"
kind_string(::Marker) = "Marker"
kind_string(::Cavity) = "Cavity"

function Base.show(io::IO, e::AbstractElement)
    s = "$(kind_string(e))($(e.name), L=$(elem_length(e)) m"
    e isa SBend && (s *= ", angle=$(e.angle) rad, k1=$(e.k1)")
    e isa Quadrupole && (s *= ", k1=$(e.k1) 1/m^2")
    e isa Sextupole && (s *= ", k2=$(e.k2) 1/m^3")
    (e isa RFCavity || e isa Cavity) && (s *= ", V=$(e.voltage_MV) MV, h=$(e.harmonic), phi_s=$(e.phi_s)")
    print(io, s, ")")
end

"A ring: ordered elements, reference particle, harmonic number used for phi."
struct Lattice
    name::String
    elements::Vector{AbstractElement}
    ref::RefParticle
    harmonic::Int
end

circumference(lat::Lattice) = sum(elem_length, lat.elements)
total_bend_angle(lat::Lattice) = sum(e -> e isa SBend ? e.angle : 0.0, lat.elements)

"Distinct elements (by value) in order of first appearance."
function unique_elements(lat::Lattice)
    seen = AbstractElement[]
    for e in lat.elements
        any(isequal(e), seen) || push!(seen, e)
    end
    return seen
end

# =============================================================================
# Lattice-editing helpers (no physics, no optimisation): they only build new
# Lattice objects.  Sextupoles default to the strengths in the lattice file
# (all zero for magnetic_2); nothing here changes strengths unless asked.
# =============================================================================

"Family names used by place_sextupole_families (D_x, beta_x, beta_y maxima)."
const SEXT_FAMILIES = ["SXD", "SXB", "SYB"]

"""
find_optics_maxima(opt; frac = 0.8, dmin = 0.5)
Local maxima of D_x, beta_x, beta_y in the optics table `opt.data` (from
LinearOptics.optics_along) that are >= frac * global maximum and lie where
D_x > dmin (a sextupole acts chromatically only through x = D delta).
"""
function find_optics_maxima(opt; frac = 0.8, dmin = 0.5)
    d = opt.data
    found = NamedTuple[]
    for (q, arr) in (("D_x", d.D), ("beta_x", d.betx), ("beta_y", d.bety))
        mask = d.D .> dmin
        g = maximum(arr[mask])
        last_s = -1e9
        for k in 2:length(arr)-1
            (mask[k] && arr[k] >= arr[k-1] && arr[k] >= arr[k+1] && arr[k] >= frac * g) || continue
            d.s[k] - last_s < 0.5 && continue
            last_s = d.s[k]
            push!(found, (quantity = q, s = d.s[k], value = arr[k], betx = d.betx[k], bety = d.bety[k], D = d.D[k]))
        end
    end
    return found
end

"""
place_sextupole_families(lat, opt) -> (new_lattice, fam_idx, rows, notes)
Puts a sextupole of family SXD / SXB / SYB (strength 0) into the nearest
short slot (existing sextupole slot or drift <= 0.5 m, within 1 m) of each
D_x / beta_x / beta_y maximum.  Priority D > beta_x > beta_y for shared slots.
Drifts longer than 0.15 m are split drift - 0.10 m sextupole - drift.
"""
function place_sextupole_families(lat::Lattice, opt; kw...)
    maxima = find_optics_maxima(opt; kw...)
    fam_of = Dict("D_x" => "SXD", "beta_x" => "SXB", "beta_y" => "SYB")
    assigned = Dict{Int,String}()
    notes = String[]; rows = NamedTuple[]
    for q in ("D_x", "beta_x", "beta_y"), m in filter(m -> m.quantity == q, maxima)
        best, bd = 0, 1e9
        for (i, (s0, s1, e)) in enumerate(opt.bounds)
            (e isa Sextupole || (e isa Drift && e.L <= 0.5)) || continue
            dist = abs((s0 + s1) / 2 - m.s)
            dist < bd && ((best, bd) = (i, dist))
        end
        if best == 0 || bd > 1.0
            push!(notes, "no slot within 1 m of $(q) max at s=$(round(m.s, digits=3))"); continue
        end
        if haskey(assigned, best)
            push!(notes, "$(q) max at s=$(round(m.s, digits=3)) m shares slot $(opt.bounds[best][3].name) with family $(assigned[best])")
            continue
        end
        assigned[best] = fam_of[q]
        push!(rows, (family = fam_of[q], quantity = q, s_max = m.s, value = m.value, slot = opt.bounds[best][3].name,
                     s_slot = (opt.bounds[best][1] + opt.bounds[best][2]) / 2, betx = m.betx, bety = m.bety, D = m.D))
    end
    elems = AbstractElement[]
    fam_idx = Dict(f => Int[] for f in SEXT_FAMILIES)
    for (i, e) in enumerate(lat.elements)
        if haskey(assigned, i)
            f = assigned[i]; L = elem_length(e)
            if L <= 0.15
                push!(elems, Sextupole(f, L, 0.0)); push!(fam_idx[f], length(elems))
            else
                push!(elems, Drift(e.name * "a", (L - 0.1) / 2)); push!(elems, Sextupole(f, 0.1, 0.0))
                push!(fam_idx[f], length(elems)); push!(elems, Drift(e.name * "b", (L - 0.1) / 2))
            end
        else
            push!(elems, e)
        end
    end
    return Lattice(lat.name * "+sext", elems, lat.ref, lat.harmonic), fam_idx, rows, notes
end

"New lattice with sextupole strengths k2[family] (families not given -> 0)."
function set_family_strengths(lat::Lattice, fam_idx, k2::Dict)
    elems = copy(lat.elements)
    for (f, idx) in fam_idx, i in idx
        elems[i] = Sextupole(f, elems[i].L, get(k2, f, 0.0))
    end
    return Lattice(lat.name, elems, lat.ref, lat.harmonic)
end

"Split drift `in_element` (first occurrence) as drift - cavity - drift."
function insert_cavity(lat::Lattice, cav::Cavity; in_element = "OBM")
    i = findfirst(e -> e.name == in_element, lat.elements)
    e = lat.elements[i]
    Ld = (e.L - cav.L) / 2
    elems = vcat(lat.elements[1:i-1], [Drift(e.name * "_u", Ld), cav, Drift(e.name * "_d", Ld)], lat.elements[i+1:end])
    return Lattice(lat.name * "+RF", elems, lat.ref, cav.harmonic), i + 1
end

"Same lattice with every finite cavity set to synchronous phase phi_s."
set_cavity_phase(lat::Lattice, phi_s) =
    Lattice(lat.name, [e isa Cavity ? Cavity(e.name, e.L, e.voltage_MV, e.harmonic, phi_s) : e for e in lat.elements],
            lat.ref, lat.harmonic)

"Same elements with a different reference particle."
with_reference(lat::Lattice, ref::RefParticle) = Lattice(lat.name, lat.elements, ref, lat.harmonic)

end # module