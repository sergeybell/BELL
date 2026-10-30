# =============================================================================
# elements.jl - beam-line element types (data only)
# =============================================================================
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# WHAT
#   Plain data types describing beam-line elements, in SI units and in the
#   normalised strengths that enter the Hamiltonians:
#     Drift(name, L)
#     SBend(name, L, angle, k1)      sector dipole, no fringe fields,
#                                    h = angle/L = 1/rho, optional gradient k1
#     Quadrupole(name, L, k1)        k1 = (dBy/dx)/(B rho)   [1/m^2], k1>0 focuses x
#     Sextupole(name, L, k2)         k2 = (d2By/dx2)/(B rho) [1/m^3]
#     RFCavity(name, V, h, phi_s)    thin cavity, voltage V [MV] (per unit charge)
#     Cavity(name, L, V, h, phi_s)   finite-length RF cavity
#     Marker(name)
#
# WHY
#   The element types carry no physics.  The Hamiltonian of each type lives
#   in hamiltonians.jl.  Sequence, named parameters and Lattice live in
#   structure.jl.
#
# =============================================================================
module Elements

export AbstractElement, Drift, SBend, Quadrupole, Sextupole, RFCavity, Cavity, Marker,
       curvature, elem_length, kind_string

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
Thin RF cavity.  Energy gain per passage
    dE = Z * e * V * sin(phi_s + phi)    [MeV for V in MV]
implemented as a kick in tracking.jl.
"""
struct RFCavity <: AbstractElement
    name::String
    voltage_MV::Float64
    harmonic::Int
    phi_s::Float64
end

"""
Finite-length RF cavity, longitudinal on-axis field only:
    E_s(s,t) = -(V/L) sin(Phi),  Phi = omega (s - L/2)/v0 - omega l/v - phi_s.
For maps / linear optics (RF off) it behaves as a drift of length L.
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

end # module
