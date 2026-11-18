# =============================================================================
# elements.jl - beam-line element types (data only)
# =============================================================================
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# WHAT
#   Plain data for beam-line elements, in metres, radians, MeV, and T·m,
#   and in the normalised strengths that enter the Hamiltonians.
#   The field L is the element length in metres. Physics lives in
#   hamiltonians.jl. physics/fields.jl is the comparison formulation.
#   Sequence and Lattice live in structure.jl.
#
# =============================================================================
module Elements

export AbstractElement, Drift, SBend, Quadrupole, Sextupole, RFCavity, Marker,
       curvature, elem_length, kind_string

"""
    AbstractElement

Supertype of all beam-line elements. Element types are plain immutable
data (no physics). Every subtype has a `name::String`; elements with a
length have a field `L` [m] (see [`elem_length`](@ref)).

Physics lives elsewhere: `hamiltonian` (Hamiltonians), `element_flow`
(ElementFlow), `element_lie_map` (Tracking), `element_da_map` (DATracking),
`exact_track_element` (FieldTrack), `linear_matrix` (Linear), `edge_map`
(FringeMaps). `element_field` (Fields) is the comparison formulation.

Subtypes: [`Drift`](@ref), [`SBend`](@ref), [`Quadrupole`](@ref),
[`Sextupole`](@ref), [`RFCavity`](@ref), [`Marker`](@ref).
"""
abstract type AbstractElement end

@doc raw"""
    Drift(name, L)

Field-free straight section, length `L` [m].
Hamiltonian ``H = \mathrm{energy\_term}(\delta) - p_s``.
Structure-file card: `NAME drift L = 0.5`.
"""
struct Drift <: AbstractElement
    name::String
    L::Float64
end

@doc raw"""
    SBend(name, L, angle, k1 = 0, e1 = 0, e2 = 0, gap = 0, fint = 0)

Sector dipole: the reference orbit is an arc of radius ``\rho = L/\mathrm{angle}``.

The short constructors `SBend(name, L, angle)` and `SBend(name, L, angle, k1)`
set the pole-face angles `e1`, `e2`, the full gap and `fint` to zero. Those
zeros are edge parameters, not the bending angle. The bending angle is the
field `angle` and may be nonzero. Pole-face angles may be set to nonzero
values by the long constructor or by a structure-file card.

# Fields
- `L::Float64` – arc length [m].
- `angle::Float64` – bending angle [rad]; curvature ``h = \mathrm{angle}/L`` [1/m] ([`curvature`](@ref)).
- `k1::Float64` – normalised gradient ``(\partial B_y/\partial x)/(B\rho)`` [1/m²]. ``k_1 > 0`` focuses in ``x``.
- `e1::Float64`, `e2::Float64` – entrance / exit pole-face rotation angles [rad].
- `gap::Float64` – full magnet gap [m].
- `fint::Float64` – fringe-field integral FINT [1].

The body Hamiltonian uses `L`, `angle`, and `k1`. The comparison field
equations use the same three. `e1`, `e2`, `gap`, and `fint` are read by the
edge model (`Fringe`, `FringeMaps`).

Structure-file card: `BM sbend L = 5.33 angle = 22.5 deg k1 = 0 e1 = 0 gap = 0.05 fint = 0.5`.
"""
struct SBend <: AbstractElement
    name::String
    L::Float64
    angle::Float64
    k1::Float64
    e1::Float64
    e2::Float64
    gap::Float64
    fint::Float64
end
SBend(name, L, angle) = SBend(name, L, angle, 0.0)
SBend(name, L, angle, k1) = SBend(name, L, angle, k1, 0.0, 0.0, 0.0, 0.0)

@doc raw"""
    Quadrupole(name, L, k1)

Straight quadrupole. `L` [m], ``k_1 = (\partial B_y/\partial x)/(B\rho)`` [1/m²], ``k_1 > 0`` focuses in ``x``.
Hamiltonian ``H = \mathrm{energy\_term}(\delta) - p_s + k_1(x^2 - y^2)/2``.
Card: `QF quadrupole L = 0.21 k1 = 0.52`.
"""
struct Quadrupole <: AbstractElement
    name::String
    L::Float64
    k1::Float64
end

@doc raw"""
    Sextupole(name, L, k2)

Straight sextupole. `L` [m], ``k_2 = (\partial^2B_y/\partial x^2)/(B\rho)`` [1/m³].
Hamiltonian ``H = \mathrm{energy\_term}(\delta) - p_s + k_2(x^3 - 3xy^2)/6``.
Card: `SD sextupole L = 0.12 k2 = -3.1`.
"""
struct Sextupole <: AbstractElement
    name::String
    L::Float64
    k2::Float64
end

@doc raw"""
    RFCavity(name, voltage_MV, harmonic, phi_s)        # thin, L = 0
    RFCavity(name, voltage_MV, harmonic, phi_s, L)     # thick, length L [m]

RF cavity. A thick cavity (`L > 0`) has s-Hamiltonian
``H = H_{\mathrm{drift}} + G/L`` with [`rf_potential`](@ref) ``G`` and a
uniform gradient ``V/L`` (transit factor 1). `L = 0` is the thin limit; its
map is [`rf_kick`](@ref).

# Fields
- `L::Float64` – length [m] (`0` = thin).
- `voltage_MV::Float64` – peak voltage ``V`` [MV]. The energy change is ``ZV`` [MeV].
- `harmonic::Int` – RF harmonic number ``h``.
- `phi_s::Float64` – synchronous phase ``\varphi_s`` [rad].

Card: `RF rfcavity V = 0.01 h = 1 phi_s = 0`.
"""
struct RFCavity <: AbstractElement
    name::String
    voltage_MV::Float64
    harmonic::Int
    phi_s::Float64
    L::Float64
end
RFCavity(name::AbstractString, V::Real, h::Integer, phi_s::Real) = RFCavity(String(name), Float64(V), Int(h), Float64(phi_s), 0.0)

"""
    Marker(name)

Zero-length marker; identity map. Card: `M1 marker`.
"""
struct Marker <: AbstractElement
    name::String
end

"""
    elem_length(e) -> Float64

Length of an element [m]: the field `L` if it exists, otherwise 0
(markers, thin `RFCavity`). Used by `path_length`, slicing and maps.
"""
elem_length(e::AbstractElement) = hasproperty(e, :L) ? e.L : 0.0
"""
    curvature(e) -> h [1/m]

Curvature of the reference orbit in the element: h = angle/L for an
[`SBend`](@ref), 0 for every other element.
"""
curvature(b::SBend) = b.angle / b.L
curvature(::AbstractElement) = 0.0

"""
    kind_string(e) -> String

Short type name used in printing and error messages.
"""
kind_string(::Drift) = "Drift"
kind_string(::SBend) = "SBend"
kind_string(::Quadrupole) = "Quadrupole"
kind_string(::Sextupole) = "Sextupole"
kind_string(::RFCavity) = "RFCavity"
kind_string(::Marker) = "Marker"

_show_open(e::AbstractElement) = "$(kind_string(e))($(e.name), L=$(elem_length(e)) m"

function Base.show(io::IO, e::Drift)
    print(io, _show_open(e), ")")
end

function Base.show(io::IO, e::Marker)
    print(io, _show_open(e), ")")
end

function Base.show(io::IO, e::SBend)
    s = _show_open(e) * ", angle=$(e.angle) rad, k1=$(e.k1)"
    any(!iszero, (e.e1, e.e2, e.gap, e.fint)) &&
        (s *= ", e1=$(e.e1), e2=$(e.e2), gap=$(e.gap), fint=$(e.fint)")
    print(io, s, ")")
end

function Base.show(io::IO, e::Quadrupole)
    print(io, _show_open(e), ", k1=$(e.k1) 1/m^2)")
end

function Base.show(io::IO, e::Sextupole)
    print(io, _show_open(e), ", k2=$(e.k2) 1/m^3)")
end

function Base.show(io::IO, e::RFCavity)
    print(io, _show_open(e), ", V=$(e.voltage_MV) MV, h=$(e.harmonic), phi_s=$(e.phi_s))")
end

function Base.show(io::IO, e::AbstractElement)
    print(io, _show_open(e), ")")
end

end # module
