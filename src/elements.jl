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
#     SBend(name, L, angle, k1)      sector dipole, h = angle/L = 1/rho,
#                                    optional gradient k1
#     SBend(name, L, angle, k1, e1, e2, gap, fint)
#                                    + pole-face angles e1/e2 [rad], full gap [m]
#                                    and fringe integral FINT; read only by the
#                                    edge maps of fringe.jl (body H unchanged)
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

"""
    AbstractElement

Supertype of all beam-line elements. Element types are plain immutable
data (no physics). Every subtype has a `name::String`; elements with a
length have a field `L` [m] (see [`elem_length`](@ref)).

Physics lives elsewhere: `hamiltonian(e, ref)` (Hamiltonians),
`linear_matrix(e, ref)` (Linear), `element_lie_map(e, ref)` (Tracking),
`exact_track_element` (ExactIntegrator), `edge_map` (Fringe).
Because elements are immutable values, equal elements share cache entries
(maps are computed once per distinct element and energy).

Subtypes: [`Drift`](@ref), [`SBend`](@ref), [`Quadrupole`](@ref),
[`Sextupole`](@ref), [`RFCavity`](@ref), [`Cavity`](@ref), [`Marker`](@ref).
"""
abstract type AbstractElement end

@doc raw"""
    Drift(name, L)

Field-free straight section, length `L` [m].
Hamiltonian ``H = e(\delta) - p_s``; linear matrix with ``R_{12} = R_{34} = L``,
``R_{56} = L/\gamma_0^2``; Dragt–Finn generators ``f_3 = -LH_3``, ``f_4 = -LH_4``.
Structure-file card: `NAME drift L = 0.5`.
"""
struct Drift <: AbstractElement
    name::String
    L::Float64
end

@doc raw"""
    SBend(name, L, angle, k1 = 0, e1 = 0, e2 = 0, gap = 0, fint = 0)

Sector dipole: the reference orbit is an arc of radius ``\rho = L/\text{angle}``.

# Fields
- `L::Float64` – arc length [m].
- `angle::Float64` – bending angle [rad]; curvature ``h = \text{angle}/L`` [1/m] ([`curvature`](@ref)).
- `k1::Float64` – normalised gradient ``(\partial B_y/\partial x)/(B\rho)`` [1/m²] (combined-function bend).
- `e1::Float64`, `e2::Float64` – entrance / exit pole-face rotation angles [rad].
- `gap::Float64` – full magnet gap [m].
- `fint::Float64` – fringe-field integral FINT [1].

The body Hamiltonian uses only `L`, `angle`, `k1`:
``H = e(\delta) - (1+hx)p_s + hx + h^2x^2/2 + k_1(x^2 - y^2)/2``.
`e1, e2, gap, fint` are read only by the edge maps of the Fringe module.
Constructors with 3 and 4 arguments set the remaining fields to 0.

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
Hamiltonian ``H = e(\delta) - p_s + k_1(x^2 - y^2)/2``.
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
Hamiltonian ``H = e(\delta) - p_s + k_2(x^3 - 3xy^2)/6``.
Card: `SD sextupole L = 0.12 k2 = -3.1`.
"""
struct Sextupole <: AbstractElement
    name::String
    L::Float64
    k2::Float64
end

@doc raw"""
    RFCavity(name, voltage_MV, harmonic, phi_s)

Thin (zero-length) RF cavity: a single energy kick per passage.

# Fields
- `voltage_MV::Float64` – peak voltage ``V`` [MV] (use ``V\cdot T`` to mimic a finite cavity).
- `harmonic::Int` – RF harmonic number ``h`` (``f_\text{RF} = hf_0``).
- `phi_s::Float64` – synchronous phase ``\varphi_s`` [rad]; ``\Delta E = ZV\sin(\varphi_s + \varphi)`` [MeV].

Tracked by `rf_kick` (numbers) and `rf_kick_taylor` (polynomial map).
Card: `RF rfcavity V = 0.01 h = 1 phi_s = 0`.
"""
struct RFCavity <: AbstractElement
    name::String
    voltage_MV::Float64
    harmonic::Int
    phi_s::Float64
end

@doc raw"""
    Cavity(name, L, voltage_MV, harmonic, phi_s)

Finite-length RF cavity with purely longitudinal on-axis field
``E_s(s, t) = -(V/L)\sin\Phi``, ``\Phi = \omega(s - L/2)/v_0 - \omega l/v - \varphi_s``, ``\omega = h\,2\pi f_0``.

# Fields
- `L::Float64` – length [m]; `voltage_MV::Float64` – ``V`` [MV];
  `harmonic::Int` – ``h``; `phi_s::Float64` – ``\varphi_s`` [rad].

For maps and linear optics (RF off) it is a drift of length ``L``. Integrated with
the exact Hamiltonian ``e(\delta) - p_s + A\cos\Phi`` by `exact_track_element`; energy gain per passage is
``ZVT\sin(\varphi_s + \varphi)`` with transit factor ``T`` (`transit_factor`).
Transverse RF fields are neglected. Card: `CAV cavity L = 0.5 V = 0.01 h = 1 phi_s = 0.3`.
"""
struct Cavity <: AbstractElement
    name::String
    L::Float64
    voltage_MV::Float64
    harmonic::Int
    phi_s::Float64
end

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
(markers, thin `RFCavity`). Used by `circumference`, slicing and maps.
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

Short type name used in printing and error messages
("Drift", "SBend", "Quadrupole", "Sextupole", "RFCavity", "Cavity", "Marker").
"""
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
    e isa SBend && any(!iszero, (e.e1, e.e2, e.gap, e.fint)) &&
        (s *= ", e1=$(e.e1), e2=$(e.e2), gap=$(e.gap), fint=$(e.fint)")
    e isa Quadrupole && (s *= ", k1=$(e.k1) 1/m^2")
    e isa Sextupole && (s *= ", k2=$(e.k2) 1/m^3")
    (e isa RFCavity || e isa Cavity) && (s *= ", V=$(e.voltage_MV) MV, h=$(e.harmonic), phi_s=$(e.phi_s)")
    print(io, s, ")")
end

end # module
