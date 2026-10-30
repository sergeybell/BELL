# =============================================================================
# LinearOptics.jl - linear transfer matrices, closed Twiss, dispersion, tunes
# =============================================================================
#
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# WHAT
#   * hessian(H): A = d2H/dz dz at z = 0  (H2 = z^T A z / 2)
#   * linear_matrix(e, ref, L): R = exp(L S A) - the exact linear map of the
#     element, obtained from its Hamiltonian (no textbook element formulas);
#     symplectic by construction since S A is a Hamiltonian matrix.
#   * ring_matrix(lat): product of all element matrices (one turn, 6x6).
#   * periodic_optics(R, ...): closed (periodic) Twiss beta, alpha in x and y,
#     fractional tunes, periodic dispersion (D, D', Dy, Dy') as the fixed
#     point of the 4D map at fixed delta, and the momentum compaction
#         alpha_c = 1/gamma0^2 - (l|delta)_closed / C ,
#     (l|delta)_closed = R56 + R5,1:4 * D4  (path term on the dispersion orbit).
#     Derivation: l = -v dt and dt = L_path/v - C/v0 give
#     dl/ddelta = -C (alpha_c - 1/gamma0^2) = -C * eta_slip.
#   * optics_along(lat): Twiss, dispersion and phase advance propagated in
#     thin slices along s (for plots and the integer part of the tune), and
#     the CLOSURE check (values after one turn = start values).
#
# WHY
#   Phase 1 of the roadmap: verify the lattice/units/conventions on the
#   simplest observables (tunes, beta, D, alpha_c) before nonlinear work.
#
# CONVENTIONS
#   beta [m], alpha [1], D [m], D' [1]; D is the orbit change per unit
#   delta = dp/p.  At linear order x' = px, so these are the same in the
#   user coordinates (x, x', y, y', phi, delta).
#
# WHAT YOU SEE
#   Nothing printed here; scripts/01 prints the results.
# =============================================================================
module Linear

using LinearAlgebra, TaylorSeries
using ..Coordinates, ..Elements, ..Structure, ..Hamiltonians, ..LieAlgebra

export hessian, linear_matrix, ring_matrix, periodic_optics, optics_along

"""
    hessian(H) -> Matrix{Float64} (6×6)

Hessian A = ∂²H/∂z_i∂z_j at z = 0, so that the quadratic part of H is
H₂ = zᵀ A z / 2.
"""
hessian(H::TaylorN) = [constant_term(differentiate(differentiate(H, i), j)) for i in 1:6, j in 1:6]

const _HCACHE = Dict{Any,Matrix{Float64}}()
function _SA(e::AbstractElement, ref::RefParticle)
    get!(_HCACHE, (e, ref.p0_MeV)) do
        SYMPLECTIC_S * hessian(hamiltonian(e, ref))
    end
end

"""
    linear_matrix(e, ref, L = elem_length(e)) -> Matrix{Float64} (6×6)

Exact linear transfer matrix of element `e` (or of a slice of length `L`).

# Algorithm
R = exp(L · S · A), where A = [`hessian`](@ref)(hamiltonian(e, ref)) and S =
[`SYMPLECTIC_S`](@ref). S·A is a Hamiltonian matrix, so R is symplectic up to
round-off. No textbook element formulas are used; the matrix follows from
the same H as the nonlinear maps. S·A is cached per (element value, p0).
Markers, thin RF and zero-length slices give the identity.
"""
function linear_matrix(e::AbstractElement, ref::RefParticle, L::Real = elem_length(e))
    (e isa Marker || e isa RFCavity || L == 0) && return Matrix{Float64}(I, 6, 6)
    return exp(L * _SA(e, ref))
end

"""
    ring_matrix(lat) -> Matrix{Float64} (6×6)

One-turn matrix R = R_n ⋯ R_1 (product of [`linear_matrix`](@ref) in element order).
"""
ring_matrix(lat::Lattice) = foldl((M, e) -> linear_matrix(e, lat.ref) * M, lat.elements; init = Matrix{Float64}(I, 6, 6))

function _plane(m)
    c = (m[1, 1] + m[2, 2]) / 2
    stable = abs(c) < 1
    stable || return (stable = false, cosmu = c, beta = NaN, alpha = NaN, mu = NaN)
    s = sign(m[1, 2]) * sqrt(1 - c^2)
    mu = mod(atan(s, c), 2pi)
    return (stable = true, cosmu = c, beta = m[1, 2] / s, alpha = (m[1, 1] - m[2, 2]) / (2s), mu = mu)
end

"""
    periodic_optics(R, ref, C) -> NamedTuple

Periodic (closed) linear optics of the one-turn matrix `R` of a ring of
circumference `C` [m].

# Returns
- `x`, `y` – per plane: `(stable, cosmu, beta [m], alpha, mu [rad])` from the
  2×2 diagonal blocks: cos μ = tr/2, sin μ = sign(R12)√(1 − cos²μ),
  β = R12/sin μ, α = (R11 − R22)/(2 sin μ).
- `D4` – periodic dispersion (D_x, D_x', D_y, D_y') = (I − R₄)⁻¹ R_{1:4,6}.
- `l_delta_closed` – (l|δ) on the dispersion orbit = R56 + R_{5,1:4}·D4 [m].
- `alpha_c` – momentum compaction 1/γ0² − (l|δ)_closed/C.
- `eta_slip` – α_c − 1/γ0²; `gamma_t` – 1/√α_c (NaN if α_c ≤ 0).
- `coupling` – max |x–y block| of R (0 for an uncoupled ring); `R`.
"""
function periodic_optics(R::AbstractMatrix, ref::RefParticle, C::Real)
    px = _plane(R[1:2, 1:2]); py = _plane(R[3:4, 3:4])
    coupling = maximum(abs, vcat(vec(R[1:2, 3:4]), vec(R[3:4, 1:2])))
    D4 = (I - R[1:4, 1:4]) \ R[1:4, 6]
    l_delta = R[5, 6] + dot(R[5, 1:4], D4)
    alpha_c = 1 / ref.gamma^2 - l_delta / C
    return (x = px, y = py, D4 = D4, l_delta_closed = l_delta, alpha_c = alpha_c,
            eta_slip = alpha_c - 1 / ref.gamma^2, gamma_t = alpha_c > 0 ? 1 / sqrt(alpha_c) : NaN,
            coupling = coupling, R = R)
end

_twiss_step(m, b, a) = begin
    g = (1 + a^2) / b
    (m[1, 1]^2 * b - 2m[1, 1] * m[1, 2] * a + m[1, 2]^2 * g,
     -m[1, 1] * m[2, 1] * b + (m[1, 1] * m[2, 2] + m[1, 2] * m[2, 1]) * a - m[1, 2] * m[2, 2] * g)
end
_dmu(m, b, a) = atan(m[1, 2], m[1, 1] * b - m[1, 2] * a)

"""
    optics_along(lat, po; ds = 0.02) -> NamedTuple

Propagate the periodic Twiss functions and dispersion along the ring.

# Arguments
- `po` – result of [`periodic_optics`](@ref) (start values at s = 0).
- `ds` – maximum slice length [m].

# Algorithm
Each element is cut into ⌈L/ds⌉ slices with matrices `linear_matrix(e, ref, L/n)`.
β, α are transported with the standard 2×2 Twiss transformation. The phase
advance is accumulated slice by slice, and the dispersion vector (D4, 0, 1) is
multiplied by the slice matrix.

# Returns
`(data, bounds, Qx, Qy, closure, C)`:
- `data` – vectors `s, betx, alfx, bety, alfy, D, Dp, mux, muy`;
- `bounds` – (s_start, s_end, element) for each element;
- `Qx`, `Qy` – total tunes, integer part included;
- `closure` – differences between the values after one turn and the start values
  (they should be at round-off level);
- `C` – total length.
"""
function optics_along(lat::Lattice, po; ds = 0.02)
    ref = lat.ref
    bx, ax, by, ay = po.x.beta, po.x.alpha, po.y.beta, po.y.alpha
    eta = vcat(po.D4, 0.0, 1.0)
    s = 0.0; mx = 0.0; my = 0.0
    out = (s = [0.0], betx = [bx], alfx = [ax], bety = [by], alfy = [ay], D = [eta[1]], Dp = [eta[2]],
           mux = [0.0], muy = [0.0])
    bounds = Tuple{Float64,Float64,AbstractElement}[]
    for e in lat.elements
        L = elem_length(e)
        s0 = s
        n = L == 0 ? 0 : max(1, ceil(Int, L / ds))
        for k in 1:n
            m = linear_matrix(e, ref, L / n)
            mx += _dmu(m[1:2, 1:2], bx, ax); my += _dmu(m[3:4, 3:4], by, ay)
            bx, ax = _twiss_step(m[1:2, 1:2], bx, ax)
            by, ay = _twiss_step(m[3:4, 3:4], by, ay)
            eta = m * eta
            s += L / n
            push!(out.s, s); push!(out.betx, bx); push!(out.alfx, ax); push!(out.bety, by); push!(out.alfy, ay)
            push!(out.D, eta[1]); push!(out.Dp, eta[2]); push!(out.mux, mx); push!(out.muy, my)
        end
        push!(bounds, (s0, s, e))
    end
    closure = (betx = abs(bx - po.x.beta) / po.x.beta, alfx = abs(ax - po.x.alpha),
               bety = abs(by - po.y.beta) / po.y.beta, alfy = abs(ay - po.y.alpha),
               D = abs(eta[1] - po.D4[1]), Dp = abs(eta[2] - po.D4[2]),
               Qx_frac = abs(mod(mx / 2pi, 1) - po.x.mu / 2pi), Qy_frac = abs(mod(my / 2pi, 1) - po.y.mu / 2pi))
    return (data = out, bounds = bounds, Qx = mx / 2pi, Qy = my / 2pi, closure = closure, C = s)
end

end # module
