# =============================================================================
# linear.jl - linear transfer matrices, periodic Twiss, optics along s
# =============================================================================
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# The linear matrix of an element is R = exp(L S Hess H). That is the
# algorithm; no textbook element formulas. Other readings of the same
# linear map already exist: the matrix R of a Dragt-Finn factor, and
# linear_part of a Lie-flow map.
#
# Three functions, three jobs:
#   periodic_optics  Twiss and linear dispersion of a one-turn matrix
#   closed_orbit     (in map_readings.jl) fixed point of a nonlinear map at δ ≠ 0
#   optics_along     the same linear functions transported along s, integer tune
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

@doc raw"""
    linear_matrix(e, ref, L = elem_length(e)) -> Matrix{Float64} (6×6)

Exact linear transfer matrix of element `e` (or of a slice of length `L`).

# Algorithm
```math
R = \exp(L\,S A),\qquad A_{ij} = \left.\frac{\partial^2 H}{\partial z_i\,\partial z_j}\right|_{z=0},
```
where ``H`` = `hamiltonian(e, ref)` (only its quadratic part enters) and ``S`` =
[`SYMPLECTIC_S`](@ref). ``SA`` is a Hamiltonian matrix, so ``R`` is symplectic up to
round-off. No textbook element formulas are used. ``SA`` is cached per (element
value, ``p_0``). Markers, thin RF and zero-length slices give the identity.

The same linear map can also be read as the matrix `R` of a Dragt-Finn factor
or as [`linear_part`](@ref Bell.DiffAlgebra.linear_part) of a Lie-flow map.
"""
function linear_matrix(e::AbstractElement, ref::RefParticle, L::Real = elem_length(e))
    (e isa Marker || e isa RFCavity || L == 0) && return Matrix{Float64}(I, 6, 6)
    return exp(L * _SA(e, ref))
end

"""
    ring_matrix(lattice) -> Matrix{Float64} (6×6)

One-turn matrix R = R_n ⋯ R_1 (product of [`linear_matrix`](@ref) in element order).
"""
ring_matrix(lattice::Lattice) = foldl((M, e) -> linear_matrix(e, lattice.ref) * M, lattice.elements; init = Matrix{Float64}(I, 6, 6))

function _plane(m)
    c = (m[1, 1] + m[2, 2]) / 2
    stable = abs(c) < 1
    stable || return (stable = false, cosmu = c, beta = NaN, alpha = NaN, mu = NaN)
    s = sign(m[1, 2]) * sqrt(1 - c^2)
    mu = mod(atan(s, c), 2pi)
    return (stable = true, cosmu = c, beta = m[1, 2] / s, alpha = (m[1, 1] - m[2, 2]) / (2s), mu = mu)
end

@doc raw"""
    periodic_optics(R, ref, path_length) -> NamedTuple

Periodic (closed) linear optics of the one-turn matrix `R`. This is Twiss and
linear dispersion of a matrix. The off-momentum fixed point of a nonlinear map
is [`closed_orbit`](@ref Bell.MapReadings.closed_orbit). Transport of the same
linear functions along ``s`` is [`optics_along`](@ref).

# Returns
- `x`, `y` – per plane `(stable, cosmu, beta [m], alpha, mu [rad])` from the 2×2
  diagonal block ``m``: ``\cos\mu = \tfrac12(m_{11}+m_{22})``,
  ``\sin\mu = \operatorname{sign}(m_{12})\sqrt{1-\cos^2\mu}``, ``\beta = m_{12}/\sin\mu``,
  ``\alpha = (m_{11}-m_{22})/(2\sin\mu)``. Fractional tune is ``\mu/2\pi``.
- `D4` – periodic dispersion ``(D_x, D_x', D_y, D_y') = (I - R_4)^{-1}R_{1:4,6}``.
- `l_delta_closed` – ``(l|\delta)`` on the dispersion orbit ``= R_{56} + R_{5,1:4}\cdot D_4`` [m].
- `alpha_c` – momentum compaction ``1/\gamma_0^2 - (l|\delta)_\text{closed}/\mathrm{path\_length}``.
- `eta_slip` – ``\alpha_c - 1/\gamma_0^2``; `gamma_t` – ``1/\sqrt{\alpha_c}`` (NaN if ``\alpha_c \le 0``).
- `coupling` – max |x–y block| of `R` (0 for an uncoupled ring); `R`.
"""
function periodic_optics(R::AbstractMatrix, ref::RefParticle, path_length::Real)
    px = _plane(R[1:2, 1:2]); py = _plane(R[3:4, 3:4])
    coupling = maximum(abs, vcat(vec(R[1:2, 3:4]), vec(R[3:4, 1:2])))
    D4 = (I - R[1:4, 1:4]) \ R[1:4, 6]
    l_delta = R[5, 6] + dot(R[5, 1:4], D4)
    alpha_c = 1 / ref.gamma^2 - l_delta / path_length
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
    optics_along(lattice, po; ds = 0.02) -> NamedTuple

Propagate the periodic Twiss functions and dispersion along the ring. This is
the linear optics as functions of ``s``. Integer parts of the tunes live here:
`Qx`, `Qy` accumulate the phase slice by slice. [`periodic_optics`](@ref) gives
only the fractional tune of the one-turn matrix. [`closed_orbit`](@ref Bell.MapReadings.closed_orbit)
is the nonlinear off-momentum fixed point, not this transport.

# Arguments
- `po` – result of [`periodic_optics`](@ref) (start values at s = 0).
- `ds` – maximum slice length [m].

# Algorithm
Each element is cut into ⌈L/ds⌉ slices with matrices `linear_matrix(e, ref, L/n)`.
β, α are transported with the standard 2×2 Twiss transformation. The phase
advance is accumulated slice by slice, and the dispersion vector (D4, 0, 1) is
multiplied by the slice matrix.

# Returns
`(data, bounds, Qx, Qy, closure, path_length)`:
- `data` – vectors `s, betx, alfx, bety, alfy, D, Dp, mux, muy`;
- `bounds` – (s_start, s_end, element) for each element;
- `Qx`, `Qy` – total tunes, integer part included;
- `closure` – differences between the values after one turn and the start values
  (they should be at round-off level);
- `path_length` – total length [m].
"""
function optics_along(lattice::Lattice, po; ds = 0.02)
    ref = lattice.ref
    bx, ax, by, ay = po.x.beta, po.x.alpha, po.y.beta, po.y.alpha
    eta = vcat(po.D4, 0.0, 1.0)
    s = 0.0; mx = 0.0; my = 0.0
    data = (s = [0.0], betx = [bx], alfx = [ax], bety = [by], alfy = [ay], D = [eta[1]], Dp = [eta[2]],
            mux = [0.0], muy = [0.0])
    bounds = Tuple{Float64,Float64,AbstractElement}[]
    for e in lattice.elements
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
            push!(data.s, s); push!(data.betx, bx); push!(data.alfx, ax); push!(data.bety, by); push!(data.alfy, ay)
            push!(data.D, eta[1]); push!(data.Dp, eta[2]); push!(data.mux, mx); push!(data.muy, my)
        end
        push!(bounds, (s0, s, e))
    end
    closure = (betx = abs(bx - po.x.beta) / po.x.beta, alfx = abs(ax - po.x.alpha),
               bety = abs(by - po.y.beta) / po.y.beta, alfy = abs(ay - po.y.alpha),
               D = abs(eta[1] - po.D4[1]), Dp = abs(eta[2] - po.D4[2]),
               Qx_frac = abs(mod(mx / 2pi, 1) - po.x.mu / 2pi), Qy_frac = abs(mod(my / 2pi, 1) - po.y.mu / 2pi))
    return (data = data, bounds = bounds, Qx = mx / 2pi, Qy = my / 2pi, closure = closure, path_length = s)
end

end # module
