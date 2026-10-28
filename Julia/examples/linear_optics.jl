#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# linear_optics.jl — shared helper functions (transfer matrices), included by
# example 06 (03 is self-contained).  Not meant to be run on its own.
#
# Coordinates: (x, x') with x' = px/p0, s = path length along the reference orbit.
# Linear equation of motion  x'' + K(s) x = 0.
#   focusing quad      K = +k > 0
#   defocusing quad    K = -k < 0
#   sector dipole (x)  K = 1/ρ²  (weak focusing of the bend), plus dispersion term δ/ρ
#   sector dipole (y)  K = 0 (pure drift)
# A transfer matrix maps (x, x')_in -> (x, x')_out. For a sequence of elements
# E1, E2, ..., En traversed in this order, M = Mn * ... * M2 * M1.

drift(L) = [1.0 L;
            0.0 1.0]

"Thick focusing quadrupole, k > 0 [1/m²], length L [m]."
function quad_f(k, L)
    ω = sqrt(k)
    return [cos(ω * L)      sin(ω * L) / ω;
            -ω * sin(ω * L) cos(ω * L)]
end

"Thick defocusing quadrupole, k > 0 [1/m²] (the plane in which it defocuses)."
function quad_d(k, L)
    ω = sqrt(k)
    return [cosh(ω * L)     sinh(ω * L) / ω;
            ω * sinh(ω * L) cosh(ω * L)]
end

"Quadrupole with signed gradient k: k > 0 focuses in this plane, k < 0 defocuses."
quad(k, L) = k > 0 ? quad_f(k, L) : (k < 0 ? quad_d(-k, L) : drift(L))

"Sector dipole, horizontal plane: radius ρ, bend angle θ (L = ρθ)."
sbend_x(ρ, θ) = [cos(θ)      ρ * sin(θ);
                 -sin(θ) / ρ cos(θ)]

"Sector dipole, vertical plane: a drift of length ρθ."
sbend_y(ρ, θ) = drift(ρ * θ)

"3x3 horizontal matrix acting on (x, x', δ), δ = Δp/p0, for a sector dipole."
sbend_x3(ρ, θ) = [cos(θ)      ρ * sin(θ) ρ * (1 - cos(θ));
                  -sin(θ) / ρ cos(θ)     sin(θ);
                  0.0         0.0        1.0]

"Embed a 2x2 matrix (no dispersion generated) into 3x3 form."
to3(M) = [M zeros(2); 0.0 0.0 1.0]

"Product of a list of matrices in beam order: elements[1] is traversed first."
compose(elements) = foldl((M, E) -> E * M, elements)

"""
Periodic Twiss parameters from a one-turn (or one-cell) 2x2 matrix
  M = [cos μ + α sin μ    β sin μ;
       -γ sin μ           cos μ - α sin μ]
Returns (μ, β, α, γ). Sign of μ chosen so that β > 0.
"""
function twiss(M)
    c = (M[1, 1] + M[2, 2]) / 2
    abs(c) < 1 || error("unstable: |tr M|/2 = $(abs(c)) ≥ 1")
    μ = acos(c)
    sign(M[1, 2]) < 0 && (μ = 2π - μ)
    β = M[1, 2] / sin(μ)
    α = (M[1, 1] - M[2, 2]) / (2 * sin(μ))
    γ = -M[2, 1] / sin(μ)
    return (μ = μ, β = β, α = α, γ = γ)
end
