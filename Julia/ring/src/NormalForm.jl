# =============================================================================
# NormalForm.jl - normal-form analysis (PHASE 1: linear part + plan/stub)
# =============================================================================
#
# WHAT (implemented now)
#   linear_normalizing_matrix(po)  A = blockdiag(A_x, A_y, I) with
#       A_u = [ sqrt(beta)          0          ;
#              -alpha/sqrt(beta)    1/sqrt(beta) ]
#     so that A^-1 R A is a pure rotation by 2 pi Q_u in each transverse
#     plane (uncoupled lattice; dispersion is not removed in phase 1).
#   rotation_residual(R, po)  how far A^-1 R A (transverse 4x4 block) is from
#     an exact rotation - a check of the closed Twiss.
#
# WHAT (planned, phase 4 - see PLAN.md)
#   nonlinear_normal_form(F::LieFactor, po): from the ring factor
#   M = R o exp(:f3:) compute, in normalised coordinates, a generator F3 such
#   that exp(:F3:) M exp(-:F3:) = R o exp(:h3:) with h3 containing only
#   terms commuting with R (for order 3 in a ring without resonances:
#   only delta-dependent terms -> chromaticity / dispersion of tunes,
#   then order 4 -> amplitude-dependent tune shifts). Resonance driving
#   terms (h_jklmp) will be printed as by-products.
#
# WHY separate
#   Normal-form theory is an analysis layer on top of the map; keeping it in
#   its own module lets us change/validate it without touching tracking.
#
# WHAT YOU SEE
#   Nothing printed here.
# =============================================================================
module NormalForm

using LinearAlgebra
using ..LieOperators

export linear_normalizing_matrix, rotation_residual, nonlinear_normal_form

function linear_normalizing_matrix(po)
    A = Matrix{Float64}(I, 6, 6)
    for (k, p) in ((1, po.x), (3, po.y))
        sb = sqrt(p.beta)
        A[k:k+1, k:k+1] = [sb 0.0; -p.alpha/sb 1/sb]
    end
    return A
end

function rotation_residual(R, po)
    A = linear_normalizing_matrix(po)
    N = A \ R * A
    res = 0.0
    for (k, p) in ((1, po.x), (3, po.y))
        rot = [cos(p.mu) sin(p.mu); -sin(p.mu) cos(p.mu)]
        res = max(res, maximum(abs, N[k:k+1, k:k+1] - rot))
    end
    return res
end

"""
Plan (v0.6): 6D Lie normal form of the one-turn map WITH RF, to order 2 in
amplitudes, giving delta_eq(Jx, Jy, Js) as the constant term of delta after
the normalising transformation, and amplitude-dependent tunes.

Steps (not yet coded):
  1. Build one-turn Taylor map M including a thin RF kick (or finite cavity
     DA map) in canonical (x,px,y,py,l,delta), truncated at MAP_ORDER>=3.
  2. Linear normalisation A that block-diagonalises R into three rotations
     (Qx, Qy, Qs); remove linear dispersion into the longitudinal block.
  3. Dragt-Finn / Lie triangularisation to order 3: find generating function
     F3 so that A^{-1} M A = R o exp(:h3:) with h3 containing only invariants
     (terms in Jx, Jy, Js, and delta-powers that commute with R).
  4. Read off:
       - delta_eq = coeff of the constant (in angles) part of delta(J)
         after the inverse normalising map, expanded as
         c00 Js + c10 Jx + c01 Jy + ...  (compare note: -eta1/eta0 * <d^2>
         with <d^2>~Js, and -Lambda/eta0 with Lambda~-2pi(Jx xi_x+...)/C);
       - Qx(J), Qy(J), Qs(J) from h3 / from the linearised rotation angles.

Until this exists, use tracking + note formulas (apps/equilibrium_energy).
"""
nonlinear_normal_form(args...) =
    error("nonlinear_normal_form: planned (v0.6 plan in NormalForm.jl docstring); not implemented yet")

end # module
