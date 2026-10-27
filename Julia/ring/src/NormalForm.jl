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

nonlinear_normal_form(args...) =
    error("nonlinear_normal_form: planned for phase 4 (see PLAN.md); not implemented in phase 1")

end # module
