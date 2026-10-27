# =============================================================================
# Optimizers.jl - OPTIONAL optimisation / fitting routines (never called by
#                 default; only scripts/04_optional_sextupole_optimization.jl)
# =============================================================================
#
# WHAT
#   ring_metrics(lat)   xi1, xi2 (x,y), D2(0), max|D2|, rms D2, alpha0/1,
#                       eta0/1, fractional tunes - from the order-3 map.
#   response(lat, fam)  d(metric)/dk2 per sextupole family (central difference;
#                       xi1, alpha1, eta1, D2 are linear in k2).
#   solve_A(base,cols)  EXAMPLE target A: xi1_x = xi1_y = 0 and alpha1 = 0.
#   solve_B(base,cols)  EXAMPLE target B: xi1_x = xi1_y = 0, rms D2 minimal.
#   Targets A and B are illustrative choices, NOT physics requirements.
#
# WHY separate
#   Optimisation is a choice of the user, not part of the model.  Nothing in
#   the model modules depends on this file.
# =============================================================================
module Optimizers

using LinearAlgebra, TaylorSeries
using ..Constants, ..Elements, ..Jets, ..TaylorMap, ..MapAnalysis, ..SecondOrder

export ring_metrics, response, solve_A, solve_B
function ring_metrics(lat::Lattice)
    M = ring_taylor_map(lat)
    C = circumference(lat)
    d = dispersion_orders(M)
    p = path_expansion(M, lat.ref, C)
    ch = chromaticity_orders(M)
    al = dispersion_along(lat, d.D1, d.D2)
    return (Qx = ch.x[1], Qy = ch.y[1], xi1x = ch.x[2], xi1y = ch.y[2], xi2x = ch.x[3], xi2y = ch.y[3],
            D1_0 = d.D1[1], D2_0 = d.D2[1], D2p_0 = d.D2[2], D2max = maximum(abs, al.D2),
            D2rms = sqrt(sum(abs2, al.D2) / length(al.D2)), alpha0 = p.alpha0, alpha1 = p.alpha1,
            eta0 = p.eta0, eta1 = p.eta1, along = al, M = M)
end

const _LINEAR_KEYS = (:xi1x, :xi1y, :alpha1, :eta1, :D2_0)

"Response of metrics to k2 of each family; also returns D2(s) response vectors."
function response(lat, fam_idx; dk = 1.0)
    base = ring_metrics(set_family_strengths(lat, fam_idx, Dict{String,Float64}()))
    cols = Dict{String,Any}()
    for f in SEXT_FAMILIES
        isempty(fam_idx[f]) && continue
        mp = ring_metrics(set_family_strengths(lat, fam_idx, Dict(f => dk)))
        mm = ring_metrics(set_family_strengths(lat, fam_idx, Dict(f => -dk)))
        cols[f] = (Dict(k => (getfield(mp, k) - getfield(mm, k)) / (2dk) for k in _LINEAR_KEYS),
                   (mp.along.D2 .- mm.along.D2) ./ (2dk))
    end
    return base, cols
end

"A: xi1x = xi1y = alpha1 = 0."
function solve_A(base, cols)
    fs = [f for f in SEXT_FAMILIES if haskey(cols, f)]
    A = [cols[f][1][k] for k in (:xi1x, :xi1y, :alpha1), f in fs]
    b = -[base.xi1x, base.xi1y, base.alpha1]
    return Dict(zip(fs, A \ b))
end

"B: xi1x = xi1y = 0, minimise sum D2(s)^2 (KKT system)."
function solve_B(base, cols)
    fs = [f for f in SEXT_FAMILIES if haskey(cols, f)]
    G = hcat([cols[f][2] for f in fs]...)            # D2(s) = d0 + G k
    d0 = base.along.D2
    Ce = [cols[f][1][k] for k in (:xi1x, :xi1y), f in fs]
    ce = -[base.xi1x, base.xi1y]
    n = length(fs)
    K = [2G'G Ce'; Ce zeros(2, 2)]
    rhs = [-2G'd0; ce]
    sol = K \ rhs
    return Dict(zip(fs, sol[1:n]))
end

end # module
