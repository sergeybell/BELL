# =============================================================================
# sextupole_families.jl - optional placement and fitting of sextupole families
# =============================================================================
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# Not part of the lattice file or of the default model. Export names are kept
# so existing callers still resolve.
# =============================================================================
module SextupoleFamilies

using LinearAlgebra, TaylorSeries
using ..Coordinates, ..Elements, ..Structure, ..LieAlgebra, ..MapReadings

export SEXT_FAMILIES, place_sextupole_families, set_family_strengths,
       ring_metrics, response, solve_A, solve_B

const SEXT_FAMILIES = ["SXD", "SXB", "SYB"]

"""
    place_sextupole_families(lattice, opt; kw...) -> (lattice, fam_idx, rows, notes)

Insert a zero-strength sextupole of family SXD / SXB / SYB into the nearest
short drift slot of every D_x / β_x / β_y maximum ([`find_optics_maxima`](@ref)).
"""
function place_sextupole_families(lattice::Lattice, opt; kw...)
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
    for (i, e) in enumerate(lattice.elements)
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
    return Lattice(lattice.name * "+sext", elems, lattice.ref, lattice.harmonic), fam_idx, rows, notes
end

"""
    set_family_strengths(lattice, fam_idx, k2::Dict) -> Lattice
"""
function set_family_strengths(lattice::Lattice, fam_idx, k2::Dict)
    elems = copy(lattice.elements)
    for (f, idx) in fam_idx, i in idx
        elems[i] = Sextupole(f, elems[i].L, get(k2, f, 0.0))
    end
    return Lattice(lattice.name, elems, lattice.ref, lattice.harmonic)
end

function ring_metrics(lattice::Lattice)
    obs = observables(lattice)
    M = obs.M
    Lpath = path_length(lattice)
    d = dispersion_orders(M)
    p = path_expansion(M, lattice.ref, Lpath)
    ch = chromaticity_orders(M)
    al = dispersion_along(lattice, d.D1, d.D2)
    return (Qx = ch.x[1], Qy = ch.y[1], xi1x = ch.x[2], xi1y = ch.y[2], xi2x = ch.x[3], xi2y = ch.y[3],
            D1_0 = d.D1[1], D2_0 = d.D2[1], D2p_0 = d.D2[2], D2max = maximum(abs, al.D2),
            D2rms = sqrt(sum(abs2, al.D2) / length(al.D2)), alpha0 = p.alpha0, alpha1 = p.alpha1,
            eta0 = p.eta0, eta1 = p.eta1, along = al, M = M)
end

const _LINEAR_KEYS = (:xi1x, :xi1y, :alpha1, :eta1, :D2_0)

function response(lattice, fam_idx; dk = 1.0)
    base = ring_metrics(set_family_strengths(lattice, fam_idx, Dict{String,Float64}()))
    cols = Dict{String,Any}()
    for f in SEXT_FAMILIES
        isempty(fam_idx[f]) && continue
        mp = ring_metrics(set_family_strengths(lattice, fam_idx, Dict(f => dk)))
        mm = ring_metrics(set_family_strengths(lattice, fam_idx, Dict(f => -dk)))
        cols[f] = (Dict(k => (getfield(mp, k) - getfield(mm, k)) / (2dk) for k in _LINEAR_KEYS),
                   (mp.along.D2 .- mm.along.D2) ./ (2dk))
    end
    return base, cols
end

function solve_A(base, cols)
    fs = [f for f in SEXT_FAMILIES if haskey(cols, f)]
    A = [cols[f][1][k] for k in (:xi1x, :xi1y, :alpha1), f in fs]
    b = -[base.xi1x, base.xi1y, base.alpha1]
    return Dict(zip(fs, A \ b))
end

function solve_B(base, cols)
    fs = [f for f in SEXT_FAMILIES if haskey(cols, f)]
    G = hcat([cols[f][2] for f in fs]...)
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
