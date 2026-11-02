# =============================================================================
# optimizers/Optimizers.jl - derivative-free fitting of a scalar objective
# =============================================================================
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# Three methods, one entry point `optimize`:
#   NelderMead   – local simplex search from a given x0; stop when the
#                  simplex spread in f is below `tol` (standard local fit).
#   RandomStart  – MultiStart: draw random x0 in a box, run NelderMead from
#                  each, keep the best (method of random initial points).
#   GridSearch   – evaluate f on a fixed list of points; optionally polish
#                  the best with NelderMead.
#
# No gradients. Objective is f(x::AbstractVector) -> Real. The package does
# not call these by default; sextupole fitting or lattice matching uses them
# explicitly.
# =============================================================================
module Optimizers

using Random

export Optimizer, OptimizationResult, NelderMead, RandomStart, GridSearch, optimize

"""
    Optimizer

Supertype of fitting methods. Dispatch of [`optimize`](@ref) is on the subtype.
"""
abstract type Optimizer end

"""
    OptimizationResult

Result of [`optimize`](@ref).

# Fields
- `x::Vector{Float64}` – best point found;
- `f::Float64` – `objective(x)`;
- `converged::Bool` – local stop criterion met (`NelderMead` / polish);
- `n_eval::Int` – number of objective calls;
- `message::String` – short status text.
"""
struct OptimizationResult
    x::Vector{Float64}
    f::Float64
    converged::Bool
    n_eval::Int
    message::String
end

"""
    NelderMead(; tol = 1e-8, max_iter = 1000, initial_step = 0.05)

Standard local derivative-free fit (Nelder–Mead simplex, 1965).

From a starting point `x0` the method keeps a simplex of `n+1` vertices in
``\\mathbb{R}^n`` and replaces the worst vertex by reflection / expansion /
contraction / shrink until
``f_\\mathrm{worst} - f_\\mathrm{best} < \\mathrm{tol}\\,(1 + |f_\\mathrm{best}|)``,
or until `max_iter` iterations. No gradient of the objective is used.

`initial_step` sets the relative size of the initial simplex
(``\\delta_i = \\mathrm{initial\\_step}\\,\\max(1, |x_{0,i}|)``).
"""
struct NelderMead <: Optimizer
    tol::Float64
    max_iter::Int
    initial_step::Float64
end
NelderMead(; tol = 1e-8, max_iter = 1000, initial_step = 0.05) =
    NelderMead(Float64(tol), Int(max_iter), Float64(initial_step))

"""
    RandomStart(n; lower, upper, local_method = NelderMead(), rng = Random.default_rng())

MultiStart / method of random initial points.

Draws `n` starting vectors uniformly in the box `[lower, upper]`, runs
[`NelderMead`](@ref) from each, and returns the best local minimum. Local
search alone depends on the unknown basin of attraction of `x0`; random
restarts raise the chance of finding a better basin when the objective is
multimodal.

# Arguments
- `n` – number of random starts;
- `lower`, `upper` – vectors of the same length (box bounds);
- `local_method` – local solver (default [`NelderMead`](@ref));
- `rng` – random number generator.
"""
struct RandomStart <: Optimizer
    n::Int
    lower::Vector{Float64}
    upper::Vector{Float64}
    local_method::NelderMead
    rng::AbstractRNG
end
function RandomStart(n::Integer; lower, upper, local_method::NelderMead = NelderMead(),
                     rng::AbstractRNG = Random.default_rng())
    lo = collect(Float64, lower); up = collect(Float64, upper)
    length(lo) == length(up) || error("RandomStart: lower and upper must have the same length")
    all(lo .<= up) || error("RandomStart: need lower .<= upper")
    RandomStart(Int(n), lo, up, local_method, rng)
end

"""
    GridSearch(points; polish = true, local_method = NelderMead())

Evaluate the objective at every vector in `points`. If `polish` is true, run
[`NelderMead`](@ref) from the best grid point. Useful when a coarse lattice of
candidates is already known (parameter scan).
"""
struct GridSearch <: Optimizer
    points::Vector{Vector{Float64}}
    polish::Bool
    local_method::NelderMead
end
GridSearch(points; polish::Bool = true, local_method::NelderMead = NelderMead()) =
    GridSearch([collect(Float64, p) for p in points], polish, local_method)

# ----- Nelder–Mead core -------------------------------------------------------

function _nelder_mead(f, x0::AbstractVector{<:Real}, method::NelderMead)
    x0 = collect(Float64, x0)
    n = length(x0)
    n >= 1 || error("NelderMead: need at least one variable")
    n_eval = 0
    function feval(x)
        n_eval += 1
        return Float64(f(x))
    end
    # initial simplex: x0 and n axis-aligned neighbours
    xs = [copy(x0) for _ in 1:(n + 1)]
    for i in 1:n
        xs[i + 1][i] += method.initial_step * max(1.0, abs(x0[i]))
    end
    fs = [feval(xs[i]) for i in 1:(n + 1)]
    order = sortperm(fs)
    xs = xs[order]; fs = fs[order]
    converged = false
    message = "max_iter"
    for _ in 1:method.max_iter
        if fs[end] - fs[1] < method.tol * (1 + abs(fs[1]))
            converged = true
            message = "ftol"
            break
        end
        # centroid of the best n vertices
        c = zeros(n)
        for i in 1:n
            c .+= xs[i]
        end
        c ./= n
        xw = xs[end]
        # reflection
        xr = c + (c - xw)
        fr = feval(xr)
        if fs[1] <= fr < fs[n]
            xs[end] = xr; fs[end] = fr
        elseif fr < fs[1]
            # expansion
            xe = c + 2(c - xw)
            fe = feval(xe)
            if fe < fr
                xs[end] = xe; fs[end] = fe
            else
                xs[end] = xr; fs[end] = fr
            end
        else
            # contraction toward the better of (xw, xr)
            if fr < fs[end]
                xc = c + 0.5(xr - c)
                fc = feval(xc)
                if fc <= fr
                    xs[end] = xc; fs[end] = fc
                else
                    # shrink toward best
                    for i in 2:(n + 1)
                        xs[i] = xs[1] + 0.5(xs[i] - xs[1])
                        fs[i] = feval(xs[i])
                    end
                end
            else
                xc = c + 0.5(xw - c)
                fc = feval(xc)
                if fc < fs[end]
                    xs[end] = xc; fs[end] = fc
                else
                    for i in 2:(n + 1)
                        xs[i] = xs[1] + 0.5(xs[i] - xs[1])
                        fs[i] = feval(xs[i])
                    end
                end
            end
        end
        order = sortperm(fs)
        xs = xs[order]; fs = fs[order]
    end
    return OptimizationResult(xs[1], fs[1], converged, n_eval, message)
end

"""
    optimize(method::NelderMead, objective, x0) -> OptimizationResult

Local Nelder–Mead fit of `objective(x)` starting at `x0`.
"""
optimize(method::NelderMead, objective, x0::AbstractVector{<:Real}) =
    _nelder_mead(objective, x0, method)

"""
    optimize(method::RandomStart, objective) -> OptimizationResult

MultiStart: `method.n` random starts in `[lower, upper]`, local Nelder–Mead
from each; return the best [`OptimizationResult`](@ref).
"""
function optimize(method::RandomStart, objective)
    dim = length(method.lower)
    best = OptimizationResult(copy(method.lower), Inf, false, 0, "no start")
    total_eval = 0
    for _ in 1:method.n
        x0 = method.lower .+ (method.upper .- method.lower) .* rand(method.rng, dim)
        r = _nelder_mead(objective, x0, method.local_method)
        total_eval += r.n_eval
        if r.f < best.f
            best = OptimizationResult(r.x, r.f, r.converged, total_eval,
                                      "random_start + " * r.message)
        else
            best = OptimizationResult(best.x, best.f, best.converged, total_eval, best.message)
        end
    end
    return best
end

"""
    optimize(method::GridSearch, objective) -> OptimizationResult

Evaluate `objective` on every grid point; optionally polish the best with
Nelder–Mead.
"""
function optimize(method::GridSearch, objective)
    isempty(method.points) && error("GridSearch: empty points")
    n_eval = 0
    best_x = copy(method.points[1])
    best_f = Inf
    for p in method.points
        fp = Float64(objective(p)); n_eval += 1
        if fp < best_f
            best_f = fp; best_x = collect(Float64, p)
        end
    end
    if method.polish
        r = _nelder_mead(objective, best_x, method.local_method)
        return OptimizationResult(r.x, r.f, r.converged, n_eval + r.n_eval,
                                  "grid + " * r.message)
    end
    return OptimizationResult(best_x, best_f, false, n_eval, "grid")
end

end # module
