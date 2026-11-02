# =============================================================================
# optimizers/Optimizers.jl - stubs for lattice fitting methods
# =============================================================================
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# Working iterations are not implemented. GridSearch is exhaustive evaluation
# of a list of points. RandomStart draws random initial points.
# =============================================================================
module Optimizers

export Optimizer, GridSearch, RandomStart

"Supertype of fitting methods."
abstract type Optimizer end

"""
    GridSearch(points)

Evaluate an objective on each vector in `points`. No iteration.
"""
struct GridSearch <: Optimizer
    points::Vector{Vector{Float64}}
end

"""
    RandomStart(n; lower, upper, rng)

Draw `n` random starting points in the box [lower, upper]. No iteration.
"""
struct RandomStart <: Optimizer
    n::Int
    lower::Vector{Float64}
    upper::Vector{Float64}
end
RandomStart(n::Integer; lower, upper) = RandomStart(Int(n), collect(Float64, lower), collect(Float64, upper))

end # module
