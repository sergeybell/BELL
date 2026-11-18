# =============================================================================
# maps/abstract.jl - thin map wrappers
# =============================================================================
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# AbstractMap is the common supertype for spin and IBS state that is only stored.
# LieFactor subtypes it (declared in LieAlgebra, which is loaded next).
# compose_maps is the Taylor composition.
# linear_part and compose below are methods of the DiffAlgebra functions.
# =============================================================================
"""
    Maps

Thin map objects: a 6×6 matrix and a vector of DA components.

[`LieFactor`](@ref) is also an `AbstractMap`. It is the Dragt–Finn object;
apply it with `apply_factor` or `F(z)`. Taylor composition
remains [`compose_maps`](@ref). DA composition remains `DiffAlgebra.compose`,
with extra methods for these wrappers.
"""
module Maps

using ..DiffAlgebra: DA, nvars
import ..DiffAlgebra: linear_part, compose

export AbstractMap, LinearMap, DAMap, linear_part, compose

"""
    AbstractMap

Supertype of the map objects a tracker can hold.

Concrete wrappers in this module are [`LinearMap`](@ref) and [`DAMap`](@ref).
[`LieFactor`](@ref) is an `AbstractMap` as well, defined with the Lie algebra
so its constructors and fields stay as they are.
"""
abstract type AbstractMap end

"""
    LinearMap(R) <: AbstractMap

Linear map given by a 6×6 real matrix. The constructor copies `R`.

`linear_part(M)` returns `M.R`. `M(z)` applies that matrix to a 6-vector.
`linear_part` is the function from `DiffAlgebra`, extended for this type.
"""
struct LinearMap <: AbstractMap
    R::Matrix{Float64}
    function LinearMap(R::AbstractMatrix{<:Real})
        size(R) == (6, 6) ||
            throw(ArgumentError("LinearMap: expected a 6x6 matrix, got $(size(R))"))
        return new(Matrix{Float64}(R))
    end
end

"""
    DAMap(components) <: AbstractMap

Map given by one [`DA`](@ref) component per variable. Every component must
belong to the same `DASpace`, and there must be `nv` of them.

`linear_part(M)` is the Jacobian at the origin of the components.
`M(z)` evaluates each component at `z`.
"""
struct DAMap{T<:DA} <: AbstractMap
    components::Vector{T}
    function DAMap(components::AbstractVector{<:DA})
        isempty(components) &&
            throw(ArgumentError("DAMap: expected one component per variable"))
        nv = nvars(first(components))
        length(components) == nv ||
            throw(ArgumentError("DAMap: expected $nv components, got $(length(components))"))
        sp = first(components).sp
        for component in components
            component.sp == sp ||
                throw(ArgumentError("DAMap: components must share one DASpace"))
        end
        comps = collect(components)
        return new{eltype(comps)}(comps)
    end
end

Base.show(io::IO, M::LinearMap) = print(io, "LinearMap(6x6)")

function Base.show(io::IO, M::DAMap)
    c = first(M.components)
    print(io, "DAMap(nv=", nvars(c), ", no=", c.sp.no, ")")
    return nothing
end

"""
    linear_part(M::LinearMap) -> Matrix{Float64}
    linear_part(M::DAMap) -> Matrix

Linear matrix of a wrapped map. For a `DAMap` this is the Jacobian at the
origin of its components.
"""
linear_part(M::LinearMap) = M.R
linear_part(M::DAMap) = linear_part(M.components)

"""
    (M::LinearMap)(z) -> Vector
    (M::DAMap)(z) -> Vector

Apply the wrapped map to a coordinate vector.
"""
function (M::LinearMap)(z::AbstractVector{<:Real})
    length(z) == 6 ||
        throw(ArgumentError("LinearMap: expected a 6-vector, got length $(length(z))"))
    return M.R * z
end

(M::DAMap)(z::AbstractVector) = [component(z) for component in M.components]

"""
    compose(A::LinearMap, B::LinearMap) -> LinearMap
    compose(A::DAMap, B::DAMap) -> DAMap

Composition ``A \\circ B``: `B` is applied first, then `A`. The same order as
`compose` on DA vectors. This does not replace [`compose_maps`](@ref).
"""
compose(A::LinearMap, B::LinearMap) = LinearMap(A.R * B.R)
compose(A::DAMap, B::DAMap) = DAMap(compose(A.components, B.components))

end # module
