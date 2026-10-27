# =============================================================================
# Jets.jl - the truncated polynomial ("jet") space shared by all modules
# =============================================================================
#
# WHAT
#   One polynomial space in the 6 canonical variables (x,px,y,py,l,delta).
#   Hamiltonians are kept to total degree JET_ORDER = 4; maps obtained from
#   the Lie series are cut at MAP_ORDER = 3.  A generator is a polynomial;
#   the Poisson bracket acts on its coefficients.
#
# WHAT YOU SEE
#   Nothing printed here; poly_string() is used by the scripts to print
#   polynomials such as "0.5*px^2 + 0.5*py^2 + 0.0386*delta^2".
# =============================================================================
module Jets

using TaylorSeries

export SPACE, Z, JET_ORDER, MAP_ORDER, truncate_to, hom, maxabs, map_linear_part,
       poly_string, map_string, coeff

const JET_ORDER = 4   # phase 2: H to order 4 -> maps exact to order 3
const MAP_ORDER = 3
const SPACE = JetSpace(order = JET_ORDER, variables = [:x, :px, :y, :py, :l, :d])
"The six canonical variables as TaylorN polynomials: Z[1] = x, ..., Z[6] = delta."
const Z = variables(SPACE)

"Keep only the homogeneous parts of order 0..n of a TaylorN polynomial."
truncate_to(p::TaylorN, n::Int) = sum(TaylorN(p[k], JET_ORDER) for k in 0:n)

"Homogeneous part of total degree k."
hom(p::TaylorN, k::Int) = TaylorN(p[k], JET_ORDER)

"Largest absolute coefficient of a polynomial (used for convergence / error checks)."
maxabs(p::TaylorN) = maximum(k -> maximum(abs, p[k].coeffs; init = 0.0), 0:JET_ORDER)
maxabs(x::Real) = abs(x)

"Coefficient of the monomial z^e, e = exponent vector of length 6."
coeff(p::TaylorN, e::Vector{Int}) = getcoeff(p, e)

"Linear part (6x6 matrix R[i,j] = d M_i / d z_j at z = 0) of a polynomial map."
map_linear_part(M::AbstractVector{<:TaylorN}) =
    [constant_term(differentiate(M[i], j)) for i in 1:6, j in 1:6]

function exponents_of_order(k::Int)
    out = Vector{Vector{Int}}()
    function rec(prefix, left, pos)
        if pos == 6
            push!(out, vcat(prefix, left)); return
        end
        for a in left:-1:0
            rec(vcat(prefix, a), left - a, pos + 1)
        end
    end
    rec(Int[], k, 1)
    return out
end

"ASCII string of a polynomial; names = variable names, terms below `tol` are skipped."
function poly_string(p::TaylorN; names = ["x", "px", "y", "py", "l", "delta"],
                     digits = 6, tol = 1e-14, maxorder = JET_ORDER)
    terms = String[]
    for k in 0:maxorder, e in exponents_of_order(k)
        c = getcoeff(p, e)
        abs(c) < tol && continue
        mono = join([e[i] == 1 ? names[i] : "$(names[i])^$(e[i])" for i in 1:6 if e[i] > 0], "*")
        cs = string(round(c, sigdigits = digits))
        push!(terms, isempty(mono) ? cs : cs * "*" * mono)
    end
    isempty(terms) && return "0"
    return replace(join(terms, " + "), "+ -" => "- ")
end

"Print all components of a map, one line each."
function map_string(M; names = ["x", "px", "y", "py", "l", "delta"], kw...)
    join(["  $(rpad(names[i], 6)) = " * poly_string(M[i]; names = names, kw...) for i in 1:6], "\n")
end

end # module
