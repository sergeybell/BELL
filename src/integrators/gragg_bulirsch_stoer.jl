# Gragg-Bulirsch-Stoer extrapolation of the modified midpoint rule.

@doc raw"""
    GraggBulirschStoer(k = 8)

Gragg-Bulirsch-Stoer extrapolation: modified midpoint rule with ``n_j = 2j`` substeps
(``j = 1..k``), whose error has an expansion in ``h^2`` only (Gragg), extrapolated to
``h \to 0`` by the Aitken-Neville scheme in ``h^2``:
```math
T_{j,m+1} = T_{j,m} + \frac{T_{j,m} - T_{j-1,m}}{(n_j/n_{j-m})^2 - 1}.
```
``T_{k,k}`` has order ``2k``; ``\|T_{k,k} - T_{k,k-1}\|`` is the embedded error estimate.
An explicit method of arbitrary even order (``k = 8`` gives order 16).
"""
struct GraggBulirschStoer <: ODEMethod
    k::Int
end
GraggBulirschStoer() = GraggBulirschStoer(8)
method_order(m::GraggBulirschStoer) = 2m.k
_errexp(m::GraggBulirschStoer) = 2m.k - 1

function _modified_midpoint(f, z, s, H, n)
    h = H / n
    zm = z; z1 = _axpy(z, h, f(z, s))
    for m in 1:n-1
        zn = _axpy(zm, 2h, f(z1, s + m * h)); zm = z1; z1 = zn
    end
    [(zm[i] + z1[i] + h * f(z1, s + H)[i]) / 2 for i in eachindex(z)]
end

function _step(m::GraggBulirschStoer, f, z, s, H)
    k = m.k; ns = [2j for j in 1:k]
    T = Vector{Any}(undef, k); nf = 0
    for j in 1:k
        Tj = Vector{Any}(undef, j)
        Tj[1] = _modified_midpoint(f, z, s, H, ns[j]); nf += ns[j] + 1
        for q in 1:j-1
            r = (ns[j] / ns[j-q])^2 - 1
            Tj[q+1] = [Tj[q][i] + (Tj[q][i] - T[j-1][q][i]) / r for i in eachindex(z)]
        end
        T[j] = Tj
    end
    zk = T[k][k]
    e = k > 1 ? [zk[i] - T[k][k-1][i] for i in eachindex(z)] : nothing
    return ([x for x in zk], e, nf)
end
