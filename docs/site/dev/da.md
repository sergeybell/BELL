# Differential algebra (`Bell.DiffAlgebra`)

`src/da/DiffAlgebra.jl` is a self-contained truncated power-series algebra
``{}_nD_v`` (Berz, *Modern Map Methods*, ch. 2): the number of variables ``v`` and
the order ``n`` are run-time parameters of a [`DASpace`](@ref Bell.DiffAlgebra.DASpace).
It is independent of the accelerator code and of `Bell.DifferentialAlgebra`
(`src/differential_algebra.jl`), which is only the TaylorN cross-check route of the Lie maps.

```julia
using Bell.DiffAlgebra
sp = DASpace(3, 8)                 # v = 3, n = 8
x, y, z = da_vars(sp)              # identity variables
a = exp(x + 0.3y) * sin(z) / (1 + x^2)
coefficient(a, [1, 0, 2])          # coefficient of x z^2
M  = [x + 0.1y^2, y - 0.2x*z, z + x^3]
Mi = inverse_map(M);  compose(M, Mi)   # identity to ~1e-16
```

## Representation
Monomials ``x^\alpha``, ``|\alpha|\le n``, are stored in graded order
(``\binom{n+v}{v}`` coefficients). Products use the additive key
``k(\alpha)=\sum_i \alpha_i (n+1)^{i-1}`` (dense table for small spaces,
`Dict` otherwise), so ``k(\alpha+\beta)=k(\alpha)+k(\beta)``; terms with
``|\alpha+\beta|>n`` are dropped (truncation).

## Intrinsics
For ``a=a_0+\bar a`` with nilpotent ``\bar a`` (``\bar a^{\,n+1}=0``),
``f(a)=\sum_{k=0}^{n} f^{(k)}(a_0)\,\bar a^k/k!`` evaluated by Horner. Implemented:
`exp log sqrt cbrt inv sin cos tan sinh cosh tanh atan abs`, `^` (integer by squaring, real by the binomial series).
`atan` uses ``\operatorname{atan}a=\operatorname{atan}a_0+\operatorname{atan}\frac{\bar a}{1+a_0 a}``.

## Calculus, composition, inversion
`deriv(a, i)` = ``\partial_i a`` (drops order by one), `antideriv(a, i)` = ``\partial_i^{-1}a``.
`compose(M, N)` substitutes the DA vector ``N`` into ``M``; evaluation `a(x)` accepts numbers or DA.
`inverse_map(M)` (``M(0)=0``, invertible ``L``) iterates
``M^{-1}=L^{-1}\circ(\mathrm{id}-\mathcal N\circ M^{-1})``, gaining one order per iteration.

## Tests
`test/runtests.jl`: series coefficients of exp, sin, log, atan, sqrt (exact to 1e-16);
identities ``\sin^2+\cos^2=1``, ``e^{\log b}=b``, ``\tan\operatorname{atan}``, ``b\cdot b^{-1}``;
``M\circ M^{-1}=M^{-1}\circ M=\mathrm{id}`` below ``10^{-14}`` (``v=3,n=10``).

## Limitations
Coefficients are `Float64` (generic `T` in the type, not exercised); no sparsity
beyond truncation; `inverse_map` requires a zero constant part.