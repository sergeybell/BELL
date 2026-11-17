# Hamiltonians, maps and Lie algebra

```@meta
CurrentModule = Bell
```

Modules: `LieAlgebra` (`src/lie_algebra.jl`), `Hamiltonians`
(`src/hamiltonians.jl`), `DifferentialAlgebra` (`src/differential_algebra.jl`),
`Symplectic` (`src/symplectic.jl`).

## Polynomial space

All polynomials live in one `TaylorSeries` jet space with six variables and
order [`JET_ORDER`](@ref) = 4. [`Z`](@ref) holds the coordinate polynomials, so
`collect(Z)` is the identity map. Maps are truncated at [`MAP_ORDER`](@ref) = 3
by [`truncate_to`](@ref).

Changing the orders is the only global switch of the nonlinear accuracy. A degree-4
Hamiltonian gives maps exact to degree 3, and a Dragt-Finn factor then has
generators f3 and f4.

## Element Hamiltonians

[`hamiltonian`](@ref) is the single source of element physics. The exact
expression is evaluated in `TaylorN` arithmetic, so the square root
``p_s`` and the energy term [`energy_term`](@ref) are expanded by the
polynomial engine, not by hand. A check that H has no linear part
guarantees that the reference orbit solves the equations of motion.

## Poisson bracket and Lie transformations

With ``[f, g] = \sum_k \partial_{q_k} f\,\partial_{p_k} g - \partial_{p_k} f\,\partial_{q_k} g``
and the Lie operator ``{:}f{:}\,g = [f, g]``, the flow of an s-independent
Hamiltonian ``H`` over length ``L`` is

```math
\mathcal M = e^{-L{:}H{:}}\, z = z - L[H, z] + \tfrac{L^2}{2}[H,[H,z]] - \dots
```

[`lie_exp`](@ref) sums this series term by term until it converges in relative
terms. Because brackets lower the degree by 2 and the space is truncated, the
series for nonlinear generators terminates exactly.

## Dragt–Finn factorisation

A symplectic map with a fixed point at the origin is written

```math
\mathcal M = R \circ e^{{:}f_3{:}} \circ e^{{:}f_4{:}}
```

with R the 6×6 linear part and ``f_3, f_4`` homogeneous generators.
[`dragt_finn_factor`](@ref) extracts them from a polynomial map. Each
order is inverted with the Euler-homogeneity identity
``f = \frac{1}{m}\sum_{k=1}^{3} (q_k\,g_{p_k} - p_k\, g_{q_k})`` for a homogeneous
generator of degree ``m``, where ``g`` is the corresponding homogeneous part of the map.

**Algorithm of [`dragt_finn_factor`](@ref)** (input: polynomial map ``M`` of degree 3,
``M(0)=0``):

1. ``R = \partial M/\partial z|_0``; ``N = R^{-1}M = z + N_2 + N_3``.
2. ``f_3 = \tfrac13\sum_k (q_k N_{2,p_k} - p_k N_{2,q_k})`` (the degree-2 part of ``N`` is the
   vector field of ``f_3``).
3. ``P = e^{-{:}f_3{:}}\circ N`` (truncated substitution), then
   ``f_4 = \tfrac14\sum_k (q_k P_{3,p_k} - p_k P_{3,q_k})``.

Applied to a point the factor acts as ``z \mapsto R\,\big(e^{:f_3:}(e^{:f_4:}z)\big)``:
first ``f_4``, then ``f_3``, then ``R``. Example (drift, see [Elements](elements.md)):
``f_3 = -LH_3``, ``f_4 = -LH_4`` exactly.

**Products.** [`concatenate`](@ref) multiplies two factors with the
Baker-Campbell-Hausdorff formula truncated at degree 4. The second generators
are pulled through the first linear part by [`pull_generator`](@ref):

```math
R = R_b R_a,\quad f_3 = f_{3a} + f_{3b}\!\circ\! R_a,\quad
f_4 = f_{4a} + f_{4b}\!\circ\! R_a - \tfrac12\,[f_{3b}\!\circ\! R_a,\; f_{3a}].
```

The ring is the product of all element factors
([`product_factors`](@ref), [`ring_lie_factor`](@ref)). Unlike repeated
truncated substitution ([`compose_maps`](@ref)), each factor stays an exactly
symplectic map, and only the generators are truncated.

**Symplecticity.** Every Lie factor ``e^{{:}f{:}}`` is symplectic for any f, and
R = exp(L·S·A) is symplectic because S·A is Hamiltonian. The truncated
polynomial `factor_to_map(F)` is symplectic only up to degree MAP_ORDER.
That is why tracking applies the factor by symplectic flows instead (see
[Tracking](tracking.md)).

## Differential-algebra (RK4) maps

An independent route to the same maps integrates ``dz/ds = S\nabla H`` with
classical RK4, where z is a vector of polynomials starting from the identity.
The integration is 4th order in the step and not symplectic. It is used only to
cross-check the Lie product and as the input of the map-based normal form.

## Symplecticity checks
