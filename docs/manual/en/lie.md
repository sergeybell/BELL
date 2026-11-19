# Lie product of the ring

The optics object of the ring is a product of Lie transformations. A generator is a polynomial in the canonical coordinates. The Poisson bracket acts on its coefficients. Substituting truncated maps into one another is a check, not the definition of the one-turn map. The figures in the [home page](index.md) show the product and the homogeneous readings of its polynomial.

## Poisson bracket and the exponential

Hamilton's equations are ``\mathrm{d}z/\mathrm{d}s=[z,H]=-{:}H{:}\,z``. The operator ``{:}f{:}`` is the Hamiltonian vector field of ``f``,

```math
{:}f{:}\,g=[f,g]
=\sum_i\Bigl(\partial_{q_i}f\,\partial_{p_i}g-\partial_{p_i}f\,\partial_{q_i}g\Bigr),
```

summed over the three pairs ``(x,p_x)``, ``(y,p_y)``, ``(l,\delta)``. Over a length ``L`` with ``H`` independent of ``s`` the map is the Lie series ``\exp(-L{:}H{:})``.

## Dragt–Finn factor

The series is written as a Dragt–Finn factor

```math
\mathcal{M}=R\circ e^{{:}f_3{:}}\circ e^{{:}f_4{:}}.
```

The current factors carry the generators ``f_3`` and ``f_4``. The Hamiltonian polynomial is kept through total degree 4, so the map is complete through degree 3. Applied to a point, the factor means: first the time-1 flow of ``f_4``, then of ``f_3``, then the linear matrix ``R``.

Each nonlinear flow is a canonical transformation for any polynomial generator. Tracking integrates those flows by a symplectic stepper, so each factor is symplectic to round-off at any amplitude. The polynomial expansion of the same factor is symplectic only up to the degree that was kept. The distinction is the subject of [Lie factors and differential algebra](lie_and_da.md).

## Product and the BCH term

Two factors are multiplied in the algebra, not by substituting maps. If ``\mathcal{M}_a`` is applied first and ``\mathcal{M}_b`` second,

```math
R=R_b R_a,\qquad
f_3=f_{3a}+f_{3b}\circ R_a,\qquad
f_4=f_{4a}+f_{4b}\circ R_a-\tfrac12[f_{3b}\circ R_a,\,f_{3a}].
```

The Poisson bracket of two cubics is quartic. Nested brackets are degree 5 and are dropped. The minus in front of the bracket is the Baker–Campbell–Hausdorff term for the coordinate maps ``z\mapsto\exp(:f:)z``. It is fixed by comparing the product of two elements with the Lie series of that pair.

The one-turn magnetic factor is this product over the lattice.

## Thin RF kick

A thin cavity whose synchronous phase fixes the origin is a canonical transformation: energy changes at fixed arrival time, then ``\delta`` and ``l`` are recomputed, as in [Element Hamiltonians](hamiltonians.md). It enters the one-turn product as one more Dragt–Finn factor. The sine of the kick is expanded to the same degree as the magnetic generators. The exact sine is kept in the kick itself and in the integrator of the unexpanded Hamiltonian.

A cavity with a net energy gain per turn shifts the origin. The product above is written for maps that fix the origin. A ramp, or a finite cavity with a closed-orbit offset, is a separate problem ([Acceleration](acceleration.md)).

The linear piece ``R=\exp(:f_2:)`` is the closed form in [Linear optics](linear.md). The normal form of ``\mathcal{F}`` is in [Normal form](normal_form.md).
