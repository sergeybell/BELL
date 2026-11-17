# Lie product of the ring

The object of the ring is a product of Lie transformations. A generator is a
polynomial in the canonical coordinates. The Poisson bracket acts on its
coefficients. Substituting truncated maps into one another is a check, not
the definition of the one-turn map. The figures in the [home page](index.md)
show the product, its three uses, and the homogeneous readings of its
polynomial.

## From a Hamiltonian to a factor

Hamilton's equations are ``\mathrm{d}z/\mathrm{d}s=[z,H]=-:H:z``. Over a
length ``L`` with ``H`` independent of ``s`` the map is the Lie series
``\exp(-L:H:)``. The series is cut at total degree 3 in ``z``. The result is
written as a Dragt–Finn factor

```math
\mathcal{M}=R\circ e^{:f_3:}\circ e^{:f_4:}.
```

Applied to a point, this means: first the time-1 flow of ``f_4``, then of
``f_3``, then the linear matrix ``R``. Each nonlinear flow is a canonical
transformation for any generator. Tracking integrates those flows by the
implicit midpoint rule, so each factor is symplectic to round-off at any
amplitude.

The polynomial expansion of the same factor is symplectic only up to the
dropped degree.

## Product of two factors

Two factors are multiplied in the algebra, not by substituting maps. If
``\mathcal{M}_a`` is applied first and ``\mathcal{M}_b`` second,

```math
R=R_b R_a,\qquad
f_3=f_{3a}+f_{3b}\circ R_a,\qquad
f_4=f_{4a}+f_{4b}\circ R_a-\tfrac12[f_{3b}\circ R_a,\,f_{3a}].
```

The Poisson bracket of two cubics is quartic. Nested brackets are degree 5
and are dropped. The minus in front of the bracket is the BCH term for the
coordinate maps ``z\mapsto\exp(:f:)z``; it is fixed by comparing the product
of two elements with the Lie series of that pair.

The one-turn magnetic factor is this product over the lattice.

## Thin RF kick

A thin cavity at ``\varphi_s=\pi`` is a canonical transformation that fixes
the origin: energy changes at fixed arrival time, then ``\delta`` and ``l``
are recomputed. It enters the one-turn product as one more Dragt–Finn
factor. The sine of the kick is expanded to the same degree as the magnetic
generators. The exact sine is kept only in the independent integrator of the
unexpanded Hamiltonian.

A cavity with a net energy gain per turn shifts the origin. The product
above is written for maps that fix the origin, so a ramp or a finite cavity
with a closed-orbit offset is a separate problem.

## Normal form

The linear part ``R`` of the one-turn factor is put into three rotations by
a symplectic matrix ``A``, ``A^\mathsf{T}SA=S``. The generators in the
normal coordinates ``w=A^{-1}z`` are the pull-backs ``f(Aw)``.

A rotation acts on a complex monomial
``h_{\mathbf{j}\mathbf{n}}=\prod_k(q_k+ip_k)^{j_k}(q_k-ip_k)^{n_k}`` by the
phase ``e^{-i(\mathbf{j}-\mathbf{n})\cdot\boldsymbol{\mu}}``. The homological
operator on a generator is therefore the divisor

```math
F-F\circ r(\mu),\qquad
F_{\mathbf{j}\mathbf{n}}=\frac{-f_{\mathbf{j}\mathbf{n}}}{1-e^{-i(\mathbf{j}-\mathbf{n})\cdot\boldsymbol{\mu}}}.
```

At degree 3, ``\mathbf{j}=\mathbf{n}`` is impossible, so in the absence of
a third-order resonance the whole cubic generator is removed. At degree 4
the kernel is the span of products of actions; those terms stay as the
normal-form Hamiltonian ``h_4=\sum_{i\le j}c_{ij}J_i J_j``. The cubic
conjugation contributes ``\tfrac12[f_3,F_3]`` to the quartic generator
before this step.

The homological operator is applied to the generators. Conjugation of a
full polynomial map is not used to build the normal form.

## The torus

After the normal form there are three amplitudes ``J_x,J_y,J_s`` and three
phases ``\theta_x,\theta_y,\theta_s``. The torus of those actions is the set
of all triples of phases at those fixed amplitudes: three circles, not a
trajectory and not a list of turns.

One turn adds constant angles to the phases and does not change the
amplitudes. Laboratory ``\delta`` is a function of the three phases. The
equilibrium momentum is the average of that function over the whole torus
at once. Pieces odd in the phases average to zero. What remains is linear
in the actions,

```math
\delta_{\mathrm{eq}}=\delta^\ast+d_x J_x+d_y J_y+d_s J_s.
```

For thin RF on the stationary separatrix, ``\delta^\ast`` is a numerical
zero. The cubic-RF correction in ``\hat\varphi^2`` is degree ``J^2`` and is
not kept.

A long tracking average equals this torus average when the tunes are
incommensurate. Tracking does not build the torus.

## Independent checks

Two constructions are not the one-turn Lie product.

1. Hamilton's equations integrated in the truncated polynomial algebra
   (a differential-algebra map) and composed by substitution. The
   coefficient-wise difference with the expanded Lie product is a
   truncation and step-size error. It is printed; it is not substituted
   into ``\delta_{\mathrm{eq}}``.
2. The unexpanded Hamiltonian (the square root is not a polynomial)
   integrated by a symplectic step, with an exact thin energy kick at the
   end of the turn. This is the amplitude-exact tracking standard.
