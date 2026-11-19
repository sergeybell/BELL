# Lie factors and differential algebra

Differential algebra and a Lie factor are two ways to carry the same ``H``. Neither is a second physical model.

## Truncated polynomial arithmetic

A differential-algebra number is a polynomial of total degree at most ``N``, with arithmetic that drops every product beyond ``N``. A map is six such polynomials, one per coordinate. The algebra knows nothing about magnets, cavities, or ``S``. It is the arithmetic used to integrate ``\mathrm{d}z/\mathrm{d}s=S\nabla H`` on the identity jet and to read coefficients afterwards.

A truncated map is symplectic only through order ``N``. The identity ``J^{\mathsf T}SJ=S`` involves products of derivatives. Those products generate degrees that the truncation has already thrown away, so the defect starts at the first omitted degree. Raising ``N`` pushes the defect out. It does not cancel it. A Taylor map has the same limitation: it is a truncated series, and it is symplectic only up to that cut.

## Lie representation

The Dragt–Finn form

```math
\mathcal{M}=R\circ e^{:f_3:}\circ e^{:f_4:}\circ\cdots
```

stores the linear map and a list of homogeneous generators. The current factors stop at ``f_4``, which matches a map of degree 3 and a Hamiltonian of degree 4. Applying ``e^{:f:}`` means integrating the Hamiltonian vector field of ``f``, not substituting a truncated polynomial. That flow is canonical for any generator, at any amplitude, up to the round-off of the stepper. Truncating the generator limits which aberrations are present. It does not make the flow of the kept generator non-symplectic.

So the same words "order ``N``" mean two different cuts. On a polynomial map they cut the symplecticity. On a Lie factor they cut the generator, and the application stays symplectic.

## Extracting generators

Dragt–Finn factorisation runs in the other direction. Any origin-fixing polynomial map, whether it came from a Lie series or from integrating a vector field in the truncated algebra, can be written as ``R`` followed by ``e^{:f_3:}``, ``e^{:f_4:}``, and so on up to the degree the polynomial still contains. The normal form in [Normal form](normal_form.md) is then a calculation on those generators. One does not conjugate the full polynomial map to build it.

## Two writings of the vector field

The production vector field is ``S\nabla H``, with ``H`` from [Element Hamiltonians](hamiltonians.md). The partial derivatives are taken either by a dual number (a first-order jet at the current point) or, on a differential-algebra state, by differentiating the polynomial once and composing it with the current jet. On a jet of order ``N`` the dual partial and the polynomial derivative agree through order ``N-1``, because differentiation itself drops one degree.

The field formulas are a second writing of the same hard-edge magnet: ``\mathrm{d}z/\mathrm{d}s`` expressed from ``h``, ``k_1`` and ``k_2`` rather than by differentiating ``H``. They are a check that the scalar ``H`` and the expanded multipole force describe one magnet. They are not the map the ring is built from.

## What a side-by-side comparison shows

Take three degree-3 maps of one element, or of one lattice, and factor each of them:

* **A.** The Lie series of the polynomial ``H``, expanded to a map and factored.
* **B.** The differential-algebra map of ``\mathrm{d}z/\mathrm{d}s=S\nabla H``, integrated with a fixed-step Runge–Kutta method, then factored.
* **C.** The differential-algebra map of the explicit field formulas, integrated with the same method and the same step count, then factored.

**B** and **C** agree. They are two writings of one vector field, stepped the same way. The residual is the difference between differentiating ``H`` and evaluating the expanded force, at the working order.

**A** stands apart from both, by the truncation of the integrator. The Lie series is not that Runge–Kutta stepper. The coefficient gap between **A** and **B** is a step-size error of the numerical flow, not a second Hamiltonian and not a failure of the factorisation. The same pattern appears on a single sextupole, on a FODO cell, and on the FODO ring: **B** with **C**, and **A** offset by the stepper.

## Target

One engine. The element is ``H``. A point and a differential-algebra jet are integrated from ``H``. Dragt–Finn turns any such jet into generators ``f_3,\ldots,f_{N+1}``. The Lie product, the normal form and the optical readings consume those generators. The explicit field formulas stay a check against ``S\nabla H``. A Taylor map stays a check. A thin cavity stays the exact kick of the integrated generator ``G``, which is already the flow ``e^{-:G:}`` and does not need a second model.
