# =============================================================================
# da/fields.jl - equations of motion from hard-edge multipole fields
# =============================================================================
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# Independent of the Lie / Hamiltonian path. Equations dz/ds are written from
# hard-edge multipole fields and relativistic kinematics in canonical
# coordinates z = (x, px, y, py, l, δ). The same f(z) evaluates on Float64 or
# on DA numbers. Lives next to DiffAlgebra under src/da/; DiffAlgebra itself
# stays pure mathematics (no elements).
# =============================================================================
module Fields

using ..Coordinates, ..Elements

export magnetic_parameters, element_field, element_rhs

"Curvature and normalised strengths (h, k1, k2) of a magnetic element."
magnetic_parameters(e::Drift) = (0.0, 0.0, 0.0)
magnetic_parameters(e::Quadrupole) = (0.0, e.k1, 0.0)
magnetic_parameters(e::Sextupole) = (0.0, 0.0, e.k2)
magnetic_parameters(e::SBend) = (curvature(e), e.k1, 0.0)
magnetic_parameters(e::AbstractElement) = (0.0, 0.0, 0.0)

@doc raw"""
    element_field(z, ref, h, k1, k2) -> dz/ds

Equations of motion for a hard-edge multipole in canonical ``z``, written from
the field strengths ``(h, k_1, k_2)`` and relativistic kinematics.

With longitudinal momentum
``p_s = \sqrt{(1+\delta)^2 - p_x^2 - p_y^2}``:

```math
\begin{aligned}
x' &= (1+hx)\,p_x/p_s,\\
p_x' &= h p_s - h - h^2 x - k_1 x - \tfrac{k_2}{2}(x^2 - y^2),\\
y' &= (1+hx)\,p_y/p_s,\\
p_y' &= k_1 y + k_2 x y,\\
l' &= v(\delta)/v_0 - (1+hx)(1+\delta)/p_s,\\
\delta' &= 0.
\end{aligned}
```

No expansion in ``p_x, p_y, \delta``. ``z`` may be numbers or DA.
"""
function element_field(z, ref::RefParticle, h, k1, k2)
    x, px, y, py, l, d = z
    ps = sqrt((1 + d)^2 - px^2 - py^2)
    hx1 = 1 + h * x
    return [
        hx1 * px / ps,
        h * ps - h - h^2 * x - k1 * x - k2 / 2 * (x^2 - y^2),
        hx1 * py / ps,
        k1 * y + k2 * x * y,
        velocity_ratio(d, ref) - hx1 * (1 + d) / ps,
        zero(x),
    ]
end

"""
    element_rhs(e, ref) -> (z, s) -> dz/ds

Right-hand side of [`element_field`](@ref) for element `e` (s-independent).
Markers and thin `RFCavity` have no continuous field; use identity / `rf_kick`.
"""
function element_rhs(e::AbstractElement, ref::RefParticle)
    (e isa Marker || e isa RFCavity) &&
        error("element_rhs: $(typeof(e)) has no continuous field")
    h, k1, k2 = magnetic_parameters(e)
    return (z, s) -> element_field(z, ref, h, k1, k2)
end

end # module
