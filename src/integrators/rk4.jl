# Classical 4-stage Runge-Kutta.

@doc raw"""
    RK4()

Classical 4-stage Runge-Kutta method, order 4. Fixed step only.
```math
z_{n+1} = z_n + \tfrac{h}{6}(k_1 + 2k_2 + 2k_3 + k_4).
```
"""
struct RK4 <: ODEMethod end
method_order(::RK4) = 4

function _step(::RK4, f, z, s, h)
    k1 = f(z, s); k2 = f(_axpy(z, h / 2, k1), s + h / 2)
    k3 = f(_axpy(z, h / 2, k2), s + h / 2); k4 = f(_axpy(z, h, k3), s + h)
    return ([z[i] + h / 6 * (k1[i] + 2k2[i] + 2k3[i] + k4[i]) for i in eachindex(z)], nothing, 4)
end
