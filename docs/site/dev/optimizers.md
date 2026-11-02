# Optimizers (`Bell.Optimizers`)

```@meta
CurrentModule = Bell
```

Derivative-free fitting of a scalar objective `f(x::Vector) -> Real`. The ring
code does not call these by default; a script builds `f` from optics / maps and
calls [`optimize`](@ref).

## Methods

| type | role |
|:--|:--|
| [`NelderMead`](@ref) | Standard local fit. Simplex search from a given `x0`; stop when the spread of `f` on the simplex is below `tol`. No gradients. |
| [`RandomStart`](@ref) | **MultiStart** / method of random initial points. Draw `n` starts in a box, run Nelder–Mead from each, keep the best. Needed when the landscape is multimodal and the basin of the true minimum is unknown. |
| [`GridSearch`](@ref) | Evaluate `f` on a fixed list of points; optionally polish the best with Nelder–Mead. |

```julia
rosen(x) = (1 - x[1])^2 + 100(x[2] - x[1]^2)^2
r1 = optimize(NelderMead(tol = 1e-12), rosen, [0.0, 0.0])
r2 = optimize(RandomStart(20; lower = [-2, -2], upper = [2, 2]), rosen)
r3 = optimize(GridSearch([[x, y] for x in -1:0.5:1.5 for y in -1:0.5:1.5]), rosen)
```

All return [`OptimizationResult`](@ref) with fields `x`, `f`, `converged`, `n_eval`, `message`.
