# optimizers

Derivative-free fitting of a scalar `f(x) -> Real`: `NelderMead`, `RandomStart`, and `GridSearch`, through one entry point `optimize`. The objective is a scalar function of a real vector.

The module depends only on `Random`. No other submodule imports it. A script calls `optimize` itself. `Bell` re-exports the names.

A new derivative-free method goes in `Optimizers.jl` as a subtype of `Optimizer`.
