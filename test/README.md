# Tests

`test/` is the automatic correctness suite. From the repository root:

```text
julia --project=. -e "using Pkg; Pkg.test()"
```

That command runs the fast suite in `test/runtests.jl`. The optics regression, high-order differential algebra, Aqua, and the optimizer tests run only when `BELL_FULL_TESTS=1`.

These files pass or fail. They do not read `benchmarks/`. Physical comparisons and performance live in `benchmarks/` and are run by hand.
