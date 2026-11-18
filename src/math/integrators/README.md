# integrators

Generic ODE methods and the driver `integrate`: RK4, Dormand–Prince 5(4), DOP853, Gragg–Bulirsch–Stoer, implicit midpoint, Yoshida compositions, and splitting methods. One file per method. The same methods accept `Float64` and DA.

The module depends on `LinearAlgebra` and `DiffAlgebra`. `FieldTrack`, Lie `Tracking`, `DATracking`, `Fringe`, and `Trackers` call `integrate`. See [src/README.md](../../README.md).

A new method is a new file in this folder, included from `Integrators.jl`. The element model stays in `src/physics/`.
