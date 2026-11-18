# ext

`BellPlotsExt.jl` is a Julia package extension. It loads only after `using Plots`, so plotting is not a dependency of the core. The methods are `plot_twiss`, `plot_dispersion2`, `plot_chromatic`, `plot_acceleration`, `separatrix!`, `gif_ramp`, and `gif_elements`.

The extension depends on `Bell` and on `Plots` (a weak dependency in `Project.toml`). A script reaches the methods with `using Plots` after `using Bell`. Stubs of the same names live in `src/core/plotting.jl` and do not import Plots.

A new plot declares its name in `src/core/plotting.jl` and implements the method in `BellPlotsExt.jl`.
