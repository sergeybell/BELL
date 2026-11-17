# Plotting and optimisation

```@meta
CurrentModule = Bell
```

## Plotting

`plot_twiss`, `plot_dispersion2`, `plot_chromatic`, `plot_acceleration`,
`separatrix!`, `gif_ramp` and `gif_elements` are exported by Bell. Their
methods live in the package extension `BellPlotsExt` and are loaded by
`using Plots`.

```@docs
plot_twiss
plot_dispersion2
plot_chromatic
plot_acceleration
separatrix!
gif_ramp
gif_elements
clear_caches!
```

## `Optimizers`

```@autodocs
Modules = [Bell.Optimizers]
```

## `SextupoleFamilies`

```@autodocs
Modules = [Bell.SextupoleFamilies]
```
