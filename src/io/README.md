# io

`optim_converter.jl` (`OptimConverter`) reads a `.opt` lattice file and writes a Bell structure file. It does not build a `Lattice`; `load_structure` does that.

It depends on `Coordinates`, `Elements`, and `Structure`. The `Bell` facade re-exports it. No other submodule imports it. Edges are listed in [src/README.md](../README.md).

Another lattice-file reader goes in this folder as its own module, included from `Bell.jl` on the `io` layer.
