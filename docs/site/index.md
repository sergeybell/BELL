# Bell

**BELL - Beam Evolution in Accelerator Lattice**

Bell is a Julia package for the beam dynamics of a ring: element Hamiltonians,
Lie-algebraic one-turn maps, linear and nonlinear optics, exact 6D tracking with
RF and acceleration, fringe-field models and 6D normal forms.

The documentation has four parts:

| part | contents |
|:--|:--|
| [Physics Reference](physics.md) | the Beam Physics Manual (PDF): derivations and assumptions |
| [User's Guide](guide/lattice.md) | how to use Bell, chapter by chapter, built on `examples/` |
| [Developer's Guide](dev/architecture.md) | internals: algorithms, matrices, pseudocode, how to extend |
| [API Reference](api/core.md) | docstrings of every module |

Quick start (from the repository root):

```
julia --project=examples examples/10_linear_optics/run.jl
```
