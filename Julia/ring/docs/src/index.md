# LieRing Beam Physics Manual

LieRing is a modular Julia program for symplectic modelling of a storage ring
(here: the Nuclotron-based lattice `magnetic_2.opt`, deuteron, 270 MeV) built
on Hamiltonians and Lie algebra. Every map in the code is derived from the
element Hamiltonians written in this manual; no external code was used for the
physics. OptiM and COSY INFINITY results are used only as external numeric checks.

This manual describes **only what exists in the code** (status: v0.4). Each
section names the source file that implements it (see [Code map](codemap.md)).

Contents:

1. [Coordinates and reference particle](coordinates.md)
2. [Coordinate transformation user <-> canonical](transform.md)
3. [Element Hamiltonians and Lie maps](elements.md)
4. [Linear optics](linear.md)
5. [Second-order dispersion, momentum compaction and slip factor](second_order.md)
6. [Chromaticity](chromaticity.md)
7. [RF cavity of finite length](rf.md)
8. [Acceleration and adiabatic damping](acceleration.md)
9. [Longitudinal motion, separatrix and synchrotron tune](longitudinal.md)
10. [Symplecticity checks and Taylor-map comparison](symplectic.md)
11. [Sextupoles and the optional optimiser](sextupoles.md)
12. [Assumptions and hypotheses](assumptions.md)
13. [Code map](codemap.md)
14. [Programmer manual (stub)](programmer.md)

All scripts are run from `Julia/ring` as `julia --project=. scripts/NN_name.jl`;
the numbers quoted below are copied from their printed output.
