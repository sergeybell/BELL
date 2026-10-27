# Assumptions

Remaining modelling assumptions (accepted design choices such as injection
energy, RF parameters, CODATA mass, and the interpretation of the lattice
file are recorded in the Programmer Manual, not here).

1. Hard-edge magnets, no fringe fields.
2. No space charge, impedances, radiation or alignment errors.
3. RF field purely longitudinal on axis (no transverse RF).
4. Magnets and RF frequency scale with the reference momentum during
   acceleration.
5. Single RF station per turn for the averaged longitudinal Hamiltonian.
6. The one-turn object is the Lie product of Dragt–Finn factors, truncated
   at degree 4 in the generators. Its polynomial expansion is symplectic
   only up to the dropped degree. Integration of the unexpanded
   Hamiltonian is symplectic to round-off.
7. The product is written for maps that fix the origin (thin RF at
   ``\varphi_s=\pi``). A net energy gain per turn is a separate problem.