# Assumptions

Accepted design choices such as an injection energy, the RF parameters, the mass constants, and the layout of the structure file belong in the programmer manual. What follows is the physics the present model leaves out, or cuts on purpose.

1. The body of a magnet is a hard-edge Hamiltonian. Fringe maps are optional and are built beside the body, not inside it. A quadrupole soft edge is not modelled. The fringe scalar is expanded to second order in ``y``.
2. No space charge, no impedance, no radiation, and no alignment error.
3. The cavity field kept in ``H`` is longitudinal and on axis. The transverse piece of ``\mathrm{TM}_{010}``, of order ``k^2 r^2``, is omitted. The voltage on the element is the effective voltage, transit time included. The intra-gap standing-wave phase is not retained.
4. During acceleration the magnets and the RF frequency scale with the reference momentum. Normalised strengths stay fixed.
5. The averaged longitudinal Hamiltonian assumes one RF station per turn and drops transverse–longitudinal coupling.
6. Lie factors carry ``f_3`` and ``f_4``. Their polynomial expansion is symplectic only through the kept degree. The flow of a generator, and the flow of the unexpanded ``H``, are symplectic up to the integrator. The explicit multipole field is a check, not a second ring.
7. The one-turn Lie product is written for maps that fix the origin. A net energy gain per turn is the moving reference of [Acceleration](acceleration.md), not another term in that product.
8. Spin is stored and not propagated ([Spin](spin.md)).
