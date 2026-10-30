# Fringe fields

```@meta
CurrentModule = Bell
```

Module: `Fringe` (`src/fringe.jl`). The physics is described in the Beam Physics
chapter "Fringe fields"; this page documents the implementation.

## Design

The lattice stays physical. An [`SBend`](@ref) carries `e1, e2, gap, fint`, and
no element has a fringe on/off flag. The fringe model is a separate object,
[`FringeModel`](@ref), passed when the ring is built:

```julia
model = FringeModel(dipole = SoftEdge(profile = Enge(gap = 0.05),
                                      method  = Splitting(n = 20, order = 4)),
                    quad = HardEdge(order = 2))
F   = ring_lie_factor(lat, model)       # method added by Fringe
obs = fringe_observables(lat, model)
```

**Expansion.** [`fringe_sequence`](@ref) replaces each magnet whose model is
not `NoEdge` with `EdgeIn → magnet → EdgeOut` ([`EdgeMap`](@ref) items).
`build_factors(lat, model)` turns every item into a Lie factor; edges are
separate factors. `FringeModel()` passes the list through unchanged, so it
reproduces the plain ring exactly.

## Edge maps

- **Hard edge.** Thin Lie maps ``e^{{:}-K{:}}`` with the linear kick plus FINT
  correction (order 1), or the nonlinear effective-angle term (order 2 or 3). For
  quadrupoles the standard cubic generator is used.
- **Soft edge.** The field ``b(s) = h_0F(s)`` is integrated over ``[-W, W]``
  around the EFB in the frame that follows the field. The edge pieces are:
  - EdgeIn = Body(−W) ∘ Soft ∘ Drift(−W);
  - EdgeOut = Drift(−W) ∘ Soft ∘ Body(−W).

  This layout keeps the lattice length unchanged. Profiles
  ([`profile_function`](@ref)) are centred on the equal-area EFB. Derivatives
  of F are averaged over a small window, so kinked profiles (LinearRamp,
  Tabulated) stay integrable.
- **Integration methods** for the s-dependent H ([`fringe_hamiltonian`](@ref)):
  - [`Splitting`](@ref): Lie slices, optionally Yoshida-4;
  - [`Magnus`](@ref): orders 2 and 4;
  - [`TaylorIntegrate`](@ref): RK4 with step halving, the reference.

```@docs
FringeModel
edge_model
EdgeModel
NoEdge
HardEdge
SoftEdge
FringeProfile
Enge
ENGE_DEFAULT
Tanh
LinearRamp
Tabulated
FringeMethod
Splitting
Magnus
TaylorIntegrate
EdgeMap
fringe_sequence
edge_map
profile_function
fringe_integral
fringe_hamiltonian
fringe_observables
```
