# Normal form

```@meta
CurrentModule = Bell
```

Module: `NormalForm` (`src/normal_form.jl`). It provides two routes.

## 1. Dragt–Finn 6D normal form (primary)

[`dragt_finn_nf6d`](@ref) works on the one-turn **Lie factor** with RF
([`one_turn_factor_rf`](@ref)), or on a polynomial map. The steps are:

1. **Linear normalisation.** [`symplectic_eigen_A`](@ref) builds A from the
   eigenvectors with ``a\cdot S b = 1``, so that ``A^{-1}RA = r(\mu)``,
   ``A^TSA = S``. The planes are assigned by eigenvector weights, and coupling is
   allowed.
2. **Normal coordinates.** The generators are pulled through A:
   ``f_n = f\circ A`` ([`pull_generator`](@ref)).
3. **Order 3.** Solve the homological equation ``F_3 - F_3\circ r = -f_3`` with
   [`homological_operator`](@ref). In the complex basis its eigenvalues are
   ``1-e^{-i(j-n)\cdot\mu}``.
4. **Order 4.** Form ``f_4' = f_4 + \tfrac12[f_3, F_3]`` and solve
   ``F_4 - F_4\circ r + h_4 = -f_4'``. The resonant part
   ``h_4 = \sum_{i\le j} c_{ij}J_iJ_j`` is solved for together with F4 by a
   pseudo-inverse. h4 gives the amplitude-dependent tunes ([`nf6d_tunes`](@ref)).
5. **Transformation.** Build ``T = e^{{:}F_3{:}}e^{{:}F_4{:}}``. The equilibrium
   momentum is ``\delta_\text{eq}(J) = z^*_6 + \langle(AT(w))_6\rangle_\theta = z^*_6 + d\cdot J``.
   It is checked against a numerical torus average (`d_grid`).

```@docs
dragt_finn_nf6d
NF6D
NormalForm._nf6d_from_factor
one_turn_factor_rf
symplectic_eigen_A
homological_operator
nf6d_delta_eq
nf6d_actions
nf6d_tunes
```

## 2. Map-based RF normal form (older route)

[`normal_form_rf`](@ref) builds the polynomial one-turn map with a thin RF
kick and finds its 6D fixed point. It normalises the map with an approximate
linear A (dispersion shear, Courant-Snyder blocks and a synchrotron β), then
fits ``\delta_\text{eq}(J)`` from the phase-slip condition on several tori.

```@docs
normal_form_rf
NormalFormRF
delta_eq_nf
tunes_nf
nonlinear_normal_form
delta_eq_phase_condition
NormalForm.linear_A_6d
```

## Helpers

```@docs
one_turn_map_rf
rf_kick_taylor
fixed_point_6d
center_map_at
evaluate_map
linear_normalizing_matrix
rotation_residual
```
