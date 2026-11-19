# RF cavity and acceleration

```@meta
CurrentModule = Bell
```

## The element

```julia
RFCavity(name, voltage_MV, harmonic, phi_s)        # thin, L = 0
RFCavity(name, voltage_MV, harmonic, phi_s, L)     # thick, length L in metres
```

A cavity has the same s-Hamiltonian interface as the magnets
(`src/physics/hamiltonians.jl`), except that a thin cavity has no continuous
``H``. For ``L > 0``,

```math
H = \mathrm{energy\_term}(\delta) - p_s + \frac{G}{L},
```

```math
G(l,\delta) = \frac{ZV}{p_0 c\, k}\cos\!\Big(\varphi_s + \frac{k\, l}{\beta_0\, v(\delta)/v_0}\Big),
\qquad k = \frac{2\pi h \beta_0}{C}.
```

`G` is [`rf_potential`](@ref)`(cav, ref, C)` and `k` is [`rf_wavenumber`](@ref).
`V` is in MV, `p0·c` is in MeV, and `C` is [`path_length`](@ref)`(lattice)`.
The gap is uniform and the transit factor is 1: `V` is the effective voltage.

**Signs.** In the pair ``(q, P) = (-c\Delta t,\ E/(p_0 c))`` the force of the gap is

```math
\frac{dP}{ds} = -\frac{\partial H}{\partial q} = \frac{ZV}{p_0 c\, L}\sin(\varphi_s + kq).
```

Bell's pair ``(l, \delta)`` is the point transformation ``\delta = \delta(P)``,
``l = q\, dP/d\delta = q\beta = -v\Delta t``. The Hamiltonian is unchanged, with
``q = l / (\beta_0\, v/v_0)``. Then ``kq = 2\pi h\, l / (C\, v/v_0)`` is the
phase of [`rf_kick`](@ref) and of [`canonical_to_laboratory`](@ref). Positive
``l`` is an early arrival, and that map sends it to ``\varphi > 0``.

A storage-ring phase (``\sin\varphi_s = 0``) leaves the reference momentum
fixed. The one-turn map does not depend on the turn index, so the Lie factor,
the 6D normal form, and ``Q_s`` are well-defined. ``\sin\varphi_s \ne 0`` puts
a linear term in ``l`` into the Taylor polynomial of ``H``: the reference
itself gains energy. The code does not cancel that term inside the element.

The voltage does not depend on ``x`` or ``y``. Every particle uses the
reference species.

## Maps

| cavity | Lie path (`src/maps/lie/tracking.jl`) | point track and DA (`element_flow`) |
|:--|:--|:--|
| thin, ``L = 0`` | ``\exp(-{:}G{:})``, which is [`rf_kick`](@ref) | [`rf_kick`](@ref) on numbers or on DA. [`element_flow`](@ref) throws |
| thick, ``L > 0`` | ``\exp(-L{:}H{:})`` | [`element_flow`](@ref) of ``H`` over ``L``. Needs `path_length > 0` |

[`exact_track_element`](@ref), [`element_da_map`](@ref), and
[`taylor_map_element`](@ref) use that split. There is no separate RF branch
and no drift/kick/drift product in the production map. A drift/kick/drift
splitting remains only as a limit checked in `test/test_element_flow_rf.jl`:
as the slices shrink, the flow of ``H`` approaches that splitting. It is not
the step the library takes.

[`magnetic_lattice`](@ref) drops a thin cavity and replaces a thick cavity by
a drift of the same length. Chromaticity, momentum compaction, and Twiss of
the magnetic ring are read from that lattice.

[`build_factors`](@ref), [`track`](@ref), [`ring_lie_factor`](@ref), and
[`one_turn_map_da`](@ref) include the cavity. Pass `magnetic_lattice` when the
reading must be the magnets only.

## Acceleration

[`accelerate`](@ref) (`src/longitudinal_models/acceleration.jl`) is a thin
cavity plus a new reference after each kick. [`exact_track_element`](@ref)
moves every particle (magnets through [`element_flow`](@ref), the cavity
through [`rf_kick`](@ref)). [`renormalize`](@ref) then makes particle 1 the
reference: ``p_0' = p_0(1+\delta_s)``, transverse momenta scale by
``1/(1+\delta_s)``, and ``l`` is rebuilt from the arrival-time offset at the
velocity before that rescaling. The formulas are in
[Tracking](tracking.md). Normalised magnet strengths stay fixed.

That rebuild is outside the cavity Hamiltonian. A continuous reference
``p_0(s)`` inside the element, with

```math
\Delta H = \frac{p_0'}{p_0}\big(x p_x + y p_y + l(1+\delta)\big)
```

cancelling the linear term in ``l`` when ``dE_0/ds = ZV\sin\varphi_s / L``, is
not implemented. ``\beta_0(s)`` and ``k(s)`` would make ``H`` depend on ``s``
explicitly, and the one-turn map would depend on the turn. The ramp that
runs today is the discrete [`renormalize`](@ref) after each thin
[`RFCavity`](@ref).

The averaged ``(\varphi, \delta)`` model in `src/longitudinal_models/smooth.jl`
(`synchrotron_tune`, `LongitudinalParams`, `bucket_height`) is a comparison.
[`accelerate`](@ref) does not call it.

## References for the Hamiltonian

- A. Wolski, Lecture 4, “RF cavities”, Cockcroft Institute (2014). The mode
  TM010, the transit factor, and the thin kick. In that lecture
  ``\delta_W = E/(p_0 c) - 1/\beta_0`` is an energy variable and ``z = q``.
  Here ``1+\delta = \sqrt{(\delta_W + 1/\beta_0)^2 - 1/(\beta_0\gamma_0)^2}``
  and ``l = \beta z``.
- A. Wolski, Lecture 7, “Longitudinal dynamics”, Cockcroft Institute (2014).
- F. Tecker, “Longitudinal beam dynamics”, CAS 2020, arXiv:2011.02932.
- J. Le Duff, “Longitudinal beam dynamics in circular accelerators”, CAS,
  CERN 94-01.
