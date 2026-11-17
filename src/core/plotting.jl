# Plotting names. Methods live in ext/BellPlotsExt.jl and load with Plots.

"""
    plot_twiss(opt, path; title = "") -> path

PNG with three panels sharing s: β_x, β_y [m]; D_x [m]; and a magnet layout strip.
In the strip, bends are grey, focusing quads are red above the axis, defocusing quads
are blue below, and sextupoles are green. `opt` is the result of [`optics_along`](@ref).

The method is the package extension `BellPlotsExt`, loaded by `using Plots`.
If `GKSwstype` is unset when Bell is loaded, it is set to `"100"` (headless GR).
"""
function plot_twiss end

"""
    plot_dispersion2(alongs, labels, sext_s, path) -> path

Two panels with D_x(s) and D2(s) for several [`dispersion_along`](@ref) tables, with
sextupole positions `sext_s` marked.

Loaded by `using Plots` (extension `BellPlotsExt`).
"""
function plot_dispersion2 end

"""
    plot_chromatic(fits, path) -> path

Fractional tunes against δ. Each fit is a NamedTuple `(title, d, qx_track, qy_track,
qx_poly, qy_poly)`. Tracked tunes are drawn as points and map polynomials
Q0 + ξ1δ + ξ2δ²/2 as lines.

Loaded by `using Plots` (extension `BellPlotsExt`).
"""
function plot_chromatic end

"""
    plot_acceleration(record, info, path) -> path

2×2 summary of an `accelerate` run (`record` is the named tuple of that call):

- E_kin(turn);
- end-of-ramp phase space with the separatrix;
- betatron amplitudes against the adiabatic-damping prediction A0√(p_inj/p0);
- Q_s from the formulas against tracking.

`info` carries the derived series. Loaded by `using Plots` (extension `BellPlotsExt`).
"""
function plot_acceleration end

"""
    separatrix!(p, lp::LongitudinalParams; phimin, phimax, dmax, label, color) -> p

Add to plot `p` the contour H(φ, δ) = H_sep of [`long_H`](@ref) and the stable and
unstable fixed points.

Loaded by `using Plots` (extension `BellPlotsExt`).
"""
function separatrix! end

"""
    gif_ramp(record, frames, lps, path; ens, philim, dlim, fps = 15) -> path

Animated GIF of the ensemble in (φ, δ) turn by turn, with the separatrix of the
current energy. A second panel shows δ scaled by the bucket height.

Loaded by `using Plots` (extension `BellPlotsExt`).
"""
function gif_ramp end

"""
    gif_elements(snaps, names, svals, lp, path; ens, philim, dlim, fps = 8) -> path

GIF of one turn, element by element.

The left panel is absolute `(φ, δ)` with the injection separatrix. Axes stay fixed:
`φ` in `philim`, `δ` in `±dlim`.

The right panel is the change since the start of the turn. Its limits are fixed by
the lattice and the RF, not by the particle cloud:
`Δφ` within `1.5 |2π h η0| d_ens` and `Δδ` within `1.5 |A|`,
where `d_ens` is the largest `|δ|` at the start of the turn and `A` is the RF kick
amplitude in [`LongitudinalParams`](@ref).

Loaded by `using Plots` (extension `BellPlotsExt`).
"""
function gif_elements end
