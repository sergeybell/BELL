# =============================================================================
# Plotting.jl - Twiss / dispersion plot with a magnet layout strip
# =============================================================================
#
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# WHAT
#   plot_twiss(opt, path; title)  three stacked panels sharing the s axis:
#     top    : beta_x (blue), beta_y (red)  [m]
#     middle : horizontal dispersion D_x [m]
#     bottom : layout strip - dipoles = grey boxes on the axis,
#              focusing quads = red boxes above, defocusing quads = blue
#              boxes below, sextupoles = small green boxes, drifts = line.
#   Uses Plots with the GR backend in headless mode (GKSwstype = 100).
#
# WHY
#   Visual check of periodicity/symmetry of the optics and of the element
#   positions read from the OptiM file.
#
# WHAT YOU SEE
#   A PNG file; nothing printed.
# =============================================================================
module Plotting

ENV["GKSwstype"] = "100"
using Plots
using ..Elements, ..Longitudinal

export plot_twiss, plot_dispersion2, plot_chromatic, plot_acceleration, separatrix!, gif_ramp, gif_elements

function _box(s0, s1, y0, y1)
    Shape([s0, s1, s1, s0], [y0, y0, y1, y1])
end

function plot_twiss(opt, path::AbstractString; title = "")
    d = opt.data
    gr()
    p1 = plot(d.s, d.betx, label = "beta_x", lw = 2, color = :blue, ylabel = "beta [m]",
              title = title, legend = :topright, xformatter = _ -> "")
    plot!(p1, d.s, d.bety, label = "beta_y", lw = 2, color = :red)
    p2 = plot(d.s, d.D, label = "D_x", lw = 2, color = :green, ylabel = "D_x [m]", legend = :topright,
              xformatter = _ -> "")
    hline!(p2, [0.0], color = :black, lw = 0.5, label = "")
    p3 = plot([0, opt.C], [0, 0], color = :black, lw = 1, label = "", ylims = (-1.2, 1.2),
              yticks = false, xlabel = "s [m]", framestyle = :box, grid = false)
    for (s0, s1, e) in opt.bounds
        if e isa SBend
            plot!(p3, _box(s0, s1, -0.45, 0.45), color = :gray60, lw = 0, label = "")
        elseif e isa Quadrupole
            plot!(p3, e.k1 > 0 ? _box(s0, s1, 0, 1) : _box(s0, s1, -1, 0),
                  color = e.k1 > 0 ? :red : :blue, lw = 0, label = "")
        elseif e isa Sextupole
            plot!(p3, _box(s0, s1, -0.6, 0.6), color = :green, lw = 0, label = "")
        end
    end
    p = plot(p1, p2, p3, layout = grid(3, 1, heights = [0.5, 0.32, 0.18]), size = (1300, 800),
             link = :x, left_margin = 8Plots.mm, xlims = (0, opt.C))
    savefig(p, path)
    return path
end


function plot_dispersion2(alongs, labels, sext_s, path)
    gr()
    p1 = plot(ylabel = "D_x [m]", legend = :topright, title = "first- and second-order dispersion")
    p2 = plot(ylabel = "D2 = (x|delta^2)_closed [m]", xlabel = "s [m]", legend = :topright)
    for (a, lab) in zip(alongs, labels)
        plot!(p1, a.s, a.D, label = lab, lw = 2)
        plot!(p2, a.s, a.D2, label = lab, lw = 2)
    end
    vline!(p1, sext_s, color = :green, ls = :dot, label = "sextupoles"); vline!(p2, sext_s, color = :green, ls = :dot, label = "")
    p = plot(p1, p2, layout = (2, 1), size = (1300, 800), left_margin = 8Plots.mm)
    savefig(p, path); return path
end

function plot_chromatic(fits, path)
    gr()
    ps = []
    for f in fits
        p = plot(title = f.title, xlabel = "delta", ylabel = "fractional tune", legend = :best)
        scatter!(p, f.d, f.qx_track, label = "Q_x tracking (exact H)", color = :blue)
        scatter!(p, f.d, f.qy_track, label = "Q_y tracking (exact H)", color = :red)
        dd = range(minimum(f.d), maximum(f.d), length = 100)
        plot!(p, dd, f.qx_poly.(dd), label = "Q_x map: Q0+xi1 d+xi2 d^2/2", color = :blue)
        plot!(p, dd, f.qy_poly.(dd), label = "Q_y map", color = :red)
        push!(ps, p)
    end
    p = plot(ps..., layout = (1, length(ps)), size = (1400, 550), left_margin = 8Plots.mm, bottom_margin = 8Plots.mm)
    savefig(p, path); return path
end

# ---- longitudinal phase space helpers (phase 3) -----------------------------
"Draw the analytic separatrix H(phi,delta) = H(phi_u, 0) and the fixed points."
function separatrix!(p, lp; phimin = -2pi, phimax = pi, dmax = 0.01, label = true, color = :black)
    phis = range(phimin, phimax, length = 500); ds = range(-dmax, dmax, length = 500)
    Hs, _ = separatrix(lp)
    Hg = [long_H(lp, f, d) for d in ds, f in phis]
    contour!(p, phis, ds, Hg, levels = [Hs], color = color, lw = 2, colorbar = false, label = "")
    fp = fixed_points(lp)
    scatter!(p, [fp.stable], [0.0], m = :circle, color = :green, ms = 5, label = label ? "stable fixed point" : "")
    scatter!(p, [fp.unstable], [0.0], m = :x, color = :red, ms = 7, msw = 3, label = label ? "unstable fixed point" : "")
    label && plot!(p, [NaN], [NaN], color = color, lw = 2, label = "separatrix H = H(phi_u, 0)")
    return p
end

const PHI_LABEL = "phi = -h w0 dt  [rad]   (phi > 0: early)"
const DELTA_LABEL = "delta = dp/p0"

function plot_acceleration(rec, info, path)
    gr()
    t = rec.turn
    p1 = plot(t, rec.Ekin, xlabel = "turn", ylabel = "E_kin [MeV]", label = "reference particle", lw = 2,
              title = "reference energy: $(round(rec.Ekin[1], digits=1)) -> $(round(rec.Ekin[end], digits=2)) MeV")
    vline!(p1, [info.flat_start], label = "start of flat top (phi_s = 180 deg)", ls = :dash)
    p2 = plot(xlabel = PHI_LABEL, ylabel = DELTA_LABEL, legend = :bottomleft,
              title = "end of ramp (turn $(info.n_end), $(round(rec.Ekin[info.n_end+1], digits=1)) MeV)")
    separatrix!(p2, info.lp_end; dmax = info.dmax_end, phimin = info.philim[1], phimax = info.philim[2])
    scatter!(p2, info.ens_phi, info.ens_delta, ms = 2, msw = 0, color = :blue, alpha = 0.6, label = "ensemble particles")
    plot!(p2, info.tr_phi, info.tr_delta, color = :orange, lw = 1, label = "test particle #4, last 300 turns")
    p3 = plot(t, info.Ax .* 1e3, label = "x amplitude of x-test particle", xlabel = "turn", ylabel = "amplitude [mm]",
              title = "adiabatic damping (betatron)", lw = 1)
    plot!(p3, t, info.Ay .* 1e3, label = "y amplitude of y-test particle", lw = 1)
    plot!(p3, t, info.Ax_pred .* 1e3, label = "A0 sqrt(p_inj/p0(n))", ls = :dash, lw = 2, color = :black)
    p4 = plot(t, info.qs_smooth, label = "formula, smooth (every turn)", lw = 2, xlabel = "turn", ylabel = "Q_s",
              title = "synchrotron tune", color = :black)
    plot!(p4, t, info.qs_kick, label = "formula, one kick per turn (every turn)", lw = 1, ls = :dash, color = :red)
    scatter!(p4, info.qs_turn, info.qs_meas, label = "tracking: fit over $(info.win)-turn window (centre)", ms = 3, color = :blue)
    p = plot(p1, p2, p3, p4, layout = (2, 2), size = (1500, 1000), left_margin = 10Plots.mm, bottom_margin = 8Plots.mm)
    savefig(p, path); return path
end

"GIF (a): ensemble turn by turn with the separatrix of the current energy and phi_s."
function gif_ramp(rec, frames, lps, path; ens = 5:length(rec.u), philim = (-2pi, pi), dlim = 0.01, fps = 15)
    gr()
    anim = @animate for (k, n) in enumerate(frames)
        phi = [rec.u[j][n+1][5] for j in ens]; dl = [rec.u[j][n+1][6] for j in ens]
        bh = bucket_height(lps[k])[1]
        pa = plot(xlabel = PHI_LABEL, ylabel = DELTA_LABEL, xlims = philim, ylims = (-dlim, dlim), legend = :bottomleft,
                  title = "turn $n, E_kin = $(round(rec.Ekin[n+1], digits=2)) MeV, phi_s = $(round(rad2deg(lps[k].phi_s), digits=0)) deg")
        separatrix!(pa, lps[k]; dmax = dlim, phimin = philim[1], phimax = philim[2])
        scatter!(pa, phi, dl, ms = 2.5, msw = 0, color = :blue, label = "particles")
        pb = plot(xlabel = PHI_LABEL, ylabel = "delta / bucket half-height", xlims = philim, ylims = (-1.3, 1.3),
                  legend = false, title = "same, delta scaled by bucket height $(round(bh, sigdigits=3))")
        separatrix!(pb, LongitudinalParams(lps[k].h, lps[k].eta0, lps[k].eta1 * bh, lps[k].A / bh^2, lps[k].phi_s);
                    dmax = 1.3, phimin = philim[1], phimax = philim[2], label = false)
        scatter!(pb, phi, dl ./ bh, ms = 2.5, msw = 0, color = :blue)
        plot(pa, pb, layout = (1, 2), size = (1400, 600), left_margin = 8Plots.mm, bottom_margin = 8Plots.mm)
    end
    gif(anim, path, fps = fps)
    return path
end

"GIF (b): one turn, element by element.
Left: absolute (phi, delta) with the analytic injection separatrix; fixed axes
  phi in philim, delta in +/- dlim (caller sets dlim from the bucket height).
Right: CHANGE relative to the turn start, with FIXED physical limits:
  dphi_lim = 2 pi |eta0| * dlim_ens  (expected one-turn phase slip at the
              largest |delta| in the ensemble, from dphi = -2 pi h eta0 delta
              with h absorbed into the definition used for eta*phi convention;
              here we use the one-turn synchrotron map scale 2 pi |eta0| * d_max
              times the harmonic already in phi = -h w0 dt, so
              |Delta phi|_expected ~ 2 pi h |eta0| |delta|  -- wait: phi already
              contains h, and dphi/dn = -2 pi h eta0 delta, so for a fraction
              f = s/C of the turn the expected slip is f * that.  Full-turn
              scale: dphi_lim = 2 pi h |eta0| * d_ens with a factor 1.5 margin.
  ddelta_lim = 1.5 * A   (A = Z V T / (beta^2 E) = per-turn RF kick amplitude
              in delta; only the cavity changes delta within a turn).
Both limits are therefore set by lattice+RF physics, not by the particle cloud."
function gif_elements(snaps, names, svals, lp, path; ens = 5:length(snaps[1]),
                      philim = (-2pi, pi), dlim = 0.01, fps = 8)
    gr()
    phi0 = [snaps[1][j][5] for j in ens]; d0 = [snaps[1][j][6] for j in ens]
    d_ens = maximum(abs, d0)
    # Full-turn phase-slip scale: |dphi| = | -2 pi h eta0 delta |; lp stores h, eta0.
    dphi_lim = 1.5 * abs(2pi * lp.h * lp.eta0) * max(d_ens, 1e-6)
    # Only the cavity changes delta; kick amplitude A is in LongitudinalParams.
    ddelta_lim = 1.5 * abs(lp.A)
    anim = @animate for k in eachindex(snaps)
        phi = [snaps[k][j][5] for j in ens]; dl = [snaps[k][j][6] for j in ens]
        frac = svals[k] / max(svals[end], eps())
        pa = plot(xlabel = PHI_LABEL, ylabel = DELTA_LABEL, xlims = philim, ylims = (-dlim, dlim),
                  legend = :bottomleft,
                  title = "after element $(k-1): $(names[k])  (s = $(round(svals[k], digits=3)) m,  $(round(100frac, digits=1))% of turn)")
        separatrix!(pa, lp; dmax = dlim, phimin = philim[1], phimax = philim[2])
        scatter!(pa, phi, dl, ms = 2.5, msw = 0, color = :blue, label = "particles")
        pb = scatter(phi .- phi0, dl .- d0, ms = 2.5, msw = 0, color = :blue, legend = :topright,
                     xlims = (-dphi_lim, dphi_lim), ylims = (-ddelta_lim, ddelta_lim),
                     xlabel = "Delta phi since turn start [rad]",
                     ylabel = "Delta delta since turn start",
                     title = "intra-turn change (fixed: +/-1.5 x full-turn slip, +/-1.5 A)")
        # reference scales as dashed lines at +/- expected full-turn values
        plot!(pb, [NaN], [NaN], color = :gray, ls = :dash, label = "axes = 1.5|2pi h eta0 d_ens|, 1.5|A|")
        plot(pa, pb, layout = (1, 2), size = (1400, 600), left_margin = 8Plots.mm, bottom_margin = 8Plots.mm)
    end
    gif(anim, path, fps = fps)
    return path
end

end # module