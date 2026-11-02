# Overlay + difference plots (code - Bell) from */raw/twiss_std.csv.
#   julia --project=examples benchmarks/twiss/plot_compare.jl
include(joinpath(@__DIR__, "common.jl"))
const HERE = @__DIR__
codes = ["bell", "optim", "cosy", "madx"]; lbl = Dict("bell"=>"Bell", "optim"=>"OptiM", "cosy"=>"COSY", "madx"=>"MAD-X")
T = Dict(c => read_std_twiss(joinpath(HERE, c, "raw", "twiss_std.csv")) for c in codes)
sty = Dict("bell"=>(:solid, 3, :black), "optim"=>(:dash, 1.5, :red), "cosy"=>(:dot, 2, :blue), "madx"=>(:dashdot, 1.5, :orange))
ps = map([(:betx, "βx [m]"), (:bety, "βy [m]"), (:D, "Dx [m]"), (:alfx, "αx")]) do (f, yl)
    p = plot(ylabel = yl, legend = :outerright)
    for c in codes
        v = getfield(T[c], f); all(isnan, v) && continue
        plot!(p, T[c].s, v, label = lbl[c], ls = sty[c][1], lw = sty[c][2], color = sty[c][3])
    end; p
end
savefig(plot(ps..., layout = (4, 1), size = (1100, 1100), xlabel = "s [m]", plot_title = "magnetic_2: Twiss overlay"), joinpath(HERE, "twiss_overlay.png"))

# differences at each code's own output points (element exits); Bell value = Bell grid point at that s
b = T["bell"]
bellat(s, f) = (i = argmin(abs.(b.s .- s)); abs(b.s[i] - s) < 1e-6 ? getfield(b, f)[i] : NaN)
spec = [(:betx, "βx/βx_Bell − 1", true), (:bety, "βy/βy_Bell − 1", true), (:alfx, "αx − αx_Bell", false),
        (:alfy, "αy − αy_Bell", false), (:D, "Dx − Dx_Bell [m]", false)]
for (fname, excl) in (("twiss_diff_vs_bell.png", String[]), ("twiss_diff_vs_bell_no_optimD.png", ["optim"]))
    qs = map(spec) do (f, yl, rel)
        p = plot(ylabel = yl, legend = :outerright)
        for c in codes[2:end]
            v = getfield(T[c], f); all(isnan, v) && continue
            (f == :D && c in excl) && continue
            d = [rel ? v[i] / bellat(T[c].s[i], f) - 1 : v[i] - bellat(T[c].s[i], f) for i in eachindex(v)]
            plot!(p, T[c].s, d, label = lbl[c], color = sty[c][3], marker = :circle, ms = 2.5, lw = 0.8)
        end; p
    end
    savefig(plot(qs..., layout = (5, 1), size = (1100, 1250), xlabel = "s [m]",
                 plot_title = "difference vs Bell at element exits" * (isempty(excl) ? "" : " (OptiM D omitted: non-periodic start)")), joinpath(HERE, fname))
end
# summary of max |diff|
open(joinpath(HERE, "twiss_diff_summary.csv"), "w") do io
    println(io, "code,max_rel_dbetx,max_rel_dbety,max_abs_dalfx,max_abs_dalfy,max_abs_dDx_m,npoints")
    for c in codes[2:end]
        t = T[c]; m(f, rel) = (v = getfield(t, f); all(isnan, v) ? NaN :
            maximum(abs, filter(!isnan, [rel ? v[i] / bellat(t.s[i], f) - 1 : v[i] - bellat(t.s[i], f) for i in eachindex(v)])))
        @printf(io, "%s,%.3e,%.3e,%.3e,%.3e,%.3e,%d\n", c, m(:betx, true), m(:bety, true), m(:alfx, false), m(:alfy, false), m(:D, false), length(t.s))
    end
end
print(read(joinpath(HERE, "twiss_diff_summary.csv"), String))
