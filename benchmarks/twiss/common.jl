# Shared helpers: standard Twiss CSV + per-code Twiss PNG (Plots.jl from the examples env).
using Plots, Printf
const STD_HEADER = "s_m,betx_m,alfx,bety_m,alfy,Dx_m"
function write_std_twiss(path, s, bx, ax, by, ay, D)
    open(path, "w") do io
        println(io, STD_HEADER)
        for i in eachindex(s)
            @printf(io, "%.9g,%.12g,%.12g,%.12g,%.12g,%.12g\n", s[i], bx[i], ax[i], by[i], ay[i], D[i])
        end
    end
end
function read_std_twiss(path)
    L = readlines(path)[2:end]; A = [parse.(Float64, split(l, ",")) for l in L]
    return (s = getindex.(A, 1), betx = getindex.(A, 2), alfx = getindex.(A, 3), bety = getindex.(A, 4), alfy = getindex.(A, 5), D = getindex.(A, 6))
end
function plot_code_twiss(csv, png, title)
    t = read_std_twiss(csv)
    p1 = plot(t.s, t.betx, label = "βx", ylabel = "β [m]", title = title, lw = 1.5)
    plot!(p1, t.s, t.bety, label = "βy", lw = 1.5)
    p2 = plot(t.s, t.D, label = "Dx", ylabel = "D [m]", color = :green, lw = 1.5)
    p3 = plot(t.s, t.alfx, label = "αx", ylabel = "α", xlabel = "s [m]", lw = 1.2)
    all(isnan, t.alfy) || plot!(p3, t.s, t.alfy, label = "αy", lw = 1.2)
    savefig(plot(p1, p2, p3, layout = (3, 1), size = (1000, 800), legend = :outerright, left_margin = 5Plots.mm), png)
    return png
end
