# OptiM CLI on magnetic_2 -> results.csv, timings.csv
#   julia --project=examples benchmarks/twiss/optim/run.jl ["C:\Program Files\OptiMX-x86-64\optimx64.exe"]
# OptiM CLI gives Twiss only (-f / -s). Chromaticity, alpha_c, I1..I5 are only in GUI
# (View | Integrals) - see README. Nothing here is copied from the GUI.
using Printf
include(joinpath(@__DIR__, "..", "common.jl"))
const HERE = @__DIR__; const RAW = joinpath(HERE, "raw"); mkpath(RAW)
const EXE = length(ARGS) >= 1 ? ARGS[1] : raw"C:\Program Files\OptiMX-x86-64\optimx64.exe"
const OPT = normpath(joinpath(HERE, "..", "..", "..", "examples", "lattices", "magnetic_2.opt"))
cp(OPT, joinpath(RAW, "magnetic_2.opt"); force = true)
runo(m, out) = @elapsed run(pipeline(`$EXE -$m $(joinpath(RAW,"magnetic_2.opt")) $(joinpath(RAW,out))`, stdout = devnull))
runo("s", "twiss_s.txt")
tf = [runo("f", "twiss_f.txt") for _ in 1:5]; ts = [runo("s", "twiss_s.txt") for _ in 1:5]
rows = []
t_parse = @elapsed for ln in eachline(joinpath(RAW, "twiss_s.txt"))
    t = split(ln); (isempty(t) || startswith(t[1], "#")) && continue
    v(i) = parse(Float64, t[i])
    push!(rows, (name=t[2], s=v(3), bx=v(4), ax=v(5), Qx=v(6), D=v(7), Dp=v(8), by=v(9), ay=v(10), Qy=v(11)))
end
open(joinpath(RAW, "twiss.csv"), "w") do io
    println(io, "name,s_m,betx_m,alfx,bety_m,alfy,Dx_m,Dpx,Qx,Qy")
    for r in rows; @printf(io, "%s,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g\n", r.name, r.s, r.bx, r.ax, r.by, r.ay, r.D, r.Dp, r.Qx, r.Qy); end
end
r0, re = rows[1], rows[end]
res = [("C_m", re.s), ("Qx", re.Qx), ("Qy", re.Qy), ("betx0_m", r0.bx), ("bety0_m", r0.by), ("alfx0", r0.ax),
       ("alfy0", r0.ay), ("Dx0_m", r0.D), ("Dpx0", r0.Dp), ("Dx_end_m", re.D),
       ("max_betx_m", maximum(r.bx for r in rows)), ("max_bety_m", maximum(r.by for r in rows)), ("max_Dx_m", maximum(r.D for r in rows))]
open(joinpath(HERE, "results.csv"), "w") do io
    println(io, "quantity,value"); for (k, v) in res; @printf(io, "%s,%.12g\n", k, v); end
end
open(joinpath(HERE, "timings.csv"), "w") do io
    println(io, "stage,seconds")
    @printf(io, "optimx64 -f process (min of 5),%.4f\noptimx64 -s process (min of 5),%.4f\nparse table,%.4f\n", minimum(tf), minimum(ts), t_parse)
end
write_std_twiss(joinpath(RAW, "twiss_std.csv"), [r.s for r in rows], [r.bx for r in rows], [r.ax for r in rows], [r.by for r in rows], [r.ay for r in rows], [r.D for r in rows])
plot_code_twiss(joinpath(RAW, "twiss_std.csv"), joinpath(HERE, "twiss_optim.png"), "OptiM CLI (-s), magnetic_2 (element exits; D from Initial, not periodic)")
print(read(joinpath(HERE, "timings.csv"), String))
