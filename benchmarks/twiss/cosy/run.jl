# COSY Infinity 10.2 on magnetic_2 (lattice_md.fox copied from SCT_2/chromaticity/COSY)
#   julia --project=examples benchmarks/twiss/cosy/run.jl [cosy.exe]
# lattice_md.fox INCLUDEs the compiled COSY system C:/GIT_REPOS/SCT_2/COSY/src/COSY (read only).
using Printf
include(joinpath(@__DIR__, "..", "common.jl"))
const HERE = @__DIR__; const RAW = joinpath(HERE, "raw"); mkpath(RAW)
const EXE = length(ARGS) >= 1 ? ARGS[1] : raw"C:\GIT_REPOS\SCT_2\COSY\src\cosy.exe"
function cosy(f)
    cd(HERE) do
        @elapsed run(pipeline(`$EXE $f.fox`, stdout = joinpath(RAW, "$f.log"), stderr = joinpath(RAW, "$f.log")))
    end
end
t_lat = cosy("lattice_md")                              # compile lattice -> lattice_md.bin
t_tw = minimum(cosy("twiss") for _ in 1:5)              # OV 1 3 0: periodic Twiss + linear 6D map
t_ch = minimum(cosy("ov321") for _ in 1:5)              # OV 3 2 1: TP tunes(delta_K) -> chromaticity
kv(f) = Dict(strip(a) => parse(Float64, replace(strip(b), r"(?<=\d)([+-]\d)" => s"E\1"))
             for (a, b) in (split(l, "=") for l in eachline(joinpath(RAW, f)) if occursin("=", l)))
lin = kv("cosy_linear.txt"); ch = kv("cosy_ov321.txt")
g = ch["gamma0"]; F = ch["factor_1p1g"]; C = lin["length"]
D, Dp = lin["D0_dK"] * F, lin["Dp0_dK"] * F         # per delta_p
eta = -(lin["M51"] * lin["D0_dK"] + lin["M52"] * lin["Dp0_dK"] + lin["M56"]) / C * F^2
tw = [split(l) for l in eachline(joinpath(RAW, "cosy_twiss.txt")) if !startswith(l, "#")]
bx = [parse(Float64, t[3]) for t in tw]; by = [parse(Float64, t[4]) for t in tw]; dx = [parse(Float64, t[5]) * F for t in tw]
# integer tune is not given by TP; integer part taken from phase advance count = 2 (cross-checked with other codes)
res = [("C_m", C), ("Qx_frac", ch["Qx_frac"]), ("Qy_frac", ch["Qy_frac"]), ("xi_x", ch["xi_x_p"]), ("xi_y", ch["xi_y_p"]),
       ("alpha_c", eta + 1 / g^2), ("eta_slip", eta), ("betx0_m", lin["betx0"]), ("bety0_m", lin["bety0"]),
       ("alfx0", lin["alfx0"]), ("alfy0", lin["alfy0"]), ("Dx0_m", D), ("Dpx0", Dp),
       ("max_betx_m", maximum(bx)), ("max_bety_m", maximum(by)), ("max_Dx_m", maximum(dx))]
open(joinpath(HERE, "results.csv"), "w") do io
    println(io, "quantity,value"); for (k, v) in res; @printf(io, "%s,%.12g\n", k, v); end
end
open(joinpath(HERE, "timings.csv"), "w") do io
    @printf(io, "stage,seconds\ncosy lattice compile (process),%.4f\ncosy twiss OV1 (process; min of 5),%.4f\ncosy chrom OV3 TP (process; min of 5),%.4f\n", t_lat, t_tw, t_ch)
end
sv = [parse(Float64, t[2]) for t in tw]; nanv = fill(NaN, length(sv))
# twiss.fox writes beta, D only (alpha not written); s of START row = 0
write_std_twiss(joinpath(RAW, "twiss_std.csv"), sv, bx, nanv, by, nanv, dx)
plot_code_twiss(joinpath(RAW, "twiss_std.csv"), joinpath(HERE, "twiss_cosy.png"), "COSY Infinity 10.2, magnetic_2 (element exits)")
print(read(joinpath(HERE, "timings.csv"), String))
