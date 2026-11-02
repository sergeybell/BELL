# MAD-X on magnetic_2: TWISS+CHROM and PTC_TWISS (exact) -> results.csv
#   julia --project=examples benchmarks/twiss/madx/run.jl [madx.exe]
# MAD-X DQ, DX are per PT = dE/(p0 c); converted here to delta_p: x beta0.
using Printf
include(joinpath(@__DIR__, "..", "common.jl"))
const HERE = @__DIR__; const RAW = joinpath(HERE, "raw"); mkpath(RAW)
const EXE = length(ARGS) >= 1 ? ARGS[1] : raw"C:\Users\Professional\Desktop\madx.exe"
mx(f) = cd(HERE) do; @elapsed run(pipeline(`$EXE $f`, stdout = joinpath(RAW, replace(f, ".madx" => ".log")))); end
t_tw = minimum(mx("magnetic_2.madx") for _ in 1:5)
t_ptc = minimum(mx("ptc_exact_true.madx") for _ in 1:3)
function tfs_summ(f)
    L = readlines(joinpath(RAW, f)); h = split(L[findfirst(startswith("*"), L)])[2:end]
    v = split(L[end]); return Dict(h[i] => parse(Float64, v[i]) for i in eachindex(h) if tryparse(Float64, v[i]) !== nothing)
end
s = tfs_summ("summ.tfs"); p = tfs_summ("ptc_summ_exact_true.tfs")
E = 1.87561294257 + 0.270; b0 = sqrt(1 - (1.87561294257 / E)^2)
L = readlines(joinpath(RAW, "twiss.tfs")); h = split(L[findfirst(startswith("*"), L)])[2:end]
rows = [split(l) for l in L if !startswith(l, "@") && !startswith(l, "*") && !startswith(l, "\$")]
col(n) = [parse(Float64, r[findfirst(==(n), h)]) for r in rows]
bx, by, ax, ay, dx, dpx = col("BETX"), col("BETY"), col("ALFX"), col("ALFY"), col("DX") .* b0, col("DPX") .* b0
res = [("C_m", s["LENGTH"]), ("Qx", s["Q1"]), ("Qy", s["Q2"]), ("xi_x", s["DQ1"] * b0), ("xi_y", s["DQ2"] * b0),
       ("xi_x_ptc_exact", p["DQ1"] * b0), ("xi_y_ptc_exact", p["DQ2"] * b0), ("alpha_c", s["ALFA"]), ("alpha_c_ptc", p["ALPHA_C"]),
       ("eta_slip_ptc", p["ETA_C"]), ("betx0_m", bx[1]), ("bety0_m", by[1]), ("alfx0", ax[1]), ("alfy0", ay[1]),
       ("Dx0_m", dx[1]), ("Dpx0", dpx[1]), ("max_betx_m", maximum(bx)), ("max_bety_m", maximum(by)), ("max_Dx_m", maximum(dx)),
       ("I1", s["SYNCH_1"]), ("I2", s["SYNCH_2"]), ("I3", s["SYNCH_3"]), ("I4", s["SYNCH_4"]), ("I5", s["SYNCH_5"])]
open(joinpath(HERE, "results.csv"), "w") do io
    println(io, "quantity,value"); for (k, v) in res; @printf(io, "%s,%.12g\n", k, v); end
end
open(joinpath(HERE, "timings.csv"), "w") do io
    @printf(io, "stage,seconds\nmadx twiss+chrom (process; min of 5),%.4f\nmadx ptc_twiss exact no=3 (process; min of 3),%.4f\n", t_tw, t_ptc)
end
write_std_twiss(joinpath(RAW, "twiss_std.csv"), col("S"), bx, ax, by, ay, dx)
plot_code_twiss(joinpath(RAW, "twiss_std.csv"), joinpath(HERE, "twiss_madx.png"), "MAD-X TWISS, magnetic_2 (element exits, D per dp/p)")
print(read(joinpath(HERE, "timings.csv"), String))
