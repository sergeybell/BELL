# Bell on magnetic_2: Twiss, tunes, chromaticity, alpha_c, with split timings.
#   julia --project=examples benchmarks/twiss/bell/run.jl
t_start = time()
using Printf
t0 = time(); using Bell; t_using = time() - t0
include(joinpath(@__DIR__, "..", "common.jl"))
const HERE = @__DIR__; const RAW = joinpath(HERE, "raw"); mkpath(RAW)
const OPT = normpath(joinpath(HERE, "..", "..", "..", "examples", "lattices", "magnetic_2.opt"))
bellfile = joinpath(RAW, "magnetic_2.bell")

function linear(lat)
    po = periodic_optics(ring_matrix(lat), lat.ref, circumference(lat))
    return po, optics_along(lat, po)
end
work(lat) = (linear(lat)..., ring_metrics(lat))

t0 = time(); convert_optim(OPT, bellfile); lat = load_structure(bellfile); t_load1 = time() - t0
t0 = time(); linear(lat); t_lin1 = time() - t0
t0 = time(); ring_metrics(lat); t_rm1 = time() - t0
t_lin2 = minimum(@elapsed(linear(lat)) for _ in 1:5)
t_rm2  = minimum(@elapsed(ring_metrics(lat)) for _ in 1:3)
po, opt = linear(lat); rm = ring_metrics(lat)

# radiation integrals I1..I5 from Bell Twiss (sector bends, k1 = 0 in bends)
I = zeros(5); d = opt.data
for k in 2:length(d.s)
    sm = (d.s[k] + d.s[k-1]) / 2
    e = opt.bounds[findfirst(b -> b[1] <= sm < b[2], opt.bounds)][3]
    e isa SBend || continue
    h = e.angle / e.L; ds = d.s[k] - d.s[k-1]; i = k
    D = (d.D[k] + d.D[k-1]) / 2; Dp = (d.Dp[k] + d.Dp[k-1]) / 2
    bx = (d.betx[k] + d.betx[k-1]) / 2; ax = (d.alfx[k] + d.alfx[k-1]) / 2
    H = (1 + ax^2) / bx * D^2 + 2ax * D * Dp + bx * Dp^2
    I .+= ds .* [h * D, h^2, abs(h)^3, D * h * (h^2 + 2e.k1), abs(h)^3 * H]
end

open(joinpath(RAW, "twiss.csv"), "w") do io
    println(io, "s_m,betx_m,alfx,bety_m,alfy,Dx_m,Dpx,mux_2pi,muy_2pi")
    for k in eachindex(d.s)
        @printf(io, "%.6f,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g\n", d.s[k], d.betx[k], d.alfx[k],
                d.bety[k], d.alfy[k], d.D[k], d.Dp[k], d.mux[k] / 2pi, d.muy[k] / 2pi)
    end
end
C = circumference(lat)
res = [("C_m", C), ("Qx", opt.Qx), ("Qy", opt.Qy), ("xi_x", rm.xi1x), ("xi_y", rm.xi1y),
       ("alpha_c", po.alpha_c), ("eta_slip", po.eta_slip), ("betx0_m", po.x.beta), ("bety0_m", po.y.beta),
       ("alfx0", po.x.alpha), ("alfy0", po.y.alpha), ("Dx0_m", po.D4[1]), ("Dpx0", po.D4[2]),
       ("max_betx_m", maximum(d.betx)), ("max_bety_m", maximum(d.bety)), ("max_Dx_m", maximum(d.D)),
       ("I1", I[1]), ("I2", I[2]), ("I3", I[3]), ("I4", I[4]), ("I5", I[5])]
open(joinpath(HERE, "results.csv"), "w") do io
    println(io, "quantity,value"); for (k, v) in res; @printf(io, "%s,%.12g\n", k, v); end
end
t_total_in = time() - t_start
open(joinpath(HERE, "timings.csv"), "w") do io
    println(io, "stage,seconds")
    for (k, v) in [("using_Bell", t_using), ("convert+load (first)", t_load1),
                   ("linear optics first call (JIT)", t_lin1), ("ring_metrics first call (JIT)", t_rm1),
                   ("linear optics warm", t_lin2), ("ring_metrics warm", t_rm2),
                   ("script total inside julia", t_total_in)]
        @printf(io, "%s,%.4f\n", k, v)
    end
end
write_std_twiss(joinpath(RAW, "twiss_std.csv"), d.s, d.betx, d.alfx, d.bety, d.alfy, d.D)
plot_code_twiss(joinpath(RAW, "twiss_std.csv"), joinpath(HERE, "twiss_bell.png"), "Bell (optics_along, ds = 2 cm), magnetic_2")
print(read(joinpath(HERE, "timings.csv"), String))
