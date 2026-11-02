# Bell: linear one-turn matrix by three internal routes, each in a FRESH julia process
# so "first call" includes its own JIT.  Usage (one method per process):
#   julia --project=examples benchmarks/twiss/bell/linear_routes.jl matrix|dragtfinn|flow1|flow2
# Driver: julia --project=examples benchmarks/twiss/bell/linear_routes.jl all
using Printf
const HERE = @__DIR__
meth = isempty(ARGS) ? "all" : ARGS[1]
if meth == "all"
    out = joinpath(HERE, "linear_routes.csv")
    open(out, "w") do io
        println(io, "route,using_s,first_s,warm_s,Qx_frac,Qy_frac,betx0,alfx0,bety0,alfy0,Dx0,Dpx0,alpha_c,max_abs_dR_vs_matrix,process_s")
    end
    for m in ("matrix", "dragtfinn", "flow1", "flow2")
        t = @elapsed run(`$(Base.julia_cmd()) --project=$(Base.active_project()) $(@__FILE__) $m`)
        open(io -> @printf(io, "%.3f\n", t), out, "a")
    end
    print(read(out, String)); exit()
end
t0 = time(); using Bell, LinearAlgebra; t_using = time() - t0
lat = load_structure(joinpath(HERE, "raw", "magnetic_2.bell")); C = circumference(lat)
R(m) = m == "matrix"    ? ring_matrix(lat) :
       m == "dragtfinn" ? ring_lie_factor(lat).R :
       m == "flow1"     ? Bell.DiffAlgebra.linear_part(one_turn_map_flow(lat, 1; cache = false)) :
                          Bell.DiffAlgebra.linear_part(one_turn_map_flow(lat, 2; cache = false))
clearcache() = (empty!(Bell.Linear._HCACHE); isdefined(Bell.LieFlow, :_FLOWCACHE) && empty!(Bell.LieFlow._FLOWCACHE))
t_first = @elapsed (Rm = R(meth))
t_warm = minimum(begin clearcache(); @elapsed R(meth) end for _ in 1:3)
po = periodic_optics(Rm, lat.ref, C)
R0 = ring_matrix(lat)
open(joinpath(HERE, "linear_routes.csv"), "a") do io
    @printf(io, "%s,%.4f,%.4f,%.5f,%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,%.3e,", meth, t_using, t_first, t_warm,
            po.x.mu / 2pi, po.y.mu / 2pi, po.x.beta, po.x.alpha, po.y.beta, po.y.alpha, po.D4[1], po.D4[2], po.alpha_c, maximum(abs, Rm - R0))
end
