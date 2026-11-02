# Bell
`julia --project=examples benchmarks/twiss/bell/run.jl` — converts .opt with convert_optim, linear optics (ring_matrix, periodic_optics, optics_along), ring_metrics (order-3 map: xi, alpha_c), I1..I5 integrated from Twiss.
timings.csv splits: `using Bell`, first (JIT) and warm calls.

linear_routes.jl all — one-turn linear matrix by ring_matrix / ring_lie_factor(lat).R / one_turn_map_flow N=1,2, each in a fresh process (first call incl. JIT, warm, values) -> linear_routes.csv.
