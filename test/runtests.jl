# Test entry point. Suites live in test_*.jl.
# Diagnostics: set BELL_TEST_VERBOSE=1.
#   julia --project=. -e "using Pkg; Pkg.test()"
using Test, Printf
using Bell

function vprintln(args...)
    get(ENV, "BELL_TEST_VERBOSE", "") == "1" || return nothing
    println(args...)
    return nothing
end

function vprintf(fmt::AbstractString, args...)
    get(ENV, "BELL_TEST_VERBOSE", "") == "1" || return nothing
    print(Printf.format(Printf.Format(fmt), args...))
    return nothing
end

@testset "Bell" begin
    include("test_diffalgebra.jl")
    include("test_integrators.jl")
    include("test_da_tracking.jl")
    include("test_lie.jl")
    include("test_structure.jl")
    include("test_coordinates.jl")
    include("test_optics_regression.jl")
    include("test_optimizers.jl")
    include("test_aqua.jl")
end
