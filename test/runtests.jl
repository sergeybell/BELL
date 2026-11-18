# Test entry point. Suites live in test_*.jl.
# Diagnostics: set BELL_TEST_VERBOSE=1.
# Default Pkg.test() is the fast suite. Optics regression, high-order DA,
# Aqua, and optimizer tests run only when BELL_FULL_TESTS=1.
#   julia --project=. -e "using Pkg; Pkg.test()"
using Test, Printf
using Bell

full_tests() = get(ENV, "BELL_FULL_TESTS", "") == "1"

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
    # Before any `using TaylorSeries`: both export `inverse_map`.
    if full_tests()
        include("test_diffalgebra.jl")
    end
    include("test_integrators.jl")
    include("test_phase_space_gradient.jl")
    include("test_src_review.jl")
    include("test_da_tracking.jl")
    include("test_element_flow.jl")
    include("test_element_flow_rf.jl")
    include("test_lie.jl")
    include("test_tracking_interfaces.jl")
    include("test_physics_consistency.jl")
    include("test_structure.jl")
    include("test_coordinates.jl")
    include("test_cache_keys.jl")
    if full_tests()
        include("test_optics_regression.jl")
        include("test_optimizers.jl")
        include("test_aqua.jl")
    end
end
