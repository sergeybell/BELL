using Aqua

# ambiguities = false: pre-existing method ambiguities in the Lie and DA stacks
# and in TaylorSeries. This task does not rewrite the package to clear them.
#
# deps_compat: the stdlibs LinearAlgebra, Printf and Random have no [compat]
# entry, and neither does the Test extra. Those bounds are left as they were.
# Plots is a weakdep with a compat bound and is checked.
Aqua.test_all(
    Bell;
    ambiguities = false,
    deps_compat = (
        ignore = [:LinearAlgebra, :Printf, :Random],
        check_extras = (ignore = [:Test],),
    ),
)

@testset "plot stubs without Plots" begin
    @test isempty(methods(plot_twiss))
    @test isempty(methods(plot_dispersion2))
    @test isempty(methods(plot_chromatic))
    @test isempty(methods(plot_acceleration))
    @test isempty(methods(separatrix!))
    @test isempty(methods(gif_ramp))
    @test isempty(methods(gif_elements))
end
