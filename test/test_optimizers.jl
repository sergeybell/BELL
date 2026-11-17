using LinearAlgebra, Random

@testset "Optimizers" begin
    # Rosenbrock banana: min at (1,1), f=0
    rosen(x) = (1 - x[1])^2 + 100(x[2] - x[1]^2)^2
    @testset "NelderMead local fit" begin
        r = optimize(NelderMead(tol = 1e-12), rosen, [0.0, 0.0])
        @test r.converged
        @test r.f < 1e-10
        @test norm(r.x .- [1.0, 1.0]) < 1e-5
        vprintln("  NelderMead: f = ", r.f, " x = ", r.x, " n_eval = ", r.n_eval)
    end
    @testset "RandomStart (MultiStart)" begin
        rng = Random.MersenneTwister(1)
        r = optimize(RandomStart(8; lower = [-2.0, -2.0], upper = [2.0, 2.0],
                                 local_method = NelderMead(tol = 1e-10), rng = rng), rosen)
        @test r.f < 1e-8
        @test norm(r.x .- [1.0, 1.0]) < 1e-4
        vprintln("  RandomStart: f = ", r.f, " x = ", r.x, " n_eval = ", r.n_eval)
    end
    @testset "GridSearch" begin
        pts = [[x, y] for x in -1.0:0.5:1.5 for y in -1.0:0.5:1.5]
        r = optimize(GridSearch(pts; polish = true, local_method = NelderMead(tol = 1e-10)), rosen)
        @test r.f < 1e-8
        vprintln("  GridSearch: f = ", r.f, " x = ", r.x)
    end
end
