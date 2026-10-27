# =============================================================================
# LieRing.jl - top-level module of the ring-modelling program
# =============================================================================
#
# WHAT
#   Loads all sub-modules in dependency order and re-exports their public
#   names, so a script only needs
#       using LieRing          # after Pkg.develop(path=...)
#   or  include(...); using .LieRing   (scripts in this repo)
#
#   Constants    coordinates, reference particle, canonical <-> user transform
#   Jets         polynomial space (degree 4 in H, degree 3 in maps)
#   Elements     element types + Lattice
#   Hamiltonians generator H(z) of each element
#   LieOperators Poisson bracket, exp(:f:), Dragt-Finn, product of factors
#   OptiMParser  .opt reader
#   LinearOptics exp(L S A) matrices, closed Twiss, dispersion, tunes, alpha_c
#   Symplectic   J^T S J - S checks
#   Tracking     flow of each factor, thin RF kick, ring product
#   MapAnalysis  closed orbit, eta, xi, observables from a polynomial map
#   TaylorMap    independent DA map (check)
#   NormalForm   one-turn Lie factor with thin RF; homological NF; torus average
#   Plotting     figures
#   SecondOrder  longitudinal H, BLonD (dt, dE) bridge
#   ExactIntegrator  unexpanded H, finite RF cavity
#   Optimizers   OPTIONAL sextupole response
#   Acceleration reference-energy ramp
#
# WHY
#   One entry point; every physics block stays in its own file.
#
# WHAT YOU SEE
#   Nothing printed.
# =============================================================================
module LieRing

include("Constants.jl")       # coordinates, reference particle, user <-> canonical transform
include("Jets.jl")            # truncated polynomial space
include("Elements.jl")        # ALL element types + lattice-editing helpers
include("Hamiltonians.jl")    # H of every element
include("LieOperators.jl")    # Poisson bracket, exp(:f:), Dragt-Finn factors
include("OptiMParser.jl")     # .opt reader
include("LinearOptics.jl")    # Twiss, dispersion, tunes, alpha_c
include("Symplectic.jl")      # J^T S J - S checks
include("Tracking.jl")        # factorised Lie-map tracking
include("MapAnalysis.jl")     # read D, eta, xi; observables from the Lie product
include("TaylorMap.jl")       # independent truncated Taylor map (RK4 in DA)
include("SecondOrder.jl")     # longitudinal H, BLonD drifts, separatrix
include("ExactIntegrator.jl") # exact-H integrator incl. finite RF cavity
include("Acceleration.jl")    # reference ramp, tune measurement
include("NormalForm.jl")      # 6D Lie normal form with RF
include("Plotting.jl")        # all plots
include("Optimizers.jl")      # OPTIONAL optimisation (not used by default)

using .Constants, .Jets, .Elements, .Hamiltonians, .LieOperators, .OptiMParser,
      .LinearOptics, .Symplectic, .Tracking, .MapAnalysis, .TaylorMap, .NormalForm,
      .SecondOrder, .ExactIntegrator, .Acceleration, .Plotting, .Optimizers

for mod in (Constants, Jets, Elements, Hamiltonians, LieOperators, OptiMParser,
            LinearOptics, Symplectic, Tracking, MapAnalysis, TaylorMap, NormalForm,
            SecondOrder, ExactIntegrator, Acceleration, Plotting, Optimizers)
    for n in names(mod)
        n === nameof(mod) && continue
        @eval export $n
    end
end

end # module
