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
#   Jets         TaylorN jet space (order 3 in 6 variables), printing helpers
#   Elements     element types + Lattice
#   Hamiltonians H(z) of each element, expanded to order 3
#   LieOperators Poisson bracket, exp(:f:), Dragt-Finn factors, products
#   OptiMParser  .opt reader
#   LinearOptics exp(L S A) matrices, closed Twiss, dispersion, tunes, alpha_c
#   Symplectic   J^T S J - S checks
#   Tracking     factor-by-factor symplectic tracking (+ thin RF kick)
#   TaylorMap    independent truncated Taylor map of the ring (RK4 in DA)
#   NormalForm   linear normalisation; nonlinear part planned (phase 4)
#   Plotting     Twiss + layout PNG, phase-2 plots
#   SecondOrder  D2, alpha0/alpha1, eta0/eta1, xi1/xi2, longitudinal H (phase 2)
#   ExactIntegrator  exact-H symplectic integrator incl. finite RF cavity (phase 2)
#   Optimizers   OPTIONAL: sextupole response + example solutions A/B (not used by default)
#   Acceleration reference-energy ramp, renormalisation, tune measurement (phase 2)
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
include("TaylorMap.jl")       # independent truncated Taylor map (RK4 in DA)
include("SecondOrder.jl")     # D2, alpha0/1, eta0/1, xi1/2, longitudinal H, separatrix
include("ExactIntegrator.jl") # exact-H integrator incl. finite RF cavity
include("Acceleration.jl")    # reference ramp, tune measurement
include("NormalForm.jl")      # 6D Lie normal form with RF
include("Plotting.jl")        # all plots
include("Optimizers.jl")      # OPTIONAL optimisation (not used by default)

using .Constants, .Jets, .Elements, .Hamiltonians, .LieOperators, .OptiMParser,
      .LinearOptics, .Symplectic, .Tracking, .TaylorMap, .NormalForm, .SecondOrder,
      .ExactIntegrator, .Acceleration, .Plotting, .Optimizers

for mod in (Constants, Jets, Elements, Hamiltonians, LieOperators, OptiMParser,
            LinearOptics, Symplectic, Tracking, TaylorMap, NormalForm, SecondOrder,
            ExactIntegrator, Acceleration, Plotting, Optimizers)
    for n in names(mod)
        n === nameof(mod) && continue
        @eval export $n
    end
end

end # module
