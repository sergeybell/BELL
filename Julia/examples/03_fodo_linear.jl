# =============================================================================
# 03_fodo_linear.jl - standard linear-optics procedure for one FODO cell
# =============================================================================
#
# WHAT THIS EXAMPLE TEACHES
#   The standard procedure of linear beam optics for a periodic lattice, done
#   "by hand" with 2x2 / 3x3 transfer matrices, for one FODO cell of a small
#   ring (16 cells, thick quadrupoles + sector dipoles):
#   periodic Twiss parameters, dispersion, phase advance / tune, chromaticity,
#   and plots of beta_x, beta_y and D_x along the cell.
#
# PHYSICS BACKGROUND
#   Coordinates: horizontal (x, x', delta) with x' = dx/ds and delta = dp/p0
#   (relative momentum deviation), vertical (y, y').  s = path length.
#   Linear equation of motion (Hill's equation):
#       x'' + K(s) x = h(s) delta,    h = 1/rho (curvature in the dipoles)
#     focusing quad   (QF):  Kx = +k,     Ky = -k
#     defocusing quad (QD):  Kx = -k,     Ky = +k
#     sector dipole   (B) :  Kx = 1/rho^2 (weak focusing), Ky = 0 (drift), h = 1/rho
#     drift           (O) :  Kx = Ky = 0
#   Each element is a matrix acting on (x, x') [or (x, x', delta)]; the cell
#   matrix is the product in beam order, M = M_n * ... * M_2 * M_1.
#
#   TWISS PARAMETERS.  In a periodic lattice a particle oscillates around the
#   design orbit (betatron oscillation):
#       x(s) = sqrt(2J beta(s)) cos(mu(s) + phi0)
#   beta(s) [m] is the envelope ("beta function"): the beam size is
#   sigma = sqrt(eps * beta); alpha = -beta'/2 is its slope; gamma =
#   (1 + alpha^2)/beta.  mu(s) = int ds/beta is the phase advance.  The
#   one-cell matrix of a stable periodic cell can always be written as
#       M = [cos mu + alpha sin mu        beta sin mu;
#            -gamma sin mu                cos mu - alpha sin mu]
#   so beta, alpha, mu at the cell start are read off from M ("closure" or
#   "periodic solution").  It exists only if |tr M| < 2 (stability).
#   Tune Q = (phase advance of the whole ring)/2pi = N mu_cell / 2pi.
#   DISPERSION D(s): the orbit shift per unit delta, x_delta = D(s) delta.
#   CHROMATICITY xi = dQ/d(delta) = -1/(4pi) int beta(s) K(s) ds
#   (here only the quadrupole part, dipole terms neglected).
#
# PROCEDURE (= SECTIONS OF THIS FILE)
#   1. Build the cell from a list of thick elements and compute the cell
#      matrices Mx (3x3, with dispersion) and My (2x2)
#   2. Closure: stability check |tr M| < 2, periodic Twiss parameters
#      (beta, alpha, gamma, mu) and periodic dispersion (D, D') from M
#   3. Propagate beta, alpha, D, D' and the phase mu along s in fine slices:
#      B1 = M B0 M^T  with  B = [beta -alpha; -alpha gamma]
#   4. Closure check: the propagated functions at the cell end must equal the
#      start values (they are periodic), and the accumulated phase must equal
#      mu from the trace
#   5. Main parameters: tunes, beta / D extrema, natural chromaticity
#   6. Plots of beta_x, beta_y, D_x and the cell layout -> output/*.png
#
# Cell (symmetric, starts and ends in the middle of the focusing quad QF):
#   QF/2 - O - B - O - QD - O - B - O - QF/2
#
# HOW TO RUN (from the Julia folder)
#   julia --project=. examples/03_fodo_linear.jl
#   Needs Plots.jl (listed in Project.toml).  The first run takes a minute or
#   more because Plots is precompiled; later runs take a few seconds.
#
# EXPECTED OUTPUT / HOW TO READ IT
#   * cell matrices; det M = 1 (symplectic) and |tr M| < 2 -> STABLE
#   * periodic solution at the cell start: mu_x = 57.5 deg, mu_y = 48.6 deg,
#     beta_x = 10.99 m, beta_y = 5.49 m, alpha = 0 (symmetry point), D_x = 3.46 m
#   * a table of beta, alpha, D and mu/2pi at the end of every element
#   * closure check: |end - start| ~ 1e-13 -> "Twiss functions CLOSED"
#   * ring tunes Qx = 2.55418, Qy = 2.16017; natural chromaticity ~ -2
#   * two PNG files in output/: one cell, and three cells (periodicity)

using LinearAlgebra, Printf

# -----------------------------------------------------------------------------
# Printing helpers (ASCII only - the Windows console garbles Unicode).
#   banner(title, lines...)   big title + "What this example does"
#   section(title, lines...)  section heading + 1-3 lines of explanation
#   print_summary(lines...)   "Summary / what you learned" at the end
#   showmat(M)                print a matrix row by row (no Unicode header)
# -----------------------------------------------------------------------------
function banner(title, lines...)
    println("="^72)
    println(" ", title)
    println("="^72)
    println("What this example does:")
    foreach(l -> println("  ", l), lines)
end

function section(title, lines...)
    println()
    println("--- ", title, " ", "-"^max(0, 66 - length(title)))
    foreach(l -> println("  > ", l), lines)
    isempty(lines) || println()
end

function print_summary(lines...)
    println()
    println("="^72)
    println(" Summary / what you learned")
    println("="^72)
    foreach(l -> println("  * ", l), lines)
end

showmat(M) = (Base.print_matrix(IOContext(stdout, :compact => true, :limit => true), M); println())

banner("Example 03: linear optics of a FODO cell (Twiss, dispersion, tunes)",
       "Builds one FODO cell (QF/2 - O - B - O - QD - O - B - O - QF/2) of a",
       "16-cell ring from thick-element matrices, finds the periodic (matched)",
       "Twiss parameters and dispersion, propagates them along the cell, checks",
       "that they close on themselves, prints tunes / beta / D / chromaticity",
       "and saves plots of beta_x, beta_y and D_x to the output/ folder.")

# =====================================================================
# 1. Parameters and element list.  Every element is a small record (name,
#    type, length, strength); the cell is just a vector of them in beam order.
#    The two dipoles per cell bend by 2pi/(2N) each, so N cells close the ring.
# =====================================================================
N     = 16                # number of cells in the ring
k     = 0.5               # quadrupole strength |k| [1/m^2]
Lq    = 0.5               # full quadrupole length [m]
Lo    = 0.5               # drift length [m]
rho   = 10.0              # dipole bending radius [m]
theta = 2pi / (2N)        # bend angle per dipole [rad] (2 dipoles per cell close the ring)
Lb    = rho * theta       # dipole arc length [m]

# An element is described by its name, type, length and strength.
struct Element
    name::String
    kind::Symbol          # :drift, :quad, :bend
    L::Float64            # length [m]
    k::Float64            # quad strength in the horizontal plane (+ focusing, - defocusing)
    rho::Float64          # bending radius (only for :bend)
end

Drift(name, L)     = Element(name, :drift, L, 0.0, Inf)
Quad(name, L, k)   = Element(name, :quad, L, k, Inf)
Bend(name, L, rho) = Element(name, :bend, L, 0.0, rho)

cell = [Quad("QF/2", Lq / 2, +k),
        Drift("O", Lo),
        Bend("B", Lb, rho),
        Drift("O", Lo),
        Quad("QD", Lq, -k),
        Drift("O", Lo),
        Bend("B", Lb, rho),
        Drift("O", Lo),
        Quad("QF/2", Lq / 2, +k)]

Lcell = sum(e.L for e in cell)

# =====================================================================
# Transfer matrices of a piece of length l of an element
# =====================================================================
# 2x2 matrix for x'' + K x = 0 (K > 0 focusing, K < 0 defocusing, K = 0 drift)
function mat2(K, l)
    if K > 0
        w = sqrt(K)
        return [cos(w * l)      sin(w * l) / w;
                -w * sin(w * l) cos(w * l)]
    elseif K < 0
        w = sqrt(-K)
        return [cosh(w * l)     sinh(w * l) / w;
                w * sinh(w * l) cosh(w * l)]
    else
        return [1.0 l;
                0.0 1.0]
    end
end

# Horizontal 3x3 matrix acting on (x, x', delta)
function matx(e::Element, l)
    if e.kind == :bend
        # sector dipole, angle phi = l / rho: weak focusing 1/rho^2 + dispersion terms
        phi = l / e.rho
        return [cos(phi)          e.rho * sin(phi)  e.rho * (1 - cos(phi));
                -sin(phi) / e.rho cos(phi)          sin(phi);
                0.0               0.0               1.0]
    else
        M = mat2(e.k, l)          # quad (Kx = k) or drift (k = 0); no dispersion generated
        return [M zeros(2); 0.0 0.0 1.0]
    end
end

# Vertical 2x2 matrix acting on (y, y')
function maty(e::Element, l)
    if e.kind == :bend
        return mat2(0.0, l)       # sector dipole is a drift vertically
    else
        return mat2(-e.k, l)      # quad: sign of k flips in the vertical plane
    end
end

# Cell matrices: product in beam order, M = M_n * ... * M_2 * M_1
# (foldl walks through the cell and multiplies each new element from the left)
Mx3 = foldl((M, e) -> matx(e, e.L) * M, cell; init = Matrix(1.0I, 3, 3))
My  = foldl((M, e) -> maty(e, e.L) * M, cell; init = Matrix(1.0I, 2, 2))
Mx  = Mx3[1:2, 1:2]

section("1. Build the cell: parameters and cell matrices",
        "Cell QF/2 - O - B - O - QD - O - B - O - QF/2, matrices multiplied in beam order.",
        "Mx acts on (x, x', delta): its 3rd column (M13, M23) is the dispersion",
        "generated by the dipoles.  My acts on (y, y').")
@printf("N = %d cells, k = %.3f 1/m^2, Lq = %.2f m, Lo = %.2f m\n", N, k, Lq, Lo)
@printf("rho = %.2f m, bend angle = %.5f rad (%.3f deg), L_dipole = %.4f m\n",
        rho, theta, rad2deg(theta), Lb)
@printf("cell length = %.4f m, ring circumference = %.4f m\n", Lcell, N * Lcell)

println("\nHorizontal cell matrix Mx (x, x', delta):")
showmat(round.(Mx3, digits = 6))
println("Vertical cell matrix My (y, y'):")
showmat(round.(My, digits = 6))

# =====================================================================
# 2. Closure: periodic Twiss parameters and dispersion at the cell start
#    We look for the beam ellipse / orbit that is reproduced by the cell
#    matrix (periodic solution).  It exists only if |tr M| < 2.
# =====================================================================
# A periodic 2x2 matrix can be written as
#   M = [cos(mu) + alpha sin(mu)    beta sin(mu);
#        -gamma sin(mu)             cos(mu) - alpha sin(mu)]
function periodic_twiss(M)
    c = tr(M) / 2
    abs(c) < 1 || error("cell is unstable: |tr M| = $(abs(tr(M))) >= 2")
    mu = acos(c) * sign(M[1, 2])          # sign of M12 fixes the sign of sin(mu)
    mu < 0 && (mu += 2pi)                 # keep 0 < mu < 2pi
    beta  = M[1, 2] / sin(mu)
    alpha = (M[1, 1] - M[2, 2]) / (2 * sin(mu))
    gamma = (1 + alpha^2) / beta
    return (mu = mu, beta = beta, alpha = alpha, gamma = gamma)
end

section("2a. Stability check",
        "A cell is stable (bounded oscillations) if |tr M| < 2.  det M = 1 is a",
        "check of the matrices (linear Hamiltonian motion preserves phase-space area).")
for (plane, M) in (("x", Mx), ("y", My))
    @printf("plane %s: det M = %.12f (must be 1),  |tr M| = %.6f -> %s\n",
            plane, det(M), abs(tr(M)), abs(tr(M)) < 2 ? "STABLE (< 2)" : "UNSTABLE")
end

tx0 = periodic_twiss(Mx)
ty0 = periodic_twiss(My)

# Periodic dispersion: (D, D') must repeat after one cell:
#   [D; D'] = M2x2 [D; D'] + [M13; M23]   =>   (I - M2x2) [D; D'] = [M13; M23]
D0 = (I - Mx) \ Mx3[1:2, 3]

section("2b. Periodic (matched) solution at the cell start (middle of QF)",
        "From M = [cos mu + alpha sin mu, beta sin mu; -gamma sin mu, cos mu - alpha sin mu]",
        "and (I - M)(D, D') = (M13, M23).  alpha = D' = 0 here because the cell is",
        "symmetric about its start (values ~1e-16 are rounding).")
@printf("mu_x    = %.6f rad = %.3f deg    mu_y    = %.6f rad = %.3f deg\n",
        tx0.mu, rad2deg(tx0.mu), ty0.mu, rad2deg(ty0.mu))
@printf("beta_x  = %9.5f m                 beta_y  = %9.5f m\n", tx0.beta, ty0.beta)
@printf("alpha_x = %9.2e                   alpha_y = %9.2e   (0 by symmetry)\n", tx0.alpha, ty0.alpha)
@printf("gamma_x = %9.5f 1/m               gamma_y = %9.5f 1/m\n", tx0.gamma, ty0.gamma)
@printf("D_x     = %9.5f m                 D'_x    = %9.2e   (0 by symmetry)\n", D0[1], D0[2])

# =====================================================================
# 3. Propagate Twiss functions and dispersion along s in fine slices
#    Starting from the periodic values, every element is cut into slices of
#    length <= ds_max; the slice matrix is applied repeatedly.  This gives
#    beta(s), alpha(s), D(s), mu(s) for the table and the plots.
# =====================================================================
# Twiss (beta) matrix B = [beta -alpha; -alpha gamma] transforms as  B1 = M B0 M'.
# Phase advance over a slice with matrix M:  dmu = atan(M12, beta0 M11 - alpha0 M12).
# Dispersion vector (D, D', 1) transforms with the 3x3 matrix.
ds_max = 0.02                                   # maximum slice length [m]

twiss_matrix(t) = [t.beta -t.alpha; -t.alpha t.gamma]
phase_step(M, B) = atan(M[1, 2], B[1, 1] * M[1, 1] + B[1, 2] * M[1, 2])  # B[1,2] = -alpha

s_list   = [0.0]
betax    = [tx0.beta];  alphax = [tx0.alpha]
betay    = [ty0.beta];  alphay = [ty0.alpha]
Dx       = [D0[1]];     Dpx    = [D0[2]]
mux      = [0.0];       muy    = [0.0]
chrom_x  = 0.0;         chrom_y = 0.0           # integrals of beta*K ds over the quads
elem_end = Int[]                                # index in s_list of each element end

let Bx = twiss_matrix(tx0), By = twiss_matrix(ty0), d = [D0; 1.0], s = 0.0, phx = 0.0, phy = 0.0
    for e in cell
        nslice = 2 * ceil(Int, e.L / (2 * ds_max))  # even number: element centre is sampled
        l = e.L / nslice
        Sx = matx(e, l)                         # the same matrix for every slice of e
        Sy = maty(e, l)
        for _ in 1:nslice
            bx_old, by_old = Bx[1, 1], By[1, 1]
            phx += phase_step(Sx[1:2, 1:2], Bx)
            phy += phase_step(Sy, By)
            Bx = Sx[1:2, 1:2] * Bx * Sx[1:2, 1:2]'
            By = Sy * By * Sy'
            d  = Sx * d
            s += l
            if e.kind == :quad              # chromaticity integral, trapezoidal rule
                global chrom_x += 0.5 * (bx_old + Bx[1, 1]) * (+e.k) * l
                global chrom_y += 0.5 * (by_old + By[1, 1]) * (-e.k) * l
            end
            push!(s_list, s)
            push!(betax, Bx[1, 1]); push!(alphax, -Bx[1, 2])
            push!(betay, By[1, 1]); push!(alphay, -By[1, 2])
            push!(Dx, d[1]);        push!(Dpx, d[2])
            push!(mux, phx);        push!(muy, phy)
        end
        push!(elem_end, length(s_list))
    end
end
section("3. Propagation of the Twiss functions along the cell",
        "Twiss matrix B = [beta -alpha; -alpha gamma] is transported as B1 = M B0 M^T,",
        "dispersion as (D, D', 1) -> Mx (D, D', 1).  Table: values at each element end;",
        "beta_x is large at QF, beta_y at QD; mu/2pi is the accumulated phase in turns.")
@printf("Propagated along s with %d slices (ds <= %.3f m)\n\n", length(s_list) - 1, ds_max)
@printf("%-6s %8s %9s %9s %9s %9s %8s %9s %9s\n",
        "elem", "s [m]", "beta_x", "alpha_x", "beta_y", "alpha_y", "D_x", "mu_x/2pi", "mu_y/2pi")
@printf("%-6s %8.4f %9.4f %9.4f %9.4f %9.4f %8.4f %9.5f %9.5f\n", "start",
        s_list[1], betax[1], alphax[1], betay[1], alphay[1], Dx[1], 0.0, 0.0)
for (e, i) in zip(cell, elem_end)
    @printf("%-6s %8.4f %9.4f %9.4f %9.4f %9.4f %8.4f %9.5f %9.5f\n", e.name,
            s_list[i], betax[i], alphax[i], betay[i], alphay[i], Dx[i], mux[i] / 2pi, muy[i] / 2pi)
end

# =====================================================================
# 4. Closure check: end of cell must equal start of cell
#    The propagated functions are periodic only if steps 2 and 3 are
#    consistent, so this is a strong test of the whole calculation.
# =====================================================================
section("4. Closure check (start vs end of the cell)",
        "The periodic solution must reproduce itself after one cell: |end - start|",
        "should be at rounding level (~1e-13), and the accumulated phase must equal",
        "mu from the trace of M.")
@printf("%-10s %14s %14s %12s\n", "quantity", "start", "end", "|end-start|")
checks = [("beta_x",  betax[1],  betax[end]),
          ("alpha_x", alphax[1], alphax[end]),
          ("beta_y",  betay[1],  betay[end]),
          ("alpha_y", alphay[1], alphay[end]),
          ("D_x",     Dx[1],     Dx[end]),
          ("D'_x",    Dpx[1],    Dpx[end])]
for (name, a, b) in checks
    @printf("%-10s %14.8f %14.8f %12.2e\n", name, a, b, abs(b - a))
end
maxdiff = maximum(abs(b - a) for (_, a, b) in checks)
@printf("phase: accumulated mu_x = %.8f rad, from trace mu_x = %.8f rad, diff = %.2e\n",
        mux[end], tx0.mu, abs(mux[end] - tx0.mu))
@printf("phase: accumulated mu_y = %.8f rad, from trace mu_y = %.8f rad, diff = %.2e\n",
        muy[end], ty0.mu, abs(muy[end] - ty0.mu))
closed = maxdiff < 1e-9
println(closed ? "Twiss functions CLOSED" : "Twiss functions NOT closed",
        @sprintf(": max |end-start| = %.2e", maxdiff))

# =====================================================================
# 5. Main parameters: tunes, extrema of beta and D, natural chromaticity
# =====================================================================
Qx_cell, Qy_cell = tx0.mu / 2pi, ty0.mu / 2pi
Qx_ring, Qy_ring = N * Qx_cell, N * Qy_cell
# Natural chromaticity: xi = -1/(4 pi) * integral beta*K ds.
# Only the quadrupole contribution is included (dipole terms neglected).
xi_x_cell = -chrom_x / 4pi
xi_y_cell = -chrom_y / 4pi

section("5. Main parameters",
        "Tune Q = N mu/2pi (betatron oscillations per turn); avoid integer/half-integer.",
        "Natural chromaticity xi = -1/(4pi) int beta K ds (quads only): how much Q",
        "changes per unit delta = dp/p; it is negative and later corrected by sextupoles.")
@printf("cell length                 : %.4f m\n", Lcell)
@printf("ring circumference (N = %d) : %.4f m\n", N, N * Lcell)
@printf("phase advance per cell      : mu_x = %.3f deg,  mu_y = %.3f deg\n",
        rad2deg(tx0.mu), rad2deg(ty0.mu))
@printf("tune per cell               : Qx = %.5f,  Qy = %.5f\n", Qx_cell, Qy_cell)
@printf("ring tune (N cells)         : Qx = %.5f,  Qy = %.5f\n", Qx_ring, Qy_ring)
@printf("beta_x max / min            : %.4f m (s = %.3f m) / %.4f m (s = %.3f m)\n",
        maximum(betax), s_list[argmax(betax)], minimum(betax), s_list[argmin(betax)])
@printf("beta_y max / min            : %.4f m (s = %.3f m) / %.4f m (s = %.3f m)\n",
        maximum(betay), s_list[argmax(betay)], minimum(betay), s_list[argmin(betay)])
@printf("D_x max / min               : %.4f m (s = %.3f m) / %.4f m (s = %.3f m)\n",
        maximum(Dx), s_list[argmax(Dx)], minimum(Dx), s_list[argmin(Dx)])
@printf("natural chromaticity (quads only, dipoles neglected):\n")
@printf("    per cell                : xi_x = %.4f,  xi_y = %.4f\n", xi_x_cell, xi_y_cell)
@printf("    ring                    : xi_x = %.4f,  xi_y = %.4f\n", N * xi_x_cell, N * xi_y_cell)
@printf("(QF centre at s = 0 and s = %.4f m, QD centre at s = %.4f m)\n", Lcell, Lcell / 2)

# =====================================================================
# 6. Plots: beta_x, beta_y (top), D_x (middle) and the element layout
#    (bottom) for one cell and for three cells.  The GR backend draws
#    off-screen and the figures are only saved as PNG files.
# =====================================================================
section("6. Plots",
        "beta_x, beta_y (top), D_x (middle) and the layout (QF red above, QD red",
        "below, dipoles blue) for 1 and 3 cells; the 3-cell plot shows periodicity.")
ENV["GKSwstype"] = "100"      # GR: draw off-screen (no window), we only save files
using Plots
gr()

outdir = joinpath(@__DIR__, "..", "output")
mkpath(outdir)

# Draw the cell layout into plot p, shifted by s0:
# dipoles = blue boxes on the axis, QF = red box above, QD = red box below, drifts = line.
function draw_layout!(p, cell, s0)
    s = s0
    for e in cell
        s1, s2 = s, s + e.L
        if e.kind == :bend
            plot!(p, Shape([s1, s2, s2, s1], [-0.3, -0.3, 0.3, 0.3]),
                  color = :royalblue, linecolor = :black, label = "")
            annotate!(p, (s1 + s2) / 2, 0.0, text("B", 8, :white))
        elseif e.kind == :quad
            h = e.k > 0 ? 1.0 : -1.0          # focusing above, defocusing below the axis
            plot!(p, Shape([s1, s2, s2, s1], [0, 0, h, h]),
                  color = :red, linecolor = :black, label = "")
            annotate!(p, (s1 + s2) / 2, h / 2, text(e.k > 0 ? "QF" : "QD", 7, :white))
        end
        s = s2
    end
    return p
end

function layout_panel(cell, ncells, Lcell)
    p = plot([0, ncells * Lcell], [0, 0], color = :black, lw = 1, label = "",
             ylims = (-1.3, 1.3), yticks = false, xlabel = "s [m]", ylabel = "layout",
             framestyle = :box)
    for n in 0:ncells-1
        draw_layout!(p, cell, n * Lcell)
    end
    return p
end

tune_text = @sprintf("mu_x = %.1f deg, mu_y = %.1f deg per cell;  Q_x = %.3f, Q_y = %.3f (ring, N = %d)",
                     rad2deg(tx0.mu), rad2deg(ty0.mu), Qx_ring, Qy_ring, N)

function twiss_figure(ncells, title)
    # repeat the one-cell solution ncells times (it is periodic)
    s  = vcat([s_list[1:end-1] .+ n * Lcell for n in 0:ncells-1]..., ncells * Lcell)
    rep(v) = vcat(repeat(v[1:end-1], ncells), v[end])
    bx, by, dx = rep(betax), rep(betay), rep(Dx)
    xl = (0, ncells * Lcell)

    p1 = plot(s, bx, lw = 2, color = :red, label = "beta_x", ylabel = "beta [m]",
              xlims = xl, ylims = (0, 1.3 * max(maximum(bx), maximum(by))),
              title = tune_text, titlefontsize = 9, legend = :top, legend_columns = 2,
              framestyle = :box, xformatter = _ -> "")
    plot!(p1, s, by, lw = 2, color = :blue, label = "beta_y")
    p2 = plot(s, dx, lw = 2, color = :green, label = "D_x", ylabel = "D_x [m]",
              xlims = xl, ylims = (0, 1.3 * maximum(dx)), legend = :top,
              framestyle = :box, xformatter = _ -> "")
    p3 = layout_panel(cell, ncells, Lcell)
    plot!(p3, xlims = xl)
    return plot(p1, p2, p3, layout = grid(3, 1, heights = [0.5, 0.3, 0.2]),
                size = (900, 750), plot_title = title, plot_titlefontsize = 13,
                left_margin = 5Plots.mm)
end

title1 = closed ? "FODO cell: periodic (closed) Twiss functions" :
                  "FODO cell: Twiss functions (NOT closed!)"
f1 = twiss_figure(1, title1)
file1 = abspath(joinpath(outdir, "fodo_twiss.png"))
savefig(f1, file1)
println("saved plot: ", file1)

f3 = twiss_figure(3, "FODO: 3 cells, the periodic Twiss functions repeat")
file3 = abspath(joinpath(outdir, "fodo_twiss_3cells.png"))
savefig(f3, file3)
println("saved plot: ", file3)

print_summary("A periodic cell has a unique matched beta, alpha, D, read off its matrix M.",
              "It exists only if |tr M| < 2; it is propagated with B1 = M B0 M^T.",
              "Closure (end = start) confirms the periodic solution; tune Q = N mu / 2pi.",
              @sprintf("This lattice: Qx = %.5f, Qy = %.5f, xi_x = %.2f, xi_y = %.2f (ring).",
                       Qx_ring, Qy_ring, N * xi_x_cell, N * xi_y_cell))
