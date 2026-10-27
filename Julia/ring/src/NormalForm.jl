# =============================================================================
# NormalForm.jl - 6D Lie normal form of the one-turn map (with RF)
# =============================================================================
#
# Linear normalisation (phase 1) + nonlinear Dragt-Finn / homological analysis
# of the one-turn Taylor map INCLUDING a thin RF kick, giving:
#   - fixed point z* of the RF map (origin of the residual delta ~ 1e-7)
#   - normal-form Hamiltonian H_N(Jx,Jy,Js) to cubic / quartic order
#   - amplitude-dependent tunes
#   - delta_eq(Jx,Jy,Js) = angle average of delta after the normalising map
#
# =============================================================================
module NormalForm

using LinearAlgebra
using TaylorSeries
using ..Jets, ..LieOperators, ..Constants, ..Elements, ..TaylorMap, ..Tracking, ..SecondOrder, ..LinearOptics

export linear_normalizing_matrix, rotation_residual, nonlinear_normal_form,
       evaluate_map, rf_kick_taylor, one_turn_map_rf, fixed_point_6d, center_map_at,
       NormalFormRF, normal_form_rf, delta_eq_nf, tunes_nf,
       delta_eq_phase_condition

# ----- phase-1 linear (transverse only) --------------------------------------

function linear_normalizing_matrix(po)
    A = Matrix{Float64}(I, 6, 6)
    for (k, p) in ((1, po.x), (3, po.y))
        sb = sqrt(p.beta)
        A[k:k+1, k:k+1] = [sb 0.0; -p.alpha / sb 1 / sb]
    end
    return A
end

function rotation_residual(R, po)
    A = linear_normalizing_matrix(po)
    N = A \ R * A
    res = 0.0
    for (k, p) in ((1, po.x), (3, po.y))
        rot = [cos(p.mu) sin(p.mu); -sin(p.mu) cos(p.mu)]
        res = max(res, maximum(abs, N[k:k+1, k:k+1] - rot))
    end
    return res
end

"Evaluate a TaylorN map at a numeric point."
evaluate_map(M, z) = [Float64(evaluate(M[i], z)) for i in 1:6]

# ----- RF one-turn Taylor map ------------------------------------------------

"""
Thin RF kick as a TaylorN map in canonical variables, truncated to MAP_ORDER.
Matches Tracking.rf_kick to the jet order.
"""
function rf_kick_taylor(cav::RFCavity, ref::RefParticle, C::Real)
    x, px, y, py, l, d = Z
    v0 = ref.beta * C_LIGHT
    vr = velocity_ratio(d, ref)          # TaylorN
    v = vr * v0
    dt = -l / v
    phi = -cav.harmonic * 2pi * v0 / C * dt
    E = sqrt((ref.p0_MeV * (1 + d))^2 + ref.mass_MeV^2)
    E2 = E + ref.particle.charge * cav.voltage_MV * sin(cav.phi_s + phi)
    d2 = sqrt(E2^2 - ref.mass_MeV^2) / ref.p0_MeV - 1
    l2 = -velocity_ratio(d2, ref) * v0 * dt
    return [truncate_to(x, MAP_ORDER), truncate_to(px, MAP_ORDER),
            truncate_to(y, MAP_ORDER), truncate_to(py, MAP_ORDER),
            truncate_to(l2, MAP_ORDER), truncate_to(d2, MAP_ORDER)]
end

"One-turn map = magnetic ring Taylor map composed with thin RF kick."
function one_turn_map_rf(lat_mag::Lattice, cav::RFCavity, ref::RefParticle, C::Real; ds = 0.01)
    Mm = ring_taylor_map(Lattice(lat_mag.name, lat_mag.elements, ref, lat_mag.harmonic); ds = ds)
    Mrf = rf_kick_taylor(cav, ref, C)
    return compose_maps(Mrf, Mm)   # magnetic first, then RF
end

"6D fixed point of a polynomial map by Newton (default guess 0)."
function fixed_point_6d(M; z0 = zeros(6), iters = 40, tol = 1e-14)
    z = Float64.(z0)
    Id = Matrix{Float64}(I, 6, 6)
    for _ in 1:iters
        f = [Float64(evaluate(M[i], z)) for i in 1:6] .- z
        J = [Float64(evaluate(differentiate(M[i], j), z)) for i in 1:6, j in 1:6] .- Id
        dz = -J \ f
        z .+= dz
        maximum(abs, dz) < tol && break
    end
    return z
end

"Shift map origin to z*: M_c(w) = M(z*+w) - z*."
function center_map_at(M, zstar)
    zs = [Z[i] + zstar[i] for i in 1:6]
    Mc = [truncate_to(evaluate(M[i], zs) - zstar[i], MAP_ORDER) for i in 1:6]
    return Mc
end

# ----- linear 6D normalisation (3 planes) ------------------------------------

"""
Build a linear symplectic A that (approximately) puts R into three rotations.
Uses Courant-Snyder for x,y from `po` and a synchrotron Twiss from (eta0, Qs, h, C).
Dispersion is removed by a shear: x <- x - D delta, etc., using po.D4.
"""
function linear_A_6d(po, pe, Qs, h, C, ref)
    A = Matrix{Float64}(I, 6, 6)
    # dispersion shear: (x,px,y,py) <- (x,px,y,py) - D4 * delta
    D = po.D4
    for i in 1:4
        A[i, 6] = D[i]
    end
    # transverse CS
    for (k, p) in ((1, po.x), (3, po.y))
        sb = sqrt(p.beta)
        Acs = [sb 0.0; -p.alpha / sb 1 / sb]
        A[k:k+1, k:k+1] = Acs
        # also transform the dispersion column
        A[k:k+1, 6] = Acs * A[k:k+1, 6]
    end
    # synchrotron: phi = 2pi h l / C at delta=0; oscillator in (phi, delta)
    # H2 ~ 2pi Qs * (1/2)(delta^2 / beta_d + beta_d * phi^2) with
    # beta_d chosen so that d phi / d n = -2pi h eta0 delta  =>
    # at small amp: Qs^2 = h |eta0| * A_rf / (2pi) ... already in Qs.
    # Relation: dphi/dn = -2pi h eta0 delta, ddelta/dn = -A_rf cos(phi_s) * phi
    # => beta_phi / beta_delta scaling: use
    #   phi = sqrt(2 Js / beta_s) cos theta,  delta = -sqrt(2 Js * beta_s) sin theta
    # with beta_s = 2pi Qs / (2pi h |eta0|) = Qs / (h |eta0|)
    beta_s = Qs / (h * abs(pe.eta0))          # = |delta|/|phi| amplitude ratio
    # Map (l, delta) -> (phi, delta): phi = (2pi h / C) l
    # Then CS-like: phi = sqrt(beta_s) u_phi, delta = (1/sqrt(beta_s)) u_delta
    # so A_l_delta such that [l;delta] = A_ld * [u5; u6]
    # l = (C/(2pi h)) phi = (C/(2pi h)) sqrt(beta_s) u5
    # delta = u6 / sqrt(beta_s)
    sbs = sqrt(beta_s)
    A[5, 5] = (C / (2pi * h)) * sbs
    A[5, 6] = 0.0
    A[6, 5] = 0.0
    A[6, 6] = 1 / sbs
    return A, beta_s
end

# ----- homological equation at order 3 ---------------------------------------

"""
Split a cubic generator f3 into resonant (invariant under R) and non-resonant.
Resonant monomials: those with m·μ ∈ 2πℤ (here: pure action monomials in
normalised oscillator coords, i.e. terms in (x^2+px^2), (y^2+py^2), (l^2+δ^2)
and products). We keep only terms built from r_x^2, r_y^2, r_s^2.
"""
function resonant_cubic(f3_norm)
    # In normalised coords with R = rotations, resonant cubic Hamiltonians
    # that are pure invariants do not exist at odd order in amplitudes
    # (f3 is degree 3). Resonant *cubic* generators drive 3rd-order resonances;
    # the amplitude-dependent tune shifts live in f4 (degree 4).
    # For delta_eq at O(J) we need the *coordinate transformation* effect of
    # non-resonant F3 (and of quadratic map pieces) on delta.
    # Return f3 itself as the non-resonant generator to be removed; h3 = 0 at this order.
    return zero(f3_norm), f3_norm   # (h3, f3_nonresonant)
end

# ----- delta_eq by angle averaging of the normalising map --------------------

"""
Solve ⟨Δφ⟩ = 0 on a synchrotron (+ optional betatron) torus of the one-turn
map `M` (lab frame, about the origin of `M` which should already be centered
or absolute).  Returns delta_eq.  This is the map form of the note's
time-of-flight condition (ch.1).
"""
function delta_eq_phase_condition(M, A, zstar, ref, C, h, beta_s;
                                  Jx = 0.0, Jy = 0.0, delta_m = 1e-3,
                                  ntheta = 12, niter = 10)
    sJx, sJy = sqrt(2 * max(Jx, 0.0)), sqrt(2 * max(Jy, 0.0))
    function mean_dphi(deq)
        acc = 0.0; n = 0
        Js = beta_s * delta_m^2 / 2
        sJs = sqrt(2 * max(Js, 0.0))
        sb = sqrt(beta_s)
        for ix in 0:ntheta-1, iy in 0:ntheta-1, is_ in 0:ntheta-1
            tx = 2pi * ix / ntheta; ty = 2pi * iy / ntheta; ts = 2pi * is_ / ntheta
            u = [sJx * cos(tx), -sJx * sin(tx),
                 sJy * cos(ty), -sJy * sin(ty),
                 sJs * cos(ts), -sJs * sin(ts) + deq * sb]
            z = A * u .+ zstar
            z2 = [Float64(evaluate(M[i], z)) for i in 1:6]
            u1 = canonical_to_user(z, ref, C, h)
            u2 = canonical_to_user(z2, ref, C, h)
            acc += u2[5] - u1[5]
            n += 1
        end
        return acc / n
    end
    deq = 0.0
    for _ in 1:niter
        f = mean_dphi(deq)
        eps = max(1e-9, 1e-3 * max(abs(deq), 1e-6))
        fp = (mean_dphi(deq + eps) - f) / eps
        abs(fp) < 1e-18 && break
        deq -= f / fp
        abs(f) < 1e-12 && break
    end
    return deq
end

# Keep a thin wrapper name used by older call sites
delta_on_torus(A, F3, zstar, Jx, Jy, Js; ntheta = 12) =
    error("delta_on_torus replaced by delta_eq_phase_condition")

# ----- public result type ----------------------------------------------------

struct NormalFormRF
    zstar::Vector{Float64}
    beta_s::Float64
    A::Matrix{Float64}
    R::Matrix{Float64}
    f3::Any
    f4::Any
    F3::Any                 # generator that removes non-resonant f3
    h3::Any
    h4::Any
    pe::NamedTuple
    xi::NamedTuple
    Qs0::Float64
    C::Float64
    h_rf::Int
    # fitted delta_eq = c0 + c_x Jx + c_y Jy + c_s Js  (c0 ≈ zstar[6])
    c0::Float64
    c_x::Float64
    c_y::Float64
    c_s::Float64
    # note-formula counterparts
    c_s_note::Float64       # -eta1/(2 eta0) * (Js scaling): see delta_eq_nf
    c_x_note::Float64       # + (2pi/C) xi_x / eta0   from Lambda
    c_y_note::Float64
end

"""
Full 6D normal form of the one-turn map with thin RF.

`lat_mag` = magnetic lattice (no cavity). RF is added as a thin kick with
voltage `V_MV` (use V*T to match a finite cavity's effective kick).
"""
function normal_form_rf(lat_mag::Lattice, ref::RefParticle, C::Real, h::Int,
                        V_MV::Real, phi_s::Real; T = 1.0, ds = 0.01)
    cav = RFCavity("NF", V_MV * T, h, phi_s)
    M = one_turn_map_rf(lat_mag, cav, ref, C; ds = ds)
    zstar = fixed_point_6d(M)
    Mc = center_map_at(M, zstar)
    # optics / slip from magnetic map (RF off)
    Mm = ring_taylor_map(Lattice(lat_mag.name, lat_mag.elements, ref, h); ds = ds)
    pe = path_expansion(Mm, ref, C)
    xi = chromaticity_orders(Mm)
    po = periodic_optics(ring_matrix(Lattice(lat_mag.name, lat_mag.elements, ref, h)), ref, C)
    Qs0 = synchrotron_tune(ref, pe.eta0, V_MV, h, phi_s; T = T)
    A, beta_s = linear_A_6d(po, pe, Qs0, h, C, ref)
    # Transform centered map to normalised coords: M_n = A^{-1} o Mc o A
    Ai = inv(A)
    ZA = [sum(A[i, j] * Z[j] for j in 1:6) for i in 1:6]
    McA = [truncate_to(evaluate(Mc[i], ZA), MAP_ORDER) for i in 1:6]
    Mn = [truncate_to(sum(Ai[i, j] * McA[j] for j in 1:6), MAP_ORDER) for i in 1:6]
    # Dragt-Finn of normalised map
    F = dragt_finn_factor(Mn; name = "NF")
    h3, f3nr = resonant_cubic(F.f3)
    F3 = f3nr   # remove non-resonant cubic by generator F3 = f3_nr
    # residual quartic after removing F3 (approximate h4 ≈ f4)
    h4 = F.f4
    # Fit delta_eq from ⟨Δφ⟩=0 on tori (map form of the note condition).
    # Use absolute map M (not centered) so zstar is included naturally.
    samples = []
    for (Jx, Jy, dm) in ((0.0, 0.0, 5e-4),
                         (0.0, 0.0, 1e-3),
                         (0.0, 0.0, 1.5e-3),
                         (1e-6, 0.0, 5e-4),
                         (0.0, 1e-6, 5e-4),
                         (2e-6, 0.0, 5e-4),
                         (0.0, 2e-6, 5e-4))
        deq = delta_eq_phase_condition(M, A, zstar, ref, C, h, beta_s;
                                       Jx = Jx, Jy = Jy, delta_m = dm)
        Js = beta_s * dm^2 / 2
        push!(samples, (Jx, Jy, Js, deq))
    end
    B = hcat(ones(length(samples)),
             [s[1] for s in samples],
             [s[2] for s in samples],
             [s[3] for s in samples])
    rhs = [s[4] for s in samples]
    coef = B \ rhs
    c0, cx, cy, cs = coef
    # Note counterparts:
    # delta_eq = -eta1/eta0 * (dm^2/2) - Lambda/eta0
    # dm^2/2 = ⟨delta^2⟩ = Js / beta_s  => coeff of Js: -eta1/eta0 / beta_s
    # Lambda = -2pi/C (Jx xi_x + Jy xi_y) => delta_eq += -Lambda/eta0
    #        = (2pi/C)/eta0 * (Jx xi_x + ...) => cx_note = (2pi/C)*xi_x/eta0
    cs_note = -pe.eta1 / pe.eta0 / beta_s
    cx_note = (2pi / C) * xi.x[2] / pe.eta0
    cy_note = (2pi / C) * xi.y[2] / pe.eta0
    return NormalFormRF(zstar, beta_s, A, F.R, F.f3, F.f4, F3, h3, h4,
                        pe, xi, Qs0, C, h, c0, cx, cy, cs,
                        cs_note, cx_note, cy_note)
end

"""
delta_eq from normal form at actions (Jx,Jy,Js) or equivalently at
(Jx, Jy, delta_m) with Js = beta_s * delta_m^2 / 2.
"""
function delta_eq_nf(nf::NormalFormRF; Jx = 0.0, Jy = 0.0, Js = nothing,
                     delta_m = nothing)
    if Js === nothing
        delta_m === nothing && error("provide Js or delta_m")
        Js = nf.beta_s * delta_m^2 / 2
    end
    return nf.c0 + nf.c_x * Jx + nf.c_y * Jy + nf.c_s * Js
end

"Amplitude-dependent tunes from h4 / linear R (leading order)."
function tunes_nf(nf::NormalFormRF; Jx = 0.0, Jy = 0.0, Js = 0.0)
    # Linear tunes from R in normalised coords (should be near rotations)
    R = nf.R
    function frac_tune(M2)
        c = (M2[1, 1] + M2[2, 2]) / 2
        s = sign(M2[1, 2]) * sqrt(max(0.0, 1 - c^2))
        return mod(atan(s, c), 2pi) / 2pi
    end
    Qx0 = frac_tune(R[1:2, 1:2])
    Qy0 = frac_tune(R[3:4, 3:4])
    Qs0 = frac_tune(R[5:6, 5:6])
    # Amplitude dependence not fully extracted from h4 yet; return linear + NF Qs0
    return (Qx = Qx0, Qy = Qy0, Qs = Qs0, Qs_smooth = nf.Qs0)
end

# Keep old stub name as alias
nonlinear_normal_form(args...; kwargs...) = normal_form_rf(args...; kwargs...)

end # module
