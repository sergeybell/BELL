# =============================================================================
# NormalForm.jl - 6D Lie normal form of the one-turn map (with RF)
# =============================================================================
#
#
# Author: Kolokolchikov Sergey Dmitrievich <sergey.bell13@gmail.com>
#
# One-turn Lie factor (magnetic product plus thin RF) and its Dragt-Finn
# normal form.  The homological equation is solved on the generators.
# delta_eq is the average of laboratory delta over the torus of actions.
#
# =============================================================================
module NormalForm

using LinearAlgebra
using TaylorSeries
using ..LieAlgebra, ..Coordinates, ..Elements, ..Structure, ..DifferentialAlgebra, ..Tracking,
      ..HigherOrder, ..Longitudinal, ..Linear

export linear_normalizing_matrix, rotation_residual, nonlinear_normal_form,
       evaluate_map, rf_kick_taylor, one_turn_map_rf, fixed_point_6d, center_map_at,
       NormalFormRF, normal_form_rf, delta_eq_nf, tunes_nf,
       delta_eq_phase_condition

# ----- phase-1 linear (transverse only) --------------------------------------

"""
    linear_normalizing_matrix(po) -> Matrix{Float64} (6×6)

Transverse Courant-Snyder normalising matrix built from `periodic_optics` output
`po`. Each of the x and y blocks is [√β 0; −α/√β 1/√β]; the longitudinal
block is the identity. A⁻¹ R A is then a rotation in x and y for an uncoupled
ring without dispersion.
"""
function linear_normalizing_matrix(po)
    A = Matrix{Float64}(I, 6, 6)
    for (k, p) in ((1, po.x), (3, po.y))
        sb = sqrt(p.beta)
        A[k:k+1, k:k+1] = [sb 0.0; -p.alpha / sb 1 / sb]
    end
    return A
end

"""
    rotation_residual(R, po) -> Float64

max |(A⁻¹ R A)_block − rotation(μ)| over the x and y blocks, where A =
[`linear_normalizing_matrix`](@ref)(po). It checks the linear normalisation.
"""
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

"""
    evaluate_map(M, z) -> Vector{Float64}

Value of the polynomial map `M` at the numeric point `z`.
"""
evaluate_map(M, z) = [Float64(evaluate(M[i], z)) for i in 1:6]

# ----- RF one-turn Taylor map ------------------------------------------------

"""
    rf_kick_taylor(cav::RFCavity, ref, C) -> Vector{TaylorN}

The thin RF kick of [`rf_kick`](@ref) as a polynomial map in canonical
variables (same formulas evaluated on `Z`), truncated to MAP_ORDER.
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

"""
    one_turn_map_rf(lat_mag, cav, ref, C; ds = 0.01) -> Vector{TaylorN}

Polynomial one-turn map with RF: the magnetic ring map (`ring_taylor_map`)
followed by [`rf_kick_taylor`](@ref) (`compose_maps(Mrf, Mm)`).
"""
function one_turn_map_rf(lat_mag::Lattice, cav::RFCavity, ref::RefParticle, C::Real; ds = 0.01)
    Mm = ring_taylor_map(Lattice(lat_mag.name, lat_mag.elements, ref, lat_mag.harmonic); ds = ds)
    Mrf = rf_kick_taylor(cav, ref, C)
    return compose_maps(Mrf, Mm)   # magnetic first, then RF
end

"""
    fixed_point_6d(M; z0 = zeros(6), iters = 40, tol = 1e-14) -> Vector{Float64}

6D fixed point M(z*) = z* by Newton's method with the exact polynomial Jacobian.
With RF the map has a non-zero fixed point (synchronous-particle offset).
"""
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

"""
    center_map_at(M, z*) -> Vector{TaylorN}

Map about its fixed point: M_c(w) = M(z* + w) − z*, truncated to MAP_ORDER.
M_c has no constant term, so it can be Dragt-Finn factorised.
"""
function center_map_at(M, zstar)
    zs = [Z[i] + zstar[i] for i in 1:6]
    Mc = [truncate_to(evaluate(M[i], zs) - zstar[i], MAP_ORDER) for i in 1:6]
    return Mc
end

# ----- linear 6D normalisation (3 planes) ------------------------------------

"""
    linear_A_6d(po, pe, Qs, h, C, ref) -> (A, β_s)

Approximate 6D linear normalisation: a dispersion shear (x, px, y, py) += D4 δ,
Courant-Snyder blocks in x and y, and a synchrotron "Twiss" for (l, δ) with
β_s = Q_s/(h|η0|), mapping l to φ = 2π h l/C. Used by [`normal_form_rf`](@ref).
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
    delta_eq_phase_condition(M, A, z*, ref, C, h, β_s; Jx = 0, Jy = 0,
                             delta_m = 1e-3, ntheta = 12, niter = 10) -> δ_eq

Equilibrium momentum offset from the condition ⟨Δφ⟩ = 0 per turn, averaged over
a torus (12³ angle grid) with betatron actions Jx, Jy and synchrotron amplitude
`delta_m`. The offset δ_eq of the torus centre is found by a secant/Newton
iteration on the mean phase slip.
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

"""
    NormalFormRF

Result of [`normal_form_rf`](@ref) (approximate, map-based route).

# Fields
- `zstar` – 6D fixed point; `beta_s` – synchrotron β; `A`, `R` – normalising and
  normalised linear matrices.
- `f3`, `f4` – generators of the normalised map; `F3` – removing generator;
  `h3`, `h4` – kept (resonant) parts.
- `pe` – [`path_expansion`](@ref) of the RF-off map; `xi` – [`chromaticity_orders`](@ref);
  `Qs0` – smooth synchrotron tune; `C` [m]; `h_rf`.
- `c0, c_x, c_y, c_s` – fitted δ_eq = c0 + c_x Jx + c_y Jy + c_s Js.
- `c_s_note, c_x_note, c_y_note` – analytic counterparts
  (−η1/(η0 β_s), (2π/C) ξ_x/η0, (2π/C) ξ_y/η0).
"""
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
    normal_form_rf(lat_mag, ref, C, h, V_MV, phi_s; T = 1.0, ds = 0.01) -> NormalFormRF

6D normal-form analysis of the one-turn map with a thin RF kick of voltage V·T.

# Algorithm
1. M = [`one_turn_map_rf`](@ref); z* = [`fixed_point_6d`](@ref); M_c = map about z*.
2. Compute slip and chromaticity from the RF-off Taylor map; take periodic optics
   from the ring matrix and Q_s from `synchrotron_tune`.
3. Use A = [`linear_A_6d`](@ref) to normalise M_c, then Dragt-Finn factorise it.
4. Fit δ_eq on 7 tori (different Jx, Jy, δ_m) with
   [`delta_eq_phase_condition`](@ref) and compare it with the analytic coefficients.

For the full Dragt-Finn normal form use [`dragt_finn_nf6d`](@ref).
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
    delta_eq_nf(nf::NormalFormRF; Jx = 0, Jy = 0, Js = nothing, delta_m = nothing) -> δ_eq

Fitted δ_eq = c0 + c_x Jx + c_y Jy + c_s Js; give either `Js` or the synchrotron
amplitude `delta_m` (Js = β_s δ_m²/2).
"""
function delta_eq_nf(nf::NormalFormRF; Jx = 0.0, Jy = 0.0, Js = nothing,
                     delta_m = nothing)
    if Js === nothing
        delta_m === nothing && error("provide Js or delta_m")
        Js = nf.beta_s * delta_m^2 / 2
    end
    return nf.c0 + nf.c_x * Jx + nf.c_y * Jy + nf.c_s * Js
end

"""
    tunes_nf(nf::NormalFormRF; Jx = 0, Jy = 0, Js = 0) -> (Qx, Qy, Qs, Qs_smooth)

Linear fractional tunes from the normalised matrix `nf.R` and the smooth Q_s.
Amplitude dependence is not included here; see [`nf6d_tunes`](@ref).
"""
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
"""
    nonlinear_normal_form(args...; kw...)

Alias of [`normal_form_rf`](@ref) (kept for older scripts).
"""
nonlinear_normal_form(args...; kwargs...) = normal_form_rf(args...; kwargs...)


# =============================================================================
# 6D Dragt-Finn normal form of the one-turn Lie factor (thin RF)
# =============================================================================
#
#   One-turn object: product of element factors, then the thin RF factor.
#   R, A: symplectic eigen-normalisation, A^{-1} R A = rotations r(μ).
#   Generators in normal coordinates: f_n(w) = f(A w).
#   Order 3: F3 − F3 ∘ r = −f3  (divisor 1 − e^{-i(j-n)·μ} on monomials).
#   Order 4: f4' = f4 + (1/2)[f3, F3];  F4 − F4 ∘ r = −f4'_nonres;
#            h4 = resonant part = sum c_ij J_i J_j.
#   Lab coordinates  z = z* + A (exp(:F3:) o exp(:F4:))(w)  on the torus J.
#   delta_eq(J) = z*_6 + <(A T(w))_6>_theta = z*_6 + d·J
# =============================================================================

export NF6D, dragt_finn_nf6d, nf6d_delta_eq, nf6d_actions, nf6d_tunes, symplectic_eigen_A,
       one_turn_factor_rf, homological_operator

_mono(e) = prod(Z[i]^e[i] for i in 1:6)
_hp_basis(k) = [_mono(e) for e in exponents_of_order(k)]
_hcoeffs(p::TaylorN, k) = [Float64(getcoeff(p, e)) for e in exponents_of_order(k)]

"""
    symplectic_eigen_A(R) -> (A, μ)

Symplectic eigen-normalisation of a stable 6×6 one-turn matrix: A⁻¹ R A =
[`rotation_matrix_6`](@ref)(μ), and Aᵀ S A = S.

# Algorithm
For each eigenvalue with Im λ > 0 and eigenvector v = a + ib, normalise so that
a·S b = 1 (flip b and μ if negative). Assign the planes as follows: s is the pair
with the largest δ-component weight, and x is the pair (of the other two) with
the larger weight in (x, px). The columns of A are (a, b) per plane, and μ = arg λ.
Throws if the matrix is not stable (fewer than 3 elliptic pairs). Coupling is
handled naturally.
"""
function symplectic_eigen_A(R)
    S = SYMPLECTIC_S
    ev = eigen(R)
    pairs = []
    for i in 1:6
        imag(ev.values[i]) > 0 || continue
        v = ev.vectors[:, i]; a = real.(v); b = imag.(v)
        n = dot(a, S * b); mu = angle(ev.values[i])
        if n < 0
            b = -b; n = -n; mu = -mu
        end
        push!(pairs, (a = a ./ sqrt(n), b = b ./ sqrt(n), mu = mu))
    end
    length(pairs) == 3 || error("one-turn matrix not stable")
    w(p, r) = sum(abs2, p.a[r]) + sum(abs2, p.b[r])
    is = argmax([abs(p.a[6]) + abs(p.b[6]) for p in pairs] ./ [sqrt(w(p, 1:6)) for p in pairs])
    rest = setdiff(1:3, is)
    ix = w(pairs[rest[1]], 1:2) > w(pairs[rest[2]], 1:2) ? rest[1] : rest[2]
    iy = setdiff(rest, ix)[1]
    A = zeros(6, 6); mu = zeros(3)
    for (k, j) in enumerate((ix, iy, is))
        A[:, 2k-1] = pairs[j].a; A[:, 2k] = pairs[j].b; mu[k] = pairs[j].mu
    end
    return A, mu
end

"""
    NF6D

Dragt-Finn 6D normal form of a one-turn map ([`dragt_finn_nf6d`](@ref)).

# Fields
- `zstar::Vector{Float64}` – fixed point.
- `A::Matrix{Float64}` – symplectic normalising matrix; `mu::Vector{Float64}` – phase advances (x, y, s) [rad].
- `Rn_residual`, `symp_residual` – |A⁻¹RA − r(μ)|, |AᵀSA − S|.
- `F3`, `F4::TaylorN` – generators removing the non-resonant cubic and quartic terms.
- `f3_before`, `f3_after`, `f4_nonres_after` – max coefficients before/after (quality check).
- `c::Matrix{Float64}` – detuning: h4 = Σ_{i≤j} c[i,j] J_i J_j.
- `d::Vector{Float64}` – δ_eq(J) = z*₆ + d·J; `d_grid` – same from a numerical torus average.
- `T`, `Tinv` – normal-form transformation w ↦ u and its inverse (polynomial maps);
  `U6` – δ-component of A·T(w).
- `sigma_delta_s`, `sigma_l_s` – δ and l amplitudes per √(2Js).
"""
struct NF6D
    zstar::Vector{Float64}
    A::Matrix{Float64}
    mu::Vector{Float64}          # linear phase advances (rad), planes x,y,s
    Rn_residual::Float64         # |A^-1 R A - rotations|
    symp_residual::Float64       # |A^T S A - S|
    F3::TaylorN{Float64}
    F4::TaylorN{Float64}
    f3_before::Float64           # max|coef| of f3 before / after
    f3_after::Float64
    f4_nonres_after::Float64     # max|f4 - h4| after order-4 step
    c::Matrix{Float64}           # h4 = sum_{i<=j} c[i,j] J_i J_j
    d::Vector{Float64}           # delta_eq = z*_6 + d . J
    d_grid::Vector{Float64}      # same, from numerical torus average (check)
    T::Vector{TaylorN{Float64}}  # w -> u
    Tinv::Vector{TaylorN{Float64}}
    U6::TaylorN{Float64}
    sigma_delta_s::Float64       # delta amplitude per sqrt(2 Js): sqrt(A65^2+A66^2)
    sigma_l_s::Float64           # l amplitude per sqrt(2 Js)
end

_lie_map_of(F) = [truncate_to(lie_exp(F, Z[i]), MAP_ORDER) for i in 1:6]
_conj(Mx, F) = compose_maps(_lie_map_of(-F), compose_maps(Mx, _lie_map_of(F)))

"""
    one_turn_factor_rf(lat_mag, cav::RFCavity, ref, C) -> LieFactor

One-turn Lie factor with RF: [`ring_lie_factor`](@ref) of the magnetic ring,
concatenated with the Dragt-Finn factor of [`rf_kick_taylor`](@ref)
(read from the polynomial of the kick itself).
"""
function one_turn_factor_rf(lat_mag::Lattice, cav::RFCavity, ref::RefParticle, C::Real)
    Fmag = ring_lie_factor(lat_mag)
    Mrf = rf_kick_taylor(cav, ref, C)
    Frf = dragt_finn_factor(Mrf; name = cav.name)
    return concatenate(Fmag, Frf; name = lat_mag.name * "+RF")
end

"""
    homological_operator(k, μ) -> (L, basis, Rrot)

Matrix of the linear operator L F = F − F ∘ r(μ) on homogeneous polynomials of
degree k in the monomial basis [`exponents_of_order`](@ref)(k). Its eigenvalues
in the complex basis are 1 − exp(−i(j − n)·μ). Non-resonant terms are invertible;
resonant ones (zero divisors) stay in the normal form.
"""
function homological_operator(k::Int, mu)
    basis = _hp_basis(k)
    Rrot = rotation_matrix_6(mu)
    L = hcat((_hcoeffs(e - pull_generator(e, Rrot), k) for e in basis)...)
    return L, basis, Rrot
end

"""
    _nf6d_from_factor(F::LieFactor, z*) -> NF6D

Core of the normal-form algorithm.

1. (A, μ) = [`symplectic_eigen_A`](@ref)(F.R); generators in normal coordinates f_n = f ∘ A.
2. Order 3: solve L₃ F3 = −f3 (least squares; no cubic resonances for generic tunes).
3. Order 4: f4' = f4 + ½[f3, F3]; solve [L₄  −B] (F4, c) = −f4' with B the basis of
   J_iJ_j products (pseudo-inverse). The resonant part h4 = Σ c_ij J_i J_j is the detuning.
4. T = exp(:F3:) ∘ exp(:F4:); U6 = (A T)₆; d_k = ½(∂²U6/∂w²_{2k−1} + ∂²U6/∂w²_{2k}) gives δ_eq(J);
   it is checked against a numerical torus average.
"""
function _nf6d_from_factor(F::LieFactor, zstar)
    R = F.R
    A, mu = symplectic_eigen_A(R)
    Rrot = rotation_matrix_6(mu)
    Rn = inv(A) * R * A
    f3n = pull_generator(F.f3, A)
    f4n = pull_generator(F.f4, A)
    L3, b3, _ = homological_operator(3, mu)
    f30 = _hcoeffs(f3n, 3)
    x3 = L3 \ (-f30)
    F3 = sum(x3[j] * b3[j] for j in eachindex(b3))
    f3_after = f3n + F3 - pull_generator(F3, Rrot)
    f4p = hom(f4n + poisson(f3n, F3) / 2, 4)
    L4, b4, _ = homological_operator(4, mu)
    f40 = _hcoeffs(f4p, 4)
    Iq = [(Z[2k-1]^2 + Z[2k]^2) / 2 for k in 1:3]
    pairs = [(i, j) for i in 1:3 for j in i:3]
    B = hcat([_hcoeffs(Iq[i] * Iq[j], 4) for (i, j) in pairs]...)
    x = pinv([L4 -B]) * (-f40)
    nb = length(b4)
    F4 = sum(x[j] * b4[j] for j in 1:nb)
    cvec = x[nb+1:end]
    cmat = zeros(3, 3)
    for (k, (i, j)) in enumerate(pairs)
        cmat[i, j] = cvec[k]
    end
    h4 = sum(cvec[k] * Iq[i] * Iq[j] for (k, (i, j)) in enumerate(pairs))
    f4_after = f4p + F4 - pull_generator(F4, Rrot)
    T = compose_maps(_lie_map_of(F3), _lie_map_of(F4))
    Tinv = compose_maps(_lie_map_of(-F4), _lie_map_of(-F3))
    U6 = sum(A[6, j] * T[j] for j in 1:6)
    d = [0.5 * (constant_term(differentiate(differentiate(U6, 2k-1), 2k-1)) +
                constant_term(differentiate(differentiate(U6, 2k), 2k))) for k in 1:3]
    dg = zeros(3)
    for k in 1:3
        J = zeros(3); J[k] = 1e-8
        dg[k] = _torus_mean(U6, J) / 1e-8
    end
    sds = hypot(A[6, 5], A[6, 6]); sls = hypot(A[5, 5], A[5, 6])
    return NF6D(zstar, A, mu, maximum(abs, Rn - Rrot),
                maximum(abs, A' * SYMPLECTIC_S * A - SYMPLECTIC_S),
                F3, F4, maximum(abs, f30), maxabs(f3_after), maxabs(f4_after - h4),
                cmat, d, dg, T, Tinv, U6, sds, sls)
end

"""
    dragt_finn_nf6d(F::LieFactor) -> NF6D
    dragt_finn_nf6d(M) -> NF6D

Normal form of a one-turn factor (primary) or of a polynomial map (the
map is first read as a Dragt-Finn factor).  The homological equation is
solved on the generators: F − F ∘ r(μ) = −f_nonres.
"""
dragt_finn_nf6d(F::LieFactor) = _nf6d_from_factor(F, zeros(6))

"""
    dragt_finn_nf6d(F::LieFactor) -> NF6D
    dragt_finn_nf6d(M::AbstractVector) -> NF6D

Dragt-Finn 6D normal form. The primary input is a one-turn factor (e.g.
[`one_turn_factor_rf`](@ref)); a polynomial map is first centred on its fixed
point and factorised. The homological equation F − F ∘ r(μ) = −f_nonres is solved
on the generators. See [`_nf6d_from_factor`](@ref Bell.NormalForm._nf6d_from_factor).
"""
function dragt_finn_nf6d(M::AbstractVector)
    zstar = fixed_point_6d(M)
    Mc = center_map_at(M, zstar)
    return _nf6d_from_factor(dragt_finn_factor(Mc), zstar)
end

function _torus_mean(P, J; nth = 8)
    acc = 0.0; n = 0
    r = sqrt.(2 .* J)
    for i in 0:nth-1, j in 0:nth-1, k in 0:nth-1
        t = 2pi .* (i, j, k) ./ nth
        w = [r[1] * cos(t[1]), -r[1] * sin(t[1]), r[2] * cos(t[2]), -r[2] * sin(t[2]),
             r[3] * cos(t[3]), -r[3] * sin(t[3])]
        acc += Float64(evaluate(P, w)); n += 1
    end
    return acc / n
end

"""
    nf6d_delta_eq(nf, J) -> δ_eq

Equilibrium momentum offset δ_eq(J) = z*₆ + d·J for actions J = (Jx, Jy, Js).
"""
nf6d_delta_eq(nf::NF6D, J) = nf.zstar[6] + dot(nf.d, J)

"""
    nf6d_actions(nf, z; linear = false) -> [Jx, Jy, Js]

Normal-form actions of a lab-frame canonical point:
u = A⁻¹(z − z*), w = T⁻¹(u) (skipped if `linear`), J_k = (w²_{2k−1} + w²_{2k})/2.
"""
function nf6d_actions(nf::NF6D, z; linear = false)
    u = nf.A \ (collect(z) .- nf.zstar)
    w = linear ? u : [Float64(evaluate(nf.Tinv[i], u)) for i in 1:6]
    return [(w[2k-1]^2 + w[2k]^2) / 2 for k in 1:3]
end

"""
    nf6d_tunes(nf, J) -> (mu, Q, dmu)

Amplitude-dependent phase advances μ_k(J) = μ_k − ∂h4/∂J_k and tunes μ/2π mod 1.
"""
function nf6d_tunes(nf::NF6D, J)
    c = nf.c
    g = [2c[k, k] * J[k] + sum(c[min(k, j), max(k, j)] * J[j] for j in 1:3 if j != k) for k in 1:3]
    muJ = nf.mu .- g
    return (mu = muJ, Q = mod.(muJ ./ (2pi), 1.0), dmu = -g)
end

end # module
