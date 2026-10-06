"""Non-inductive current diagnostics: neutral-beam current drive (NBCD) and
bootstrap current (2026-10-06).

Post-processing only -- nothing here feeds back on the operating point. Both
need the radial profiles that solve_operating_point stores when
profile_averaging is on (OperatingPoint.rho_profile, ne/Te/Ti profiles, the
per-beam birth rates S_j h_j(rho) and birth pitch xi0_j(rho)); with
profile_averaging off there are no gradients and current_drive returns None.

NBCD -- local 1-D Cordey / Start-Cordey model
---------------------------------------------
Steady-state slowing-down distribution of a mono-energetic beam with
electron + ion drag and ion pitch-angle scattering (no energy diffusion,
no orbit effects beyond the deposition smoothing already in h_j):

    j_f = Z_b e S h tau_s v_b xi0 J(u_c, beta)
    J   = int_0^1 u^3/(u^3+u_c^3) [u^3 (1+u_c^3)/(u^3+u_c^3)]^(beta/3) du

with u = v/v_b, u_c^3 = (E_c/E_b)^(3/2), tau_s the Spitzer slowing-down time
(dv/dt = -(v/tau_s)(1 + v_c^3/v^3)), E_c = 14.8 A_b T_e Zh^(2/3),
Zh = sum_j n_j Z_j^2 / (n_e A_j), and beta = Z_eff / Zbar,
Zbar = sum_j n_j Z_j^2 (m_b/m_j) / n_e (ratio of the ion pitch-angle
scattering rate to the ion speed drag; the l = 1 Legendre exponent, written
l(l+1) beta'/3 with beta' = Z_eff/(2 Zbar) elsewhere in the literature).

Electron shielding with the trapped-electron correction (Start & Cordey
1980; Wesson, Tokamaks):

    j_NB = j_f [1 - (Z_b/Z_eff)(1 - G)],
    G    = (1.55 + 0.85/Z_eff) sqrt(eps) - (0.20 + 1.55/Z_eff) eps,

eps = rho a / R0 (G clipped to [0, 1]). Total I = (V / 2 pi R0) <j>_vol.

Bootstrap -- Sauter, Angioni & Lin-Liu, Phys. Plasmas 6 (1999) 2834,
erratum 9 (2002) 5140, in the cylindrical form

    j_bs = -(1/B_p) [L31 dp/dr + L32 n_e dT_e/dr + L34 alpha n_i dT_i/dr]

with p the thermal pressure (fast ions excluded), f_t from Lin-Liu & Miller,
collisionalities nu*_e, nu*_i from the local profiles and an ASSUMED current
profile j ~ (1 - rho^2)^nu for q(rho) and B_p(rho): nu = q_cyl/q0 - 1 with
q_cyl = physics.safety_factor_cyl_edge and q0 an input,
B_p = mu0 I(rho) / L_p(rho), I(rho) = I_p [1 - (1-rho^2)^(nu+1)],
L_p = 2 pi rho a sqrt((1+kappa^2)/2); q(rho) = q_cyl(rho) [1 + (F-1) rho^2],
F = q95/q_cyl (q95 = physics.safety_factor_cyl_edge_arbitrary_A).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from . import physics

MU0 = 4.0e-7 * math.pi
_MASS_NUMBER_C = 12.0


@dataclass
class CurrentDrive:
    rho: list[float]
    j_nb_per_beam: list[list[float]]    # shielded NBCD per beam [A/m^2]
    j_fast_per_beam: list[list[float]]  # unshielded fast-ion current per beam [A/m^2]
    j_nb: list[float]                   # sum over beams [A/m^2]
    j_bs: list[float]                   # bootstrap [A/m^2]
    shielding: list[float]              # 1 - (Z_b/Z_eff)(1 - G)  (Z_b = 1)
    q: list[float]                      # assumed q(rho)
    nu_star_e: list[float]
    I_nb_per_beam_A: list[float]
    I_nb_A: float
    I_bs_A: float
    I_p_A: float
    I_ind_A: float                      # I_p - I_nb - I_bs (inductive remainder, can be < 0)
    f_ni: float                         # (I_nb + I_bs) / I_p
    f_bs: float
    eta_cd: float                       # I_nb R0 n_line,20 / P_NB,inj  [1e20 A W^-1 m^-2]
    eta_cd_per_beam: list[float]
    n_line_m3: float
    q0: float
    q_edge: float                       # q95 (arbitrary-A formula)
    q_cyl: float                        # cylindrical edge q (sets nu_j)
    nu_j: float                         # current-profile exponent
    beta_scatter: list[float] = field(default_factory=list)   # beta = Z_eff/Zbar of beam 1, on rho


def _pitch_integral(uc3: np.ndarray, beta: np.ndarray, n: int = 200) -> np.ndarray:
    """J(u_c, beta) on arrays (same shape), midpoint rule in u."""
    u = (np.arange(n) + 0.5) / n
    u3 = u[None, :] ** 3
    c = uc3.reshape(-1, 1)
    b = beta.reshape(-1, 1)
    xi = (u3 * (1.0 + c) / (u3 + c)) ** (b / 3.0)
    return (np.sum(u3 / (u3 + c) * xi, axis=1) / n).reshape(uc3.shape)


def _sauter_coefficients(ft, nue, nui, Z):
    """L31, L32, L34, alpha (Sauter 1999 + 2002 erratum)."""
    sqe = np.sqrt(nue)
    x31 = ft / (1.0 + (1.0 - 0.1 * ft) * sqe + 0.5 * (1.0 - ft) * nue / Z)

    def l31(X):
        return ((1.0 + 1.4 / (Z + 1.0)) * X - 1.9 / (Z + 1.0) * X ** 2
                + 0.3 / (Z + 1.0) * X ** 3 + 0.2 / (Z + 1.0) * X ** 4)

    L31 = l31(x31)
    xee = ft / (1.0 + 0.26 * (1.0 - ft) * sqe + 0.18 * (1.0 - 0.37 * ft) * nue / math.sqrt(Z))
    xei = ft / (1.0 + (1.0 + 0.6 * ft) * sqe + 0.85 * (1.0 - 0.37 * ft) * nue * (1.0 + Z))
    f_ee = ((0.05 + 0.62 * Z) / (Z * (1.0 + 0.44 * Z)) * (xee - xee ** 4)
            + 1.0 / (1.0 + 0.22 * Z) * (xee ** 2 - xee ** 4 - 1.2 * (xee ** 3 - xee ** 4))
            + 1.2 / (1.0 + 0.5 * Z) * xee ** 4)
    f_ei = (-(0.56 + 1.93 * Z) / (Z * (1.0 + 0.44 * Z)) * (xei - xei ** 4)
            + 4.95 / (1.0 + 2.48 * Z) * (xei ** 2 - xei ** 4 - 0.55 * (xei ** 3 - xei ** 4))
            - 1.2 / (1.0 + 0.5 * Z) * xei ** 4)
    L32 = f_ee + f_ei
    x34 = ft / (1.0 + (1.0 - 0.1 * ft) * sqe + 0.5 * (1.0 - 0.5 * ft) * nue / Z)
    L34 = l31(x34)
    sqi = np.sqrt(nui)
    a0 = -1.17 * (1.0 - ft) / (1.0 - 0.22 * ft - 0.19 * ft ** 2)
    alpha = (((a0 + 0.25 * (1.0 - ft ** 2) * sqi) / (1.0 + 0.5 * sqi) + 0.315 * nui ** 2 * ft ** 6)
             / (1.0 + 0.15 * nui ** 2 * ft ** 6))
    return L31, L32, L34, alpha


def current_drive(config, beams, op, q0: float = 1.0) -> CurrentDrive | None:
    """NBCD + bootstrap for a solved OperatingPoint (see module docstring).
    `beams` are the input BeamSpecs in the order passed to the solver.
    Returns None without profiles (profile_averaging off or infeasible)."""
    if not getattr(op, "feasible", False) or not op.rho_profile or not op.ne_profile_m3:
        return None
    geo = config.geometry
    a, R0, kappa = geo.minor_radius, geo.major_radius, geo.elongation
    V = geo.volume()
    Zeff = max(float(config.Zeff), 1.0)
    rho = np.asarray(op.rho_profile, dtype=float)
    ne = np.maximum(np.asarray(op.ne_profile_m3, dtype=float), 1.0e10)
    te = np.maximum(np.asarray(op.Te_profile_keV, dtype=float), 0.01)
    ti = np.maximum(np.asarray(op.Ti_profile_keV, dtype=float), 0.01)
    n_fuel = np.asarray(op.n_fuel_profile_m3, dtype=float)
    _ne1, nc1 = physics.compute_charge_neutrality(1.0, Zeff)
    n_c = (nc1 / _ne1) * ne                                  # carbon-like impurity
    mD, mT = float(config.mix_D), float(config.mix_T)
    nD, nT = mD * n_fuel, mT * n_fuel
    Zc = physics.Z_IMPURITY
    eps = np.clip(rho * a / R0, 0.0, 0.999)

    # ---------------- NBCD ----------------
    zh = (nD / 2.0 + nT / 3.0 + n_c * Zc ** 2 / _MASS_NUMBER_C) / ne        # sum n Z^2 / (n_e A)
    G = np.clip((1.55 + 0.85 / Zeff) * np.sqrt(eps) - (0.20 + 1.55 / Zeff) * eps, 0.0, 1.0)
    shield = 1.0 - (1.0 / Zeff) * (1.0 - G)
    j_fast, j_nb, I_beam, eta_beam = [], [], [], []
    vol_fac = V / (2.0 * math.pi * R0)
    nline = float(np.trapezoid(ne, rho)) if hasattr(np, "trapezoid") else float(np.trapz(ne, rho))
    beta_out = np.zeros_like(rho)
    for k, b in enumerate(beams):
        if k >= len(op.beam_birth_profiles_m3s):
            break
        src = np.asarray(op.beam_birth_profiles_m3s[k], dtype=float)
        xi0 = np.asarray(op.beam_pitch_profiles[k], dtype=float)
        Ab = physics.beam_mass_number(b.species)
        vb = float(physics.beam_velocity(b.Eb_keV, b.species))
        tau_s = np.array([physics.slowing_down_time(float(n), float(t), b.species) for n, t in zip(ne, te)])
        Ec = 14.8 * Ab * te * zh ** (2.0 / 3.0)
        uc3 = (Ec / b.Eb_keV) ** 1.5
        zbar = (nD * Ab / 2.0 + nT * Ab / 3.0 + n_c * Zc ** 2 * Ab / _MASS_NUMBER_C) / ne
        beta = Zeff / np.maximum(zbar, 1e-6)
        if k == 0:
            beta_out = beta
        jf = physics.E_CHARGE * src * tau_s * vb * xi0 * _pitch_integral(uc3, beta)
        jn = jf * shield
        I_k = vol_fac * physics.profile_volume_average(jn, rho)
        j_fast.append(jf)
        j_nb.append(jn)
        I_beam.append(I_k)
        eta_beam.append(I_k * R0 * nline * 1e-20 / b.P_NB_W if b.P_NB_W > 0.0 else 0.0)
    jnb_tot = np.sum(j_nb, axis=0) if j_nb else np.zeros_like(rho)
    I_nb = float(sum(I_beam))
    P_inj = sum(b.P_NB_W for b in beams)
    eta = I_nb * R0 * nline * 1e-20 / P_inj if P_inj > 0.0 else 0.0

    # ---------------- bootstrap (Sauter) ----------------
    Ip = config.Ip_MA * 1e6
    # Assumed current profile j ~ (1-rho^2)^nu. Its peaking comes from the
    # CYLINDRICAL edge q (consistent with B_p = mu0 I(rho)/L_p): nu = q_cyl/q0 - 1.
    # The low-aspect-ratio toroidal enhancement of q (q95 / q_cyl, large on an
    # ST) is not a current-peaking effect, so it is applied to q(rho) only --
    # blended as 1 + (F - 1) rho^2 so that q(0) = q0 and q(1) = q95.
    q_cyl = physics.safety_factor_cyl_edge(config.Ip_MA, config.Bt0, R0, a, kappa)
    q_edge = physics.safety_factor_cyl_edge_arbitrary_A(
        config.Ip_MA, config.Bt0, R0, a, kappa, geo.triangularity)
    q0 = max(float(q0), 0.05)
    nu = max(q_cyl / q0 - 1.0, 0.0)
    F = q_edge / q_cyl
    Lp_fac = 2.0 * math.pi * a * math.sqrt((1.0 + kappa ** 2) / 2.0)
    keV = 1.0e3 * physics.E_CHARGE

    # Analytic profile shapes of the solve, X(rho) = X(0) (1-rho^2)^(2p):
    # the gradients of the low-peaking shapes (2p < 1) diverge at rho = 1,
    # so a finite-difference gradient on the solver's grid would make the
    # edge bootstrap grid-dependent. Gradients are taken analytically and
    # the current integrated on a fine grid clustered at the edge (the
    # integrand stays finite there: nu* -> infinity as T -> 0).
    p_n = max(config.density_peaking, 0.0)
    p_te = max(config.temperature_peaking, 0.0)
    p_ti = p_te if config.temperature_peaking_i < 0.0 else max(config.temperature_peaking_i, 0.0)
    ne0_, te0_, ti0_ = ne[0], te[0], ti[0]
    nsum_per_ne = 1.0 / _ne1
    nc_per_ne = nc1 / _ne1
    nb_c = np.maximum(nsum_per_ne * ne - n_fuel, 0.0)      # fast ions on the solver grid
    dnb_c = np.gradient(nb_c, rho)

    def shape(r, pk):
        b = np.maximum(1.0 - r ** 2, 1e-12)
        return b ** (2.0 * pk), -4.0 * pk * r * b ** (2.0 * pk - 1.0)

    def bootstrap(r):
        sn, dsn = shape(r, p_n)
        ste, dste = shape(r, p_te)
        sti, dsti = shape(r, p_ti)
        ne_r, dne = ne0_ * sn, ne0_ * dsn
        te_r, dte = np.maximum(te0_ * ste, 0.01), te0_ * dste
        ti_r, dti = np.maximum(ti0_ * sti, 0.01), ti0_ * dsti
        nb_r, dnb = np.interp(r, rho, nb_c), np.interp(r, rho, dnb_c)
        n_main = np.maximum(nsum_per_ne * ne_r - nb_r, 1.0e10)     # thermal D + T
        ni = n_main + nc_per_ne * ne_r                             # all thermal ions
        dni = (nsum_per_ne + nc_per_ne) * dne - dnb
        rr = np.maximum(r, 1e-4)
        frac = 1.0 - np.maximum(1.0 - rr ** 2, 0.0) ** (nu + 1.0)
        q = q_cyl * rr ** 2 / np.maximum(frac, 1e-12) * (1.0 + (F - 1.0) * rr ** 2)
        Bp = MU0 * Ip * frac / (Lp_fac * rr)
        eps_r = np.clip(rr * a / R0, 1e-4, 0.999)
        ft = np.array([physics.trapped_particle_fraction(float(e)) for e in eps_r])
        te_ev, ti_ev = te_r * 1e3, ti_r * 1e3
        ne_s = np.maximum(ne_r, 1.0e10)
        lnle = 31.3 - np.log(np.sqrt(ne_s) / te_ev)
        lnli = 30.0 - np.log(np.sqrt(n_main) / ti_ev ** 1.5)
        nue = 6.921e-18 * q * R0 * ne_s * Zeff * lnle / (te_ev ** 2 * eps_r ** 1.5)
        nui = 4.90e-18 * q * R0 * n_main * lnli / (ti_ev ** 2 * eps_r ** 1.5)
        L31, L32, L34, alpha = _sauter_coefficients(ft, nue, nui, Zeff)
        dp = (dne * te_r + ne_r * dte + dni * ti_r + ni * dti) * keV / a
        j = -(L31 * dp + L32 * ne_r * dte * keV / a + L34 * alpha * ni * dti * keV / a) / np.maximum(Bp, 1e-9)
        j = np.where(np.isfinite(j) & (r < 1.0), j, 0.0)
        return j, q, nue

    r_f = np.sin(np.linspace(0.0, 0.5 * math.pi, 4001))       # clustered at rho = 1
    jbs_f, _q, _n = bootstrap(r_f)
    I_bs = vol_fac * 2.0 * float(np.sum(0.5 * (r_f[1:] * jbs_f[1:] + r_f[:-1] * jbs_f[:-1]) * np.diff(r_f)))
    jbs, q, nue = bootstrap(rho)
    jbs[0] = jbs[1] if len(jbs) > 1 else 0.0
    q[0] = q_cyl / (nu + 1.0)

    return CurrentDrive(
        rho=rho.tolist(),
        j_nb_per_beam=[x.tolist() for x in j_nb], j_fast_per_beam=[x.tolist() for x in j_fast],
        j_nb=jnb_tot.tolist(), j_bs=jbs.tolist(), shielding=shield.tolist(), q=q.tolist(),
        nu_star_e=nue.tolist(), I_nb_per_beam_A=I_beam, I_nb_A=I_nb, I_bs_A=I_bs, I_p_A=Ip,
        I_ind_A=Ip - I_nb - I_bs, f_ni=(I_nb + I_bs) / Ip, f_bs=I_bs / Ip,
        eta_cd=eta, eta_cd_per_beam=eta_beam, n_line_m3=nline, q0=q0, q_edge=q_edge, q_cyl=q_cyl, nu_j=nu,
        beta_scatter=beta_out.tolist(),
    )
