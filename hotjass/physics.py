"""Physics primitives for the HotJass beam-plasma model.

Most of this module is ported near-verbatim from the predecessor project
(Jassby/Code, src/beam_plasma_optimizer.py) -- the classical Coulomb
slowing-down physics, Bosch-Hale DT reactivity, and charge-neutrality
algebra are unchanged by the switch from a max-P_f search to a direct
fixed-input solve (see docs/model.md). Only single-energy-component beams
are supported here (the old project's multi-component beam-spectrum
machinery is dropped -- not needed for this model's scope).

New in this module (not present in the predecessor): electron_heating_fraction()
and ion_heating_fraction(), the L_e/L_i beam-power partition -- see
docs/model.md for the derivation.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

MU0 = 4.0 * math.pi * 1.0e-7
E_CHARGE = 1.602176634e-19
E_FUSION_J = 17.6e6 * E_CHARGE
M_D = 3.344494e-27  # neutral deuterium ATOM mass (kg) -- matches the beam as a neutral before ionization
M_T = 5.008268e-27  # neutral tritium ATOM mass (kg), same atomic-mass convention as M_D
M_E = 9.10938356e-31
M_PROTON = 1.67262192369e-27
LOG_LAMBDA = 17.0
Z_IMPURITY = 6.0
E_ALPHA_MEV = 3.5
E_FUSION_MEV = 17.6
F_ALPHA = 1.0


@dataclass
class TokamakGeometry:
    major_radius: float
    minor_radius: float
    elongation: float = 1.0
    triangularity: float = 0.0

    def volume(self) -> float:
        return 2.0 * math.pi**2 * self.major_radius * self.minor_radius**2 * self.elongation

    def aspect_ratio(self) -> float:
        return self.major_radius / max(self.minor_radius, 1e-6)


def _trapezoidal_integral(x: np.ndarray, y: np.ndarray) -> float:
    dx = np.diff(x)
    return float(np.sum(0.5 * (y[:-1] + y[1:]) * dx))


def bosch_hale_dt_reactivity(temperature_keV: float | np.ndarray) -> float | np.ndarray:
    """Bosch-Hale (1992, Nucl. Fusion 32, 611) DT thermal reactivity fit.

    Valid for 0.2-100 keV. Returns <sigma v> in m^3/s. This is a function of
    ION temperature (relative D-T velocity distribution) -- callers must
    pass T_i, not T_e.
    """
    t = np.asarray(temperature_keV, dtype=float)
    t_safe = np.maximum(t, 1e-6)
    B_G = 34.3827
    MC2 = 1124656.0
    C1, C2, C3 = 1.17302e-9, 1.51361e-2, 7.51886e-2
    C4, C5, C6, C7 = 4.60643e-3, 1.35000e-2, -1.06750e-4, 1.36600e-5

    theta = t_safe / (1.0 - (t_safe * (C2 + t_safe * (C4 + t_safe * C6))) / (1.0 + t_safe * (C3 + t_safe * (C5 + t_safe * C7))))
    xi = (B_G**2 / (4.0 * theta)) ** (1.0 / 3.0)
    sigma_v_cm3_s = C1 * theta * np.sqrt(xi / (MC2 * t_safe**3)) * np.exp(-3.0 * xi)
    sigma_v_m3_s = np.where(t > 0.0, sigma_v_cm3_s * 1e-6, 0.0)
    return sigma_v_m3_s if np.ndim(temperature_keV) else float(sigma_v_m3_s)


def bosch_hale_dt_cross_section(energy_keV: float | np.ndarray) -> float | np.ndarray:
    """Bosch-Hale (1992, Nucl. Fusion 32, 611) DT cross-section fit.

    Valid for 0.5-550 keV. Returns sigma in m^2.
    """
    e = np.maximum(np.asarray(energy_keV, dtype=float), 0.5)
    B_G = 34.3827
    A1, A2, A3, A4, A5 = 6.927e4, 7.454e8, 2.050e6, 5.2002e4, 0.0
    B1, B2, B3, B4 = 6.38e1, -9.95e-1, 6.981e-5, 1.728e-4

    numerator = A1 + e * (A2 + e * (A3 + e * (A4 + e * A5)))
    denominator = 1.0 + e * (B1 + e * (B2 + e * (B3 + e * B4)))
    s_of_e = numerator / denominator
    sigma_mb = s_of_e / (e * np.exp(B_G / np.sqrt(e)))
    sigma_m2 = sigma_mb * 1e-31
    return sigma_m2 if np.ndim(energy_keV) else float(sigma_m2)


def beam_velocity(energy_keV: float | np.ndarray, species: str = "D") -> float | np.ndarray:
    energy = np.asarray(energy_keV, dtype=float)
    energy_j = energy * 1e3 * E_CHARGE
    return np.sqrt(2.0 * energy_j / beam_mass_kg(species))


def beam_target_reactivity_spectrum(energies_keV: np.ndarray, species: str = "D") -> np.ndarray:
    """sigma(E)*v(E) for a beam ion of the given species hitting a stationary
    target. v(E) uses the ACTUAL beam mass (the real relative velocity, since
    the target is at rest). sigma(E), however, is Bosch-Hale's D-T fit
    (bosch_hale_dt_cross_section, Table VII "T(d,n)4He"), which by that
    table's own convention takes E = the DEUTERON's lab kinetic energy with
    the TRITON at rest -- i.e. its argument is not just "beam energy" but
    specifically a deuteron-frame energy. For a D beam this is already
    exactly what energies_keV means, so nothing changes (species="D" is an
    identity conversion below). For a T beam, the physical cross section
    depends on the same relative velocity v_rel = v_T, but expressing that
    v_rel as "the deuteron energy that would produce it" requires converting
    via E_D_equiv = (1/2)*M_D*v_rel^2 = E_T_lab*(M_D/M_T) -- i.e. scale the
    triton's own lab energy down by M_D/M_T before handing it to the D-frame
    fit. Getting this wrong (e.g. plugging E_T_lab directly into the D-frame
    fit) would silently use the wrong v_rel <-> cross-section mapping for a
    T beam.
    """
    mass_kg = beam_mass_kg(species)
    deuteron_equivalent_energy_keV = np.asarray(energies_keV, dtype=float) * (M_D / mass_kg)
    return bosch_hale_dt_cross_section(deuteron_equivalent_energy_keV) * beam_velocity(energies_keV, species)


def compute_charge_neutrality(n_sum: float, Zeff: float) -> tuple[float, float]:
    """(n_electron, n_impurity) implied by an ion sum density n_sum=n_D+n_T+n_b,
    at charge Z_eff, assuming a single carbon-like (Z_IMPURITY=6) impurity
    species. Pure density algebra -- no temperature dependence.
    """
    if Zeff <= 1.0 or Z_IMPURITY <= Zeff:
        return n_sum, 0.0
    n_impurity = ((Zeff - 1.0) * n_sum) / (Z_IMPURITY * (Z_IMPURITY - Zeff))
    n_electron = n_sum + Z_IMPURITY * n_impurity
    return n_electron, n_impurity


def critical_energy_keV(Te_keV: float, species: str = "D") -> float:
    """E_c = (m_beam/m_e)^(1/3) * T_e -- where electron and bulk-ion drag on a
    fast ion are equal (PPPL-1280 Eq. 2.5 simplified single-species form,
    generalized here from its original D-only (m_D/m_e)^(1/3) to use the
    ACTUAL beam species' mass -- PPPL-1280's full Eq. 2.5 has the beam mass
    A_b entering explicitly, so a heavier beam ion having a higher E_c is the
    right qualitative direction; this project's reduced single-species form
    only ever used m_D because only a D beam existed until now). Depends on
    T_e and the fast ion's own species -- NOT on the background D:T mix
    (that part of the simplification is unchanged).
    """
    return (beam_mass_kg(species) / M_E) ** (1.0 / 3.0) * Te_keV


def slowing_down_time(ne0: float, Te_keV: float, species: str = "D") -> float:
    """Raw Spitzer (electron-drag) slowing-down time tau_se (PPPL-1280 Eq. 2.9
    normalization) -- proportional to the beam ion's mass number A_b. NOT the
    full thermalization time -- use thermalization_time() for that.
    """
    ne_cm3 = max(ne0 * 1e-6, 1e3)
    Te_eV = max(Te_keV * 1e3, 1.0)
    return 6.27e8 * beam_mass_number(species) * Te_eV**1.5 / (ne_cm3 * LOG_LAMBDA)


def thermalization_time(ne0: float, Te_keV: float, energy_keV: float, species: str = "D") -> float:
    """Full fast-ion thermalization time (PPPL-1280 Eq. 2.9):
        tau_s = (tau_se/3) * ln[1 + (E/Ec)^1.5]
    """
    tau_se = slowing_down_time(ne0, Te_keV, species)
    Ec = critical_energy_keV(Te_keV, species)
    ratio_cubed = (energy_keV / max(Ec, 1e-6)) ** 1.5
    return (tau_se / 3.0) * math.log(1.0 + ratio_cubed)


def equipartition_power_w(
    n_i_m3: float, Te_keV: float, Ti_keV: float, ne0_m3: float, species: str, volume_m3: float,
) -> float:
    """Electron-ion Coulomb equipartition power [W] flowing from electrons
    to a single THERMAL ion species (density n_i_m3 -- pass nD0 or nT0, not
    n_thermal) -- positive when Te>Ti (electrons heating ions), negative
    when Ti>Te (ions heating electrons back). NRL Plasma Formulary
    "Thermal Equilibration" rate (2023 ed., "Collisions and Transport" --
    web-verified 2026-08-24 directly against the primary source PDF, not a
    secondary citation):

        nu_bar_ei = 1.8e-19 * (m_i*m_e)^0.5 * Z_i^2*Z_e^2 * n_e * lnLambda
                    / (m_i*Te + m_e*Ti)^1.5   [sec^-1]

    -- masses in GRAMS, temperatures in eV, n_e in cm^-3 (NRL practical-unit
    convention, same style already used for slowing_down_time/PPPL-1280
    here); dT_i/dt = nu_bar_ei*(Te-Ti). Converted to a power for one ion
    species via P_ei = (3/2)*n_i*(Te-Ti)*(1e3*e)*V*nu_bar_ei -- using the
    ION-perspective rate (this species' own mass) together with its own
    density, so D and T are two independent equilibration channels (this
    model already tracks nD0/nT0 separately -- summing two mass-correct
    channels is more faithful than lumping D+T into one "effective ion
    mass"). Z=1 for both D and T (hydrogenic); impurities are NOT included
    (this model never gives them a thermal population of their own -- see
    compute_charge_neutrality) -- so this only covers D/T <-> electron
    exchange, not the (smaller) impurity channel.
    """
    if n_i_m3 <= 0.0:
        return 0.0
    m_i_g = beam_mass_kg(species) * 1e3
    m_e_g = M_E * 1e3
    ne_cm3 = ne0_m3 * 1e-6
    Te_eV = max(Te_keV, 1e-6) * 1e3
    Ti_eV = max(Ti_keV, 1e-9) * 1e3
    Z_i = 1.0
    denom = (m_i_g * Te_eV + m_e_g * Ti_eV) ** 1.5
    nu_bar_ei = 1.8e-19 * (m_i_g * m_e_g) ** 0.5 * Z_i**2 * ne_cm3 * LOG_LAMBDA / denom
    return 1.5 * n_i_m3 * (Te_keV - Ti_keV) * 1e3 * E_CHARGE * volume_m3 * nu_bar_ei


def slowing_down_distribution(
    Te_keV: float, nb0: float, Eb_keV: float, energies_keV: np.ndarray, species: str = "D",
) -> np.ndarray:
    """Steady-state slowing-down distribution f(E) (m^-3 keV^-1) for a
    single-energy-component beam injected at Eb_keV:
        f(E) ~ sqrt(E) / (E^1.5 + Ec^1.5),  E <= Eb
    normalized so integrating f(E) over (0,Eb] recovers nb0.
    """
    energies = np.maximum(np.asarray(energies_keV, dtype=float), 1e-3)
    Ec = critical_energy_keV(Te_keV, species)
    kernel = np.sqrt(energies) / (energies**1.5 + Ec**1.5)
    kernel = np.where(energies <= Eb_keV, kernel, 0.0)
    norm_grid = np.linspace(1.0e-3, Eb_keV, 501)
    norm_kernel = np.sqrt(norm_grid) / (norm_grid**1.5 + Ec**1.5)
    integral = _trapezoidal_integral(norm_grid, norm_kernel)
    if integral <= 0.0:
        return np.zeros_like(energies)
    return nb0 * kernel / integral


def ion_heating_fraction(x: float) -> float:
    """L_i(x), x=Eb/Ec(Te): fraction of injected beam power delivered to
    thermal IONS as a fast ion slows from Eb down to zero, under the same
    slowing-down kernel as slowing_down_distribution().

        L_i(x) = (1/x) * integral_0^x du / (1 + u^1.5)

    Derivation (docs/model.md): in steady-state one-directional energy-space
    slowing with no losses before E=0, the particle flux in energy space
    f(E)*|dE/dt| is constant across (0,Eb) -- so power dissipation is
    *uniform* in energy space, and the ion-fraction integral is just the
    local ion-drag fraction Ec^1.5/(E^1.5+Ec^1.5) averaged uniformly over
    E in (0,Eb). Verified against Jassby (1975) Eq. 15's W_bar_b formula.
    Independent of density and of nb0 -- depends only on x=Eb/Ec(Te).
    Limits: L_i->1 as x->0 (cold/ion-dominated), L_i->0 as x->infinity
    (fast-beam/electron-dominated).
    """
    if x <= 0.0:
        return 1.0
    u = np.linspace(0.0, x, 2001)
    integrand = 1.0 / (1.0 + u**1.5)
    return float(_trapezoidal_integral(u, integrand) / x)


def electron_heating_fraction(x: float) -> float:
    """L_e(x) = 1 - L_i(x). See ion_heating_fraction()."""
    return 1.0 - ion_heating_fraction(x)


def compute_pressure(
    ne0: float, n_thermal: float, n_fast: float, Te_keV: float, Ti_keV: float,
    fast_ion_mean_pitch2: float, average_fast_energy_keV_value: float,
) -> float:
    """Total plasma pressure p_e + p_i + p_fast [Pa]. p_fast uses the
    anisotropic pitch correction p_fast=(1-<zeta^2>)*u_fast (isotropic
    default <zeta^2>=1/3 recovers p_fast=(2/3)*u_fast).
    """
    p_e = ne0 * Te_keV * 1e3 * E_CHARGE
    p_i = n_thermal * Ti_keV * 1e3 * E_CHARGE
    u_fast = n_fast * average_fast_energy_keV_value * 1e3 * E_CHARGE
    p_fast = (1.0 - fast_ion_mean_pitch2) * u_fast
    return p_e + p_i + p_fast


def compute_beta_t(pressure_pa: float, Bt0: float) -> float:
    magnetic_pressure = Bt0**2 / (2.0 * MU0)
    return pressure_pa / max(magnetic_pressure, 1e-12)


def thermal_energy_density(ne0: float, n_thermal: float, Te_keV: float, Ti_keV: float) -> float:
    """U_t = (3/2)*(ne*Te + n_thermal*Ti) [J/m^3]."""
    return 1.5 * (ne0 * Te_keV + n_thermal * Ti_keV) * 1e3 * E_CHARGE


def fast_ion_energy_density(nb0: float, average_fast_energy_keV_value: float) -> float:
    return nb0 * average_fast_energy_keV_value * 1e3 * E_CHARGE


def average_fast_energy_keV(Te_keV: float, Eb_keV: float, species: str = "D") -> float:
    """Mean energy of the steady-state slowing-down distribution (Jassby
    1975 Eq. 15's W_bar_b), always < Eb.
    """
    energies = np.linspace(1.0e-3, Eb_keV, 1001)
    Ec = critical_energy_keV(Te_keV, species)
    kernel = np.sqrt(energies) / (energies**1.5 + Ec**1.5)
    integral = _trapezoidal_integral(energies, kernel)
    if integral <= 0.0:
        return Eb_keV
    return float(_trapezoidal_integral(energies, energies * kernel) / integral)


def alpha_heating_power(pf_total_w: float, F_alpha: float = F_ALPHA) -> float:
    return F_alpha * (E_ALPHA_MEV / E_FUSION_MEV) * pf_total_w


def thermal_fusion_power(nD0: float, nT0: float, Ti_keV: float, volume_m3: float) -> float:
    """Thermal D-T fusion power [W]. IMPORTANT: the temperature argument is
    T_i (ion temperature -- Bosch-Hale reactivity depends on relative ion
    velocities), not T_e. The predecessor project's equivalent method took
    a parameter literally named Te0 here, which was only ever safe because
    every call site there enforced Ti=Te; this model solves Te and Ti
    independently (no e-i equipartition, see docs/model.md), so the
    distinction is load-bearing here. Flat radial profile assumed (no
    profile shaping in this model, unlike the predecessor project).
    """
    if nD0 <= 0.0 or nT0 <= 0.0:
        return 0.0
    reactivity = bosch_hale_dt_reactivity(Ti_keV)
    return nD0 * nT0 * reactivity * volume_m3 * E_FUSION_J


_BEAM_MASS_NUMBER = {"D": 2.0, "T": 3.0}
_BEAM_MASS_KG = {"D": M_D, "T": M_T}


def beam_mass_number(species: str) -> float:
    """Integer mass number A for the beam species. Simplification consistent
    with this project's other fixed-constant conventions (e.g. LOG_LAMBDA=17,
    Z_IMPURITY=6).
    """
    return _BEAM_MASS_NUMBER[species]


def beam_mass_kg(species: str) -> float:
    """Actual (neutral-atom) mass [kg] for the beam species -- M_D or M_T.
    Unlike beam_mass_number() (the unitless A used for the "per amu" stopping
    cross section), this is the real mass used everywhere a velocity,
    gyroradius, or critical-energy calculation needs actual kg.
    """
    return _BEAM_MASS_KG[species]


def beam_stopping_cross_section_m2(energy_per_amu_keV: float) -> float:
    """Effective (electron/ion impact ionization + charge-exchange) stopping
    cross section for a hydrogenic neutral beam in a hydrogenic plasma, as a
    function of beam energy per nucleon [keV/amu].

    PPPL-1280 Sect. 4.1 (this project's own beam-physics reference) defines
    exactly this quantity via the attenuation law dI_b/dx = -I_b/lambda_t
    (its Eq. 4.1) and cites Riviere (1971, Nucl. Fusion 11, 363) for its
    energy dependence -- but only as a plot (its Fig. 7a), not a closed-form
    fit. The coefficients below are a widely-reproduced practical fit
    commonly attributed to Riviere; they have NOT been independently
    verified against the original paper here. Treat as structurally-correct
    / order-of-magnitude pending that check (see docs/model.md "Possible
    extensions" -- this was an explicit, informed choice: use a standard
    literature fit now to test the model's qualitative behavior, flagged for
    later verification, rather than blocking on digitizing the original
    figure).
    """
    E = max(float(energy_per_amu_keV), 1e-6)
    sigma_cm2 = 0.6937e-14 * (1.0 - 0.155 * math.log10(E)) ** 2 / (1.0 + 0.1112e-14 * E**3.3)
    return sigma_cm2 * 1e-4  # cm^2 -> m^2


def tangential_path_length(minor_radius_m: float) -> float:
    """Chord length through the plasma for tangential injection aimed at
    tangency radius R_t = R_0 (through the magnetic axis) -- PPPL-1280
    Sect. 4.1: "the path length for tangential injection is approximately
    2a_p". Circular-cross-section approximation (no elongation correction),
    consistent with this project's other geometric simplifications. Only
    valid for R_t = R_0; a general R_t != R_0 chord would need a different
    (asymmetric) formula -- not implemented here.
    """
    return 2.0 * minor_radius_m


def captured_power_fraction(ne0: float, Eb_keV: float, species: str, path_length_m: float) -> float:
    """Fraction of injected neutral-beam power actually absorbed by the
    plasma (1 - shine-through), under Beer-Lambert attenuation of the
    neutral beam flux (PPPL-1280 Eq. 4.1: I(x)=I(0)*exp(-x/lambda_t)) along
    a straight chord of length path_length_m. Depends only on ne0 and the
    beam's own energy/species -- NOT on T_e -- so it can be evaluated before
    the T_e solve.
    """
    A = beam_mass_number(species)
    sigma_m2 = beam_stopping_cross_section_m2(Eb_keV / A)
    optical_depth = ne0 * sigma_m2 * path_length_m
    return 1.0 - math.exp(-optical_depth)


def iter98y2_confinement_time(
    Ip_MA: float, Bt_T: float, ne0_m3: float, P_heat_W: float,
    R0_m: float, minor_radius_m: float, elongation: float, M_eff_amu: float,
) -> float:
    """IPB98(y,2) ELMy H-mode energy confinement scaling (ITER Physics Basis
    1999; standard reference form, web-verified against multiple sources
    2026-08-23 -- not independently re-derived from the original ITER
    Physics Basis table here):

        tau_E [s] = 0.0562 * Ip[MA]^0.93 * Bt[T]^0.15 * n19^0.41
                    * P[MW]^-0.69 * R[m]^1.97 * kappa_a^0.78 * eps^0.58
                    * M[amu]^0.19

    where n19 = ne0_m3/1e19 and eps = minor_radius/R0.

    KNOWN BIAS FOR THIS PROJECT'S GEOMETRY: this scaling was fit to
    conventional-aspect-ratio devices (A~2.5-4); it is documented in the ST
    literature to mis-predict confinement at low aspect ratio (A<2 here,
    ~1.86) -- real spherical-tokamak fits (e.g. Globus-M2's
    Ip^0.43*Bt^1.19, vs this formula's Ip^0.93*Bt^0.15 -- Kaye, Gusev et
    al., arXiv:2509.02214) show much weaker I_p and much stronger B_t
    dependence. This matters only when I_p or B_t are themselves being
    varied; for a scan at fixed I_p/B_t (as in this project's n_e/P_NB
    scans), that bias is just a fixed multiplicative offset on tau_E,
    not a shape distortion -- see docs/model.md "Possible extensions" for
    the full reasoning behind using this formula anyway.
    """
    n19 = ne0_m3 / 1e19
    P_MW = P_heat_W / 1e6
    eps = minor_radius_m / R0_m
    return (
        0.0562
        * Ip_MA**0.93
        * Bt_T**0.15
        * n19**0.41
        * max(P_MW, 1e-9) ** -0.69
        * R0_m**1.97
        * elongation**0.78
        * eps**0.58
        * M_eff_amu**0.19
    )


def kaye_nstx_confinement_time(Ip_A: float, Bt_T: float, ne0_m3: float, P_heat_W: float) -> float:
    """Kaye (2006, Nucl. Fusion 46, 848) NSTX H-mode energy confinement
    scaling -- OLSR "Case 1" fit (its Table 1, all 85 H-mode points,
    RMSE=0.145), a genuine low-aspect-ratio (spherical tokamak) fit, unlike
    IPB98(y,2) above:

        tau_E,th [s] = 4.69e-9 * Ip[A]^0.57 * Bt[T]^1.08 * ne[m^-3]^0.44
                       * P[W]^-0.73

    IMPORTANT UNITS (differ from iter98y2_confinement_time() above): Ip in
    AMPERES (not MA), ne in m^-3 (not 1e19 m^-3), P in WATTS (not MW) --
    exactly as stated in the paper's Table 1 caption ("units for Ip, BT, ne,
    PL,th and tau_E,th are A, T, m^-3, W and sec respectively"). Verified
    directly against the primary source (references/Kaye_2006_NF_46_848_NSTX.pdf,
    Table 1, Case 1), including the leading coefficient -- unlike the
    IPB98(y,2) formula above, this one was read from the actual paper, not a
    secondary citation, so there is no "not yet verified" caveat on the
    numbers themselves. The aspect-ratio-appropriateness of the underlying
    NSTX dataset relative to this project's own geometry (A~1.86 vs NSTX's
    A~1.3-1.5) is a separate, still-open question -- see docs/model.md.

    Weaker I_p exponent (0.57 vs IPB98(y,2)'s 0.93) and much stronger B_t
    exponent (1.08 vs 0.15) confirm the qualitative ST trend already
    discussed there, now with an equally strong P_heat/n_e dependence
    (-0.73/0.44, similar order to IPB98(y,2)'s -0.69/0.41).
    """
    return 4.69e-9 * Ip_A**0.57 * Bt_T**1.08 * ne0_m3**0.44 * max(P_heat_W, 1e-3) ** -0.73


def kaye_nstx_lmode_confinement_time(Ip_A: float, Bt_T: float, ne0_m3: float, P_loss_W: float) -> float:
    """Kaye (2006, Nucl. Fusion 46, 848) NSTX L-mode GLOBAL energy
    confinement scaling -- an OLSR fit to NSTX L-mode discharges (its Eq. 5
    and Fig. 9), the L-mode counterpart to kaye_nstx_confinement_time()'s
    H-mode fit above, same primary source, same low-aspect-ratio dataset:

        tau_E [s] = 4.73e-4 * Ip[A]^1.01 * Bt[T]^0.70 * ne[m^-3]^-0.07
                    * P_loss[W]^-0.37

    UNITS: same convention as kaye_nstx_confinement_time() -- Ip in AMPERES,
    ne in m^-3, P in WATTS. P_loss = P_NBI + P_OH - dW/dt - P_shine-thru
    (paper's own definition); in this project's steady-state, no-Ohmic-
    heating model that reduces to P_useful_w (Step 0's captured, post-
    orbit-loss, post-CX-loss beam power) exactly.

    THE COEFFICIENT IS NOT IN THE PAPER'S OWN TEXT: Eq. 5 is stated there
    only as a proportionality ("tau_E ~ Ip^1.01 Bt^0.70 ne^0.07 Ploss^-0.37",
    no leading constant -- the paper says L-mode data quality precluded even
    a thermal-energy fit, so this global-tau_E fit was only shown
    graphically). The 4.73e-4 coefficient used here was recovered from
    Figure 9's own axis label (the regression tool's literal fit-equation
    printout under the plot, not retyped prose) -- web-verified 2026-08-29
    directly against the primary source PDF page image.

    A GENUINE PRIMARY-SOURCE DISCREPANCY, resolved by testing rather than by
    preference: the body text's Eq. 5 states the density exponent as
    ne^+0.07, but Figure 9's own fit label -- which is also the only place
    the 4.73e-4 coefficient appears at all -- states ne^-0.07. These are NOT
    independently swappable: the coefficient and every exponent came out of
    one joint regression, so pairing 4.73e-4 with the text's +0.07 instead
    of the figure's own -0.07 doesn't give "weakly increasing instead of
    weakly decreasing with density" -- it gives tau_E~=29s at a
    representative NSTX point (Ip=1 MA, Bt=0.5 T, ne=5e19 m^-3, P_loss=4 MW)
    versus ~0.05s with -0.07 (Fig. 9's own 0-90 ms axis range) -- ~600x off,
    not a sane alternative. So -0.07 (Figure 9's own value, the only
    self-consistent pairing with 4.73e-4) is used here; the density
    dependence is weak either way (|exponent|=0.07), so over this model's
    actual n_e scan range (0.5-15e19, a 30x span) this is only a +/-27%
    effect, not a large physics change.
    """
    return 4.73e-4 * Ip_A**1.01 * Bt_T**0.70 * ne0_m3**-0.07 * max(P_loss_W, 1e-3) ** -0.37


def safety_factor_cyl_edge(Ip_MA: float, Bt_T: float, R0_m: float, a_m: float, elongation: float) -> float:
    """Approximate edge (cylindrical + elongation-corrected) safety factor:

        q_a ~= 5 * a^2 * Bt * (1+kappa^2)/2 / (R0 * Ip[MA])

    Standard large-aspect-ratio circular-cylinder q_cyl = 2*pi*a^2*Bt /
    (mu0*R0*Ip), rearranged into the common practical-units form (a, R0 in
    m, Bt in T, Ip in MA) with the usual (1+kappa^2)/2 elongation correction
    (e.g. Wesson "Tokamaks", Freidberg "Plasma Physics and Fusion Energy")
    -- an approximation, not a real MHD equilibrium q, same level of rigor
    as this project's other geometry-only estimates (e.g.
    tangential_path_length). Used only inside neoclassical_tau_Ei_s() below,
    since no other part of this model has needed q before now.
    """
    return 5.0 * a_m**2 * Bt_T * (1.0 + elongation**2) / 2.0 / (R0_m * max(Ip_MA, 1e-9))


def ion_ion_collision_frequency_hz(Ti_keV: float, ni_m3: float, species: str = "D") -> float:
    """Ion-ion collision frequency [s^-1], NRL Plasma Formulary "Collision
    Rates" (Z_i=1 hydrogenic ions, same Coulomb logarithm convention as
    equipartition_power_w's LOG_LAMBDA=17):

        nu_ii = 4.80e-8 * Z_i^4 * n_i[cm^-3] * lnLambda / (mu^0.5 * Ti[eV]^1.5)

    where mu = m_i/m_proton. Practical-unit NRL form, same style as
    equipartition_power_w's nu_bar_ei above.
    """
    if ni_m3 <= 0.0:
        return 0.0
    ni_cm3 = ni_m3 * 1e-6
    Ti_eV = max(Ti_keV, 1e-9) * 1e3
    mu = beam_mass_kg(species) / M_PROTON
    Z_i = 1.0
    return 4.80e-8 * Z_i**4 * ni_cm3 * LOG_LAMBDA / (mu**0.5 * Ti_eV**1.5)


def neoclassical_tau_Ei_s(
    Ti_keV: float, n_thermal_m3: float, Bt_T: float, Ip_MA: float,
    R0_m: float, a_m: float, elongation: float, mix_D: float, mix_T: float,
) -> float:
    """Neoclassical (banana-regime) ion energy confinement time -- the
    textbook order-unity-coefficient estimate (Wesson "Tokamaks" /
    Chang-Hinton low-collisionality limit, NOT the full Chang-Hinton
    interpolation across all collisionality regimes):

        chi_i,neo ~= q^2 * rho_i^2 * nu_ii / eps^1.5
        tau_Ei,neo = a^2 / chi_i,neo = a^2 * eps^1.5 / (q^2 * rho_i^2 * nu_ii)

    where eps=a/R0, q from safety_factor_cyl_edge() above, rho_i the
    THERMAL ion Larmor radius (reuses larmor_radius_m()'s exact formula,
    just evaluated at Ti_keV instead of a beam's Eb_keV -- same
    m_i*v/(e*Bt) convention, v=sqrt(2E/m)), and nu_ii from
    ion_ion_collision_frequency_hz() above.

    SIMPLIFICATION: this model tracks D and T as two separate thermal
    populations (mix_D/mix_T), but Chang-Hinton itself is already only an
    order-unity estimate -- so rather than a full two-species neoclassical
    combination, this uses ONE effective ion mass
    (mix_D*M_D + mix_T*M_T, weighted by density fraction) and the total
    n_thermal, treating the ion channel as a single effective species for
    rho_i/nu_ii. Consistent with the base formula's own level of rigor;
    flagged here as a modeling choice, not asserted as exact.

    tau_Ei,neo depends on Ti itself (~Ti^0.5, since rho_i^2~Ti and
    nu_ii~Ti^-1.5) -- callers must treat this as part of a fixed-point/
    root-find for Ti, not a closed-form substitution (see
    solve._solve_ti_neoclassical_keV).
    """
    eps = a_m / R0_m
    q = safety_factor_cyl_edge(Ip_MA, Bt_T, R0_m, a_m, elongation)
    m_eff_kg = mix_D * M_D + mix_T * M_T
    v_thi = math.sqrt(2.0 * max(Ti_keV, 1e-9) * 1e3 * E_CHARGE / m_eff_kg)
    rho_i = m_eff_kg * v_thi / (E_CHARGE * Bt_T)
    nu_ii_D = ion_ion_collision_frequency_hz(Ti_keV, mix_D * n_thermal_m3, "D")
    nu_ii_T = ion_ion_collision_frequency_hz(Ti_keV, mix_T * n_thermal_m3, "T")
    nu_ii = nu_ii_D + nu_ii_T
    chi_i = q**2 * rho_i**2 * nu_ii / eps**1.5
    if chi_i <= 0.0:
        return float("inf")
    return a_m**2 / chi_i


def larmor_radius_m(Eb_keV: float, Bt_T: float, species: str = "D") -> float:
    """Bare fast-ion Larmor (gyro-)radius [m] at the full injection energy in
    field Bt_T: rho_Li = m_beam * v_b / (e * Bt) -- the un-shape-corrected
    factor inside passing_orbit_width()'s q*-multiplied drift-orbit width,
    exposed standalone since it also drives the separate, direction-
    independent Larmor prompt-loss criterion in first_orbit_loss_fraction()
    (an ion gyrates through +/-rho_Li of its guiding center every cyclotron
    period, regardless of which way the much-slower guiding-center drift
    orbit is shifted). Species matters twice over: v_b itself is smaller for
    a heavier ion at the same Eb_keV, but m_beam is bigger -- net effect,
    rho_Li ~ sqrt(m_beam*Eb) grows with mass at fixed energy, so a T ion is
    NOT simply "slower, same gyroradius" -- it has a larger gyroradius too.
    """
    v_b = beam_velocity(Eb_keV, species)
    return beam_mass_kg(species) * v_b / (E_CHARGE * Bt_T)


def passing_orbit_width(Eb_keV: float, Bt_T: float, Ip_MA: float, R0_m: float,
                          minor_radius_m: float, elongation: float, species: str = "D") -> float:
    """Radial guiding-center orbit width [m] for a strongly co-passing fast
    ion born at the full injection energy -- the standard large-aspect-
    ratio estimate

        Delta_r ~ q* * rho_Li,  q* = (2 pi a^2 Bt)/(mu0 R0 Ip) * (1+kappa^2)/2

    (cylindrical/shape-corrected safety factor times the toroidal-field
    gyroradius at birth energy -- textbook approximation, e.g. Wesson
    "Tokamaks"; circular-flux-surface, large-aspect-ratio assumption, same
    level of rigor as this project's other geometric simplifications like
    tangential_path_length()). Applies to R_t=R_0 tangential injection
    (either direction), where fast ions are born launched nearly purely
    parallel to B (strongly passing, not trapped) -- only the *magnitude*
    of the shift, not its sign, so the same width feeds
    first_orbit_loss_fraction()'s co_current parameter for either
    direction. species is passed straight through to larmor_radius_m() --
    see its docstring for why a heavier beam ion isn't just "slower."
    """
    q_star = (2.0 * math.pi * minor_radius_m**2 * Bt_T) / (MU0 * R0_m * Ip_MA * 1e6) * (1.0 + elongation**2) / 2.0
    return q_star * larmor_radius_m(Eb_keV, Bt_T, species)


def first_orbit_loss_fraction(
    ne0: float, Eb_keV: float, species: str, path_length_m: float,
    minor_radius_m: float, orbit_width_m: float, co_current: bool = True,
    include_larmor_loss: bool = False, larmor_radius_m_value: float = 0.0,
) -> float:
    """Fraction of CAPTURED (post shine-through) beam ions promptly lost to
    first-orbit loss -- born close enough to the LCFS that their passing-
    particle orbit (width orbit_width_m, from passing_orbit_width()) carries
    them outside it before any slowing-down or charge-exchange can occur.

    DIRECTION MATTERS: a passing ion's guiding center is displaced from its
    birth flux surface by ~orbit_width, INWARD for co-current injection
    (favorable) and OUTWARD for counter-current (unfavorable) -- confirmed
    against this project's own PPPL-1280 reference ("the stagnation axis is
    located ... outside (inside) the magnetic axis for ions injected
    parallel (antiparallel) to the plasma current", i.e. co=inward,
    counter=outward) and cross-checked against the literature. Modeled here
    as a simple mean-shift of the effective birth radius:
    rho_eff = rho -/+ orbit_width/a (co/counter); the ion is lost if
    rho_eff > 1. For co_current=True this can never trigger (rho<=1 and the
    shift is negative), so this simple shift-only model gives EXACTLY zero
    first-orbit loss for co-current -- the correct qualitative limit (real
    devices see a small but nonzero residual from finite orbit width beyond
    this mean-shift treatment, e.g. ~6% co- vs ~26-34% counter- in one
    published comparison), just not something this simplified model
    captures quantitatively. co_current=False (counter) reproduces the
    original threshold rho > 1 - orbit_width/a.

    OPTIONAL second loss channel -- LARMOR (gyro-orbit) PROMPT LOSS
    (include_larmor_loss=True, with larmor_radius_m_value =
    physics.larmor_radius_m(Eb_keV, Bt_T), the *bare* gyroradius, NOT the
    q*-multiplied orbit_width_m): even under co_current=True, where the
    guiding-center drift orbit shifts inward and this shift-only model
    gives exactly zero loss, an ion born within one gyro-diameter of the
    LCFS still gyrates -- radius rho_Li, independent of injection direction
    -- through positions outside the LCFS within a single cyclotron period,
    much faster than the guiding-center drift shift develops over a full
    poloidal transit. Modeled as a second, direction-independent cutoff:
    lost if rho > 1 - 2*rho_Li/a, OR-ed with the drift-shift criterion above
    (not summed, to avoid double-counting ions already caught by it).
    Disabled by default: existing first_orbit_loss_fraction() results
    (including the exact co_current=True zero) are unaffected unless
    explicitly turned on.

    Reuses the same Beer-Lambert deposition-density-along-the-chord already
    used for shine-through (Step 0a) -- ionization rate density
    n_e*sigma*exp(-n_e*sigma*x) -- mapped to normalized minor radius via the
    tangential (R_t=R_0) chord geometry rho(x)=|x-a|/a for x in [0,2a]
    (circular-cross-section approximation, consistent with
    tangential_path_length()), then integrates the deposited fraction born
    in the loss layer. A single 1-D numerical integral, the same cost/style
    as L_i(x) or the beam-target fusion integral already in this module.
    """
    A = beam_mass_number(species)
    sigma_m2 = beam_stopping_cross_section_m2(Eb_keV / A)
    x = np.linspace(0.0, path_length_m, 501)
    density = np.exp(-ne0 * sigma_m2 * x)
    rho = np.abs(x - minor_radius_m) / minor_radius_m
    shift = -orbit_width_m / minor_radius_m if co_current else orbit_width_m / minor_radius_m
    rho_eff = rho + shift
    lost_mask = rho_eff > 1.0
    if include_larmor_loss:
        lost_mask = lost_mask | (rho > 1.0 - 2.0 * larmor_radius_m_value / minor_radius_m)
    lost_density = np.where(lost_mask, density, 0.0)
    total = _trapezoidal_integral(x, density)
    if total <= 0.0:
        return 0.0
    lost = _trapezoidal_integral(x, lost_density)
    return float(lost / total)


def beam_target_fusion_power(
    nb0: float, n_target0: float, Te_keV: float, Eb_keV: float, volume_m3: float, species: str = "D",
) -> float:
    """Beam-target fusion power [W]: a fast beam ion (given species, D or T)
    slowing down through a stationary (T_i=0) target-ion population of the
    OTHER D-T species -- matching PPPL-1280 Fig. 20(a)'s stationary-target
    convention, an inherited simplification, so this still ignores the
    newly-nonzero T_i that this model otherwise solves for (see
    docs/model.md). Te_keV here is genuinely T_e (shapes the slowing-down
    kernel via Ec(Te)), not a stand-in for Ti. n_target0 is whichever
    density the CALLER decides is the target (nT0 for a D beam, nD0 for a T
    beam -- solve.py picks this per-beam).

    Generalized from "fast-D on thermal-T" to either beam species: species
    is threaded through to slowing_down_distribution() (species-correct
    Ec(Te)) and beam_target_reactivity_spectrum() (species-correct v_rel and
    the deuteron-equivalent-energy conversion its cross-section fit needs --
    see that function's docstring).

    D-D beam-target reactions (relevant once mix_D>0) are NOT included --
    an inherited gap from the predecessor project, live from day one here.
    """
    energies = np.linspace(1.0e-3, Eb_keV, 601)
    distribution = slowing_down_distribution(Te_keV, nb0, Eb_keV, energies, species)
    sigma_v = beam_target_reactivity_spectrum(energies, species)
    reaction_density = n_target0 * distribution * sigma_v
    reaction_rate = volume_m3 * _trapezoidal_integral(energies, reaction_density)
    return reaction_rate * E_FUSION_J
