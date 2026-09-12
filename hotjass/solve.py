"""The HotJass direct solve: given fixed beam power(s), beam energy(ies),
D:T mix, and energy confinement time tau_E, find the self-consistent
operating point (T_e, T_i, densities, fusion power) at a given electron
density n_e -- no search, no optimization. See docs/model.md for the full
derivation.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from . import physics


@dataclass(frozen=True)
class BeamSpec:
    """One neutral-beam source: injection energy and injected power.

    Phase 1 uses a single D-beam (one-element beams list everywhere below).
    Phase 2 (not yet implemented) adds a second T-beam -- the electron/ion
    balance below already sums over an arbitrary list of beams, so that
    extension needs no change to solve_operating_point()'s structure, only
    a second BeamSpec and (separately, unimplemented) a beam-beam fusion
    reaction channel.
    """

    Eb_keV: float
    P_NB_W: float
    species: str = "D"  # "D" or "T" -- which ion species this beam injects
    tangent_R_m: float | None = None
    tangent_Z_m: float = 0.0
    shine_through_model: str = "manual"
    manual_shine_through_fraction: float = 0.01
    co_current: bool = True  # per-beamline injection direction (co- vs counter-Ip);
    #   used for first-orbit loss instead of the (now-legacy) config.orbit_loss_co_current


@dataclass
class TokamakConfig:
    geometry: physics.TokamakGeometry
    Bt0: float
    Ip_MA: float
    Zeff: float = 1.0
    mix_D: float = 0.5
    mix_T: float = 0.5
    density_peaking: float = 0.0
    temperature_peaking: float = 0.0        # electron T profile exponent
    temperature_peaking_i: float = -1.0     # ion T profile exponent; < 0 -> same as temperature_peaking
    centrepost_radius_m: float = -1.0       # central-column outer radius; < 0 -> R0 - a
    profile_averaging: bool = False
    # profile_averaging=False (default): the 0-D balance treats ne0_m3 and the
    # solved T as spatially UNIFORM -- the historical behaviour, byte-identical.
    #   The reported T is then a flat-plasma effective temperature, and comparing
    #   it to a measured on-axis T0 needs a peaking factor supplied by the user.
    # profile_averaging=True: ne0_m3 is taken to be the ON-AXIS density; the
    #   energy balance and the tau_E scalings run on the volume-averaged
    #   <n_e> = ne0_m3 / (1 + 2*density_peaking) and return VOLUME-AVERAGED
    #   T_e, T_i (OperatingPoint.Te_keV / Ti_keV). The on-axis values
    #   Te0_keV = Te_keV*(1+2*temperature_peaking), Ti0_keV similarly, are
    #   reported alongside and are what the fusion / pressure integrals use as
    #   their central values (with the (1-rho^2) shape from the peaking params).
    fast_ion_mean_pitch2: float = 1.0 / 3.0
    enable_orbit_loss: bool = False
    orbit_loss_co_current: bool = True
    include_larmor_loss: bool = False
    orbit_model: str = "large_aspect"
    # First-orbit-loss orbit-width model (only used when enable_orbit_loss):
    #   "large_aspect" (default): physics.passing_orbit_width (q* rho_Li) +
    #     the direction-only mean-shift in first_orbit_loss_fraction -- gives
    #     EXACTLY zero loss for co-current. include_larmor_loss adds the
    #     optional bare-gyroradius channel. Unchanged from before.
    #   "st_meanshift": physics.first_orbit_loss_fraction_st(variant=
    #     "meanshift") -- arbitrary-A q_a, banana channel, gyro channel;
    #     co-current loss is NOT zero on a spherical tokamak.
    #   "st_pitch": physics.first_orbit_loss_fraction_st(variant="pitch") --
    #     a pitch-angle-resolved 2-D integral (co-passing / counter-passing /
    #     trapped classes). Both ST variants are offered so they can be
    #     compared; neither touches the large-aspect result.
    # First-orbit loss (physics.first_orbit_loss_fraction/passing_orbit_width):
    # a fraction of CAPTURED beam power (post shine-through) is lost
    # promptly -- born close enough to the LCFS that a fast ion's passing-
    # orbit is shifted back outside before any slowing-down or CX can occur
    # -- for R_t=R_0 tangential injection (this project's assumed geometry).
    # DIRECTION MATTERS: co-current shifts the orbit INWARD (favorable --
    # confirmed against this project's own PPPL-1280 reference), so with
    # orbit_loss_co_current=True (this project's actual stated injection
    # direction, see e.g. Figures 8-13's titles) this simple shift-only
    # model gives EXACTLY ZERO first-orbit loss -- turning enable_orbit_loss
    # on has no effect unless orbit_loss_co_current is also set False
    # (counter-current, the unfavorable/outward-shifted case). Unlike
    # cx_loss_fraction, f_orbit_loss itself is NOT a free knob: it's
    # computed from (ne0, Eb, geometry, direction) with no adjustable
    # magnitude, so there's no default *value* to set beyond these two
    # on/off switches. Defaults (False, True) leave existing results
    # unaffected unless explicitly enabled -- same reasoning as
    # cx_loss_fraction's default. See docs/model.md.
    # include_larmor_loss (physics.larmor_radius_m, folded into
    # first_orbit_loss_fraction): a SECOND, direction-independent loss
    # channel -- an ion born within one gyro-diameter (2*rho_Li) of the LCFS
    # is lost to its own gyro-motion regardless of co/counter, since
    # gyration happens far faster than the guiding-center drift shift can
    # develop. Only takes effect when enable_orbit_loss is also True.
    # Default False: leaves the existing (co_current=True -> exactly zero)
    # first-orbit-loss result unaffected unless explicitly enabled.
    cx_loss_fraction: float = 0.0
    # Charge-exchange model selector (docs/CX_model.tex Sec.4):
    #   "manual_fraction" (default): the flat cx_loss_fraction knob above,
    #     UNCHANGED -- byte-identical to every existing result.
    #   "manual_n0": user supplies n0/ne directly (cx_n0_over_ne); the real
    #     Sec.2 survival integral (physics.cx_loss_fractions) replaces the
    #     flat knob. n0 is spatially uniform (0-D), applied via
    #     solve_operating_point's cx_loss_fraction_override -- see that
    #     function's docstring for why it stays a passthrough rather than
    #     computing this itself (Te-dependence -> a fixed point, closed by
    #     an outer Picard loop in hotjass_core, same pattern as alpha
    #     heating below).
    #   "penetration": n0 is instead derived from cx_n0_lcfs_over_ne via the
    #     Sec.3.2 edge-penetration profile (physics.neutral_penetration_profile),
    #     volume-averaged, then fed through the same Sec.2 integral.
    # Only engaged when cx_model != "manual_fraction"; default leaves this
    # whole feature inert.
    cx_model: str = "manual_fraction"
    cx_n0_over_ne: float = 1.0e-5
    cx_n0_lcfs_over_ne: float = 0.02
    # Bulk toroidal rotation (physics.nbi_torque_Nm / toroidal_rotation_velocity_ms):
    #   "off" (default): v_phi=0 -- byte-identical to every existing result
    #     (beam_target_fusion_power/beam_target_dd_fusion_power's v_phi_m_s=0
    #     default reproduces the stationary-target result exactly).
    #   "manual": v_phi = manual_v_phi_m_s directly (positive = co-current).
    #   "momentum_balance": v_phi from a 0-D angular-momentum balance -- net
    #     NBI torque (summed over beams, signed by each beam's own
    #     co_current) against a momentum-confinement loss L/tau_phi, with
    #     tau_phi = tau_phi_over_tauEi * tau_Ei_s (no validated tau_phi
    #     scaling exists, so this ties it to whatever ion energy confinement
    #     time -- fixed or physics-scaling-based -- is already in use).
    # v_phi then corrects the beam-target reactivity's relative velocity for
    # each beam (its own co/counter direction against the bulk rotation);
    # first-orbit loss, shine-through and CX loss are NOT touched -- v_phi is
    # negligible next to the beam velocity for those channels (see the
    # session that added this feature for the estimate).
    rotation_model: str = "off"  # "off" | "manual" | "momentum_balance"
    manual_v_phi_m_s: float = 0.0
    tau_phi_over_tauEi: float = 1.0
    # Beam-beam fusion (physics.beam_beam_fusion_power): a reduced,
    # monoenergetic-population estimate for the ONE pairwise beam-beam
    # reaction between two distinct NBI sources (D-T or D-D; see that
    # function's docstring for why counter-injecting the second beam
    # matters here -- it maximizes the beam-beam relative velocity). Off by
    # default: leaves every existing pf_dt_w/pf_dd_w/neutron_rate_s result
    # unaffected unless explicitly enabled. Only the first two useful beams
    # are paired (this project's tested scope is exactly two NBI sources);
    # a third beam, if ever added, would need generalizing to all pairs.
    enable_beam_beam: bool = False
    enable_equipartition: bool = False
    # Fusion alpha self-heating. When enable_alpha_heating is True the caller
    # (hotjass_core) closes a fixed-point P_alpha = f_alpha*(3.5/17.6)*P_fusion,
    # feeding P_alpha back into the T_e/T_i balance via solve_operating_point's
    # extra_e_source_w/extra_i_source_w (split electron/ion by the alpha
    # slowing-down fraction, physics.alpha_electron_heating_fraction). f_alpha is
    # the fraction of alpha power confined and thermalised (1.0 = fully confined).
    # solve_operating_point itself never iterates -- it just accepts the two extra
    # source terms. Default False -> byte-identical to the alpha-free result.
    enable_alpha_heating: bool = False
    f_alpha: float = 1.0
    tau_Ee_mode: str = "fixed"
    tau_Ei_mode: str = "fixed"
    # Physics-based confinement-time scalings (docs/model.md), replacing the
    # externally-fixed tau_E_s/tau_Ei_s arguments to solve_operating_point()
    # with a computed value -- OPT IN, per-channel, independent of each
    # other:
    #   tau_Ee_mode="kaye_nstx_lmode": electron balance (Step 1) uses
    #     physics.kaye_nstx_lmode_confinement_time(Ip_A, Bt0, ne0, P_useful_w)
    #     instead of the passed-in tau_E_s. Kaye (2006, NF 46, 848) NSTX
    #     L-mode fit (its Fig. 9) -- an ST-appropriate scaling, unlike
    #     IPB98(y,2). No fixed-point issue: P_useful_w (Step 0) is known
    #     before Te is solved, so tau_Ee itself doesn't depend on Te.
    #   tau_Ei_mode="neoclassical": ion balance (Step 5) uses
    #     physics.neoclassical_tau_Ei_s(Ti, n_thermal, Bt0, Ip_MA, geometry,
    #     mix_D, mix_T) instead of the passed-in tau_Ei_s. UNLIKE the
    #     electron case, this DOES depend on the quantity being solved for
    #     (tau_Ei,neo ~ Ti^0.5) -- Ti's closed form becomes a genuine 1-D
    #     fixed-point root-find (_solve_ti_neoclassical_keV) rather than a
    #     one-line substitution.
    #   tau_Ei_mode="kaye_nstx_lmode": ion balance uses the SAME
    #     kaye_nstx_lmode_confinement_time() value as tau_Ee_mode above
    #     (computed once and shared -- it doesn't depend on Te or Ti, so
    #     there's no reason to recompute it per channel). Setting BOTH
    #     tau_Ee_mode and tau_Ei_mode to "kaye_nstx_lmode" recovers the
    #     common-tau_E model (tau_Ee=tau_Ei), just with a physics-based
    #     tau_E instead of a fixed input -- no fixed-point issue on this
    #     side either, so this combination works with enable_equipartition
    #     too (unlike "neoclassical").
    #   NOT YET SUPPORTED: tau_Ei_mode="neoclassical" together with
    #     enable_equipartition=True (would need neoclassical Ti nested
    #     inside the coupled Te/Ti bisection's inner solve, on top of its
    #     existing Te-Ti coupling -- raises ValueError instead of silently
    #     giving a wrong answer; not yet requested, so not yet implemented).
    # Defaults ("fixed", "fixed") leave every existing call site and result
    # completely unaffected -- both new modes are purely additive.
    # Electron-ion Coulomb equipartition (physics.equipartition_power_w,
    # NRL Plasma Formulary "Thermal Equilibration" rate): couples the T_e
    # and T_i balances via a P_ei(T_e,T_i) exchange term (subtracted from
    # the electron balance, added to the ion balance), replacing the
    # default's two INDEPENDENT solves with a genuinely coupled pair --
    # see docs/model.md's roadmap section for the derivation and the
    # "structural cost" this implies (a nested bisection, see
    # _solve_te_ti_coupled_keV, rather than the default's single 1-D
    # root-find). D and T thermal populations equilibrate with electrons as
    # two independent mass-correct channels; impurities are not included
    # (this model gives them no thermal population of their own). Default
    # False: leaves the existing decoupled Step-1/Step-5 cascade (and every
    # existing result) completely unaffected unless explicitly enabled.
    # Fraction of CAPTURED beam power (post shine-through, Step 0) lost to
    # charge-exchange with background neutrals during slowing-down, rather
    # than delivered to electron/ion heating -- a fast ion that
    # charge-exchanges becomes a fast NEUTRAL and escapes the plasma before
    # fully thermalizing. PPPL-1280 Sect. 6.2.2/6.2.3 quotes this loss as
    # "up to 10%" of injected power for Zeff=1 (a rough qualitative figure,
    # not a formula). Modeled here as a single fixed efficiency knob
    # (same pattern as physics.F_ALPHA) -- a flat user-supplied constant,
    # NOT a function of Eb, ne, or any background neutral density (confirmed
    # explicitly, not assumed: the physically correct version needs a new
    # n0 parameter this model doesn't otherwise carry, plus a CX-specific
    # cross section, feeding a survival-probability integral over the
    # slowing-down process -- see docs/model.md Step 0b for that formula).
    # Default 0.0 (off) so existing results are unaffected
    # unless explicitly enabled. See docs/model.md Step 0.
    # No beta_N / pressure_limit() here -- the predecessor project's Troyon
    # limit only ever existed to bound a max-Pf search; this model computes
    # beta_t as a reported diagnostic (see OperatingPoint.beta_t), never as
    # a constraint, since nothing is being optimized or searched here.


@dataclass
class OperatingPoint:
    ne0_m3: float
    feasible: bool
    Te_keV: float | None = None
    Ti_keV: float | None = None
    nb0_m3: float = 0.0
    nD0_m3: float = 0.0
    nT0_m3: float = 0.0
    n_thermal_m3: float = 0.0
    n_thermal_fraction: float = 0.0  # n_thermal / (n_thermal + n_b0): how close to the
    # n_thermal<=0 feasibility edge this point sits. Near 0, T_i blows up (a fixed P_i
    # heats an almost-empty thermal population) -- an expected consequence of dropping
    # e-i equipartition (docs/model.md), not a solver bug, but worth watching.
    P_e_w: float = 0.0
    P_i_w: float = 0.0
    P_ei_w: float = 0.0  # equipartition power, electrons->ions (can be negative); 0 unless config.enable_equipartition
    P_alpha_w: float = 0.0  # fusion alpha heating power actually folded into the T_e/T_i balance
    #   (0 unless the caller ran the alpha fixed-point -- solve_operating_point only sees it via
    #   extra_e_source_w/extra_i_source_w; the Picard loop that closes P_alpha ~ 0.2*P_fus lives
    #   in hotjass_core so solve_operating_point itself stays a single non-iterative solve)
    P_aux_e_w: float = 0.0  # prescribed auxiliary heating to electrons (ECRH + ICRH), folded into Step 1
    P_aux_i_w: float = 0.0  # prescribed auxiliary heating to ions (ECRH + ICRH), folded into Step 5
    Le: list[float] = field(default_factory=list)  # per-beam electron-heating fraction
    Li: list[float] = field(default_factory=list)  # per-beam ion-heating fraction
    P_NB_total_w: float = 0.0  # sum of injected (pre-shine-through) beam power, echoed for reference
    P_shine_w: float = 0.0  # total shine-through loss (not absorbed by the plasma)
    P_capt_w: float = 0.0  # P_NB_total_w - P_shine_w
    f_capture: list[float] = field(default_factory=list)  # per-beam captured fraction (1 - shine-through)
    P_orbit_loss_w: float = 0.0  # first-orbit loss (fraction of P_capt_w), only if config.enable_orbit_loss
    f_orbit_loss: list[float] = field(default_factory=list)  # per-beam first-orbit-loss fraction
    P_cx_loss_w: float = 0.0  # charge-exchange loss during slowing-down (fraction of post-orbit-loss power)
    f_cx_loss: list[float] = field(default_factory=list)  # per-beam CX power-loss fraction actually applied
    cx_f_particle: list[float] = field(default_factory=list)  # per-beam f_cx,N (0 in "manual_fraction" mode)
    cx_gamma: list[float] = field(default_factory=list)  # per-beam gamma_cx figure of merit (0 in "manual_fraction" mode)
    cx_escape_probability: list[float] = field(default_factory=list)  # per-beam informational escape factor (1.0 in "manual_fraction" mode)
    cx_n0_m3: list[float] = field(default_factory=list)  # per-beam background-neutral density actually used (0 in "manual_fraction" mode)
    P_useful_w: float = 0.0  # P_capt_w - P_orbit_loss_w - P_cx_loss_w == P_e_w + P_i_w -- what actually heats the plasma
    pf_thermal_w: float = 0.0   # D-T thermal fusion power
    pf_beam_w: float = 0.0      # D-T beam-target fusion power
    pf_dt_w: float = 0.0        # pf_thermal_w + pf_beam_w -- the D-T total (drives the alpha closure)
    pf_dd_thermal_w: float = 0.0  # D-D thermal fusion power (both branches)
    pf_dd_beam_w: float = 0.0     # D-D beam-target fusion power (fast D on thermal D)
    pf_dd_w: float = 0.0          # pf_dd_thermal_w + pf_dd_beam_w
    pf_total_w: float = 0.0       # pf_dt_w + pf_dd_w -- all fusion power
    neutron_rate_s: float = 0.0  # 14 MeV (D-T) + 2.45 MeV (D-D n-branch) neutrons per second
    neutron_rate_thermal_s: float = 0.0  # thermal-thermal contribution (D-T + D-D) to neutron_rate_s
    neutron_rate_beam_s: float = 0.0     # beam-target contribution (D-T + D-D) to neutron_rate_s
    neutron_rate_bb_s: float = 0.0       # beam-beam contribution (D-T + D-D) to neutron_rate_s; 0 unless config.enable_beam_beam
    pf_bb_dt_w: float = 0.0  # beam-beam D-T fusion power (already folded into pf_dt_w)
    pf_bb_dd_w: float = 0.0  # beam-beam D-D fusion power (already folded into pf_dd_w)
    # neutron_rate_thermal_s + neutron_rate_beam_s + neutron_rate_bb_s == neutron_rate_s exactly.
    v_phi_m_s: float = 0.0  # bulk toroidal rotation velocity (positive = co-current); 0 unless config.rotation_model != "off"
    torque_total_Nm: float = 0.0  # net NBI torque actually used to close v_phi ("momentum_balance" only; 0 otherwise)
    Te0_keV: float | None = None  # on-axis T_e (== Te_keV unless config.profile_averaging)
    Ti0_keV: float | None = None  # on-axis T_i (== Ti_keV unless config.profile_averaging)
    pressure_pa: float = 0.0  # isotropic fast-ion pressure assumption (config.fast_ion_mean_pitch2, default 1/3 -> p_fast=(2/3)u_fast)
    beta_t: float = 0.0  # beta_t from pressure_pa (isotropic, per config.fast_ion_mean_pitch2)
    pressure_anisotropic_pa: float = 0.0  # SAME point, fast-ion pressure recomputed at pitch2=1.0 (purely
    # parallel/passing -- this project's tangential R_t=R_0 injection geometry at birth, before any
    # pitch-angle scattering has isotropized the population) -- p_fast=0 in this limit, so this is the
    # OTHER bounding case from the isotropic default, not a separate free parameter. Always computed
    # (not gated by config.fast_ion_mean_pitch2), so both bounds are available on every OperatingPoint.
    beta_t_anisotropic: float = 0.0  # beta_t from pressure_anisotropic_pa
    avg_fast_energy_keV: float = 0.0  # mean fast-ion energy (physics.average_fast_energy_keV, power-weighted
    # across beams) -- exposed so pressure_pa/pressure_anisotropic_pa can both be recomputed from stored
    # OperatingPoint fields alone, without re-solving Te/Ti (fast-ion pitch anisotropy doesn't affect the solve).
    R_fast_thermal: float = 0.0
    tau_E_s: float = 0.0  # electron-side confinement time actually used (Step 1's Te balance)
    tau_Ei_s: float = 0.0  # ion-side confinement time actually used (Step 5's Ti balance) -- equals
    # tau_E_s unless solve_operating_point's tau_Ei_s argument was explicitly given a different
    # value (split tau_E, see TokamakConfig-adjacent docs/model.md section); NOT part of
    # TokamakConfig, since (like tau_E_s itself) it's a per-call scanned quantity, not a static
    # device/config parameter.
    infeasible_reason: str | None = None


def _n_sum_for_ne(ne0_m3: float, Zeff: float) -> float:
    """Exact inverse of physics.compute_charge_neutrality(): the (D+T+fast)
    ion sum density that produces a given target electron density ne0.
    """
    if Zeff <= 1.0 or physics.Z_IMPURITY <= Zeff:
        return ne0_m3
    alpha_i = (Zeff - 1.0) / (physics.Z_IMPURITY * (physics.Z_IMPURITY - Zeff))
    alpha_e = 1.0 + physics.Z_IMPURITY * alpha_i
    return ne0_m3 / alpha_e


def _solve_te_keV(ne0_m3: float, beams: list[BeamSpec], tau_E_s: float, volume_m3: float,
                  extra_e_source_w: float = 0.0) -> float:
    """Step 1: solve the electron power balance
        (3/2)*ne0*Te*(1e3*e)*V/tau_E = sum_i Le(Eb_i/Ec_i(Te)) * P_NB_i + P_extra_e
    for Te alone (proof of a unique root: docs/model.md). extra_e_source_w is a
    Te-independent additive electron source (fusion alpha heating, closed by the
    caller's fixed-point) -- adding a positive constant to the RHS only shifts the
    unique root up, so the monotonicity argument is unaffected. Ec is now PER BEAM
    (critical_energy_keV depends on the beam's own species/mass -- Phase 2,
    D+T beams) rather than shared -- the monotonicity proof is unaffected,
    since each Ec_i(Te) is still individually strictly increasing in Te, so
    x_i=Eb_i/Ec_i(Te) is still strictly decreasing and L_e(x_i) strictly
    increasing, term by term. Self-expanding bracket + bisection, no scipy
    dependency.
    """

    def residual(Te_keV: float) -> float:
        lhs = 1.5 * ne0_m3 * Te_keV * 1e3 * physics.E_CHARGE * volume_m3 / tau_E_s
        rhs = extra_e_source_w
        for beam in beams:
            Ec = physics.critical_energy_keV(Te_keV, beam.species)
            x = beam.Eb_keV / max(Ec, 1e-12)
            rhs += physics.electron_heating_fraction(x) * beam.P_NB_W
        return lhs - rhs

    lo, hi = 1.0e-3, 1.0
    while residual(hi) < 0.0:
        hi *= 2.0
        if hi > 1.0e4:
            raise RuntimeError(
                f"Te root-find did not bracket below 1e4 keV (ne0={ne0_m3:.3e}, tau_E={tau_E_s:.3e}, "
                f"beams={beams}) -- check inputs, this should be unreachable per the monotonicity proof."
            )

    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if residual(mid) < 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _solve_te_ti_coupled_keV(
    ne0_m3: float, useful_beams: list[BeamSpec], tau_Ee_s: float, tau_Ei_s: float, volume_m3: float,
    config: TokamakConfig, extra_e_source_w: float = 0.0, extra_i_source_w: float = 0.0,
) -> tuple[float, float, float, float, list[float], list[float], float, float, float]:
    """Coupled electron/ion power balance with equipartition
    (config.enable_equipartition):
        (3/2)*ne0*Te*(1e3e)*V/tau_Ee = sum_i Le(x_i)*P_NB_i - P_ei(Te,Ti)
        (3/2)*n_thermal*Ti*(1e3e)*V/tau_Ei = sum_i Li(x_i)*P_NB_i + P_ei(Te,Ti)

    tau_Ee_s and tau_Ei_s may differ (split-tau_E model, docs/model.md) --
    each balance uses its own confinement time; passing tau_Ee_s==tau_Ei_s
    exactly reproduces the original common-tau_E behavior.

    docs/model.md's roadmap flags this as a genuine 2-D root-find (tau_ei
    depends on both temperatures) with monotonicity "not yet checked."
    Implemented here as a NESTED pair of 1-D bisections instead -- outer
    bisection on Te (same self-expanding-bracket style as _solve_te_keV),
    with an INNER bisection on Ti at each trial Te -- rather than a
    hand-rolled multivariate Newton solver or a scipy dependency. This
    works because n_thermal(Te), the L_e/L_i split, and nb0(Te) all depend
    on Te ALONE (equipartition doesn't touch beam slowing-down physics, only
    adds a second, separate thermal-thermal energy channel) -- so at each
    trial Te the inner Ti-solve is a genuine 1-D problem, and the outer
    residual is well-defined once that inner solve returns.

    n_thermal(Te) can be <=0 for some Te values tried DURING the search
    (before convergence) -- clamped to a small positive floor for use
    inside this search only, exactly mirroring how the original decoupled
    solve never even evaluates n_thermal until Te is already final. The
    REAL (unclamped) n_thermal at the converged Te is what the caller uses
    to decide feasibility, same as before.
    """

    def n_thermal_raw(Te_keV: float) -> tuple[float, float]:
        nb0 = 0.0
        for beam in useful_beams:
            tau_s = physics.thermalization_time(ne0_m3, Te_keV, beam.Eb_keV, beam.species)
            nb0 += beam.P_NB_W * tau_s / (beam.Eb_keV * 1e3 * physics.E_CHARGE * volume_m3)
        n_sum = _n_sum_for_ne(ne0_m3, config.Zeff)
        return n_sum - nb0, nb0

    def le_li_split(Te_keV: float):
        P_e_beam, P_i_beam = 0.0, 0.0
        Le_list, Li_list = [], []
        for beam in useful_beams:
            Ec = physics.critical_energy_keV(Te_keV, beam.species)
            x = beam.Eb_keV / max(Ec, 1e-12)
            li = physics.ion_heating_fraction(x)
            le = 1.0 - li
            Le_list.append(le)
            Li_list.append(li)
            P_e_beam += le * beam.P_NB_W
            P_i_beam += li * beam.P_NB_W
        return P_e_beam, P_i_beam, Le_list, Li_list

    def p_ei(Te_keV: float, Ti_keV: float, n_thermal_eff: float) -> float:
        nD0 = config.mix_D * n_thermal_eff
        nT0 = config.mix_T * n_thermal_eff
        return (physics.equipartition_power_w(nD0, Te_keV, Ti_keV, ne0_m3, "D", volume_m3)
                + physics.equipartition_power_w(nT0, Te_keV, Ti_keV, ne0_m3, "T", volume_m3))

    def solve_ti(Te_keV: float, n_thermal_eff: float, P_i_beam: float) -> float:
        def residual(Ti_keV: float) -> float:
            lhs = 1.5 * n_thermal_eff * Ti_keV * 1e3 * physics.E_CHARGE * volume_m3 / tau_Ei_s
            rhs = P_i_beam + p_ei(Te_keV, Ti_keV, n_thermal_eff) + extra_i_source_w
            return lhs - rhs

        lo, hi = 1.0e-4, max(Te_keV * 2.0, 1.0)
        for _ in range(80):
            if residual(hi) >= 0.0:
                break
            hi *= 2.0
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if residual(mid) < 0.0:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    def residual_te(Te_keV: float):
        n_thermal, nb0 = n_thermal_raw(Te_keV)
        n_thermal_eff = max(n_thermal, 1.0e10)  # search-only floor, see docstring
        P_e_beam, P_i_beam, Le_list, Li_list = le_li_split(Te_keV)
        Ti_keV = solve_ti(Te_keV, n_thermal_eff, P_i_beam)
        P_ei = p_ei(Te_keV, Ti_keV, n_thermal_eff)
        lhs = 1.5 * ne0_m3 * Te_keV * 1e3 * physics.E_CHARGE * volume_m3 / tau_Ee_s
        rhs = P_e_beam - P_ei + extra_e_source_w
        return (lhs - rhs), Ti_keV, n_thermal, nb0, P_e_beam, P_i_beam, Le_list, Li_list, P_ei

    lo, hi = 1.0e-3, 1.0
    for _ in range(80):
        r_hi, *_ = residual_te(hi)
        if r_hi >= 0.0:
            break
        hi *= 2.0
        if hi > 1.0e4:
            raise RuntimeError(
                f"Coupled Te/Ti root-find did not bracket below 1e4 keV (ne0={ne0_m3:.3e}, "
                f"tau_Ee={tau_Ee_s:.3e}, tau_Ei={tau_Ei_s:.3e}, beams={useful_beams}) -- monotonicity of "
                "the coupled system has not been proven (docs/model.md); check inputs."
            )

    for _ in range(60):
        mid = 0.5 * (lo + hi)
        r_mid, *_ = residual_te(mid)
        if r_mid < 0.0:
            lo = mid
        else:
            hi = mid
    Te_keV = 0.5 * (lo + hi)
    _, Ti_keV, n_thermal, nb0, P_e_beam, P_i_beam, Le_list, Li_list, P_ei = residual_te(Te_keV)
    return Te_keV, Ti_keV, n_thermal, nb0, Le_list, Li_list, P_e_beam, P_i_beam, P_ei


def _solve_ti_neoclassical_keV(
    P_i_w: float, n_thermal_m3: float, volume_m3: float, config: "TokamakConfig",
) -> tuple[float, float]:
    """Step 5 (config.tau_Ei_mode="neoclassical"): solve
        (3/2)*n_thermal*Ti*(1e3e)*V / tau_Ei,neo(Ti) = P_i_w
    for Ti, where tau_Ei,neo(Ti)=physics.neoclassical_tau_Ei_s(Ti,...) is
    ITSELF a function of Ti (~Ti^0.5, see that function's docstring) --
    unlike the default closed form, this is a genuine fixed-point 1-D
    root-find, not a substitution. LHS ~ Ti/tau_Ei,neo(Ti) ~ Ti^0.5, still
    monotonically increasing in Ti, so a unique root exists -- same
    self-expanding-bracket + bisection style as _solve_te_keV.

    Returns (Ti_keV, tau_Ei_neo_s_at_root) -- the second value is what
    OperatingPoint.tau_Ei_s reports: the confinement time ACTUALLY realized
    at the converged Ti, not an externally supplied input (there is none in
    this mode).
    """
    geom = config.geometry

    def tau_ei(Ti_keV: float) -> float:
        return physics.neoclassical_tau_Ei_s(
            Ti_keV, n_thermal_m3, config.Bt0, config.Ip_MA,
            geom.major_radius, geom.minor_radius, geom.elongation,
            config.mix_D, config.mix_T,
            arbitrary_A=(config.tau_Ei_mode == "neoclassical_arbA"),
            triangularity=geom.triangularity,
        )

    def residual(Ti_keV: float) -> float:
        lhs = 1.5 * n_thermal_m3 * Ti_keV * 1e3 * physics.E_CHARGE * volume_m3 / tau_ei(Ti_keV)
        return lhs - P_i_w

    # Bracket ceiling is deliberately huge (1e9 keV, not the usual 1e4 used
    # elsewhere in this file): tau_Ei,neo(Ti) ~ Ti^0.5 grows so slowly that
    # LHS(Ti) ~ Ti/tau_Ei,neo(Ti) ~ sqrt(Ti) itself only grows slowly, so
    # genuinely converged roots at extreme Ti (thousands of keV and up) are
    # possible and expected in this mode -- see docs/model.md: pure
    # neoclassical ion transport, with no anomalous channel, is known to
    # radically UNDERESTIMATE real ion transport (that's precisely why real
    # tokamaks are anomalous-transport-dominated), so it can predict a very
    # weak ion heat sink and correspondingly very hot, often physically
    # implausible Ti -- a genuine property of this model choice, not a
    # solver bug. Doubling from 1e-4 to 1e9 is only ~44 iterations either
    # way, so the generous ceiling costs nothing computationally.
    lo, hi = 1.0e-4, 1.0
    while residual(hi) < 0.0:
        hi *= 2.0
        if hi > 1.0e9:
            raise RuntimeError(
                f"Neoclassical Ti root-find did not bracket below 1e9 keV "
                f"(n_thermal={n_thermal_m3:.3e}, P_i={P_i_w:.3e}) -- check inputs."
            )
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if residual(mid) < 0.0:
            lo = mid
        else:
            hi = mid
    Ti_keV = 0.5 * (lo + hi)
    return Ti_keV, tau_ei(Ti_keV)


def solve_operating_point(
    ne0_m3: float,
    beams: list[BeamSpec],
    tau_E_s: float,
    config: TokamakConfig,
    tau_Ei_s: float | None = None,
    extra_e_source_w: float = 0.0,
    extra_i_source_w: float = 0.0,
    aux_e_source_w: float = 0.0,
    aux_i_source_w: float = 0.0,
    cx_loss_fraction_override: list[dict] | None = None,
) -> OperatingPoint:
    """The full steps 0-7 cascade (docs/model.md) at one scanned n_e0.

    tau_E_s governs the electron balance (Step 1); tau_Ei_s, if given,
    governs the ion balance (Step 5) instead -- a split-tau_E model
    (e.g. tau_E_s=0.01 for electrons, tau_Ei_s=0.02 for ions), motivated by
    the same physical picture as this model's hot-ion-mode drivers
    (docs/model.md): electrons and ions need not actually share one
    transport time. Defaulting tau_Ei_s to tau_E_s when omitted (the
    common-tau_E case) reproduces every existing call site and result
    exactly -- this parameter is purely additive.

    cx_loss_fraction_override: one dict per beam (physics.cx_loss_fractions()'s
    return value, plus "escape_probability" and "n0_m3"), used in place of
    the flat config.cx_loss_fraction when config.cx_model != "manual_fraction".
    This function stays a pure passthrough for that dict -- the physics-based
    CX modes need Te (via Ec(Te)) to compute it, which Te itself is being
    solved for here, so the Te-dependence is a genuine fixed point closed by
    an outer Picard loop in hotjass_core (same architectural pattern as the
    alpha-heating fixed point: this function itself stays a single
    non-iterative solve). None (default) reproduces every existing result
    exactly -- purely additive.
    """
    tau_Ei_s = tau_E_s if tau_Ei_s is None else tau_Ei_s
    if config.tau_Ee_mode not in ("fixed", "kaye_nstx_lmode", "kaye_nstx_hmode", "iter98y2"):
        raise ValueError(f"Unknown tau_Ee_mode {config.tau_Ee_mode!r}")
    if config.tau_Ei_mode not in (
        "fixed", "neoclassical", "neoclassical_arbA", "kaye_nstx_lmode", "kaye_nstx_hmode", "iter98y2",
    ):
        raise ValueError(f"Unknown tau_Ei_mode {config.tau_Ei_mode!r}")
    if config.tau_Ei_mode in ("neoclassical", "neoclassical_arbA") and config.enable_equipartition:
        raise ValueError(
            "tau_Ei_mode='neoclassical' combined with enable_equipartition=True is not yet "
            "implemented (would need neoclassical Ti nested inside the coupled Te/Ti bisection's "
            "inner solve) -- see TokamakConfig docstring."
        )
    V = config.geometry.volume()
    path_length_m = physics.tangential_path_length(config.geometry.minor_radius)

    # Profile averaging (config.profile_averaging). pk_n / pk_te / pk_ti convert
    # between the on-axis and the volume-averaged value of a (1-rho^2)^(2p)
    # profile: <X>/X0 = 1/(1 + 2p). The electron and ion temperature profiles
    # can have different peaking (temperature_peaking vs temperature_peaking_i)
    # -- important in a hot-ion regime where T_i is far more peaked than T_e.
    # All three are 1.0 (a no-op) when profile_averaging is off, so every
    # existing result is byte-identical.
    p_te = max(config.temperature_peaking, 0.0)
    p_ti = p_te if config.temperature_peaking_i < 0.0 else max(config.temperature_peaking_i, 0.0)
    if config.profile_averaging:
        pk_n = 1.0 + 2.0 * max(config.density_peaking, 0.0)
        pk_te = 1.0 + 2.0 * p_te
        pk_ti = 1.0 + 2.0 * p_ti
    else:
        pk_n = pk_te = pk_ti = 1.0
    ne_bal = ne0_m3 / pk_n          # density seen by the energy balance & tau_E scalings
    ne_axis = ne0_m3               # on-axis density (== ne0_m3; the input, in profile mode)

    # Step 0: beam shine-through / capture fraction, then (optionally)
    # first-orbit loss, then charge-exchange loss. All three depend only on
    # ne0 and each beam's own (Eb, species) -- not on Te -- so all are
    # computed once, up front, and every later step uses the USEFUL power
    # (useful_beams) in place of the raw injected P_NB. Ordering matters:
    # first-orbit loss removes ions PROMPTLY at birth, before slowing-down
    # or CX has a chance to occur, so it's applied to the captured
    # population before CX-loss (which represents losses accumulated over
    # the whole slowing-down lifetime of the ions that survive birth). See
    # docs/model.md.
    f_capture_list: list[float] = []
    f_orbit_loss_list: list[float] = []
    f_cx_loss_list: list[float] = []
    cx_f_particle_list: list[float] = []
    cx_gamma_list: list[float] = []
    cx_escape_list: list[float] = []
    cx_n0_list: list[float] = []
    useful_beams: list[BeamSpec] = []
    P_NB_total_w = 0.0
    P_capt_w = 0.0
    P_orbit_loss_w = 0.0
    P_useful_w = 0.0
    for beam_index, beam in enumerate(beams):
        f_capt = physics.captured_power_fraction(
            ne0_m3, beam.Eb_keV, beam.species, path_length_m, config.density_peaking,
            geometry=config.geometry, tangent_R_m=beam.tangent_R_m, tangent_Z_m=beam.tangent_Z_m,
            model=beam.shine_through_model,
            manual_shine_through_fraction=beam.manual_shine_through_fraction,
            Zeff=config.Zeff,
        )
        f_capture_list.append(f_capt)
        P_NB_total_w += beam.P_NB_W
        P_capt_i = f_capt * beam.P_NB_W
        P_capt_w += P_capt_i

        if config.enable_orbit_loss:
            # first_orbit_loss_fraction needs an attenuation SHAPE for the birth
            # profile; "manual" shine-through only fixes a fraction, not a shape,
            # so fall back to the standard Riviere shape in that case.
            orbit_stopping_model = (
                "riviere" if beam.shine_through_model == "manual" else beam.shine_through_model
            )
            co = getattr(beam, "co_current", config.orbit_loss_co_current)
            # Off-axis / vertically-shifted injection: sample the real
            # (R_t, Z_t) chord for the birth-radius profile. For the standard
            # on-axis tangential case (tangent_Z_m == 0 and R_t == R0 or
            # unset) keep chord=None so the crude linear rho(x) mapping -- and
            # every existing orbit-loss result -- is byte-for-byte unchanged.
            _R0 = config.geometry.major_radius
            _tR = beam.tangent_R_m
            orbit_chord = None
            _off_axis = beam.tangent_Z_m != 0.0 or (_tR is not None and abs(_tR - _R0) > 1e-6)
            if _off_axis:
                _rcp = config.centrepost_radius_m if config.centrepost_radius_m > 0.0 else None
                orbit_chord = physics.tangential_chord(
                    _R0, config.geometry.minor_radius, config.geometry.elongation,
                    tangent_R_m=_tR, tangent_Z_m=beam.tangent_Z_m, R_centrepost_m=_rcp,
                )
            if _off_axis and orbit_chord is None:
                # aim point outside the plasma -> nothing born -> no orbit loss
                f_orbit = 0.0
            elif config.orbit_model in ("st_meanshift", "st_pitch"):
                widths = physics.st_orbit_widths(
                    Eb_keV=beam.Eb_keV, Bt_T=config.Bt0, Ip_MA=config.Ip_MA,
                    R0_m=config.geometry.major_radius, minor_radius_m=config.geometry.minor_radius,
                    elongation=config.geometry.elongation, triangularity=config.geometry.triangularity,
                    species=beam.species,
                )
                f_orbit = physics.first_orbit_loss_fraction_st(
                    ne0_m3, beam.Eb_keV, beam.species, path_length_m,
                    config.geometry.minor_radius, config.geometry.major_radius, widths,
                    co_current=co,
                    variant="pitch" if config.orbit_model == "st_pitch" else "meanshift",
                    stopping_model=orbit_stopping_model, Zeff=config.Zeff,
                    tangent_R_m=beam.tangent_R_m, chord=orbit_chord,
                )
            else:
                orbit_width_m = physics.passing_orbit_width(
                    Eb_keV=beam.Eb_keV, Bt_T=config.Bt0, Ip_MA=config.Ip_MA,
                    R0_m=config.geometry.major_radius, minor_radius_m=config.geometry.minor_radius,
                    elongation=config.geometry.elongation, species=beam.species,
                )
                larmor_radius_m_value = physics.larmor_radius_m(Eb_keV=beam.Eb_keV, Bt_T=config.Bt0, species=beam.species)
                f_orbit = physics.first_orbit_loss_fraction(
                    ne0_m3, beam.Eb_keV, beam.species, path_length_m, config.geometry.minor_radius, orbit_width_m,
                    co_current=co, chord=orbit_chord,
                    include_larmor_loss=config.include_larmor_loss, larmor_radius_m_value=larmor_radius_m_value,
                    stopping_model=orbit_stopping_model,
                    Zeff=config.Zeff,
                )
        else:
            f_orbit = 0.0
        f_orbit_loss_list.append(f_orbit)
        P_after_orbit_i = P_capt_i * (1.0 - f_orbit)
        P_orbit_loss_w += P_capt_i - P_after_orbit_i

        if cx_loss_fraction_override is not None:
            cx_d = cx_loss_fraction_override[beam_index]
            f_cx = cx_d["f_cx_P"]
            cx_f_particle_list.append(cx_d.get("f_cx_N", 0.0))
            cx_gamma_list.append(cx_d.get("gamma_cx", 0.0))
            cx_escape_list.append(cx_d.get("escape_probability", 1.0))
            cx_n0_list.append(cx_d.get("n0_m3", 0.0))
        else:
            f_cx = config.cx_loss_fraction
            cx_f_particle_list.append(f_cx)
            cx_gamma_list.append(0.0)
            cx_escape_list.append(1.0)
            cx_n0_list.append(0.0)
        f_cx_loss_list.append(f_cx)
        P_useful_i = P_after_orbit_i * (1.0 - f_cx)
        P_useful_w += P_useful_i
        useful_beams.append(
            BeamSpec(
                Eb_keV=beam.Eb_keV, P_NB_W=P_useful_i, species=beam.species,
                tangent_R_m=beam.tangent_R_m, tangent_Z_m=beam.tangent_Z_m,
                shine_through_model=beam.shine_through_model,
                manual_shine_through_fraction=beam.manual_shine_through_fraction,
                co_current=beam.co_current,
            )
        )
    P_shine_w = P_NB_total_w - P_capt_w
    P_cx_loss_w = P_capt_w - P_orbit_loss_w - P_useful_w

    # Physics-based tau_E scalings (Kaye NSTX L-mode, IPB98(y,2)): both depend
    # only on Ip/Bt/ne0/geometry and the loss power P_loss (known after Step 0),
    # never on Te or Ti, so each is computed once, up front, and reused for
    # whichever channel(s) request it -- no fixed-point issue on either side
    # (unlike neoclassical tau_Ei). P_loss here is the captured, post-orbit,
    # post-CX beam power plus prescribed auxiliary heating (ECRH/ICRH); fusion
    # alpha heating is not included (it is closed by an outer fixed-point and is
    # not yet known at this point).
    p_loss_scaling_w = P_useful_w + aux_e_source_w + aux_i_source_w
    if config.tau_Ee_mode == "kaye_nstx_lmode" or config.tau_Ei_mode == "kaye_nstx_lmode":
        kaye_lmode_tau_s = physics.kaye_nstx_lmode_confinement_time(
            Ip_A=config.Ip_MA * 1e6, Bt_T=config.Bt0, ne0_m3=ne_bal, P_loss_W=p_loss_scaling_w,
        )
        if config.tau_Ee_mode == "kaye_nstx_lmode":
            tau_E_s = kaye_lmode_tau_s
        if config.tau_Ei_mode == "kaye_nstx_lmode":
            tau_Ei_s = kaye_lmode_tau_s
    if config.tau_Ee_mode == "kaye_nstx_hmode" or config.tau_Ei_mode == "kaye_nstx_hmode":
        kaye_hmode_tau_s = physics.kaye_nstx_confinement_time(
            Ip_A=config.Ip_MA * 1e6, Bt_T=config.Bt0, ne0_m3=ne_bal, P_heat_W=p_loss_scaling_w,
        )
        if config.tau_Ee_mode == "kaye_nstx_hmode":
            tau_E_s = kaye_hmode_tau_s
        if config.tau_Ei_mode == "kaye_nstx_hmode":
            tau_Ei_s = kaye_hmode_tau_s
    if config.tau_Ee_mode == "iter98y2" or config.tau_Ei_mode == "iter98y2":
        iter98y2_tau_s = physics.iter98y2_confinement_time(
            Ip_MA=config.Ip_MA, Bt_T=config.Bt0, ne0_m3=ne_bal, P_heat_W=p_loss_scaling_w,
            R0_m=config.geometry.major_radius, minor_radius_m=config.geometry.minor_radius,
            elongation=config.geometry.elongation,
            M_eff_amu=config.mix_D * 2.0 + config.mix_T * 3.0,
        )
        if config.tau_Ee_mode == "iter98y2":
            tau_E_s = iter98y2_tau_s
        if config.tau_Ei_mode == "iter98y2":
            tau_Ei_s = iter98y2_tau_s

    # Te-independent additive sources folded into the electron/ion balances:
    # extra_*_source_w is the fusion alpha term (closed by the caller's
    # fixed-point), aux_*_source_w is prescribed auxiliary heating (ECRH/ICRH).
    # The solvers only see the sum; the OperatingPoint reports them separately.
    ee_source_w = extra_e_source_w + aux_e_source_w
    ei_source_w = extra_i_source_w + aux_i_source_w

    P_ei_w = 0.0
    if config.enable_equipartition:
        # Coupled solve: Te and Ti (and hence P_e_w/P_i_w, which stay the
        # pure BEAM heating split -- P_ei_w is the separate equipartition
        # exchange) all come out of one nested-bisection call.
        Te_keV, Ti_keV, n_thermal, nb0_total, Le_list, Li_list, P_e_w, P_i_w, P_ei_w = \
            _solve_te_ti_coupled_keV(ne_bal, useful_beams, tau_E_s, tau_Ei_s, V, config,
                                     ee_source_w, ei_source_w)
    else:
        Te_keV = _solve_te_keV(ne_bal, useful_beams, tau_E_s, V, ee_source_w)

        Le_list = []
        Li_list = []
        P_e_w = 0.0
        P_i_w = 0.0
        nb0_total = 0.0
        for beam in useful_beams:
            Ec = physics.critical_energy_keV(Te_keV, beam.species)
            x = beam.Eb_keV / max(Ec, 1e-12)
            li = physics.ion_heating_fraction(x)
            le = 1.0 - li
            Li_list.append(li)
            Le_list.append(le)
            P_e_w += le * beam.P_NB_W
            P_i_w += li * beam.P_NB_W
            tau_s = physics.thermalization_time(ne_bal, Te_keV, beam.Eb_keV, beam.species)
            nb0_total += beam.P_NB_W * tau_s / (beam.Eb_keV * 1e3 * physics.E_CHARGE * V)

        n_thermal = _n_sum_for_ne(ne_bal, config.Zeff) - nb0_total

    n_sum = _n_sum_for_ne(ne_bal, config.Zeff)
    if n_thermal <= 0.0:
        return OperatingPoint(
            ne0_m3=ne0_m3, feasible=False, Te_keV=Te_keV, nb0_m3=nb0_total,
            P_e_w=P_e_w, P_i_w=P_i_w, P_ei_w=P_ei_w, P_alpha_w=extra_e_source_w + extra_i_source_w,
            P_aux_e_w=aux_e_source_w, P_aux_i_w=aux_i_source_w,
            Le=Le_list, Li=Li_list,
            tau_E_s=tau_E_s, tau_Ei_s=tau_Ei_s,
            P_NB_total_w=P_NB_total_w, P_shine_w=P_shine_w, P_capt_w=P_capt_w, f_capture=f_capture_list,
            P_orbit_loss_w=P_orbit_loss_w, f_orbit_loss=f_orbit_loss_list,
            P_cx_loss_w=P_cx_loss_w, f_cx_loss=f_cx_loss_list, cx_f_particle=cx_f_particle_list,
            cx_gamma=cx_gamma_list, cx_escape_probability=cx_escape_list, cx_n0_m3=cx_n0_list,
            P_useful_w=P_useful_w,
            infeasible_reason=(
                f"n_thermal={n_thermal:.3e} <= 0: fast-ion density {nb0_total:.3e} m^-3 required to "
                f"sustain P_useful={P_useful_w/1e6:.3f} MW (of P_NB={P_NB_total_w/1e6:.3f} MW injected) "
                f"exceeds what charge neutrality allows at this n_e (n_sum={n_sum:.3e})"
            ),
        )

    if not config.enable_equipartition:
        P_i_total_w = P_i_w + ei_source_w
        if config.tau_Ei_mode in ("neoclassical", "neoclassical_arbA"):
            Ti_keV, tau_Ei_s = _solve_ti_neoclassical_keV(P_i_total_w, n_thermal, V, config)
        else:
            Ti_keV = P_i_total_w * tau_Ei_s / (1.5 * n_thermal * 1e3 * physics.E_CHARGE * V)

    nD0 = config.mix_D * n_thermal
    nT0 = config.mix_T * n_thermal

    # On-axis values for the fusion / pressure integrals. With profile_averaging
    # off, pk_n == pk_te == pk_ti == 1 so these equal the balance quantities
    # exactly and every result below is byte-identical to before.
    Te0_keV = Te_keV * pk_te
    Ti0_keV = Ti_keV * pk_ti
    nD0_axis = nD0 * pk_n
    nT0_axis = nT0 * pk_n

    pf_thermal = physics.thermal_fusion_power(
        nD0_axis, nT0_axis, Ti0_keV, V, config.density_peaking, p_ti
    )

    # Bulk toroidal rotation (config.rotation_model): computed here, after
    # Te/Ti/n_thermal are known, since neither mode feeds back into the
    # energy balance -- only the beam-target reactivity below (via each
    # beam's own co/counter direction) sees v_phi, so no fixed-point/Picard
    # loop is needed (unlike the CX-loss physics modes, whose n0 depends on
    # Te *before* Te is solved).
    torque_total_Nm = 0.0
    if config.rotation_model == "manual":
        v_phi_m_s = config.manual_v_phi_m_s
    elif config.rotation_model == "momentum_balance":
        R0 = config.geometry.major_radius
        for beam in useful_beams:
            tR = beam.tangent_R_m if beam.tangent_R_m is not None else R0
            torque_total_Nm += physics.nbi_torque_Nm(
                beam.P_NB_W, beam.Eb_keV, beam.species, tR, beam.co_current)
        rho_i_kg_m3 = nD0_axis * physics.M_D + nT0_axis * physics.M_T
        tau_phi_s = tau_Ei_s * config.tau_phi_over_tauEi
        v_phi_m_s = physics.toroidal_rotation_velocity_ms(
            torque_total_Nm, rho_i_kg_m3, R0, V, tau_phi_s)
    else:
        v_phi_m_s = 0.0

    pf_beam = 0.0
    pf_dd_beam = 0.0
    dd_beam_neutrons = 0.0
    nb0_per_beam: list[float] = []
    for beam in useful_beams:
        # Per-beam fast-ion density, recomputed (cheap) rather than stored,
        # to keep beam_target_fusion_power's existing per-beam signature.
        # beam.P_NB_W here is already the USEFUL power (Step 0: captured,
        # minus CX loss) -- only particles that survive both shine-through
        # and charge-exchange become fast ions.
        tau_s = physics.thermalization_time(ne_axis, Te0_keV, beam.Eb_keV, beam.species)
        nb0_i = beam.P_NB_W * tau_s / (beam.Eb_keV * 1e3 * physics.E_CHARGE * V)
        nb0_per_beam.append(nb0_i)
        target_n = nT0_axis if beam.species == "D" else nD0_axis
        pf_beam += physics.beam_target_fusion_power(
            nb0_i, target_n, Te0_keV, beam.Eb_keV, V, beam.species,
            ne_axis, config.density_peaking, p_te,
            v_phi_m_s, beam.co_current,
        )
        # D-D beam-target: a fast D ion on the thermal-D population.
        if beam.species == "D" and nD0_axis > 0.0:
            p_dd, r_dd = physics.beam_target_dd_fusion_power(
                nb0_i, nD0_axis, Te0_keV, beam.Eb_keV, V,
                ne_axis, config.density_peaking, p_te,
                v_phi_m_s, beam.co_current,
            )
            pf_dd_beam += p_dd
            dd_beam_neutrons += r_dd
    pf_dd_thermal, dd_thermal_neutrons = physics.thermal_dd_fusion_power(
        nD0_axis, Ti0_keV, V, config.density_peaking, p_ti
    )

    # Beam-beam fusion (config.enable_beam_beam): the ONE pairwise reaction
    # between the first two useful beams (this project's tested scope).
    # Off by default -- pf_bb_dt/pf_bb_dd/neutron_rate_bb all stay 0.0.
    pf_bb_dt = 0.0
    pf_bb_dd = 0.0
    neutron_rate_bb = 0.0
    if config.enable_beam_beam and len(useful_beams) >= 2:
        b1, b2 = useful_beams[0], useful_beams[1]
        e1 = physics.average_fast_energy_keV(Te0_keV, b1.Eb_keV, b1.species)
        e2 = physics.average_fast_energy_keV(Te0_keV, b2.Eb_keV, b2.species)
        bb = physics.beam_beam_fusion_power(
            nb0_per_beam[0], e1, b1.species, b1.co_current,
            nb0_per_beam[1], e2, b2.species, b2.co_current,
            V,
        )
        pf_bb_dt = bb["pf_dt_w"]
        pf_bb_dd = bb["pf_dd_w"]
        neutron_rate_bb = bb["neutron_rate_s"]

    pf_dt = pf_thermal + pf_beam + pf_bb_dt
    pf_dd = pf_dd_thermal + pf_dd_beam + pf_bb_dd
    pf_total = pf_dt + pf_dd
    # One 14 MeV neutron per D-T reaction, plus the D-D n-branch's 2.45 MeV
    # neutron -- split by thermal-thermal / beam-target / beam-beam so the
    # three can be reported (and plotted) separately; their sum is
    # neutron_rate exactly.
    neutron_rate_thermal = pf_thermal / physics.E_FUSION_J + dd_thermal_neutrons
    neutron_rate_beam = pf_beam / physics.E_FUSION_J + dd_beam_neutrons
    neutron_rate = neutron_rate_thermal + neutron_rate_beam + neutron_rate_bb

    # Diagnostics -- reported, never enforced (see TokamakConfig docstring).
    # P_useful_w==0 (e.g. cx_loss_fraction=1.0, or ne0 low enough that
    # shine-through captures ~nothing) means nb0_total==0 too, so u_fast
    # will be 0 regardless of avg_fast_energy -- guard the division rather
    # than let a real (if extreme) input raise.
    if P_useful_w > 0.0:
        avg_fast_energy = sum(
            physics.average_fast_energy_keV(Te_keV, beam.Eb_keV, beam.species) * (beam.P_NB_W / P_useful_w)
            for beam in useful_beams
        )
    else:
        avg_fast_energy = 0.0
    n_thermal_axis = n_thermal * pk_n
    pressure_pa = physics.compute_pressure(
        ne_axis, n_thermal_axis, nb0_total, Te0_keV, Ti0_keV, config.fast_ion_mean_pitch2, avg_fast_energy,
        config.density_peaking, p_te, temperature_peaking_i=p_ti,
    )
    beta_t = physics.compute_beta_t(pressure_pa, config.Bt0)
    # Second bounding case, always computed (not config-gated): pitch2=1.0, a purely parallel/passing
    # fast-ion population -- p_fast=(1-1)*u_fast=0, so the fast-ion term drops out of the pressure sum
    # entirely. Reuses compute_pressure/compute_beta_t exactly as-is, just called again with a different
    # pitch2 -- no new physics, and Te/Ti/n_thermal etc. don't depend on this choice at all.
    pressure_anisotropic_pa = physics.compute_pressure(
        ne_axis, n_thermal_axis, nb0_total, Te0_keV, Ti0_keV, 1.0, avg_fast_energy,
        config.density_peaking, p_te, temperature_peaking_i=p_ti,
    )
    beta_t_anisotropic = physics.compute_beta_t(pressure_anisotropic_pa, config.Bt0)
    u_thermal = physics.thermal_energy_density(
        ne_axis, n_thermal_axis, Te0_keV, Ti0_keV, config.density_peaking, p_te,
        temperature_peaking_i=p_ti,
    )
    u_fast = physics.fast_ion_energy_density(nb0_total, avg_fast_energy)
    R = u_fast / u_thermal if u_thermal > 0.0 else float("inf")

    return OperatingPoint(
        ne0_m3=ne0_m3, feasible=True, Te_keV=Te_keV, Ti_keV=Ti_keV,
        nb0_m3=nb0_total, nD0_m3=nD0, nT0_m3=nT0, n_thermal_m3=n_thermal,
        n_thermal_fraction=n_thermal / n_sum,
        P_e_w=P_e_w, P_i_w=P_i_w, P_ei_w=P_ei_w, P_alpha_w=extra_e_source_w + extra_i_source_w,
            P_aux_e_w=aux_e_source_w, P_aux_i_w=aux_i_source_w,
        Le=Le_list, Li=Li_list,
        P_NB_total_w=P_NB_total_w, P_shine_w=P_shine_w, P_capt_w=P_capt_w, f_capture=f_capture_list,
        P_orbit_loss_w=P_orbit_loss_w, f_orbit_loss=f_orbit_loss_list,
        P_cx_loss_w=P_cx_loss_w, f_cx_loss=f_cx_loss_list, cx_f_particle=cx_f_particle_list,
        cx_gamma=cx_gamma_list, cx_escape_probability=cx_escape_list, cx_n0_m3=cx_n0_list,
        P_useful_w=P_useful_w,
        pf_thermal_w=pf_thermal, pf_beam_w=pf_beam, pf_dt_w=pf_dt,
        pf_dd_thermal_w=pf_dd_thermal, pf_dd_beam_w=pf_dd_beam, pf_dd_w=pf_dd,
        pf_total_w=pf_total, neutron_rate_s=neutron_rate,
        neutron_rate_thermal_s=neutron_rate_thermal, neutron_rate_beam_s=neutron_rate_beam,
        neutron_rate_bb_s=neutron_rate_bb, pf_bb_dt_w=pf_bb_dt, pf_bb_dd_w=pf_bb_dd,
        v_phi_m_s=v_phi_m_s, torque_total_Nm=torque_total_Nm,
        Te0_keV=Te0_keV, Ti0_keV=Ti0_keV,
        pressure_pa=pressure_pa, beta_t=beta_t,
        pressure_anisotropic_pa=pressure_anisotropic_pa, beta_t_anisotropic=beta_t_anisotropic,
        avg_fast_energy_keV=avg_fast_energy, R_fast_thermal=R, tau_E_s=tau_E_s, tau_Ei_s=tau_Ei_s,
    )
