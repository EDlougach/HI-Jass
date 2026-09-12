from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List

import numpy as np

from hotjass import physics
from hotjass.physics import TokamakGeometry, captured_power_fraction, equipartition_power_w, thermalization_time
from hotjass.solve import BeamSpec, TokamakConfig, solve_operating_point


@dataclass
class PlasmaParams:
    major_radius: float = 0.65
    minor_radius: float = 0.35
    elongation: float = 2.2
    triangularity: float = -0.35
    central_density: float = 1.0e20
    n_e_min: float = 1.0e19
    n_e_max: float = 1.5e20
    density_peaking: float = 0.1
    temp_peaking: float = 1.0          # electron temperature profile exponent
    temp_peaking_i: float = -1.0       # ion T profile exponent; < 0 -> same as temp_peaking
    centrepost_radius: float = -1.0    # central-column outer radius [m]; < 0 -> R0 - a
    effective_charge: float = 2.0
    toroidal_field: float = 1.5
    plasma_current: float = 1.5e6
    deuterium_fraction: float = 0.2
    tritium_fraction: float = 0.8
    tauE_e: float = 0.15
    tauE_i: float = 0.15
    alpha_heating: bool = False  # fusion alpha self-heating fed back into the T_e/T_i balance
    f_alpha: float = 1.0  # fraction of alpha power confined & thermalised
    p_ecrh_MW: float = 2.0  # ECRH power (mostly to electrons)
    ecrh_f_e: float = 1.0  # electron fraction of ECRH power (rest to ions)
    p_icrh_MW: float = 0.0  # ICRH power
    icrh_f_e: float = 0.5  # electron fraction of ICRH power
    icrh_f_i: float = 0.5  # ion fraction of ICRH power (f_e + f_i need not sum to 1)
    tau_Ee_mode: str = "fixed"
    tau_Ei_mode: str = "fixed"
    enable_orbit_loss: bool = True
    orbit_loss_co_current: bool = True
    orbit_model: str = "st_pitch"  # "large_aspect" | "st_meanshift" | "st_pitch"
    profile_averaging: bool = False  # treat central_density as ON-AXIS; balance runs on <n_e>
    cx_loss_fraction: float = 0.1
    cx_model: str = "manual_fraction"  # "manual_fraction" | "manual_n0" | "penetration"
    cx_n0_over_ne: float = 1.0e-5      # "manual_n0": uniform background-neutral n0 = ratio * ne0
    cx_n0_lcfs_over_ne: float = 0.02   # "penetration": edge n0_LCFS/ne_LCFS boundary ratio
    rotation_model: str = "off"        # "off" | "manual" | "momentum_balance"
    manual_v_phi_m_s: float = 0.0      # "manual": bulk toroidal rotation velocity (positive = co-current)
    tau_phi_over_tauEi: float = 1.0    # "momentum_balance": tau_phi = this * tau_Ei (no validated tau_phi scaling exists)
    enable_beam_beam: bool = False     # reduced monoenergetic beam-beam fusion between the first two NBI sources
    enable_equipartition: bool = True


@dataclass
class BeamParams:
    species: str = "D"
    power_MW: float = 6.0
    beam_energy_keV: float = 80.0
    beam_width: float = 0.35
    beam_shift: float = 0.12
    injection_angle_deg: float = 25.0
    shine_through_fraction: float = 0.08
    tangent_R_m: float | None = None
    tangent_Z_m: float = 0.0
    shine_through_model: str = "manual"
    manual_shine_through_fraction: float = 0.01
    co_current: bool = True  # this beamline's injection direction (co- vs counter-Ip)


class HotJassModel:
    """A compact HotJass-inspired 0D model with a low-aspect plasma baseline."""

    def __init__(self, plasma: PlasmaParams | None = None, beams: List[BeamParams] | None = None):
        self.plasma = plasma or PlasmaParams()
        self.beams = beams or [
            BeamParams(species="D", power_MW=10.0, beam_energy_keV=120.0, tangent_R_m=0.5,
                       shine_through_model="janev_suzuki", co_current=True),
            BeamParams(species="T", power_MW=0.1, beam_energy_keV=180.0, tangent_R_m=0.6,
                       shine_through_model="janev_suzuki", co_current=False),
        ]

    def rho_grid(self, n_points: int = 200) -> np.ndarray:
        return np.linspace(0.0, 1.0, n_points)

    def _tokamak_config(self) -> TokamakConfig:
        return TokamakConfig(
            geometry=TokamakGeometry(
                self.plasma.major_radius, self.plasma.minor_radius,
                self.plasma.elongation, self.plasma.triangularity,
            ),
            Bt0=self.plasma.toroidal_field,
            Ip_MA=self.plasma.plasma_current / 1.0e6,
            Zeff=self.plasma.effective_charge,
            mix_D=self.plasma.deuterium_fraction,
            mix_T=self.plasma.tritium_fraction,
            density_peaking=self.plasma.density_peaking,
            temperature_peaking=self.plasma.temp_peaking,
            temperature_peaking_i=self.plasma.temp_peaking_i,
            centrepost_radius_m=self.plasma.centrepost_radius,
            tau_Ee_mode=self.plasma.tau_Ee_mode,
            tau_Ei_mode=self.plasma.tau_Ei_mode,
            enable_orbit_loss=self.plasma.enable_orbit_loss,
            orbit_loss_co_current=self.plasma.orbit_loss_co_current,
            orbit_model=self.plasma.orbit_model,
            profile_averaging=self.plasma.profile_averaging,
            cx_loss_fraction=self.plasma.cx_loss_fraction,
            cx_model=self.plasma.cx_model,
            cx_n0_over_ne=self.plasma.cx_n0_over_ne,
            cx_n0_lcfs_over_ne=self.plasma.cx_n0_lcfs_over_ne,
            rotation_model=self.plasma.rotation_model,
            manual_v_phi_m_s=self.plasma.manual_v_phi_m_s,
            tau_phi_over_tauEi=self.plasma.tau_phi_over_tauEi,
            enable_beam_beam=self.plasma.enable_beam_beam,
            enable_equipartition=self.plasma.enable_equipartition,
            enable_alpha_heating=self.plasma.alpha_heating,
            f_alpha=self.plasma.f_alpha,
        )

    def _aux_powers_w(self) -> tuple[float, float]:
        """Prescribed auxiliary heating (ECRH + ICRH) split into (electron, ion)
        power [W]. ECRH is mostly to electrons (ecrh_f_e); ICRH uses independent
        icrh_f_e / icrh_f_i (their sum need not be 1 -- any remainder is taken as
        power that does not thermalise, e.g. a fast-ion tail or direct loss).
        """
        p = self.plasma
        ecrh = max(p.p_ecrh_MW, 0.0) * 1.0e6
        icrh = max(p.p_icrh_MW, 0.0) * 1.0e6
        aux_e = ecrh * p.ecrh_f_e + icrh * p.icrh_f_e
        aux_i = ecrh * (1.0 - p.ecrh_f_e) + icrh * p.icrh_f_i
        return aux_e, aux_i

    def _cx_fraction_dicts(self, ne0_m3: float, Te_keV: float, Ti_keV: float,
                            config: TokamakConfig, beams) -> List[dict]:
        """Per-beam CX-loss diagnostic dicts (physics.cx_loss_fractions(), plus
        escape_probability/n0_m3) for config.cx_model in ("manual_n0",
        "penetration") -- docs/CX_model.tex Sec.4's other two selector modes
        ("manual_fraction" never calls this; solve_operating_point's own
        config.cx_loss_fraction path handles it, untouched).
        """
        a = config.geometry.minor_radius
        if config.cx_model == "manual_n0":
            n0 = max(config.cx_n0_over_ne, 0.0) * ne0_m3
        elif config.cx_model == "penetration":
            rho = self.rho_grid()
            ne_profile = self.density_profile(rho)
            # ne_profile hits exactly 0 at rho=1 whenever density_peaking>0 --
            # use a near-edge shell (rho~0.95) as the physically meaningful
            # "LCFS" density instead of that zero. Also used by the GUI's
            # Profiles-tab n0(rho) panel (_render_profiles) -- keep both in
            # sync if this fraction ever changes.
            edge_idx = max(int(0.95 * (len(rho) - 1)), 0)
            n0_lcfs = max(config.cx_n0_lcfs_over_ne, 0.0) * ne_profile[edge_idx]
            n0_profile = physics.neutral_penetration_profile(
                rho, ne_profile, Te_keV, Ti_keV, a, n0_lcfs)
            n0 = physics.profile_volume_average(n0_profile, rho)
        else:
            n0 = 0.0
        out = []
        for beam in beams:
            d = physics.cx_loss_fractions(n0, ne0_m3, Te_keV, beam.Eb_keV, beam.species)
            d["escape_probability"] = physics.cx_escape_probability(
                beam.Eb_keV, beam.species, ne0_m3, a,
                stopping_model="riviere", Te_keV=Te_keV, Zeff=config.Zeff)
            d["n0_m3"] = n0
            out.append(d)
        return out

    def _solve_point(self, ne_m3: float, beams, tau_e: float, tau_i: float, config: TokamakConfig):
        """One operating point. Prescribed auxiliary heating (ECRH/ICRH) is added
        to the electron/ion balances via solve_operating_point's aux source terms.
        When config.enable_alpha_heating, an under-relaxed Picard closes the fusion
        alpha term P_alpha = f_alpha*(E_alpha/E_fus)*P_fusion on top of that; the
        alpha-off path is a single plain solve_operating_point call.

        When config.cx_model != "manual_fraction" (physics-based CX loss), an
        OUTER under-relaxed Picard loop closes the second fixed point this
        introduces: the CX-loss integral needs Te (via critical_energy_keV) --
        and, for "penetration", Ti too -- while Te/Ti themselves depend on the
        CX-reduced P_useful. solve_operating_point() stays a single
        non-iterative solve throughout (same principle as the alpha loop);
        this iterates around it instead.
        """
        aux_e, aux_i = self._aux_powers_w()

        def _run(cx_override):
            op = solve_operating_point(ne_m3, beams, tau_e, config, tau_i,
                                       aux_e_source_w=aux_e, aux_i_source_w=aux_i,
                                       cx_loss_fraction_override=cx_override)
            if not (config.enable_alpha_heating and op.feasible):
                return op
            ratio = config.f_alpha * (physics.E_ALPHA_MEV / physics.E_FUSION_MEV)
            p_alpha = 0.0
            for _ in range(30):
                # alphas come from D-T only -- not the D-D channel
                target = ratio * op.pf_dt_w
                step = 0.5 * (target - p_alpha)
                if abs(step) <= 1.0e-3 * max(target, 1.0e-9):
                    break
                p_alpha += step
                le = physics.alpha_electron_heating_fraction(max(op.Te_keV, 1.0e-3))
                trial = solve_operating_point(
                    ne_m3, beams, tau_e, config, tau_i,
                    extra_e_source_w=le * p_alpha, extra_i_source_w=(1.0 - le) * p_alpha,
                    aux_e_source_w=aux_e, aux_i_source_w=aux_i,
                    cx_loss_fraction_override=cx_override,
                )
                if not trial.feasible:
                    return trial
                op = trial
            return op

        if config.cx_model == "manual_fraction" or not beams:
            return _run(None)

        op = _run(None)  # seed pass: Te/Ti with no CX loss yet
        if not op.feasible:
            return op
        for _ in range(8):
            cx = self._cx_fraction_dicts(ne_m3, op.Te_keV, op.Ti_keV, config, beams)
            trial = _run(cx)
            if not trial.feasible:
                return trial
            converged = (
                op.Te_keV is not None and trial.Te_keV is not None
                and abs(trial.Te_keV - op.Te_keV) <= 1.0e-4 * max(abs(op.Te_keV), 1e-9)
            )
            op = trial
            if converged:
                break
        return op

    def _beam_specs(self) -> List[BeamSpec]:
        return [
            BeamSpec(
                beam.beam_energy_keV, beam.power_MW * 1.0e6, beam.species.upper(),
                beam.tangent_R_m, beam.tangent_Z_m,
                beam.shine_through_model, beam.manual_shine_through_fraction,
                bool(beam.co_current),
            )
            for beam in self.beams
        ]

    def operating_point(self, n_e: float | None = None):
        """Solve one self-consistent operating point (default: at central_density)."""
        density = self.plasma.central_density if n_e is None else float(n_e)
        return self._solve_point(
            density, self._beam_specs(),
            max(self.plasma.tauE_e, 1.0e-6), max(self.plasma.tauE_i, 1.0e-6),
            self._tokamak_config(),
        )

    def density_profile(self, rho: np.ndarray) -> np.ndarray:
        exponent = 2.0 * self.plasma.density_peaking
        return self.plasma.central_density * (1.0 - rho ** 2) ** exponent

    def density_scan(self, n_e: np.ndarray | None = None) -> Dict[str, np.ndarray]:
        """Solve the compact electron/ion 0D balance over central density."""
        densities = np.asarray(
            n_e if n_e is not None else np.linspace(self.plasma.n_e_min, self.plasma.n_e_max, 21), dtype=float
        )
        volume = self.plasma_volume()
        config = self._tokamak_config()
        beams = self._beam_specs()
        tau_e = max(self.plasma.tauE_e, 1.0e-6)
        tau_i = max(self.plasma.tauE_i, 1.0e-6)
        points = [self._solve_point(float(density), beams, tau_e, tau_i, config) for density in densities]

        def values(attribute: str, scale: float = 1.0, infeasible: float = np.nan) -> np.ndarray:
            return np.asarray([getattr(point, attribute) * scale if point.feasible else infeasible for point in points], dtype=float)

        te = values("Te_keV")
        ti = values("Ti_keV")
        p_e = values("P_e_w", 1.0e-6)
        p_i = values("P_i_w", 1.0e-6)
        p_nb = values("P_NB_total_w", 1.0e-6)
        p_useful = values("P_useful_w", 1.0e-6)
        pf_tot_mw = values("pf_total_w", 1.0e-6)
        with np.errstate(divide="ignore", invalid="ignore"):
            q_values = np.where(p_nb > 0.0, pf_tot_mw / np.where(p_nb > 0.0, p_nb, 1.0), np.nan)
        tau_s_values = []
        tau_ie_values = []
        p_ie_values = []
        for density, point in zip(densities, points):
            if not point.feasible or point.P_useful_w <= 0.0:
                tau_s_values.append(np.nan)
                tau_ie_values.append(np.nan)
                p_ie_values.append(np.nan)
                continue
            useful_power = 0.0
            weighted_tau = 0.0
            for beam in beams:
                useful = beam.P_NB_W * captured_power_fraction(
                    float(density), beam.Eb_keV, beam.species, 2.0 * self.plasma.minor_radius,
                    self.plasma.density_peaking, geometry=TokamakGeometry(
                        self.plasma.major_radius, self.plasma.minor_radius, self.plasma.elongation,
                        self.plasma.triangularity,
                    ), tangent_R_m=beam.tangent_R_m, tangent_Z_m=beam.tangent_Z_m,
                    model=beam.shine_through_model,
                    manual_shine_through_fraction=beam.manual_shine_through_fraction,
                    Zeff=self.plasma.effective_charge,
                )
                tau = thermalization_time(float(density), float(point.Te_keV), beam.Eb_keV, beam.species)
                weighted_tau += useful * tau
                useful_power += useful
            tau_s_values.append(weighted_tau / max(useful_power, 1.0e-30))
            p_ei = sum(
                equipartition_power_w(
                    fraction * point.n_thermal_m3, point.Te_keV, point.Ti_keV,
                    point.ne0_m3, species, volume,
                )
                for fraction, species in (
                    (self.plasma.deuterium_fraction, "D"),
                    (self.plasma.tritium_fraction, "T"),
                )
            )
            p_ie = -p_ei
            p_ie_values.append(p_ie)
            numerator = 1.5 * point.n_thermal_m3 * abs(point.Ti_keV - point.Te_keV) * 1.0e3 * 1.602176634e-19 * volume
            tau_ie_values.append(numerator / abs(p_ie) if abs(p_ie) > 1.0e-30 else np.nan)
        tau_s_values = np.asarray(tau_s_values)
        tau_ie_values = np.asarray(tau_ie_values)
        p_ie_values = np.asarray(p_ie_values)
        feasible_mask = np.asarray([point.feasible for point in points], dtype=bool)
        avg_energy = values("avg_fast_energy_keV")
        thermal_energy = 1.5 * (values("ne0_m3") * te + values("n_thermal_m3") * ti) * 1e3 * 1.602176634e-19
        fast_energy = values("nb0_m3") * avg_energy * 1e3 * 1.602176634e-19

        result = {
            "n_e": densities, "Te": te, "Ti": ti, "P_e": p_e, "P_i": p_i,
                "Pi_e": p_ie_values * 1.0e-6, "P_shine-through": values("P_shine_w", 1.0e-6),
                "P_alpha": values("P_alpha_w", 1.0e-6),
                "P_aux_e": values("P_aux_e_w", 1.0e-6), "P_aux_i": values("P_aux_i_w", 1.0e-6),
                "n_D": values("nD0_m3"), "n_T": values("nT0_m3"), "n_b": values("nb0_m3"),
                "Pf_tot": pf_tot_mw, "Pf_th": values("pf_thermal_w", 1.0e-6), "Pf_b": values("pf_beam_w", 1.0e-6),
                "Pf_bb": values("pf_bb_dt_w", 1.0e-6) + values("pf_bb_dd_w", 1.0e-6),
                "Pf_DT": values("pf_dt_w", 1.0e-6), "Pf_DD": values("pf_dd_w", 1.0e-6),
                # Neutron rates (Y, not R -- R is already the u_fast/u_thermal ratio elsewhere).
                "Y_neutron": values("neutron_rate_s"),
                "Y_neutron_th": values("neutron_rate_thermal_s"), "Y_neutron_b": values("neutron_rate_beam_s"),
                "Y_neutron_bb": values("neutron_rate_bb_s"),  # 0 unless plasma.enable_beam_beam
                "P_orbit": values("P_orbit_loss_w", 1.0e-6), "P_cx": values("P_cx_loss_w", 1.0e-6),
                "P_useful": p_useful, "P_lost": p_nb - p_useful, "Q": q_values,
                "Te0": values("Te0_keV"), "Ti0": values("Ti0_keV"),
                "E_fast": avg_energy, "tau_S": tau_s_values,
                "tauE_e": values("tau_E_s"), "tauE_i": values("tau_Ei_s"), "tau_IE": tau_ie_values, "R": fast_energy / np.maximum(thermal_energy, 1.0e-30),
                "Pr_th": values("pressure_pa") - (1.0 - 1.0 / 3.0) * fast_energy, "Pr_fast": (1.0 - 1.0 / 3.0) * fast_energy, "beta_T": values("beta_t") * 100.0,
        }
        result["n_e_requested"] = densities.copy()
        result["valid_mask"] = feasible_mask
        valid_densities = densities[feasible_mask]
        result["n_e_valid_min"] = np.asarray([valid_densities.min() if valid_densities.size else np.nan])
        result["n_e_valid_max"] = np.asarray([valid_densities.max() if valid_densities.size else np.nan])
        return result

    def temperature_profile(self, rho: np.ndarray, central_temperature: float,
                            ion: bool = False) -> np.ndarray:
        p = self.plasma.temp_peaking
        if ion and self.plasma.temp_peaking_i >= 0.0:
            p = self.plasma.temp_peaking_i
        return central_temperature * (1.0 - rho ** 2) ** (2.0 * p)

    def fast_ion_profile(self, rho: np.ndarray) -> np.ndarray:
        return np.full_like(rho, sum(beam.power_MW for beam in self.beams), dtype=float)

    def plasma_volume(self) -> float:
        return TokamakGeometry(
            self.plasma.major_radius,
            self.plasma.minor_radius,
            self.plasma.elongation,
            self.plasma.triangularity,
        ).volume()

    def power_balance(self) -> Dict[str, float]:
        scan = self.density_scan(np.array([self.plasma.central_density]))

        return {
            "electron_power_MW": float(scan["P_e"][0]),
            "ion_power_MW": float(scan["P_i"][0]),
            "fusion_power_MW": float(scan["Pf_tot"][0]),
            "total_beam_power_MW": float(sum(b.power_MW for b in self.beams)),
        }

    def compute_all(self, n_points: int = 200) -> Dict[str, np.ndarray]:
        rho = self.rho_grid(n_points)
        return {
            "rho": rho,
            "density": self.density_profile(rho),
            "temperature": self.temperature_profile(rho, self.density_scan(np.array([self.plasma.central_density]))["Te"][0]),
            "fast_ions": self.fast_ion_profile(rho),
        }

    def summary(self) -> Dict[str, float]:
        power = self.power_balance()
        scan = self.density_scan(np.array([self.plasma.central_density]))
        effective_power = power["electron_power_MW"] + power["ion_power_MW"]
        return {
            "line_average_density": self.plasma.central_density * 0.65,
            "central_temperature_e_keV": float(scan["Te"][0]),
            "central_temperature_i_keV": float(scan["Ti"][0]),
            "density_peaking": self.plasma.density_peaking,
            "temperature_peaking": self.plasma.temp_peaking,
            "deuterium_fraction": self.plasma.deuterium_fraction,
            "tritium_fraction": self.plasma.tritium_fraction,
            "elongation": self.plasma.elongation,
            "triangularity": self.plasma.triangularity,
            "tauE_e": self.plasma.tauE_e,
            "tauE_i": self.plasma.tauE_i,
            "beam_power_to_plasma": effective_power,
            "beam_1_power_MW": self.beams[0].power_MW,
            "beam_2_power_MW": self.beams[1].power_MW,
            "beam_1_energy_keV": self.beams[0].beam_energy_keV,
            "beam_2_energy_keV": self.beams[1].beam_energy_keV,
            "beam_1_species": self.beams[0].species,
            "beam_2_species": self.beams[1].species,
            "fusion_power_MW": power["fusion_power_MW"],
            "electron_power_loss_MW": power["electron_power_MW"],
            "ion_power_loss_MW": power["ion_power_MW"],
            "toroidal_field_T": self.plasma.toroidal_field,
            "plasma_current_MA": self.plasma.plasma_current / 1e6,
            "plasma_volume_m3": self.plasma_volume(),
        }
