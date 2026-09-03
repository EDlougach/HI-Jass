from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List

import numpy as np

from hotjass.physics import TokamakGeometry, captured_power_fraction, equipartition_power_w, thermalization_time
from hotjass.solve import BeamSpec, TokamakConfig, solve_operating_point


@dataclass
class PlasmaParams:
    major_radius: float = 0.65
    minor_radius: float = 0.35
    elongation: float = 2.2
    triangularity: float = -0.35
    central_density: float = 1.5e20
    n_e_min: float = 1.0e19
    n_e_max: float = 1.0e20
    density_peaking: float = 0.1
    temp_peaking: float = 1.0
    effective_charge: float = 2.0
    toroidal_field: float = 1.5
    plasma_current: float = 1.5e6
    deuterium_fraction: float = 0.5
    tritium_fraction: float = 0.5
    tauE_e: float = 0.02
    tauE_i: float = 0.05
    alpha_heating: bool = False


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


class HotJassModel:
    """A compact HotJass-inspired 0D model with a low-aspect plasma baseline."""

    def __init__(self, plasma: PlasmaParams | None = None, beams: List[BeamParams] | None = None):
        self.plasma = plasma or PlasmaParams()
        self.beams = beams or [
            BeamParams(species="D", power_MW=5.0, beam_energy_keV=120.0),
            BeamParams(species="T", power_MW=5.0, beam_energy_keV=180.0),
        ]

    def rho_grid(self, n_points: int = 200) -> np.ndarray:
        return np.linspace(0.0, 1.0, n_points)

    def density_profile(self, rho: np.ndarray) -> np.ndarray:
        exponent = 2.0 * self.plasma.density_peaking
        return self.plasma.central_density * (1.0 - rho ** 2) ** exponent

    def density_scan(self, n_e: np.ndarray | None = None) -> Dict[str, np.ndarray]:
        """Solve the compact electron/ion 0D balance over central density."""
        densities = np.asarray(
            n_e if n_e is not None else np.linspace(self.plasma.n_e_min, self.plasma.n_e_max, 21), dtype=float
        )
        volume = self.plasma_volume()
        config = TokamakConfig(
            geometry=TokamakGeometry(self.plasma.major_radius, self.plasma.minor_radius, self.plasma.elongation, self.plasma.triangularity),
            Bt0=self.plasma.toroidal_field,
            Ip_MA=self.plasma.plasma_current / 1.0e6,
            Zeff=self.plasma.effective_charge,
            mix_D=self.plasma.deuterium_fraction,
            mix_T=self.plasma.tritium_fraction,
            density_peaking=self.plasma.density_peaking,
            temperature_peaking=self.plasma.temp_peaking,
        )
        beams = [
            BeamSpec(
                beam.beam_energy_keV, beam.power_MW * 1.0e6, beam.species.upper(),
                beam.tangent_R_m, beam.tangent_Z_m,
                beam.shine_through_model, beam.manual_shine_through_fraction,
            )
            for beam in self.beams
        ]
        points = [solve_operating_point(float(density), beams, max(self.plasma.tauE_e, 1.0e-6), config, max(self.plasma.tauE_i, 1.0e-6)) for density in densities]

        def values(attribute: str, scale: float = 1.0, infeasible: float = np.nan) -> np.ndarray:
            return np.asarray([getattr(point, attribute) * scale if point.feasible else infeasible for point in points], dtype=float)

        te = values("Te_keV")
        ti = values("Ti_keV")
        p_e = values("P_e_w", 1.0e-6)
        p_i = values("P_i_w", 1.0e-6)
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
                "n_D": values("nD0_m3"), "n_T": values("nT0_m3"), "n_b": values("nb0_m3"),
                "Pf_tot": values("pf_total_w", 1.0e-6), "Pf_th": values("pf_thermal_w", 1.0e-6), "Pf_b": values("pf_beam_w", 1.0e-6),
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

    def temperature_profile(self, rho: np.ndarray, central_temperature: float) -> np.ndarray:
        exponent = 2.0 * self.plasma.temp_peaking
        return central_temperature * (1.0 - rho ** 2) ** exponent

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
