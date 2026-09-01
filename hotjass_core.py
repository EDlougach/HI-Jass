from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List

import numpy as np


@dataclass
class PlasmaParams:
    major_radius: float = 6.2
    minor_radius: float = 1.8
    central_density: float = 2.5e20
    density_peaking: float = 1.2
    central_temperature: float = 8.0
    temp_peaking: float = 1.4
    toroidal_field: float = 2.5
    plasma_current: float = 1.2e6


@dataclass
class NBIParams:
    injected_power: float = 5.0
    beam_energy: float = 80.0
    beam_species: str = "D"
    beam_width: float = 0.45
    beam_shift: float = 0.15
    injection_angle_deg: float = 25.0
    shine_through_fraction: float = 0.10
    slowing_down_time: float = 0.15


class HotJassModel:
    """A compact model layer that produces representative plasma and NBI profiles.

    This is intentionally lightweight and keeps the GUI logic separate from the
    numerical calculations. It is suitable as a scaffold for later replacement with
    a more complete HotJass implementation.
    """

    def __init__(self, plasma: PlasmaParams | None = None, nbi: NBIParams | None = None):
        self.plasma = plasma or PlasmaParams()
        self.nbi = nbi or NBIParams()

    def rho_grid(self, n_points: int = 200) -> np.ndarray:
        return np.linspace(0.0, 1.0, n_points)

    def density_profile(self, rho: np.ndarray) -> np.ndarray:
        exponent = 2.0 * self.plasma.density_peaking
        return self.plasma.central_density * (1.0 - rho ** 2) ** exponent

    def temperature_profile(self, rho: np.ndarray) -> np.ndarray:
        exponent = 2.0 * self.plasma.temp_peaking
        return self.plasma.central_temperature * (1.0 - rho ** 2) ** exponent

    def fast_ion_profile(self, rho: np.ndarray) -> np.ndarray:
        width = self.nbi.beam_width
        shift = self.nbi.beam_shift
        gaussian = np.exp(-0.5 * ((rho - shift) / width) ** 2)
        return self.nbi.injected_power * gaussian / (np.sqrt(2.0 * np.pi) * width)

    def compute_all(self, n_points: int = 200) -> Dict[str, np.ndarray]:
        rho = self.rho_grid(n_points)
        return {
            "rho": rho,
            "density": self.density_profile(rho),
            "temperature": self.temperature_profile(rho),
            "fast_ions": self.fast_ion_profile(rho),
        }

    def summary(self) -> Dict[str, float]:
        effective_power = self.nbi.injected_power * (1.0 - self.nbi.shine_through_fraction)
        return {
            "line_average_density": self.plasma.central_density * 0.65,
            "central_temperature_keV": self.plasma.central_temperature,
            "beam_power_to_plasma": effective_power,
            "injected_energy_keV": self.nbi.beam_energy,
            "toroidal_field_T": self.plasma.toroidal_field,
            "plasma_current_MA": self.plasma.plasma_current / 1e6,
        }
