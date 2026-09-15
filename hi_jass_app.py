"""HI-Jass desktop GUI -- v1 layout (sidebar inputs + tabbed results + run modes).

Same 0-D physics as before (hotjass_core / hotjass.solve); the interface is
reorganised per docs/HI-Jass-GUI.pdf section 8: a collapsible input rail on the
left, a run-mode toggle (Operating point / Scan), and a tabbed results area on
the right that surfaces quantities the solver already computes (power balance,
fast-ion f(E), assumptions/validity) alongside the legacy density scans.
"""

from __future__ import annotations

import dataclasses
import json
import os
import queue
import subprocess
import threading
import traceback
import webbrowser
from datetime import datetime
from pathlib import Path
from tkinter import filedialog, messagebox

import customtkinter as ctk
import matplotlib
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from matplotlib.sankey import Sankey
from PIL import Image

from hotjass import physics
from hotjass_core import HotJassModel

matplotlib.use("TkAgg")

SETTINGS_PATH = Path.home() / ".hi_jass" / "settings.json"

CONFINEMENT_MODES = {
    # v1 scope: 4 choices only. Fixed allows tau_Ee != tau_Ei (manual per-
    # channel input); every scaling below is "natural" and always sets
    # tau_Ee = tau_Ei to that scaling's single computed value -- the
    # "X e / neoclassical i" hybrids and the arbitrary-A neoclassical
    # variants are deliberately not offered here for now (memory note
    # "confinement-scaling-p-loss-definition" territory: the neoclassical
    # ion channel is a separate, still-open piece of the picture).
    "Fixed tauE (input)": ("fixed", "fixed"),
    "Kaye NSTX L-mode": ("kaye_nstx_lmode", "kaye_nstx_lmode"),
    "Kaye NSTX H-mode": ("kaye_nstx_hmode", "kaye_nstx_hmode"),
    "IPB98(y,2) ELMy H-mode": ("iter98y2", "iter98y2"),
}

ORBIT_MODELS = {
    "Large-aspect (q* rho_Li)": "large_aspect",
    "ST orbits - mean-shift (arbitrary A)": "st_meanshift",
    "ST orbits - pitch-resolved": "st_pitch",
}
CX_MODELS = {
    "Manual fraction": "manual_fraction",
    "Manual n0/ne": "manual_n0",
    "Penetration (n0_LCFS/ne)": "penetration",
}
CX_MODELS_INV = {v: k for k, v in CX_MODELS.items()}
ROTATION_MODELS = {
    "Off": "off",
    "Manual v_phi": "manual",
    "Momentum balance": "momentum_balance",
}
SHINE_LABEL_TO_MODEL = {"Riviere": "riviere", "Janev": "janev_suzuki", "Manual": "manual"}
SHINE_MODEL_TO_LABEL = {v: k for k, v in SHINE_LABEL_TO_MODEL.items()}

# Accessibility: an "Aa" button cycles through these UI-scale presets,
# rescaling both CustomTkinter widgets (labels/entries/buttons -- the
# built-in ctk.set_widget_scaling mechanism) and matplotlib plot text
# (which doesn't follow ctk scaling on its own, so font.size is scaled by
# the same factor and the visible tabs are re-rendered on each change).
UI_SCALES = [1.0, 1.15, 1.3, 1.5]
_BASE_MPL_FONT_SIZE = matplotlib.rcParams["font.size"]


def _scholar(query: str) -> str:
    """A Google Scholar search URL for a citation string -- used instead of
    hand-typed DOIs so every link resolves to the real paper (top hit) and
    nothing here is a fabricated identifier."""
    import urllib.parse
    return "https://scholar.google.com/scholar?q=" + urllib.parse.quote(query)


# key -> (citation, url).  Grouped/filtered by active model in _render_references.
REFERENCES = {
    # --- beam stopping / deposition ---
    "riviere": ("Rivière, Nucl. Fusion 11 (1971) 363 - attenuation of fast neutral hydrogen beams",
                _scholar("Riviere 1971 Nuclear Fusion 11 363 penetration neutral beam")),
    "janev_suzuki": ("Janev, Boley & Post, Nucl. Fusion 29 (1989) 2125 - beam stopping cross-sections",
                     _scholar("Janev Boley Post 1989 Nuclear Fusion 29 2125 beam stopping")),
    "suzuki": ("Suzuki et al., Plasma Phys. Control. Fusion 40 (1998) 2097 - beam-stopping fit incl. excited states",
               _scholar("Suzuki 1998 Plasma Physics Controlled Fusion beam stopping")),
    # --- confinement scalings ---
    "iter98y2": ("ITER Physics Basis, Ch. 2, Nucl. Fusion 39 (1999) 2175 - IPB98(y,2) ELMy H-mode scaling",
                 _scholar("ITER Physics Basis 1999 Nuclear Fusion 39 2175 confinement IPB98(y,2)")),
    "kaye_nstx": ("Kaye et al., Nucl. Fusion 46 (2006) 848 - NSTX energy confinement scaling (H- and L-mode)",
                  _scholar("Kaye 2006 Nuclear Fusion 46 848 NSTX confinement scaling")),
    "neoclassical_ch": ("Chang & Hinton, Phys. Fluids 25 (1982) 1493 - neoclassical ion thermal transport",
                        _scholar("Chang Hinton 1982 Physics of Fluids 25 1493 neoclassical")),
    "helander_arbA": ("Helander, Phys. Plasmas 7 (2000) 3999 - neoclassical transport, arbitrary aspect ratio & collisionality",
                      _scholar("Helander 2000 Physics of Plasmas 3999 neoclassical arbitrary aspect ratio")),
    "helander_sigmar": ("Helander & Sigmar, Collisional Transport in Magnetized Plasmas (CUP, 2002)",
                        _scholar("Helander Sigmar Collisional Transport in Magnetized Plasmas")),
    "linliu_miller": ("Lin-Liu & Miller, Phys. Plasmas 2 (1995) 1666 - trapped-particle fraction at arbitrary A",
                      _scholar("Lin-Liu Miller 1995 Physics of Plasmas 2 1666 trapped particle fraction")),
    "hinton_wiley": ("Hinton, Wiley, Düchs, Furth & Rutherford, Phys. Rev. Lett. 29 (1972) 698 - finite-orbit neoclassical ion transport",
                     "https://doi.org/10.1103/PhysRevLett.29.698"),
    "satake_fow": ("Satake et al., Phys. Plasmas 9 (2002) 3946 - finite-orbit-width neoclassical transport",
                   _scholar("Satake 2002 Physics of Plasmas finite orbit width neoclassical")),
    # --- first-orbit / fast-ion losses ---
    "akers_start": ("Akers et al., Nucl. Fusion 42 (2002) 122 - NBI heating & fast-ion losses in the START spherical tokamak",
                    _scholar("Akers 2002 Nuclear Fusion 42 122 START neutral beam spherical tokamak")),
    "gwb1981": ("Goldston, White & Boozer, Phys. Rev. Lett. 47 (1981) 1004 - confinement of high-energy trapped particles",
                "https://doi.org/10.1103/PhysRevLett.47.1004"),
    "goldston_rutherford": ("Goldston & Rutherford, Introduction to Plasma Physics (IOP, 1995) - Ch. 12, guiding-centre orbits",
                            _scholar("Goldston Rutherford Introduction to Plasma Physics 1995")),
    "wesson": ("Wesson, Tokamaks (4th ed., Oxford, 2011) - Sect. 3.10 particle orbits; neoclassical estimates",
               _scholar("Wesson Tokamaks 4th edition Oxford")),
    "uckan": ("Uckan & ITER Physics Group, ITER-TN-PH-8-6 (1988) / Uckan Fusion Technol. 14 (1988) 299 - q95 engineering formula",
              _scholar("Uckan 1988 ITER q95 engineering safety factor formula")),
    # --- fusion / plasma primitives ---
    "bosch_hale": ("Bosch & Hale, Nucl. Fusion 32 (1992) 611 - improved D-T fusion reactivity & cross-section",
                   "https://doi.org/10.1088/0029-5515/32/4/I07"),
    "nrl": ("NRL Plasma Formulary (2019 rev.) - collision rates, thermal equilibration",
            "https://www.nrl.navy.mil/News-Media/Publications/nrl-plasma-formulary/"),
    "stix": ("Stix, Plasma Phys. 14 (1972) 367 - heating of toroidal plasmas by neutral injection (L_e/L_i split)",
             _scholar("Stix 1972 Plasma Physics 14 367 heating toroidal plasmas neutral injection")),
    # --- charge-exchange loss (docs/CX_model.tex) ---
    "janev_smith_1993": ("Janev & Smith, Nucl. Fusion Suppl. 4 (1993) - CX cross-section analytic form (Sec. 2.3.1)",
                         _scholar("Janev Smith 1993 Nuclear Fusion Suppl 4 atomic plasma-material interaction data")),
    "swaczyna_2024": ("Swaczyna, Bzowski & Kubiak, arXiv:2411.13174 (2024) - refit CX cross-section (Eq. A1) used here",
                      "https://arxiv.org/abs/2411.13174"),
    "freeman_jones_1974": ("Freeman & Jones, Culham CLM-R137 (1974) - analytic atomic cross-section fits",
                          _scholar("Freeman Jones 1974 CLM-R137 analytic cross sections rate coefficients")),
    "pppl1280_cx": ("PPPL-1280, Sect. 6.2.2-6.2.3 - CX loss quoted as up to 10% of injected power (Zeff=1)",
                    _scholar("PPPL-1280 neutral beam injection heating tokamak")),
}

# Active-machine geometry references (keyed by preset name). Each device
# maps to a *list* of (citation, url) rows, so a device can carry more than
# one reference (e.g. TCV: the device paper plus its NBI upgrade paper).
MACHINE_REFERENCES = {
    "DANTE": [("DANTE - internal low-aspect design point; no external publication.", "")],
    "ITER": [("ITER Physics Basis, Ch. 1, Nucl. Fusion 39 (1999) 2137 - device description & parameters",
              _scholar("ITER Physics Basis 1999 Nuclear Fusion 39 2137 overview"))],
    "JET": [("Rebut, Bickerton & Keen, Nucl. Fusion 25 (1985) 1011 - the JET project & its prospects",
             _scholar("Rebut Bickerton Keen 1985 Nuclear Fusion 25 1011 JET project")),
            ("Ciric et al., Fusion Eng. Des. 82 (2007) 610 - JET neutral beam enhancement "
             "(2 injector boxes, up to 8 PINIs each, >34 MW design)",
             _scholar("Ciric 2007 Fusion Engineering Design 82 610 JET neutral beam enhancement")),
            ("Maggi et al., Nucl. Fusion 64 (2024) 112012 - JET DTE2 T/D-T overview; "
             "independently confirms pulse #99971's 26.5 MW NBI / 4 MW ICRH (N=1 D minority, "
             "29 MHz) / 15:85 D:T / 59 MJ record used by this preset",
             _scholar("Maggi 2024 Nuclear Fusion 64 112012 JET tritium deuterium-tritium overview")),
            ("Villari et al., Fusion Eng. Des. 217 (2025) 115133 - JET D-T nuclear operations "
             "overview (DTE2/DTE3 neutron yields, 14 MeV calibration to +/-6%)",
             _scholar("Villari 2025 Fusion Engineering Design 217 115133 JET deuterium tritium nuclear operations"))],
    "ST40": [("Gryaznevich et al., Nucl. Fusion 62 (2022) 042008 - ST40 compact high-field spherical tokamak",
              _scholar("Gryaznevich 2022 Nuclear Fusion ST40 spherical tokamak"))],
    "T-15MD": [("Khvostenko et al., Fusion Eng. Des. 146 (2019) 1108 - T-15MD tokamak construction",
                _scholar("Khvostenko 2019 Fusion Engineering Design T-15MD tokamak"))],
    "TCV": [("Hofmann et al., Plasma Phys. Control. Fusion 36 (1994) B277 - the TCV tokamak",
             _scholar("Hofmann 1994 Plasma Physics Controlled Fusion TCV tokamak")),
            ("Karpushov et al., Fusion Eng. Des. 187 (2023) 113384 - TCV second high-energy "
             "NBI (NBI-2, 1.0 MW/55 keV, counter-injected vs. NBI-1's 1.3 MW/28 keV, co) upgrade",
             _scholar("Karpushov 2023 Fusion Engineering Design 187 113384 TCV second neutral beam"))],
}

REPO_URL = "https://github.com/EDlougach/HI-Jass"
_FEEDBACK_USER = "eugenia.dlougach"
_FEEDBACK_HOST = "real-nbi.com"
_LOGO_PATH = Path(__file__).resolve().parent / "assets" / "NBR_logo_color.png"


@dataclasses.dataclass
class Result:
    ok: bool
    error: str = ""
    op: object = None
    scan: dict = None
    summary: dict = None
    volume_m3: float = 0.0


class CollapsibleSection(ctk.CTkFrame):
    """A titled frame whose body can be folded away by clicking the header."""

    def __init__(self, master, title: str, expanded: bool = True):
        super().__init__(master, fg_color=("gray92", "gray17"))
        self.grid_columnconfigure(0, weight=1)
        self._expanded = expanded
        self._title = title
        self.header = ctk.CTkButton(
            self, text=self._label(), anchor="w", fg_color="transparent",
            text_color=("gray10", "gray90"), hover_color=("gray85", "gray25"),
            command=self.toggle,
        )
        self.header.grid(row=0, column=0, sticky="ew", padx=2, pady=(2, 0))
        self.body = ctk.CTkFrame(self, fg_color="transparent")
        self.body.grid_columnconfigure(1, weight=1)
        if self._expanded:
            self.body.grid(row=1, column=0, sticky="ew", padx=6, pady=6)

    def _label(self) -> str:
        return ("▾  " if self._expanded else "▸  ") + self._title

    def toggle(self):
        # grid_forget() (not grid_remove()) -- see the matching comment in
        # HIJassApp._set_mode: CTkBaseClass replays a widget's last grid()
        # call on every UI-scale change, and grid_remove() (unlike
        # grid_forget()) doesn't clear that memory, so a collapsed section
        # would silently re-expand the next time the zoom button is pressed.
        self._expanded = not self._expanded
        self.header.configure(text=self._label())
        if self._expanded:
            self.body.grid(row=1, column=0, sticky="ew", padx=6, pady=6)
        else:
            self.body.grid_forget()


class HIJassApp(ctk.CTk):
    # Per-device Plasma-field overrides applied by _apply_preset(). Two
    # optional keys aren't Plasma text-entry fields and are special-cased in
    # _apply_preset(): "profile_averaging" (bool) sets the profile-corrected
    # 0-D checkbox, and "confinement_mode" (a CONFINEMENT_MODES label) sets
    # the Confinement dropdown to whichever scaling actually fits the
    # device's aspect ratio -- so switching devices doesn't leave e.g. an
    # ST-fit Kaye scaling active on a conventional-aspect machine.
    #
    # All 6 devices share the same peaking defaults (density_peaking=0.1,
    # temp_peaking(_i)=1.0) -- chosen after direct comparison across presets.
    # profile_averaging itself is per-device (user-tested, not aspect-ratio
    # derived): ON for ITER/JET/T-15MD, OFF for DANTE/ST40/TCV.
    PLASMA_PRESETS = {
        "DANTE": {
            "major_radius": 0.65, "minor_radius": 0.35, "elongation": 2.2,
            "triangularity": -0.35, "effective_charge": 2.0, "toroidal_field": 1.5,
            "plasma_current_MA": 1.5, "central_density": 1.0e20,
            "n_e_min": 1.0e19, "n_e_max": 1.5e20, "density_peaking": 0.1, "temp_peaking": 1.0,
            "temp_peaking_i": 1.0, "profile_averaging": False,
            "deuterium_fraction": 0.2, "tritium_fraction": 0.8, "tauE_e": 0.15, "tauE_i": 0.15,
            "confinement_mode": "Kaye NSTX H-mode",
        },
        "ITER": {
            "major_radius": 6.2, "minor_radius": 2.0, "elongation": 1.85,
            "triangularity": 0.33, "effective_charge": 1.7, "toroidal_field": 5.3,
            "plasma_current_MA": 15.0, "central_density": 1.0e20,
            "n_e_min": 5.0e19, "n_e_max": 1.5e20, "density_peaking": 0.1, "temp_peaking": 1.0,
            "temp_peaking_i": 1.0, "profile_averaging": True,
            "deuterium_fraction": 0.5, "tritium_fraction": 0.5, "tauE_e": 3.7, "tauE_i": 3.7,
            "confinement_mode": "IPB98(y,2) ELMy H-mode",
        },
        "JET": {
            "major_radius": 2.96, "minor_radius": 1.25, "elongation": 1.7,
            "triangularity": 0.32, "effective_charge": 1.5, "toroidal_field": 3.45,
            "plasma_current_MA": 4.0, "central_density": 6.0e19,
            "n_e_min": 2.0e19, "n_e_max": 1.0e20, "density_peaking": 0.1, "temp_peaking": 1.0,
            "temp_peaking_i": 1.0, "profile_averaging": True,
            "deuterium_fraction": 0.5, "tritium_fraction": 0.5, "tauE_e": 1.5, "tauE_i": 1.5,
            "confinement_mode": "IPB98(y,2) ELMy H-mode",
        },
        "ST40": {
            "major_radius": 0.45, "minor_radius": 0.30, "elongation": 1.8,
            "triangularity": 0.4, "effective_charge": 1.5, "toroidal_field": 3.0,
            "plasma_current_MA": 2.0, "central_density": 5.0e19,
            "n_e_min": 1.0e19, "n_e_max": 1.0e20, "density_peaking": 0.1, "temp_peaking": 1.0,
            "temp_peaking_i": 1.0, "profile_averaging": False,
            "deuterium_fraction": 0.5, "tritium_fraction": 0.5, "tauE_e": 0.01, "tauE_i": 0.01,
            "confinement_mode": "Kaye NSTX H-mode",
        },
        "T-15MD": {
            "major_radius": 1.5, "minor_radius": 0.67, "elongation": 1.8,
            "triangularity": 0.3, "effective_charge": 1.5, "toroidal_field": 2.0,
            "plasma_current_MA": 2.0, "central_density": 5.0e19,
            "n_e_min": 1.0e19, "n_e_max": 1.0e20, "density_peaking": 0.1, "temp_peaking": 1.0,
            "temp_peaking_i": 1.0, "profile_averaging": True,
            "deuterium_fraction": 0.5, "tritium_fraction": 0.5, "tauE_e": 0.1, "tauE_i": 0.1,
            "confinement_mode": "IPB98(y,2) ELMy H-mode",
        },
        "TCV": {
            "major_radius": 0.88, "minor_radius": 0.25, "elongation": 1.8,
            "triangularity": 0.5, "effective_charge": 2.0, "toroidal_field": 1.43,
            "plasma_current_MA": 0.4, "central_density": 5.0e19,
            "n_e_min": 5.0e18, "n_e_max": 1.0e20, "density_peaking": 0.1, "temp_peaking": 1.0,
            "temp_peaking_i": 1.0, "profile_averaging": False,
            "deuterium_fraction": 1.0, "tritium_fraction": 0.0, "tauE_e": 0.005, "tauE_i": 0.005,
            "confinement_mode": "IPB98(y,2) ELMy H-mode",
        },
    }

    # Typical/representative heating mix for each device (approximate,
    # rounded design or common-operating-point figures, not exact specs).
    # "beams" gives (species, power_MW, energy_keV, co_current) explicitly
    # for NBI-1 and NBI-2 -- explicit rather than inferred, because how a
    # device's real injected power is split across two beamlines (species,
    # energy, and direction alike) is itself a device fact (e.g. JET and
    # ITER each split one species across two real co-current injector
    # boxes; DANTE deliberately uses two different species; TCV's two real
    # beams are deliberately counter-injected relative to each other, sec.
    # E2). Species is always D or T -- this model has no mass entry for
    # hydrogen (_BEAM_MASS_NUMBER in hotjass/physics.py) -- and both beams'
    # tangent points are reset to (major_radius, 0) so they don't carry over
    # a differently-scaled device's geometry.
    # JET/TCV/DANTE are cross-checked against the sourced shots in
    # docs/Validation/HI-Jass_validation.md (sec. E); ITER/ST40/T-15MD are
    # not validated anywhere in this repo -- see that section's caveats.
    MACHINE_HEATING_PRESETS = {
        "DANTE": {"beams": [("D", 10.0, 120.0, True), ("T", 0.1, 180.0, False)],
                  "ecrh_MW": 2.0, "icrh_MW": 0.0},
        "ITER": {"beams": [("D", 16.5, 1000.0, True), ("D", 16.5, 1000.0, True)],
                 "ecrh_MW": 20.0, "icrh_MW": 20.0},
        "JET": {"beams": [("D", 13.25, 108.0, True), ("D", 13.25, 108.0, True)],
                "ecrh_MW": 0.0, "icrh_MW": 4.0},
        "ST40": {"beams": [("D", 2.7, 25.0, True), ("D", 0.0, 25.0, True)],
                 "ecrh_MW": 0.0, "icrh_MW": 0.0},
        "T-15MD": {"beams": [("D", 0.0, 100.0, True), ("D", 0.0, 100.0, True)],
                   "ecrh_MW": 10.0, "icrh_MW": 0.0},
        # NBI-1: original 2015 beam (Karpushov et al. 2017), co-current.
        # NBI-2: second, high-energy beam, counter-current -- TCV's two real
        # beams are deliberately tangentially opposed for low/no net torque.
        # Power/energy/tangent radius per Karpushov et al., Fusion Eng. Des.
        # 187 (2023) 113384 (both beams share the same 0.736 m tangent
        # radius, off-axis relative to R0=0.88 m).
        "TCV": {"beams": [("D", 1.3, 28.0, True), ("D", 1.0, 55.0, False)],
                "ecrh_MW": 4.5, "icrh_MW": 0.0, "tangent_R_m": 0.736},
    }

    PLASMA_FIELDS = [
        ("R0 [m]", "major_radius", 0.65),
        ("a [m]", "minor_radius", 0.35),
        ("Elongation kappa", "elongation", 2.2),
        ("Triangularity delta", "triangularity", -0.35),
        ("Z_eff", "effective_charge", 2.0),
        ("B0 [T]", "toroidal_field", 1.5),
        ("Ip [MA]", "plasma_current_MA", 1.5),
        ("Density peaking", "density_peaking", 0.1),
        ("Temp peaking (electron)", "temp_peaking", 1.0),
        ("Temp peaking (ion, -1=same)", "temp_peaking_i", -1.0),
        ("Centre-post R [m] (-1=R0-a)", "centrepost_radius", -1.0),
        ("Central n_e [m^-3]", "central_density", 1.0e20),
        ("Scan n_e min [m^-3]", "n_e_min", 1.0e19),
        ("Scan n_e max [m^-3]", "n_e_max", 1.5e20),
        ("D fraction", "deuterium_fraction", 0.2),
        ("T fraction", "tritium_fraction", 0.8),
    ]

    OBSERVABLES = {
        "Te, Ti": ["Te", "Ti"],
        "P_e, P_i, Pi_e, Q": ["P_e", "P_i", "Pi_e", "Q"],
        "n_D, n_T, n_b": ["n_D", "n_T", "n_b"],
        "Pf_tot, Pf_th, Pf_b, Pf_bb": ["Pf_tot", "Pf_th", "Pf_b", "Pf_bb"],
        "P_shine, P_orbit, P_cx, P_lost": ["P_shine-through", "P_orbit", "P_cx", "P_lost"],
        "Y_beam, Y_thermal, Y_bb, Y_full": ["Y_neutron_b", "Y_neutron_th", "Y_neutron_bb", "Y_neutron"],
        "<E_fast>": ["E_fast"],
        "tau_S, tauE_e, tauE_i, tau_IE": ["tau_S", "tauE_e", "tauE_i", "tau_IE"],
        "R = U_fast / U_th": ["R"],
        "Pr_th, Pr_fast (isotropic)": ["Pr_th", "Pr_fast"],
        "beta_T": ["beta_T"],
    }
    EQUIP_SENSITIVE = {"Tₑ, Tᵢ", "Pₑ, Pᵢ, Pᵢₑ, Q", "τS, τE,e, τE,i, τIE"}
    ALPHA_SENSITIVE = EQUIP_SENSITIVE | {
        "Pƒ,tot, Pƒ,th, Pƒ,b, Pƒ,bb", "Y_beam, Y_thermal, Y_bb, Y_full", "pₜₕ, pfast", "βt", "R = ufast / Uₜₕ"}
    DISPLAY_GROUPS = {
        "Tₑ, Tᵢ": "Te, Ti",
        "Pₑ, Pᵢ, Pᵢₑ, Q": "P_e, P_i, Pi_e, Q",
        "nᴅ, nₜ, nᵦ": "n_D, n_T, n_b",
        "Pƒ,tot, Pƒ,th, Pƒ,b, Pƒ,bb": "Pf_tot, Pf_th, Pf_b, Pf_bb",
        "P_shine, P_orbit, P_cx, P_lost": "P_shine, P_orbit, P_cx, P_lost",
        "Y_beam, Y_thermal, Y_bb, Y_full": "Y_beam, Y_thermal, Y_bb, Y_full",
        "⟨Efast⟩": "<E_fast>",
        "τS, τE,e, τE,i, τIE": "tau_S, tauE_e, tauE_i, tau_IE",
        "R = ufast / Uₜₕ": "R = U_fast / U_th",
        "pₜₕ, pfast": "Pr_th, Pr_fast (isotropic)",
        "βt": "beta_T",
    }
    UNITS = {
        "Te": "keV", "Ti": "keV", "P_e": "MW", "P_i": "MW", "Pi_e": "MW",
        "P_shine-through": "MW", "n_D": "m^-3", "n_T": "m^-3", "n_b": "m^-3",
        "Pf_tot": "MW", "Pf_th": "MW", "Pf_b": "MW", "Pf_bb": "MW", "P_useful": "MW", "Q": "1",
        "P_orbit": "MW", "P_cx": "MW", "P_lost": "MW",
        "Y_neutron_b": "1/s", "Y_neutron_th": "1/s", "Y_neutron_bb": "1/s", "Y_neutron": "1/s",
        "E_fast": "keV",
        "tau_S": "s", "tauE_e": "s", "tauE_i": "s", "tau_IE": "s", "R": "1",
        "Pr_th": "Pa", "Pr_fast": "Pa", "beta_T": "%",
    }
    LATEX_UNITS = {
        "keV": r"\mathrm{keV}", "MW": r"\mathrm{MW}", "m^-3": r"\mathrm{m}^{-3}",
        "s": r"\mathrm{s}", "1/s": r"\mathrm{s}^{-1}", "Pa": r"\mathrm{Pa}", "%": r"\%", "1": "1",
    }
    LATEX_NAMES = {
        "Te": r"$T_e$", "Ti": r"$T_i$", "P_e": r"$P_e$", "P_i": r"$P_i$",
        "Pi_e": r"$P_{ie}$", "P_shine-through": r"$P_{shine}$",
        "n_D": r"$n_D$", "n_T": r"$n_T$", "n_b": r"$n_{b0}$",
        "Pf_tot": r"$P_{f,tot}$", "Pf_th": r"$P_{f,th}$", "Pf_b": r"$P_{f,b}$", "Pf_bb": r"$P_{f,bb}$",
        "P_useful": r"$P_{useful}$", "Q": r"$Q$",
        "P_orbit": r"$P_{orbit}$", "P_cx": r"$P_{cx}$", "P_lost": r"$P_{lost}$",
        "Y_neutron_b": r"$Y_{n,beam}$", "Y_neutron_th": r"$Y_{n,th}$",
        "Y_neutron_bb": r"$Y_{n,bb}$", "Y_neutron": r"$Y_{n,full}$",
        "E_fast": r"$\langle E_{fast}\rangle$", "tau_S": r"$\tau_S$",
        "tauE_e": r"$\tau_{E,e}$", "tauE_i": r"$\tau_{E,i}$", "tau_IE": r"$\tau_{IE}$",
        "R": r"$R = u_{fast}/U_t$", "Pr_th": r"$p_{th}$", "Pr_fast": r"$p_{fast}$",
        "beta_T": r"$\beta_t$",
    }

    DASH_ROWS = [
        ("Feasibility", "feasibility"),
        ("tau_E,e / tau_E,i [s]", "tauE"),
        ("T_e [keV]", "Te"), ("T_i [keV]", "Ti"),
        ("  on-axis T_e0 / T_i0 [keV]", "T0"),
        ("n_e0 [m^-3]", "ne0"), ("n_b0 [m^-3]", "nb0"),
        ("P_NB injected [MW]", "P_NB"), ("P shine-through [MW]", "P_shine"),
        ("P captured [MW]", "P_capt"), ("P first-orbit loss [MW]", "P_orbit"),
        ("  NBI-1  orbit widths / f_orbit", "orb1"),
        ("  NBI-2  orbit widths / f_orbit", "orb2"),
        ("P charge-exchange loss [MW]", "P_cx"),
        ("  NBI-1  CX: f_cx,P / f_cx,N / gamma_cx", "cx1"),
        ("  NBI-2  CX: f_cx,P / f_cx,N / gamma_cx", "cx2"),
        ("P useful (to plasma) [MW]", "P_useful"),
        ("  -> electrons P_e [MW]", "P_e"), ("  -> ions P_i [MW]", "P_i"),
        ("P_ei equipartition (e->i) [MW]", "P_ei"),
        ("P_alpha self-heating [MW]", "P_alpha"),
        ("P_ECRH -> e / i [MW]", "p_ecrh"),
        ("P_ICRH -> e / i [MW]", "p_icrh"),
        ("P_heat total (NBI + aux + alpha) [MW]", "P_heat"),
        ("P_fusion total [MW]", "Pf_tot"),
        ("  D-T (thermal / beam-target / beam-beam) [MW]", "Pf_dt"),
        ("  D-D (thermal / beam-target / beam-beam) [MW]", "Pf_dd"),
        ("neutron rate [n/s]", "Y_n"), ("Q = P_fus / P_NB", "Q"),
        ("<E_fast> [keV]", "E_fast"), ("beta_t [%]", "beta_t"),
        ("q* (edge safety factor)", "q_star"),
        ("Toroidal rotation v_phi [km/s] (Mach)", "v_phi"),
        ("<n_e>/n_GW  (Greenwald, at ne_c)", "f_gw"),
        ("Dominant loss", "dominant"),
    ]

    def __init__(self):
        super().__init__()
        self.title("HI-Jass")
        self.geometry("1360x880")
        self.minsize(1100, 720)

        self.model = HotJassModel()
        self.mode = "Operating point"
        self.last_result: Result | None = None
        self._result_queue: queue.Queue[Result] = queue.Queue()
        self._running = False
        self._last_infeasible_key: str | None = None
        self._applied_snapshot: dict[str, str] = {}
        self._saved = self._load_settings()
        self.active_device = self._saved.get("_active_device", "DANTE")

        saved_scale = self._saved.get("_ui_scale", 1.0)
        self._ui_scale_idx = min(
            range(len(UI_SCALES)), key=lambda i: abs(UI_SCALES[i] - saved_scale))
        self._apply_ui_scale(rerender=False)

        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(0, weight=1)

        self._build_rail()
        self._build_results()

        saved_mode = self._saved.get("_mode")
        if saved_mode in ("Operating point", "Scan"):
            self.mode_toggle.set(saved_mode)
        self._set_mode(self.mode_toggle.get())
        self._snapshot_inputs()
        self._highlight_active_preset()

        self.protocol("WM_DELETE_WINDOW", self._on_close)
        self.after(150, self._run)

    # ------------------------------------------------------------------ rail
    def _build_rail(self):
        rail = ctk.CTkScrollableFrame(self, width=350, label_text="Inputs")
        rail.grid(row=0, column=0, sticky="nsw", padx=(10, 6), pady=10)
        rail.grid_columnconfigure(0, weight=1)
        self.rail = rail
        self.entries: dict[str, ctk.CTkEntry] = {}

        run_bar = ctk.CTkFrame(rail, fg_color="transparent")
        run_bar.grid(row=0, column=0, sticky="ew", pady=(0, 8))
        run_bar.grid_columnconfigure(0, weight=1)
        self.mode_toggle = ctk.CTkSegmentedButton(
            run_bar, values=["Operating point", "Scan"], command=self._set_mode,
        )
        self.mode_toggle.set(self.mode)
        self.mode_toggle.grid(row=0, column=0, sticky="ew")
        self.zoom_btn = ctk.CTkButton(
            run_bar, text=self._zoom_btn_text(), width=64,
            command=self._cycle_ui_scale)
        self.zoom_btn.grid(row=0, column=1, sticky="ns", padx=(6, 0))
        self.run_btn = ctk.CTkButton(run_bar, text="▶  Run", command=self._run)
        self.run_btn.grid(row=1, column=0, columnspan=2, sticky="ew", pady=(6, 0))
        self._run_btn_fg = self.run_btn.cget("fg_color")
        self._run_btn_hover = self.run_btn.cget("hover_color")
        self.status = ctk.CTkLabel(run_bar, text="Ready", anchor="w", text_color="gray")
        self.status.grid(row=2, column=0, columnspan=2, sticky="ew", pady=(4, 0))

        row = 1
        presets = CollapsibleSection(rail, "Presets", expanded=True)
        presets.grid(row=row, column=0, sticky="ew", pady=4)
        self.preset_buttons: dict[str, ctk.CTkButton] = {}
        for i, name in enumerate(self.PLASMA_PRESETS):
            btn = ctk.CTkButton(
                presets.body, text=name, width=70,
                command=lambda n=name: self._apply_preset(n),
            )
            btn.grid(row=i // 3, column=i % 3, padx=3, pady=3, sticky="ew")
            self.preset_buttons[name] = btn
        self._preset_fg_default = next(iter(self.preset_buttons.values())).cget("fg_color")
        self._preset_hover_default = next(iter(self.preset_buttons.values())).cget("hover_color")
        row += 1

        plasma = CollapsibleSection(rail, "Plasma", expanded=True)
        plasma.grid(row=row, column=0, sticky="ew", pady=4)
        self._add_entries(plasma.body, "plasma", self._plasma_defaults())
        n_plasma_fields = len(self.PLASMA_FIELDS)
        self.profile_var = ctk.BooleanVar(value=self._saved.get("_profile_averaging", False))
        ctk.CTkCheckBox(
            plasma.body, text="profile-corrected 0-D (central n_e in; <T> + T0 out)",
            variable=self.profile_var,
        ).grid(row=n_plasma_fields, column=0, columnspan=2, padx=8, pady=4, sticky="w")
        self.equip_var = ctk.BooleanVar(value=self._saved.get("_equipartition", True))
        ctk.CTkCheckBox(
            plasma.body, text="e-i equipartition (couple Te, Ti)", variable=self.equip_var,
        ).grid(row=n_plasma_fields + 1, column=0, columnspan=2, padx=8, pady=4, sticky="w")
        self.alpha_var = ctk.BooleanVar(value=self._saved.get("_alpha_heating", False))
        ctk.CTkCheckBox(
            plasma.body, text="alpha self-heating (P_a into Te, Ti)", variable=self.alpha_var,
        ).grid(row=n_plasma_fields + 2, column=0, columnspan=2, padx=8, pady=4, sticky="w")
        self._add_entries(plasma.body, "models", [
            ("f_alpha (confined fraction 0-1)", "f_alpha", self._saved.get("models.f_alpha", 1.0)),
        ], start_row=n_plasma_fields + 3)
        row += 1

        for beam_index in (0, 1):
            beam = self.model.beams[beam_index]
            section = CollapsibleSection(rail, f"NBI-{beam_index + 1}", expanded=beam_index == 0)
            section.grid(row=row, column=0, sticky="ew", pady=4)
            self._add_entries(section.body, f"beam{beam_index}", [
                ("Species", "species", beam.species),
                ("P_NB [MW]", "power_MW", beam.power_MW),
                ("E_b [keV]", "beam_energy_keV", beam.beam_energy_keV),
                ("Tangent R_t [m]", "tangent_R_m",
                 beam.tangent_R_m if beam.tangent_R_m is not None else self.model.plasma.major_radius),
                ("Tangent Z_t [m]", "tangent_Z_m", beam.tangent_Z_m),
                ("Manual shine-through frac", "manual_shine_through_fraction",
                 beam.manual_shine_through_fraction),
            ])
            var = ctk.StringVar(value=self._saved.get(f"beam{beam_index}._shine",
                                                      SHINE_MODEL_TO_LABEL.get(beam.shine_through_model, "Riviere")))
            setattr(self, f"shine_var_{beam_index}", var)
            r = len(section.body.grid_slaves()) // 2
            ctk.CTkLabel(section.body, text="Shine-through model", anchor="w").grid(
                row=r, column=0, padx=8, pady=(8, 2), sticky="w")
            ctk.CTkSegmentedButton(
                section.body, values=list(SHINE_LABEL_TO_MODEL), variable=var,
            ).grid(row=r, column=1, padx=8, pady=(8, 2), sticky="ew")
            dir_var = ctk.StringVar(value=self._saved.get(
                f"beam{beam_index}._dir", "co" if beam.co_current else "counter"))
            setattr(self, f"beam_dir_var_{beam_index}", dir_var)
            ctk.CTkLabel(section.body, text="Injection vs Ip", anchor="w").grid(
                row=r + 1, column=0, padx=8, pady=(6, 2), sticky="w")
            ctk.CTkSegmentedButton(
                section.body, values=["co", "counter"], variable=dir_var,
            ).grid(row=r + 1, column=1, padx=8, pady=(6, 2), sticky="ew")
            row += 1

        # "Losses" groups every "how is this computed" selector (confinement
        # time, first-orbit loss, charge-exchange loss) in one place, right
        # below the two NBI sections whose power these losses act on --
        # separated from the plain physics toggles now living in "Plasma"
        # (profile-corrected 0-D, e-i equipartition, alpha self-heating),
        # which aren't loss mechanisms. Shine-through's own selector stays
        # inline in each NBI-n section (already compact there).
        losses = CollapsibleSection(rail, "Losses", expanded=False)
        losses.grid(row=row, column=0, sticky="ew", pady=4)
        ctk.CTkLabel(losses.body, text="Confinement", anchor="w").grid(
            row=0, column=0, padx=8, pady=4, sticky="w")
        saved_confinement = self._saved.get("_confinement", "Fixed tauE (input)")
        if saved_confinement not in CONFINEMENT_MODES:
            # A settings.json (or an imported run record) from before a mode
            # was retired -- fall back rather than KeyError on first _run().
            saved_confinement = "Fixed tauE (input)"
        self.confinement_var = ctk.StringVar(value=saved_confinement)
        ctk.CTkOptionMenu(
            losses.body, values=list(CONFINEMENT_MODES), variable=self.confinement_var,
        ).grid(row=0, column=1, padx=8, pady=4, sticky="ew")
        self._add_entries(losses.body, "plasma", [
            ("tauE,e [s]", "tauE_e", self._saved.get("plasma.tauE_e", 0.15)),
            ("tauE,i [s]", "tauE_i", self._saved.get("plasma.tauE_i", 0.15)),
        ], start_row=1)
        self.orbit_var = ctk.BooleanVar(value=self._saved.get("_orbit_loss", True))
        ctk.CTkCheckBox(
            losses.body, text="First-orbit loss (direction set per NBI)", variable=self.orbit_var,
        ).grid(row=3, column=0, columnspan=2, padx=8, pady=4, sticky="w")
        ctk.CTkLabel(losses.body, text="Orbit model", anchor="w").grid(
            row=4, column=0, padx=8, pady=4, sticky="w")
        self.orbit_model_var = ctk.StringVar(
            value=self._saved.get("_orbit_model", "ST orbits - pitch-resolved"))
        ctk.CTkOptionMenu(
            losses.body, values=list(ORBIT_MODELS), variable=self.orbit_model_var,
        ).grid(row=4, column=1, padx=8, pady=4, sticky="ew")
        ctk.CTkLabel(losses.body, text="CX-loss model", anchor="w").grid(
            row=5, column=0, padx=8, pady=4, sticky="w")
        self.cx_model_var = ctk.StringVar(
            value=self._saved.get("_cx_model", "Manual fraction"))
        ctk.CTkOptionMenu(
            losses.body, values=list(CX_MODELS), variable=self.cx_model_var,
        ).grid(row=5, column=1, padx=8, pady=4, sticky="ew")
        self._add_entries(losses.body, "models", [
            ("CX loss fraction (0-1)", "cx_loss_fraction",
             self._saved.get("models.cx_loss_fraction", 0.1)),
            ("CX: n0/ne (manual n0)", "cx_n0_over_ne",
             self._saved.get("models.cx_n0_over_ne", 1.0e-5)),
            ("CX: n0_LCFS/ne (penetration)", "cx_n0_lcfs_over_ne",
             self._saved.get("models.cx_n0_lcfs_over_ne", 0.02)),
        ], start_row=6)
        ctk.CTkLabel(losses.body, text="Rotation model", anchor="w").grid(
            row=9, column=0, padx=8, pady=4, sticky="w")
        self.rotation_model_var = ctk.StringVar(
            value=self._saved.get("_rotation_model", "Off"))
        ctk.CTkOptionMenu(
            losses.body, values=list(ROTATION_MODELS), variable=self.rotation_model_var,
        ).grid(row=9, column=1, padx=8, pady=4, sticky="ew")
        self._add_entries(losses.body, "models", [
            ("v_phi manual [m/s] (+=co)", "manual_v_phi_m_s",
             self._saved.get("models.manual_v_phi_m_s", 0.0)),
            ("tau_phi / tauE,i (mom. balance)", "tau_phi_over_tauEi",
             self._saved.get("models.tau_phi_over_tauEi", 1.0)),
        ], start_row=10)
        self.beam_beam_var = ctk.BooleanVar(value=self._saved.get("_beam_beam", False))
        ctk.CTkCheckBox(
            losses.body, text="Beam-beam fusion (reduced, NBI-1 x NBI-2)", variable=self.beam_beam_var,
        ).grid(row=12, column=0, columnspan=2, padx=8, pady=4, sticky="w")
        row += 1

        aux = CollapsibleSection(rail, "Aux heating (ECRH / ICRH)", expanded=False)
        aux.grid(row=row, column=0, sticky="ew", pady=4)
        self._add_entries(aux.body, "aux", [
            ("P_ECRH [MW]", "p_ecrh_MW", self.model.plasma.p_ecrh_MW),
            ("ECRH f_e (rest to ions)", "ecrh_f_e", self.model.plasma.ecrh_f_e),
            ("P_ICRH [MW]", "p_icrh_MW", self.model.plasma.p_icrh_MW),
            ("ICRH f_e", "icrh_f_e", self.model.plasma.icrh_f_e),
            ("ICRH f_i", "icrh_f_i", self.model.plasma.icrh_f_i),
        ])
        row += 1

        # Second Run button at the foot of the input rail, so the user does not
        # have to scroll back up after changing a model. Kept in lock-step with
        # the top button by _configure_run_buttons().
        self.run_btn_bottom = ctk.CTkButton(rail, text="▶  Run", command=self._run)
        self.run_btn_bottom.grid(row=row, column=0, sticky="ew", pady=(10, 4))
        self._run_btns = [self.run_btn, self.run_btn_bottom]

    def _configure_run_buttons(self, **kwargs):
        for btn in getattr(self, "_run_btns", [self.run_btn]):
            btn.configure(**kwargs)

    def _plasma_defaults(self):
        preset = {**{f: d for _, f, d in self.PLASMA_FIELDS}}
        return [(label, field, preset[field]) for label, field, _ in self.PLASMA_FIELDS]

    def _add_entries(self, parent, prefix, fields, start_row=0):
        for offset, (label, field, default) in enumerate(fields):
            row = start_row + offset
            key = f"{prefix}.{field}"
            ctk.CTkLabel(parent, text=label, anchor="w").grid(
                row=row, column=0, padx=8, pady=3, sticky="w")
            var = ctk.StringVar(value=str(self._saved.get(key, default)))
            entry = ctk.CTkEntry(parent, width=120, textvariable=var)
            entry.grid(row=row, column=1, padx=8, pady=3, sticky="ew")
            self.entries[key] = entry
            if not hasattr(self, "_entry_border_default"):
                self._entry_border_default = entry.cget("border_color")
            var.trace_add("write", lambda *_a, k=key: self._mark_dirty(k))

    # --------------------------------------------------------------- results
    def _build_results(self):
        container = ctk.CTkFrame(self, fg_color="transparent")
        container.grid(row=0, column=1, sticky="nsew", padx=(6, 10), pady=10)
        container.grid_columnconfigure(0, weight=1)
        container.grid_rowconfigure(0, weight=1)

        self.tv_op = ctk.CTkTabview(container)
        self.tv_scan = ctk.CTkTabview(container)
        for tv in (self.tv_op, self.tv_scan):
            tv.grid(row=0, column=0, sticky="nsew")

        for name in ("Dashboard", "Deposition", "Power flow", "Profiles", "Fusion", "Assumptions", "References"):
            self.tv_op.add(name)
        for name in ("Scan", "Summary"):
            self.tv_scan.add(name)

        # Plot-heavy tabs (each a matplotlib Figure + FigureCanvasTkAgg +
        # NavigationToolbar2Tk) are the dominant cost of building this
        # window in CustomTkinter. Only the tab shown by default (Dashboard)
        # plus the always-referenced Assumptions/References tabs are built
        # eagerly; the rest are built on first visit via _ensure_tab_built,
        # wired to each CTkTabview's command= callback below -- so startup
        # (and the very first solve's render pass) no longer pays for tabs
        # the user hasn't opened yet.
        self._built_tabs: set[str] = set()
        self._zoom_dirty_tabs: set[str] = set()
        self._tab_builders = {
            "Deposition": self._build_deposition_tab,
            "Power flow": self._build_powerflow_tab,
            "Profiles": self._build_profiles_tab,
            "Fusion": self._build_fusion_tab,
            "Scan": self._build_scan_tab,
            "Summary": self._build_summary_tab,
        }
        self._tab_renderers = {
            "Deposition": lambda r: self._render_deposition(r),
            "Power flow": lambda r: self._render_powerflow(r),
            "Profiles": lambda r: self._render_profiles(r),
            "Fusion": lambda r: self._render_fusion(r),
            "Scan": lambda r: self._render_scan_tab(r),
            "Summary": lambda r: self._render_summary(r),
        }
        self.tv_op.configure(command=lambda: self._on_result_tab_shown(self.tv_op.get()))
        self.tv_scan.configure(command=lambda: self._on_result_tab_shown(self.tv_scan.get()))

        self._build_dashboard(self.tv_op.tab("Dashboard"))
        self._build_assumptions(self.tv_op.tab("Assumptions"))
        self._build_references(self.tv_op.tab("References"))

        # Two visually distinct rows -- "Export" (send results out) and
        # "Import" (bring inputs back in) are opposite directions of data
        # flow, so a beginner shouldn't have to spot "Load" hiding among a
        # row of otherwise-all-export buttons. The Import button also gets
        # its own accent colour (the same blue this app already uses for
        # "incoming"/injected quantities, e.g. the power-flow waterfall's
        # P_inj bar) so it reads as a different kind of action at a glance.
        export = ctk.CTkFrame(container, fg_color="transparent")
        export.grid(row=1, column=0, sticky="ew", pady=(6, 0))
        ctk.CTkLabel(export, text="Export:").pack(side="left", padx=(4, 8))
        ctk.CTkButton(export, text="Summary PNG", width=110,
                      command=lambda: self._export_summary("png")).pack(side="left", padx=4)
        ctk.CTkButton(export, text="Summary PDF", width=110,
                      command=lambda: self._export_summary("pdf")).pack(side="left", padx=4)
        ctk.CTkButton(export, text="Run record (JSON)", width=150,
                      command=self._export_json).pack(side="left", padx=4)

        run_import = ctk.CTkFrame(container, fg_color="transparent")
        run_import.grid(row=2, column=0, sticky="ew", pady=(4, 0))
        ctk.CTkLabel(run_import, text="Import:").pack(side="left", padx=(4, 8))
        ctk.CTkButton(run_import, text="Load run record (JSON)...", width=190,
                      fg_color="#3b6fb0", hover_color="#2f5a8f",
                      command=self._import_json).pack(side="left", padx=4)

        self._render_references()

    def _ensure_tab_built(self, name: str):
        """Lazily construct a plot-heavy result tab the first time it is shown."""
        builder = self._tab_builders.get(name)
        if builder is None or name in self._built_tabs:
            return
        self._built_tabs.add(name)
        builder()
        if self.last_result is not None:
            self._tab_renderers[name](self.last_result)

    def _on_result_tab_shown(self, name: str):
        """CTkTabview command= callback: build a tab on first visit, or --
        if a UI-scale change happened while it was hidden -- redraw it now
        instead of at zoom-click time, so a zoom click only ever pays for
        the one tab currently on screen."""
        if name not in self._tab_builders:
            return
        if name not in self._built_tabs:
            self._ensure_tab_built(name)
        elif name in self._zoom_dirty_tabs and self.last_result is not None:
            self._zoom_dirty_tabs.discard(name)
            self._tab_renderers[name](self.last_result)

    def _current_visible_tab(self) -> str:
        return self.tv_scan.get() if self.mode == "Scan" else self.tv_op.get()

    def _build_deposition_tab(self):
        self.dep_fig, self.dep_canvas, holder = self._plot_area(
            self.tv_op.tab("Deposition"), figsize=(9.5, 7.2))
        holder.pack(fill="both", expand=True)

    def _build_powerflow_tab(self):
        pf_tab = self.tv_op.tab("Power flow")
        pf_bar = ctk.CTkFrame(pf_tab, fg_color="transparent")
        pf_bar.pack(fill="x", pady=(2, 0))
        ctk.CTkLabel(pf_bar, text="View:").pack(side="left", padx=(4, 6))
        self.pf_view_var = ctk.StringVar(value="all")
        ctk.CTkSegmentedButton(
            pf_bar, values=["all", "waterfall", "sankey", "pie"], variable=self.pf_view_var,
            command=lambda _: self._render_powerflow(self.last_result) if self.last_result else None,
        ).pack(side="left")
        self.pf_fig, self.pf_canvas, holder = self._plot_area(pf_tab, figsize=(9.5, 6.2))
        holder.pack(fill="both", expand=True)

    def _build_profiles_tab(self):
        self.prof_fig, self.prof_canvas, holder = self._plot_area(
            self.tv_op.tab("Profiles"), figsize=(13.5, 6.2))
        holder.pack(fill="both", expand=True)

    def _build_fusion_tab(self):
        self.fusion_fig, self.fusion_canvas, holder = self._plot_area(
            self.tv_op.tab("Fusion"), figsize=(7.5, 7.0))
        holder.pack(fill="both", expand=True)

    def _build_scan_tab(self):
        self._build_scan(self.tv_scan.tab("Scan"))

    def _build_summary_tab(self):
        self.sum_fig, self.sum_canvas, holder = self._plot_area(self.tv_scan.tab("Summary"), figsize=(11, 8))
        holder.pack(fill="both", expand=True)

    def _plot_area(self, parent, figsize=(7.5, 5.5)):
        """Build a figure + canvas + toolbar inside a holder frame.

        The holder is returned unplaced; the caller pack()s or grid()s it.
        """
        holder = ctk.CTkFrame(parent, fg_color="transparent")
        fig = Figure(figsize=figsize, dpi=100)
        canvas = FigureCanvasTkAgg(fig, master=holder)
        toolbar = NavigationToolbar2Tk(canvas, holder, pack_toolbar=False)
        toolbar.update()
        toolbar.pack(side="bottom", fill="x")
        canvas.get_tk_widget().pack(side="top", fill="both", expand=True)
        return fig, canvas, holder

    def _build_dashboard(self, parent):
        frame = ctk.CTkScrollableFrame(parent)
        frame.pack(fill="both", expand=True)
        frame.grid_columnconfigure(1, weight=1)
        self.dash_values: dict[str, ctk.CTkLabel] = {}
        for row, (label, key) in enumerate(self.DASH_ROWS):
            ctk.CTkLabel(frame, text=label, anchor="w").grid(
                row=row, column=0, padx=10, pady=3, sticky="w")
            value = ctk.CTkLabel(frame, text="—", anchor="w",
                                 font=ctk.CTkFont(size=13, weight="bold"))
            value.grid(row=row, column=1, padx=10, pady=3, sticky="w")
            self.dash_values[key] = value

    def _build_assumptions(self, parent):
        self.assump_box = ctk.CTkTextbox(parent, wrap="word", font=ctk.CTkFont(family="monospace", size=12))
        self.assump_box.pack(fill="both", expand=True, padx=6, pady=6)

    def _build_scan(self, parent):
        parent.grid_columnconfigure(1, weight=1)
        parent.grid_rowconfigure(0, weight=1)
        left = ctk.CTkFrame(parent, fg_color="transparent")
        left.grid(row=0, column=0, sticky="ns", padx=6, pady=6)
        ctk.CTkLabel(left, text="n_e scan", font=ctk.CTkFont(size=15, weight="bold")).pack(pady=6)
        self.observable_var = ctk.StringVar(value=list(self.DISPLAY_GROUPS)[0])
        ctk.CTkOptionMenu(left, variable=self.observable_var, values=list(self.DISPLAY_GROUPS),
                          command=lambda _: self._render_scan_plot()).pack(pady=6)
        self.scan_equip_label = ctk.CTkLabel(
            left, text="", anchor="w", justify="left", wraplength=270,
            font=ctk.CTkFont(size=12, weight="bold"))
        self.scan_equip_label.pack(pady=(0, 0), fill="x", padx=6)
        self.scan_alpha_label = ctk.CTkLabel(
            left, text="", anchor="w", justify="left", wraplength=270,
            font=ctk.CTkFont(size=12, weight="bold"))
        self.scan_alpha_label.pack(pady=(0, 4), fill="x", padx=6)
        self.scan_text = ctk.CTkTextbox(left, width=280, height=410, wrap="word")
        self.scan_text.pack(pady=8, fill="y")
        self.scan_fig, self.scan_canvas, holder = self._plot_area(parent, figsize=(8, 6))
        holder.grid(row=0, column=1, sticky="nsew", padx=6, pady=6)

    # ----------------------------------------------------------------- modes
    def _set_mode(self, mode: str):
        # CTkBaseClass remembers each widget's last grid(**kwargs) call and
        # replays it whenever ctk.set_widget_scaling() runs (to rescale
        # padding/size in those kwargs) -- but grid_remove() doesn't clear
        # that memory the way grid_forget() does. So a tabview hidden with
        # grid_remove() silently reappears the next time the zoom button is
        # pressed. grid_forget() (paired with an explicit re-grid() here,
        # since forget -- unlike remove -- doesn't remember placement) keeps
        # a hidden tabview hidden across a scale change.
        self.mode = mode
        if mode == "Scan":
            self.tv_op.grid_forget()
            self.tv_scan.grid(row=0, column=0, sticky="nsew")
            self._on_result_tab_shown(self.tv_scan.get())
        else:
            self.tv_scan.grid_forget()
            self.tv_op.grid(row=0, column=0, sticky="nsew")
            self._on_result_tab_shown(self.tv_op.get())

    # ------------------------------------------------------------ ui scale
    def _zoom_btn_text(self) -> str:
        return f"Aa {int(round(UI_SCALES[self._ui_scale_idx] * 100))}%"

    def _apply_ui_scale(self, rerender: bool = True):
        scale = UI_SCALES[self._ui_scale_idx]
        ctk.set_widget_scaling(scale)
        matplotlib.rcParams["font.size"] = _BASE_MPL_FONT_SIZE * scale
        if hasattr(self, "zoom_btn"):
            self.zoom_btn.configure(text=self._zoom_btn_text())
        if not rerender or self.last_result is None:
            return
        # Redrawing every built plot tab on each zoom click was the cause of
        # the window appearing to "freeze" after zooming (each matplotlib
        # figure's text has to be re-laid-out at the new font size). Redraw
        # only the tab actually on screen right now; the rest are marked
        # dirty and get redrawn lazily, one at a time, the next time each is
        # actually shown (see _on_result_tab_shown).
        self._render_dashboard(self.last_result)
        visible = self._current_visible_tab()
        self._zoom_dirty_tabs |= self._built_tabs
        if visible in self._built_tabs:
            self._zoom_dirty_tabs.discard(visible)
            self._tab_renderers[visible](self.last_result)

    def _cycle_ui_scale(self):
        self._ui_scale_idx = (self._ui_scale_idx + 1) % len(UI_SCALES)
        self._apply_ui_scale()

    # ------------------------------------------------------------- run/solve
    def _apply_inputs(self) -> list[str]:
        errors: list[str] = []
        plasma = self.model.plasma
        for label, field, _ in self.PLASMA_FIELDS:
            raw = self.entries[f"plasma.{field}"].get()
            try:
                value = float(raw)
            except ValueError:
                errors.append(label)
                continue
            if field == "plasma_current_MA":
                plasma.plasma_current = value * 1.0e6
            else:
                setattr(plasma, field, value)
        if plasma.n_e_min <= 0 or plasma.n_e_max <= plasma.n_e_min:
            errors.append("scan n_e range")
        if plasma.central_density <= 0:
            errors.append("central n_e")

        # tauE,e / tauE,i live in the Losses section (Confinement subsection),
        # not the PLASMA_FIELDS loop above, since they're only meaningful
        # under "Fixed tauE (input)" -- parsed the same way, just separately.
        for label, field in (("tauE,e", "tauE_e"), ("tauE,i", "tauE_i")):
            try:
                setattr(plasma, field, float(self.entries[f"plasma.{field}"].get()))
            except ValueError:
                errors.append(label)

        ee_mode, ei_mode = CONFINEMENT_MODES[self.confinement_var.get()]
        plasma.tau_Ee_mode, plasma.tau_Ei_mode = ee_mode, ei_mode
        plasma.enable_orbit_loss = bool(self.orbit_var.get())
        plasma.orbit_model = ORBIT_MODELS[self.orbit_model_var.get()]
        plasma.profile_averaging = bool(self.profile_var.get())
        plasma.enable_equipartition = bool(self.equip_var.get())
        if plasma.enable_equipartition and plasma.tau_Ei_mode in ("neoclassical", "neoclassical_arbA"):
            errors.append("equipartition + neoclassical ion (unsupported together)")
        plasma.alpha_heating = bool(self.alpha_var.get())
        try:
            plasma.f_alpha = min(max(float(self.entries["models.f_alpha"].get()), 0.0), 1.0)
        except ValueError:
            errors.append("f_alpha")
        try:
            cx = float(self.entries["models.cx_loss_fraction"].get())
            plasma.cx_loss_fraction = min(max(cx, 0.0), 1.0)
        except ValueError:
            errors.append("CX loss fraction")
        plasma.cx_model = CX_MODELS[self.cx_model_var.get()]
        try:
            plasma.cx_n0_over_ne = max(float(self.entries["models.cx_n0_over_ne"].get()), 0.0)
        except ValueError:
            errors.append("CX n0/ne (manual)")
        try:
            plasma.cx_n0_lcfs_over_ne = max(float(self.entries["models.cx_n0_lcfs_over_ne"].get()), 0.0)
        except ValueError:
            errors.append("CX n0_LCFS/ne (penetration)")
        plasma.rotation_model = ROTATION_MODELS[self.rotation_model_var.get()]
        try:
            plasma.manual_v_phi_m_s = float(self.entries["models.manual_v_phi_m_s"].get())
        except ValueError:
            errors.append("v_phi manual")
        try:
            plasma.tau_phi_over_tauEi = max(float(self.entries["models.tau_phi_over_tauEi"].get()), 0.0)
        except ValueError:
            errors.append("tau_phi / tauE,i")
        plasma.enable_beam_beam = bool(self.beam_beam_var.get())

        for label, field, clamp01 in (
            ("P_ECRH", "p_ecrh_MW", False), ("ECRH f_e", "ecrh_f_e", True),
            ("P_ICRH", "p_icrh_MW", False), ("ICRH f_e", "icrh_f_e", True),
            ("ICRH f_i", "icrh_f_i", True),
        ):
            try:
                val = float(self.entries[f"aux.{field}"].get())
                setattr(plasma, field, min(max(val, 0.0), 1.0) if clamp01 else max(val, 0.0))
            except ValueError:
                errors.append(label)

        total = plasma.deuterium_fraction + plasma.tritium_fraction
        if total > 0:
            plasma.deuterium_fraction /= total
            plasma.tritium_fraction /= total

        for beam_index, beam in enumerate(self.model.beams):
            prefix = f"beam{beam_index}"
            beam.species = self.entries[f"{prefix}.species"].get().strip().upper() or "D"
            for field in ("power_MW", "beam_energy_keV", "tangent_R_m", "tangent_Z_m",
                          "manual_shine_through_fraction"):
                try:
                    setattr(beam, field, float(self.entries[f"{prefix}.{field}"].get()))
                except ValueError:
                    errors.append(f"NBI-{beam_index + 1} {field}")
            beam.shine_through_model = SHINE_LABEL_TO_MODEL[
                getattr(self, f"shine_var_{beam_index}").get()]
            beam.co_current = getattr(self, f"beam_dir_var_{beam_index}").get() == "co"
        return errors

    # -------------------------------------------------------- edit indicators
    def _snapshot_inputs(self):
        self._applied_snapshot = {k: e.get() for k, e in self.entries.items()}
        self._refresh_dirty()

    def _mark_dirty(self, key: str):
        entry = self.entries[key]
        dirty = entry.get() != self._applied_snapshot.get(key)
        entry.configure(border_color="#c0504d" if dirty else self._entry_border_default)
        self._update_dirty_hint()

    def _refresh_dirty(self):
        for key in self.entries:
            self._mark_dirty(key)

    def _dirty_keys(self):
        return [k for k, e in self.entries.items() if e.get() != self._applied_snapshot.get(k)]

    def _update_dirty_hint(self):
        if self._running:
            return
        n = len(self._dirty_keys())
        if n:
            self.status.configure(text=f"{n} edited input(s) — press Run", text_color="#d98b00")
        elif self.status.cget("text").endswith("press Run"):
            self.status.configure(text="Ready", text_color="gray")

    def _highlight_active_preset(self):
        for name, btn in self.preset_buttons.items():
            active = name == self.active_device
            btn.configure(fg_color="#2f7d32" if active else self._preset_fg_default,
                          hover_color="#276b2a" if active else self._preset_hover_default)

    def _run(self):
        errors = self._apply_inputs()
        if errors:
            self.status.configure(text="Invalid: " + ", ".join(errors), text_color="#c0504d")
            self._refresh_dirty()
            return
        if self._running:
            return
        self._snapshot_inputs()
        self._running = True
        self.status.configure(text="Running…", text_color="gray")
        self._configure_run_buttons(state="disabled", text="⏳  Running…",
                                    fg_color="#c0504d", hover_color="#c0504d")
        threading.Thread(target=self._worker, daemon=True).start()
        self.after(80, self._poll_result)

    def _worker(self):
        try:
            model = self.model
            scan = model.density_scan()
            op = model.operating_point()
            summary = model.summary()
            result = Result(
                ok=True, op=op, scan=scan, summary=summary,
                volume_m3=model.plasma_volume(),
            )
        except Exception as exc:  # surface, don't crash the UI thread
            result = Result(ok=False, error=f"{type(exc).__name__}: {exc}\n\n{traceback.format_exc()}")
        self._result_queue.put(result)

    def _poll_result(self):
        try:
            result = self._result_queue.get_nowait()
        except queue.Empty:
            self.after(80, self._poll_result)
            return
        self._running = False
        self._render(result)

    def _render(self, result: Result):
        self._configure_run_buttons(state="normal", text="▶  Run",
                                    fg_color=self._run_btn_fg, hover_color=self._run_btn_hover)
        if not result.ok:
            self.status.configure(text="Error — see Assumptions tab", text_color="#c0504d")
            self.assump_box.delete("1.0", "end")
            self.assump_box.insert("end", result.error)
            return
        self.last_result = result
        self.scan = result.scan
        op = result.op
        self.status.configure(
            text="Done" if op.feasible else "Done — infeasible operating point (see details)",
            text_color="gray" if op.feasible else "#c0504d",
        )
        if not op.feasible:
            reason = op.infeasible_reason or (
                "The solver returned an infeasible operating point but gave no detail "
                "(e.g. an inner temperature root-find failed to converge). Try a larger "
                "tau_E, lower beam power, or a higher central n_e.")
            if reason != self._last_infeasible_key:
                self._last_infeasible_key = reason
                messagebox.showwarning(
                    "HI-Jass — infeasible operating point",
                    "No self-consistent operating point at the central density:\n\n"
                    + reason + "\n\n(Full text is also on the Assumptions tab.)",
                    parent=self)
        else:
            self._last_infeasible_key = None
        self._render_dashboard(result)
        for name in ("Deposition", "Power flow", "Profiles", "Fusion", "Scan", "Summary"):
            if name in self._built_tabs:
                self._tab_renderers[name](result)
        self.assump_box.delete("1.0", "end")
        self.assump_box.insert("end", self._assess(self.model, result.op))
        self._render_references()

    def _render_scan_tab(self, result: Result):
        eq_on = self.model.plasma.enable_equipartition
        self.scan_equip_label.configure(
            text=("●  e-i equipartition: ON  (Te, Ti coupled)" if eq_on
                  else "○  e-i equipartition: off  (Te, Ti independent)"),
            text_color=("#2f7d32" if eq_on else "gray"))
        al_on = self.model.plasma.alpha_heating
        self.scan_alpha_label.configure(
            text=(f"●  alpha self-heating: ON  (f_alpha={self.model.plasma.f_alpha:.2g})" if al_on
                  else "○  alpha self-heating: off"),
            text_color=("#2f7d32" if al_on else "gray"))
        self._render_scan_plot()

    # --------------------------------------------------------------- renders
    @staticmethod
    def _fmt(value, spec="{:.3g}"):
        if value is None or (isinstance(value, float) and not np.isfinite(value)):
            return "—"
        return spec.format(value)

    def _render_dashboard(self, result: Result):
        op = result.op
        mw = 1.0e-6
        P_NB = op.P_NB_total_w * mw
        losses = {
            "shine-through": op.P_shine_w * mw,
            "first-orbit": op.P_orbit_loss_w * mw,
            "charge-exchange": op.P_cx_loss_w * mw,
        }
        dominant_name = max(losses, key=losses.get)
        dominant_mw = losses[dominant_name]
        pct = 100.0 * dominant_mw / P_NB if P_NB > 0 else 0.0

        plasma = self.model.plasma
        a = plasma.minor_radius
        st_orbit = getattr(plasma, "orbit_model", "large_aspect") in ("st_meanshift", "st_pitch")
        cx_mode = getattr(plasma, "cx_model", "manual_fraction")
        orb = {}
        cx = {}
        for i, beam in enumerate(self.model.beams):
            sp = beam.species.upper()
            rho_li = physics.larmor_radius_m(beam.beam_energy_keV, plasma.toroidal_field, sp)
            f_orb = op.f_orbit_loss[i] if op.f_orbit_loss else 0.0
            if st_orbit:
                w = physics.st_orbit_widths(
                    beam.beam_energy_keV, plasma.toroidal_field, plasma.plasma_current / 1e6,
                    plasma.major_radius, a, plasma.elongation, plasma.triangularity, sp)
                orb[f"orb{i + 1}"] = (
                    f"rho_th={w['rho_theta'] / a:.2f} a,  w_pass={w['w_pass'] / a:.2f} a,  "
                    f"w_ban={w['w_ban'] / a:.2f} a,  f_c/f_t={1.0 - w['f_trap']:.2f}/{w['f_trap']:.2f},  "
                    f"f_orbit={f_orb:.3f}")
            else:
                dr = physics.passing_orbit_width(
                    beam.beam_energy_keV, plasma.toroidal_field, plasma.plasma_current / 1e6,
                    plasma.major_radius, a, plasma.elongation, sp)
                orb[f"orb{i + 1}"] = (
                    f"rho_Li={rho_li * 100:.1f} cm ({rho_li / a:.2f} a),  "
                    f"dr={dr * 100:.1f} cm ({dr / a:.2f} a),  f_orbit={f_orb:.3f}")
            f_p = op.f_cx_loss[i] if i < len(op.f_cx_loss) else 0.0
            if cx_mode == "manual_fraction":
                cx[f"cx{i + 1}"] = f"{f_p:.3f}  (flat fraction -- manual mode, no gamma_cx)"
            else:
                f_n = op.cx_f_particle[i] if i < len(op.cx_f_particle) else 0.0
                gam = op.cx_gamma[i] if i < len(op.cx_gamma) else 0.0
                esc = op.cx_escape_probability[i] if i < len(op.cx_escape_probability) else 1.0
                n0 = op.cx_n0_m3[i] if i < len(op.cx_n0_m3) else 0.0
                warn = ""
                if gam > 1.0:
                    warn += "  ** gamma_cx>1: single-pass treatment breaking down"
                if f_p > 0.10:
                    warn += "  ** f_cx,P>10%: above PPPL-1280's Zeff=1 rough ceiling"
                cx[f"cx{i + 1}"] = (
                    f"{f_p:.3f} / {f_n:.3f} / {gam:.3f}   "
                    f"(n0={n0:.2e} m^-3, escape~{esc:.1e}){warn}")
        n_line, n_gw, f_gw = self._greenwald()
        Ip_MA = plasma.plasma_current / 1e6
        if st_orbit:
            q_star = physics.safety_factor_cyl_edge_arbitrary_A(
                Ip_MA, plasma.toroidal_field, plasma.major_radius, a, plasma.elongation, plasma.triangularity)
            q_star_label = f"{q_star:.2f}  (arbitrary-A / q_a form)"
        else:
            q_star = physics.safety_factor_cyl_edge(
                Ip_MA, plasma.toroidal_field, plasma.major_radius, a, plasma.elongation)
            q_star_label = f"{q_star:.2f}  (large-aspect cylindrical form)"

        prof_on = getattr(self.model.plasma, "profile_averaging", False)
        te0 = op.Te0_keV if op.Te0_keV is not None else op.Te_keV
        ti0 = op.Ti0_keV if op.Ti0_keV is not None else op.Ti_keV

        rotation_mode = getattr(plasma, "rotation_model", "off")
        if rotation_mode == "off":
            v_phi_label = "off"
        else:
            m_i_eff = plasma.deuterium_fraction * physics.M_D + plasma.tritium_fraction * physics.M_T
            cs = np.sqrt(max((te0 or 0.0) + (ti0 or 0.0), 1e-6) * 1e3 * physics.E_CHARGE / max(m_i_eff, 1e-30))
            mach = op.v_phi_m_s / cs if cs > 0 else 0.0
            extra = f",  T_NBI={op.torque_total_Nm:.3f} N*m" if rotation_mode == "momentum_balance" else ""
            v_phi_label = f"{op.v_phi_m_s / 1e3:+.2f} km/s  (M={mach:+.2f}, c_s~{cs / 1e3:.1f} km/s){extra}"

        values = {
            "tauE": f"{self._fmt(op.tau_E_s)} / {self._fmt(op.tau_Ei_s)}"
                    + ("  (input)" if self.confinement_var.get() == "Fixed tauE (input)"
                       else f"  ({self.confinement_var.get()})"),
            "Te": self._fmt(op.Te_keV) + ("  (<T>)" if prof_on else ""),
            "Ti": self._fmt(op.Ti_keV) + ("  (<T>)" if prof_on else ""),
            "T0": (f"{self._fmt(te0)} / {self._fmt(ti0)}" if prof_on
                   else "(profile-corrected 0-D off)"),
            "ne0": self._fmt(op.ne0_m3, "{:.3e}") + ("  (on-axis)" if prof_on else ""),
            "nb0": self._fmt(op.nb0_m3, "{:.3e}"),
            "P_NB": self._fmt(P_NB), "P_shine": self._fmt(op.P_shine_w * mw),
            "P_capt": self._fmt(op.P_capt_w * mw), "P_orbit": self._fmt(op.P_orbit_loss_w * mw),
            "P_cx": self._fmt(op.P_cx_loss_w * mw), "P_useful": self._fmt(op.P_useful_w * mw),
            "cx1": cx.get("cx1", "-"), "cx2": cx.get("cx2", "-"),
            "P_e": self._fmt(op.P_e_w * mw), "P_i": self._fmt(op.P_i_w * mw),
            "P_ei": self._fmt(op.P_ei_w * mw) if self.model.plasma.enable_equipartition else "off",
            "P_alpha": self._fmt(op.P_alpha_w * mw) if self.model.plasma.alpha_heating else "off",
            "p_ecrh": (f"{plasma.p_ecrh_MW * plasma.ecrh_f_e:.2f} / "
                       f"{plasma.p_ecrh_MW * (1.0 - plasma.ecrh_f_e):.2f}   (of {plasma.p_ecrh_MW:.2f})"
                       if plasma.p_ecrh_MW > 0 else "off"),
            "p_icrh": (f"{plasma.p_icrh_MW * plasma.icrh_f_e:.2f} / "
                       f"{plasma.p_icrh_MW * plasma.icrh_f_i:.2f}   (of {plasma.p_icrh_MW:.2f})"
                       if plasma.p_icrh_MW > 0 else "off"),
            "P_heat": self._fmt((op.P_useful_w + op.P_alpha_w + op.P_aux_e_w + op.P_aux_i_w) * mw),
            "Pf_tot": self._fmt(op.pf_total_w * mw),
            "Pf_dt": f"{op.pf_thermal_w * mw:.3g} / {op.pf_beam_w * mw:.3g} / {op.pf_bb_dt_w * mw:.3g}   (= {op.pf_dt_w * mw:.3g})",
            "Pf_dd": f"{op.pf_dd_thermal_w * mw:.3g} / {op.pf_dd_beam_w * mw:.3g} / {op.pf_bb_dd_w * mw:.3g}   (= {op.pf_dd_w * mw:.3g})",
            "Y_n": self._fmt(op.neutron_rate_s, "{:.3e}"),
            "Q": self._fmt(op.pf_total_w / op.P_NB_total_w if op.P_NB_total_w else None, "{:.3g}"),
            "E_fast": self._fmt(op.avg_fast_energy_keV),
            "beta_t": self._fmt(op.beta_t * 100.0),
            "q_star": q_star_label,
            "v_phi": v_phi_label,
            "orb1": orb["orb1"], "orb2": orb["orb2"],
            "f_gw": f"{f_gw:.3f}   (<n_e>={n_line:.2e}, n_GW={n_gw:.2e} m^-3)",
            "dominant": f"{dominant_name}: {dominant_mw:.2f} MW ({pct:.0f}% of P_NB)",
        }
        for key, text in values.items():
            self.dash_values[key].configure(text=text, text_color=("gray10", "gray90"))
        if f_gw > 1.0:
            self.dash_values["f_gw"].configure(text_color="#c0504d")
        badge = self.dash_values["feasibility"]
        if op.feasible:
            badge.configure(text="FEASIBLE", text_color="#2f7d32")
        else:
            badge.configure(text="INFEASIBLE", text_color="#c0504d")

    def _render_powerflow(self, result: Result):
        op = result.op
        mw = 1.0e-6
        v = {
            "P_inj": op.P_NB_total_w * mw, "shine": op.P_shine_w * mw,
            "orbit": op.P_orbit_loss_w * mw, "cx": op.P_cx_loss_w * mw,
            "P_use": op.P_useful_w * mw, "P_e": op.P_e_w * mw, "P_i": op.P_i_w * mw,
        }
        view = self.pf_view_var.get()
        fig = self.pf_fig
        fig.clear()
        drawers = {"waterfall": self._pf_waterfall, "sankey": self._pf_sankey, "pie": self._pf_pie}
        if view == "all":
            gs = fig.add_gridspec(2, 2)
            self._pf_waterfall(fig.add_subplot(gs[0, 0]), v)
            self._pf_pie(fig.add_subplot(gs[0, 1]), v)
            self._pf_sankey(fig.add_subplot(gs[1, :]), v)
        else:
            drawers.get(view, self._pf_waterfall)(fig.add_subplot(111), v)
        closure = v["P_inj"] - (v["shine"] + v["orbit"] + v["cx"] + v["P_e"] + v["P_i"])
        fig.suptitle(f"NBI power-flow audit   (closure error {closure:+.3f} MW)", fontsize=11)
        fig.tight_layout(rect=(0, 0, 1, 0.93))
        self.pf_canvas.draw_idle()

    @staticmethod
    def _pf_waterfall(ax, v):
        P_inj, shine, orbit, cx = v["P_inj"], v["shine"], v["orbit"], v["cx"]
        P_use, P_e, P_i = v["P_use"], v["P_e"], v["P_i"]
        bars = [
            (0.0, P_inj, "#3b6fb0"),
            (P_inj - shine, shine, "#c0504d"),
            (P_inj - shine - orbit, orbit, "#c0504d"),
            (P_use, cx, "#c0504d"),
            (0.0, P_use, "#4f9d5d"),
            (0.0, P_e, "#3fa7a7"),
            (0.0, P_i, "#e0913a"),
        ]
        labels = ["$P_{inj}$", "shine-\nthrough", "first-\norbit", "CX\nloss",
                  "$P_{useful}$", "electrons", "ions"]
        for i, (bottom, height, color) in enumerate(bars):
            ax.bar(i, max(height, 0.0), bottom=bottom, width=0.62, color=color,
                   edgecolor="black", linewidth=0.5)
            ax.text(i, bottom + max(height, 0.0) + 0.02 * max(P_inj, 1e-9),
                    f"{height:.2f}", ha="center", va="bottom", fontsize=8)
        ax.set_xticks(range(len(bars)))
        ax.set_xticklabels(labels, fontsize=8)
        ax.set_ylabel("power [MW]")
        ax.set_title("Waterfall", fontsize=10)
        ax.grid(axis="y", alpha=0.3)

    @staticmethod
    def _pf_pie(ax, v):
        P_inj = v["P_inj"]
        wedges = [v["shine"], v["orbit"], v["cx"], v["P_e"], v["P_i"]]
        names = ["shine-through", "first-orbit", "charge-exchange", "electron heating", "ion heating"]
        colors = ["#c0504d", "#b03a37", "#d98b88", "#3fa7a7", "#e0913a"]
        keep = [(w, n, c) for w, n, c in zip(wedges, names, colors) if w > 1e-9]
        if keep:
            ws, ns, cs = zip(*keep)
            ax.pie(ws, labels=[f"{n}\n{w:.2f} MW" for n, w in zip(ns, ws)],
                   colors=cs, autopct="%1.0f%%", textprops={"fontsize": 8}, startangle=90)
        ax.set_title(f"Split of $P_{{inj}}$ = {P_inj:.2f} MW", fontsize=10)

    @staticmethod
    def _pf_sankey(ax, v):
        ax.axis("off")
        ax.set_title("Sankey", fontsize=10, pad=2)
        P_inj = v["P_inj"]
        if P_inj <= 1e-9:
            ax.text(0.5, 0.5, "no injected power", ha="center", va="center", fontsize=9)
            return
        raw = [("shine-through", v["shine"], 1), ("first-orbit", v["orbit"], 1),
               ("CX loss", v["cx"], 1), ("electrons", v["P_e"], 0), ("ions", v["P_i"], -1)]
        outs = [(n, w, o) for n, w, o in raw if w > 1e-3 * P_inj]
        resid = P_inj - sum(w for _, w, _ in outs)
        if abs(resid) > 1e-3 * P_inj:
            outs.append(("closure", resid, -1))
        flows = [P_inj] + [-w for _, w, _ in outs]
        labels = ["$P_{inj}$"] + [f"{n}\n{w:.2f}" for n, w, _ in outs]
        orientations = [0] + [o for _, _, o in outs]
        pathlengths = [0.35]
        up = 0
        for _, _, o in outs:
            pathlengths.append(0.08 + 0.33 * up if o == 1 else 0.35)
            up += o == 1
        try:
            sankey = Sankey(ax=ax, unit=" MW", format="%.2f", scale=1.0 / P_inj,
                            gap=0.5, radius=0.08, shoulder=0.03)
            sankey.add(flows=flows, labels=labels, orientations=orientations,
                       pathlengths=pathlengths, facecolor="#3b6fb0", edgecolor="black", lw=0.6)
            for bunch in sankey.finish():
                for text in bunch.texts + [bunch.text]:
                    text.set_fontsize(7)
        except Exception:
            ax.text(0.5, 0.5, "Sankey unavailable\nfor this operating point",
                    ha="center", va="center", fontsize=9)

    # ------------------------------------------------------------ fusion tab
    def _render_fusion(self, result: Result):
        op = result.op
        mw = 1.0e-6
        fig = self.fusion_fig
        fig.clear()
        wedges = [
            (op.pf_thermal_w * mw, "D-T thermal", "#4f9d5d"),
            (op.pf_beam_w * mw, "D-T beam-target", "#3fa7a7"),
            (op.pf_bb_dt_w * mw, "D-T beam-beam", "#8064a2"),
            (op.pf_dd_thermal_w * mw, "D-D thermal", "#e0913a"),
            (op.pf_dd_beam_w * mw, "D-D beam-target", "#c0504d"),
            (op.pf_bb_dd_w * mw, "D-D beam-beam", "#4bacc6"),
        ]
        pf_tot_mw = op.pf_total_w * mw
        # A relative (not absolute) threshold: a slice under 0.5% of the total
        # would just clutter the pie with an unreadable sliver and an
        # overlapping label -- still fully available in the Dashboard's
        # D-T/D-D breakdown rows, just not worth a wedge here.
        keep = [(w, n, c) for w, n, c in wedges if w > 0.005 * max(pf_tot_mw, 1.0e-9)]
        ax = fig.add_subplot(111)
        if keep:
            ws, ns, cs = zip(*keep)
            ax.pie(ws, labels=[f"{n}\n{w:.3g} MW" for n, w in zip(ns, ws)],
                   colors=cs, autopct="%1.1f%%", pctdistance=0.75,
                   textprops={"fontsize": 8}, startangle=90)
        else:
            ax.text(0.5, 0.5, "no fusion power\nat this operating point",
                    ha="center", va="center", fontsize=10, transform=ax.transAxes)
        q_val = op.pf_total_w / op.P_NB_total_w if op.P_NB_total_w else float("nan")
        ax.set_title(
            f"$P_{{fus}}$ = {pf_tot_mw:.3g} MW   (Q={q_val:.3g},  "
            f"$Y_n$={op.neutron_rate_s:.3g} s$^{{-1}}$)", fontsize=10)
        fig.tight_layout()
        self.fusion_canvas.draw_idle()

    # --------------------------------------------------------- deposition tab
    def _beam_chord_samples(self, beam, plasma, Te_keV, n=3000):
        """Sample the beam's real (R_t, Z_t) tangent chord through the plasma:
        distance-from-entry s, major radius R, normalised flux label rho, local
        n_e, neutral survival I(s)/I0 and the (unnormalised) fast-ion birth
        rate n_e*sigma*I(s)/I0. Uses physics.tangential_chord() so off-axis /
        vertically-shifted injection is handled consistently with the solver.
        """
        R0, a = plasma.major_radius, plasma.minor_radius
        rcp = plasma.centrepost_radius if getattr(plasma, "centrepost_radius", -1.0) > 0.0 else None
        ch = physics.tangential_chord(
            R0, a, plasma.elongation, tangent_R_m=beam.tangent_R_m,
            tangent_Z_m=beam.tangent_Z_m, R_centrepost_m=rcp, n_samples=n)
        if ch is None:
            return None
        s, rho, R = ch
        y_out = 0.5 * s[-1]
        Rt = float(beam.tangent_R_m) if beam.tangent_R_m else R0
        pn = max(plasma.density_peaking, 0.0)
        ne = plasma.central_density * np.maximum(1.0 - rho ** 2, 0.0) ** (2.0 * pn)
        A = physics.beam_mass_number(beam.species.upper())
        model = beam.shine_through_model
        model = "riviere" if model == "manual" else model
        sigma = physics.stopping_cross_section_m2(
            beam.beam_energy_keV / A, model, plasma.central_density * 1e-6,
            max(Te_keV, 1.0), plasma.effective_charge)
        ds = np.gradient(s)
        tau = np.cumsum(ne * sigma * ds)
        tau = tau - tau[0]
        survival = np.exp(-tau)
        birth = ne * sigma * survival
        i_tan = int(np.argmin(R))
        return dict(s=s, R=R, rho=rho, ne=ne, survival=survival, birth=birth,
                    ds=ds, Rt=Rt, y_out=y_out, s_tan=float(s[i_tan]),
                    blocked=bool(i_tan >= len(s) - 2),
                    f_capt=float(1.0 - np.exp(-tau[-1])))

    def _orbit_cutoff_rho(self, beam, plasma):
        """Lower rho of the prompt first-orbit-loss zone for this beam
        (1.0 -> no loss zone). Mirrors the criteria in
        physics.first_orbit_loss_fraction[_st]: gyro (1 - rho_Li/a),
        passing drift (co: none, counter: 1 - w_pass/a) and trapped-tip
        (1 - 0.5 w_ban/a -/+ 0.25 w_pass/a for co/counter)."""
        if not plasma.enable_orbit_loss:
            return 1.0
        a = plasma.minor_radius
        sp = beam.species.upper()
        Ip_MA = plasma.plasma_current / 1e6
        rho_li = physics.larmor_radius_m(beam.beam_energy_keV, plasma.toroidal_field, sp)
        co = bool(beam.co_current)
        if getattr(plasma, "orbit_model", "large_aspect") in ("st_meanshift", "st_pitch"):
            w = physics.st_orbit_widths(
                beam.beam_energy_keV, plasma.toroidal_field, Ip_MA, plasma.major_radius,
                a, plasma.elongation, plasma.triangularity, sp)
            wp, wb = w["w_pass"] / a, w["w_ban"] / a
            s = -1.0 if co else 1.0
            gyro = 1.0 - rho_li / a
            passing = 1.0 if co else 1.0 - wp
            trapped = 1.0 - 0.5 * wb - s * 0.25 * wp
            return float(max(0.0, min(gyro, passing, trapped)))
        if co:
            return 1.0
        dr = physics.passing_orbit_width(
            beam.beam_energy_keV, plasma.toroidal_field, Ip_MA, plasma.major_radius,
            a, plasma.elongation, sp)
        return float(max(0.0, 1.0 - dr / a))

    def _render_deposition(self, result: Result):
        op = result.op
        model = self.model
        plasma = model.plasma
        fig = self.dep_fig
        fig.clear()
        R0, a = plasma.major_radius, plasma.minor_radius
        Te = op.Te_keV or 1.0
        colors = ["tab:blue", "tab:orange", "tab:green", "tab:red"]
        chords = [self._beam_chord_samples(b, plasma, Te) for b in model.beams]

        # (1) beam geometry, torus top view --------------------------------
        ax = fig.add_subplot(221)
        th = np.linspace(0.0, 2.0 * np.pi, 240)
        for r, st in ((R0 - a, "-"), (R0 + a, "-"), (R0, "--")):
            ax.plot(r * np.cos(th), r * np.sin(th), st, color="0.6", lw=1.0)
        r_cp = plasma.centrepost_radius if getattr(plasma, "centrepost_radius", -1.0) > 0.0 else (R0 - a)
        ax.fill(r_cp * np.cos(th), r_cp * np.sin(th), color="0.75", zorder=0)
        ax.text(0.0, 0.0, "post", ha="center", va="center", fontsize=6.5, color="0.4")
        ax.annotate("", xy=((R0 + a) * np.cos(0.55), (R0 + a) * np.sin(0.55)),
                    xytext=((R0 + a) * np.cos(0.25), (R0 + a) * np.sin(0.25)),
                    arrowprops=dict(arrowstyle="-|>", color="0.55", lw=1.4))
        ax.text((R0 + a) * 1.03, 0.0, r"$I_p$", color="0.5", fontsize=8, va="center")
        nb = len(model.beams)
        for i, (beam, ch) in enumerate(zip(model.beams, chords)):
            c = colors[i % len(colors)]
            phi = np.deg2rad(38.0 * (i - (nb - 1) / 2.0))
            Rt = min(float(beam.tangent_R_m or R0), R0 + a)
            p = np.array([Rt * np.cos(phi), Rt * np.sin(phi)])
            d = np.array([-np.sin(phi), np.cos(phi)])
            if ch:
                s_end = ch["s"][-1]
                t_entry, t_end = -ch["s_tan"], s_end - ch["s_tan"]   # offsets from tangency
            else:
                t_entry, t_end = -0.6 * a, 0.6 * a
            # co-current: beam enters from the -d (clockwise) side and travels
            # +d (with Ip, drawn CCW); counter: mirror. Flip the whole chord so
            # p1 is always the physical entry and p1->p2 the travel direction.
            if not beam.co_current:
                t_entry, t_end = -t_end, -t_entry
                d = -d
            p1, p2 = p + d * t_entry, p + d * t_end
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=c, lw=1.5, alpha=0.9)
            ax.annotate("", xy=(p1 + (p2 - p1) * 0.55), xytext=p1,
                        arrowprops=dict(arrowstyle="-|>", color=c, lw=1.7))
            ax.plot([p[0]], [p[1]], "o", color=c, ms=4)
            if ch and ch.get("blocked"):
                ax.plot([p2[0]], [p2[1]], "x", color=c, ms=7, mew=2)
            ax.text(p2[0] * 1.06, p2[1] * 1.06,
                    f"NBI-{i + 1} ({'co' if beam.co_current else 'ctr'})"
                    + ("  blocked" if (ch and ch.get('blocked')) else ""),
                    color=c, fontsize=7, ha="center", va="center")
        lim = (R0 + a) * 1.32
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.set_title(f"{self.active_device} - beam targeting geometry (top view)",
                     fontsize=10, fontweight="bold")
        ax.set_xlabel("X [m]")
        ax.set_ylabel("Y [m]")
        ax.grid(alpha=0.2)

        # (2) attenuation + birth rate along the beam ---------------------
        ax = fig.add_subplot(222)
        ax2 = ax.twinx()
        shine_lines = []
        for i, ch in enumerate(chords):
            if ch is None:
                continue
            c = colors[i % len(colors)]
            ax.plot(ch["s"], ch["survival"], color=c, lw=1.6,
                    label=f"NBI-{i + 1} survival")
            b = ch["birth"]
            ax2.plot(ch["s"], b / (b.max() if b.max() > 0 else 1.0), color=c,
                     lw=1.2, ls=":", label=f"NBI-{i + 1} birth rate")
            f_capt = op.f_capture[i] if op.f_capture else ch["f_capt"]
            shine_lines.append(f"NBI-{i + 1}: shine-through {100.0 * (1.0 - f_capt):.1f}%  "
                               f"(capture {100.0 * f_capt:.1f}%)")
        ax.set_xlabel("distance along beam from plasma entry [m]")
        ax.set_ylabel(r"neutral survival $I(s)/I_0$")
        ax.set_ylim(0.0, 1.03)
        ax2.set_ylabel("fast-ion birth rate (norm.)")
        ax2.set_ylim(0.0, 1.05)
        ax.set_title("Beam stopping & fast-ion birth rate", fontsize=10)
        ax.grid(alpha=0.2)
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, fontsize=7, loc="upper right")
        if shine_lines:
            ax.text(0.97, 0.5, "\n".join(shine_lines), transform=ax.transAxes,
                    fontsize=7.5, va="center", ha="right",
                    bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9))

        # (3) fast-ion birth vs normalised radius -----------------------
        ax = fig.add_subplot(223)
        edges = np.linspace(0.0, 1.0, 26)
        ctr = 0.5 * (edges[:-1] + edges[1:])
        dr = edges[1] - edges[0]
        smooth = np.array([0.25, 0.5, 0.25])
        total = np.zeros_like(ctr)
        cutoffs = []
        for i, (beam, ch) in enumerate(zip(model.beams, chords)):
            if ch is None:
                continue
            c = colors[i % len(colors)]
            f_capt = op.f_capture[i] if op.f_capture else ch["f_capt"]
            f_orb = op.f_orbit_loss[i] if op.f_orbit_loss else 0.0
            f_cx = op.f_cx_loss[i] if op.f_cx_loss else plasma.cx_loss_fraction
            w_mw = beam.power_MW * f_capt * (1.0 - f_orb) * (1.0 - f_cx)
            hist, _ = np.histogram(ch["rho"], bins=edges, weights=ch["birth"] * ch["ds"])
            if hist.sum() > 0:
                hist = hist / hist.sum() * w_mw / dr
            hist = np.convolve(hist, smooth, mode="same")
            ax.plot(ctr, hist, color=c, lw=1.3, label=f"NBI-{i + 1}")
            total += hist
            rc = self._orbit_cutoff_rho(beam, plasma)
            if rc < 0.999:
                cutoffs.append(rc)
                ax.axvline(rc, color=c, ls="--", lw=1.0, alpha=0.75)
        if nb > 1:
            ax.plot(ctr, total, color="k", lw=2.0, label="total")
        if cutoffs:
            rc_min = min(cutoffs)
            ax.axvspan(rc_min, 1.0, color="0.5", alpha=0.12,
                       label=r"first-orbit-loss zone ($\rho>\rho_{cut}$)")
        ax.set_xlabel(r"$\rho = |R-R_0|/a$  (midplane)")
        ax.set_ylabel(r"deposited power  [MW per unit $\rho$]")
        ax.set_title("Fast-ion birth vs normalised radius (chord-sampled)", fontsize=10)
        ax.set_xlim(0.0, 1.0)
        ax.grid(alpha=0.2)
        ax.legend(fontsize=7)

        # (4) slowing-down distribution --------------------------------
        ax = fig.add_subplot(224)
        te = op.Te_keV or 1e-3
        vol = result.volume_m3
        for i, beam in enumerate(model.beams):
            sp = beam.species.upper()
            eb = beam.beam_energy_keV
            f_capt = op.f_capture[i] if op.f_capture else 1.0
            f_orb = op.f_orbit_loss[i] if op.f_orbit_loss else 0.0
            f_cx = op.f_cx_loss[i] if op.f_cx_loss else plasma.cx_loss_fraction
            p_use = beam.power_MW * 1e6 * f_capt * (1.0 - f_orb) * (1.0 - f_cx)
            tau_s = physics.thermalization_time(op.ne0_m3, te, eb, sp)
            nb0_i = p_use * tau_s / (eb * 1e3 * physics.E_CHARGE * max(vol, 1e-9))
            grid = np.linspace(1e-3, eb, 300)
            fE = physics.slowing_down_distribution(te, nb0_i, eb, grid, sp)
            mean_e = physics.average_fast_energy_keV(te, eb, sp)
            line, = ax.plot(grid, fE, color=colors[i % len(colors)],
                            label=fr"NBI-{i + 1} {sp} {eb:.0f} keV, $\langle E\rangle$={mean_e:.0f}, $\tau_s$={tau_s:.2g}s")
            ax.axvline(eb, color=line.get_color(), linestyle=(0, (4, 3)), linewidth=1.0, alpha=0.8)
        ax.set_xlabel("E [keV]")
        ax.set_ylabel(r"$f(E)$ [$\mathrm{m}^{-3}\,\mathrm{keV}^{-1}$]")
        ax.set_ylim(bottom=0.0)
        ax.set_title("Steady-state slowing-down distribution", fontsize=10)
        ax.grid(alpha=0.2)
        ax.legend(fontsize=7)

        fig.tight_layout()
        self.dep_canvas.draw_idle()

    # --------------------------------------------------------- references tab
    def _build_references(self, parent):
        self.ref_frame = ctk.CTkScrollableFrame(parent, label_text="References for the active models")
        self.ref_frame.pack(fill="both", expand=True, padx=4, pady=(4, 2))
        self.ref_frame.grid_columnconfigure(0, weight=1)

        bar = ctk.CTkFrame(parent)
        bar.pack(fill="x", padx=4, pady=(2, 4))
        self._make_logo(bar).pack(side="left", padx=(6, 8), pady=6)
        ctk.CTkLabel(bar, text="Spotted something strange?  Please report it.",
                     anchor="w").pack(side="left", padx=(0, 8))
        ctk.CTkButton(bar, text="Send feedback", width=110,
                      command=self._open_feedback_mail).pack(side="left", padx=4)
        ctk.CTkButton(bar, text="HI-Jass on GitHub", width=140,
                      command=lambda: webbrowser.open(REPO_URL)).pack(side="left", padx=4)
        # human-readable, scraper-unfriendly (no literal user@host anywhere)
        ctk.CTkLabel(bar, text=f"{_FEEDBACK_USER}  [at]  {_FEEDBACK_HOST}",
                     text_color="gray", font=ctk.CTkFont(size=11)).pack(side="left", padx=8)

    def _open_feedback_mail(self):
        addr = _FEEDBACK_USER + "@" + _FEEDBACK_HOST
        webbrowser.open("mailto:" + addr + "?subject=" + "HI-Jass%20feedback")

    def _make_logo(self, parent, px: int = 44):
        try:
            img = Image.open(_LOGO_PATH)
        except OSError:
            return self._draw_default_logo(parent, px)
        logo_img = ctk.CTkImage(light_image=img, dark_image=img, size=(px, px))
        label = ctk.CTkLabel(parent, image=logo_img, text="")
        label._logo_image = logo_img
        return label

    def _draw_default_logo(self, parent, px: int = 44):
        fig = Figure(figsize=(px / 100.0, px / 100.0), dpi=100)
        fig.patch.set_alpha(0.0)
        ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
        ax.set_xlim(-1.35, 1.35)
        ax.set_ylim(-1.35, 1.35)
        ax.set_aspect("equal")
        ax.axis("off")
        cols = ["#e6194B", "#3cb44b", "#4363d8"]
        for k in range(3):
            a0 = np.deg2rad(90.0 + k * 120.0 + 6.0)
            a1 = np.deg2rad(90.0 + (k + 1) * 120.0 - 30.0)
            t = np.linspace(a0, a1, 24)
            ax.plot(np.cos(t), np.sin(t), color=cols[k], lw=3.0, solid_capstyle="round")
            p = np.array([np.cos(a1), np.sin(a1)])
            tang = np.array([-np.sin(a1), np.cos(a1)])
            ax.annotate("", xy=(p + tang * 0.02), xytext=(p - tang * 0.30),
                        arrowprops=dict(arrowstyle="-|>", color=cols[k], lw=2.4))
        canvas = FigureCanvasTkAgg(fig, master=parent)
        canvas.draw()
        return canvas.get_tk_widget()

    def _active_reference_keys(self):
        """(group title -> [ref keys]) for the currently selected models."""
        ee, ei = CONFINEMENT_MODES[self.confinement_var.get()]
        orbit = ORBIT_MODELS[self.orbit_model_var.get()]
        shine = {getattr(self, f"shine_var_{i}").get() for i in range(len(self.model.beams))}
        groups: dict[str, list[str]] = {}

        groups["Machine geometry"] = ["__machine__"]

        stop = []
        if "Riviere" in shine or "Manual" in shine:
            stop.append("riviere")
        if "Janev" in shine:
            stop += ["janev_suzuki", "suzuki"]
        stop.append("stix")
        groups["Beam stopping / deposition"] = stop

        conf = []
        modes = {ee, ei}
        if "iter98y2" in modes:
            conf.append("iter98y2")
        if "kaye_nstx_lmode" in modes or "kaye_nstx_hmode" in modes:
            conf.append("kaye_nstx")
        if "neoclassical" in modes:
            conf += ["neoclassical_ch", "wesson"]
        if "neoclassical_arbA" in modes:
            conf += ["helander_arbA", "helander_sigmar", "linliu_miller",
                     "hinton_wiley", "satake_fow", "uckan"]
        if not conf:
            conf = ["iter98y2"]
        groups["Confinement model"] = conf

        if self.orbit_var.get():
            if orbit in ("st_meanshift", "st_pitch"):
                groups["First-orbit loss model"] = [
                    "akers_start", "gwb1981", "linliu_miller", "uckan", "goldston_rutherford"]
            else:
                groups["First-orbit loss model"] = ["wesson", "goldston_rutherford"]

        prim = ["bosch_hale", "nrl"]
        if self.model.plasma.enable_equipartition:
            prim.append("nrl")
        groups["Fusion & plasma primitives"] = list(dict.fromkeys(prim))

        cx_mode = CX_MODELS[self.cx_model_var.get()]
        if cx_mode == "manual_fraction":
            groups["Charge-exchange loss"] = ["pppl1280_cx"]
        else:
            groups["Charge-exchange loss"] = [
                "janev_smith_1993", "swaczyna_2024", "freeman_jones_1974", "pppl1280_cx"]

        if ROTATION_MODELS[self.rotation_model_var.get()] != "off":
            groups["Bulk toroidal rotation"] = ["wesson", "stix"]
        return groups

    def _render_references(self):
        for child in self.ref_frame.winfo_children():
            child.destroy()
        row = 0
        for title, keys in self._active_reference_keys().items():
            ctk.CTkLabel(self.ref_frame, text=title, anchor="w",
                         font=ctk.CTkFont(size=13, weight="bold")).grid(
                row=row, column=0, sticky="w", padx=6, pady=(10, 2))
            row += 1
            seen = set()
            for key in keys:
                if key == "__machine__":
                    rows = MACHINE_REFERENCES.get(
                        self.active_device, [(f"{self.active_device}: no reference on file.", "")])
                else:
                    if key in seen:
                        continue
                    seen.add(key)
                    rows = [REFERENCES.get(key, (key, ""))]
                for cite, url in rows:
                    line = ctk.CTkFrame(self.ref_frame, fg_color="transparent")
                    line.grid(row=row, column=0, sticky="ew", padx=6, pady=1)
                    line.grid_columnconfigure(0, weight=1)
                    ctk.CTkLabel(line, text="- " + cite, anchor="w", justify="left",
                                 wraplength=760).grid(row=0, column=0, sticky="w")
                    if url:
                        ctk.CTkButton(line, text="open", width=54,
                                      command=lambda u=url: webbrowser.open(u)).grid(
                            row=0, column=1, padx=(8, 0))
                    row += 1
        ctk.CTkLabel(self.ref_frame, anchor="w", text_color="gray", wraplength=760,
                     text=("Links are Google Scholar searches (they resolve to the paper) except "
                           "where a DOI is certain. Please flag anything that looks wrong via "
                           "Send feedback.")).grid(row=row, column=0, sticky="w", padx=6, pady=(12, 4))

    def _render_profiles(self, result: Result):
        model = self.model
        op = result.op
        plasma = model.plasma
        rho = model.rho_grid()
        density = model.density_profile(rho)
        prof_on = getattr(plasma, "profile_averaging", False)
        te_c = (op.Te0_keV if prof_on and op.Te0_keV is not None else op.Te_keV) or 0.0
        ti_c = (op.Ti0_keV if prof_on and op.Ti0_keV is not None else op.Ti_keV) or 0.0
        te_profile = model.temperature_profile(rho, te_c)
        ti_profile = model.temperature_profile(rho, ti_c, ion=True)

        fig = self.prof_fig
        fig.clear()
        ax_n = fig.add_subplot(241)
        ax_n.plot(rho, density / 1e20, label=r"$n_e$")
        ax_n.set(title=r"$n_e(\rho)$, $p_n=%.2f$" % plasma.density_peaking,
                 xlabel=r"$\rho$", ylabel=r"$n_e$ [$10^{20}\,\mathrm{m}^{-3}$]")
        ax_n.grid(alpha=0.3)
        ax_n.legend()

        ax_t = fig.add_subplot(242)
        ax_t.plot(rho, te_profile, label=r"$T_e$")
        ax_t.plot(rho, ti_profile, label=r"$T_i$")
        p_ti_show = plasma.temp_peaking if plasma.temp_peaking_i < 0 else plasma.temp_peaking_i
        ax_t.set(title=r"$T_e(\rho),\ T_i(\rho)$; $p_{Te}=%.2f$, $p_{Ti}=%.2f$"
                 % (plasma.temp_peaking, p_ti_show),
                 xlabel=r"$\rho$", ylabel=r"$T$ [keV]")
        ax_t.grid(alpha=0.3)
        ax_t.legend()

        # tau_S(rho): local thermalization time from the local n_e(rho), T_e(rho).
        # n_e and T_e can have very different peaking, so they don't reach their
        # (unphysical, edge-of-grid) zero at the same rho -- floor each one
        # independently as a safety net against div-by-zero, but also drop the
        # last few grid points, where whichever profile floors first would
        # otherwise show up as an artificial kink right at rho=1.
        ax_tau = fig.add_subplot(243)
        edge_cut = max(len(rho) - 4, 1)
        rho_tau = rho[:edge_cut]
        ne_floor = np.maximum(density[:edge_cut], 1.0e17)
        te_floor = np.maximum(te_profile[:edge_cut], 0.05)
        for i, beam in enumerate(model.beams):
            sp = beam.species.upper()
            eb = beam.beam_energy_keV
            tau_prof = np.array([
                physics.thermalization_time(float(ne_i), float(te_i), eb, sp)
                for ne_i, te_i in zip(ne_floor, te_floor)
            ])
            ax_tau.plot(rho_tau, tau_prof * 1e3, label=f"NBI-{i + 1} {sp} {eb:.0f} keV")
        ax_tau.set(title=r"$\tau_S(\rho)$ (thermalization time)",
                   xlabel=r"$\rho$", ylabel=r"$\tau_S$ [ms]")
        ax_tau.set_yscale("log")
        ax_tau.grid(alpha=0.3, which="both")
        ax_tau.legend(fontsize=7)

        # n0(rho): background-neutral density behind the active CX-loss model
        # (off / flat in "manual fraction"; a uniform value in "manual n0/ne";
        # the real Sec.3.2 edge-penetration profile in "penetration" -- same
        # rho~0.95 edge-shell convention as hotjass_core._cx_fraction_dicts,
        # kept in sync with it by hand).
        ax_n0 = fig.add_subplot(244)
        cx_mode = getattr(plasma, "cx_model", "manual_fraction")
        if cx_mode == "manual_fraction":
            ax_n0.text(0.5, 0.5, "CX-loss model: manual fraction\n(no n0 profile)",
                       ha="center", va="center", transform=ax_n0.transAxes,
                       fontsize=9, color="0.45")
            ax_n0.set_xticks([])
            ax_n0.set_yticks([])
        elif cx_mode == "manual_n0":
            n0_flat = max(plasma.cx_n0_over_ne, 0.0) * op.ne0_m3
            ax_n0.axhline(n0_flat, color="tab:blue", lw=1.6)
            ax_n0.set(title=r"$n_0$ (manual, uniform)",
                      xlabel=r"$\rho$", ylabel=r"$n_0$ [m$^{-3}$]")
            ax_n0.set_xlim(0.0, 1.0)
            ax_n0.grid(alpha=0.3)
        else:  # "penetration"
            edge_idx = max(int(0.95 * (len(rho) - 1)), 0)
            n0_lcfs = max(plasma.cx_n0_lcfs_over_ne, 0.0) * density[edge_idx]
            n0_profile = physics.neutral_penetration_profile(
                rho, density, op.Te_keV or 0.0, op.Ti_keV or 0.0, plasma.minor_radius, n0_lcfs)
            ax_n0.plot(rho, n0_profile, color="tab:blue")
            ax_n0.axvline(rho[edge_idx], color="0.6", ls="--", lw=1.0)
            ax_n0.text(rho[edge_idx], n0_lcfs, "  LCFS ref.\n  ($\\rho$=%.2f)" % rho[edge_idx],
                       fontsize=7, color="0.4", va="bottom")
            ax_n0.set(title=r"$n_0(\rho)$ (penetration)",
                      xlabel=r"$\rho$", ylabel=r"$n_0$ [m$^{-3}$]")
            ax_n0.set_yscale("log")
            ax_n0.grid(alpha=0.3, which="both")

        ax_s = fig.add_subplot(245)
        theta = np.linspace(0, 2 * np.pi, 400)
        delta = np.clip(plasma.triangularity, -0.999, 0.999)
        ax_s.plot(plasma.major_radius + plasma.minor_radius * np.cos(theta + np.arcsin(delta) * np.sin(theta)),
                  plasma.elongation * plasma.minor_radius * np.sin(theta))
        ax_s.plot(plasma.major_radius, 0.0, marker="+", color="tab:red", markersize=9, markeredgewidth=1.5)
        ax_s.set(title=r"Shape: $R_0=%.2f$, $a=%.2f$, $\kappa=%.2f$, $\delta=%.2f$"
                 % (plasma.major_radius, plasma.minor_radius, plasma.elongation, plasma.triangularity),
                 xlabel=r"$R$ [m]", ylabel=r"$Z$ [m]")
        ax_s.set_aspect("equal")
        ax_s.grid(alpha=0.3)

        # P_fus(rho): local D-T + D-D fusion power density, thermal vs beam-plasma
        # (beam-target). Densities are converted from the OperatingPoint's
        # volume-balance nD0/nT0 to their on-axis value (pk_n) exactly as
        # solve.py does internally for the volume-integrated totals -- this
        # panel just exposes the same local density per rho instead.
        ax_pf = fig.add_subplot(246)
        pk_n = 1.0 + 2.0 * max(plasma.density_peaking, 0.0) if prof_on else 1.0
        nD0_axis = op.nD0_m3 * pk_n
        nT0_axis = op.nT0_m3 * pk_n
        ne_axis = op.ne0_m3
        p_te = plasma.temp_peaking
        th_total = (
            physics.thermal_fusion_power_density_profile(rho, nD0_axis, nT0_axis, ti_c, plasma.density_peaking, p_ti_show)
            + physics.thermal_dd_power_density_profile(rho, nD0_axis, ti_c, plasma.density_peaking, p_ti_show)
        )
        bt_total = np.zeros_like(rho)
        vol = result.volume_m3
        for i, beam in enumerate(model.beams):
            sp = beam.species.upper()
            eb = beam.beam_energy_keV
            f_capt = op.f_capture[i] if op.f_capture else 1.0
            f_orb = op.f_orbit_loss[i] if op.f_orbit_loss else 0.0
            f_cx = op.f_cx_loss[i] if op.f_cx_loss else plasma.cx_loss_fraction
            p_use = beam.power_MW * 1e6 * f_capt * (1.0 - f_orb) * (1.0 - f_cx)
            tau_s0 = physics.thermalization_time(ne_axis, te_c, eb, sp)
            nb0_axis = p_use * tau_s0 / (eb * 1e3 * physics.E_CHARGE * max(vol, 1e-9))
            target_n = nT0_axis if sp == "D" else nD0_axis
            bt_total += physics.beam_target_power_density_profile(
                rho, nb0_axis, target_n, te_c, eb, sp, ne_axis, plasma.density_peaking, p_te)
            if sp == "D" and nD0_axis > 0.0:
                bt_total += physics.beam_target_dd_power_density_profile(
                    rho, nb0_axis, nD0_axis, te_c, eb, ne_axis, plasma.density_peaking, p_te)
        ax_pf.plot(rho, th_total / 1e3, label="thermal")
        ax_pf.plot(rho, bt_total / 1e3, label="beam-plasma")
        ax_pf.set(title=r"$P_{fus}(\rho)$  (D-T + D-D)",
                  xlabel=r"$\rho$", ylabel=r"$P_{fus}$ [kW/m$^3$]")
        ax_pf.grid(alpha=0.3)
        ax_pf.legend(fontsize=7)

        ax_f = fig.add_subplot(247)
        ax_f.axis("off")
        ax_f.text(0, 0.85, "Profiles / 0-D balance", fontsize=11, weight="bold")
        ax_f.text(0, 0.60, r"$n_e(\rho)=n_{e0}(1-\rho^2)^{2p_n}$", fontsize=9)
        ax_f.text(0, 0.40, r"$T_{e,i}(\rho)=T_{e0,i0}(1-\rho^2)^{2p_T}$", fontsize=9)
        ax_f.text(0, 0.20, r"$T_e$, $T_i$ solved at the central $n_e$", fontsize=9)

        fig.tight_layout()
        self.prof_canvas.draw_idle()

    def _render_scan_plot(self):
        if not getattr(self, "scan", None):
            return
        selected = self.observable_var.get()
        keys = self.OBSERVABLES[self.DISPLAY_GROUPS[selected]]
        fig = self.scan_fig
        fig.clear()
        density_axis = self.scan["n_e"] / 1e20
        if len(keys) > 2:
            axes_list = list(fig.subplots(2, 2, squeeze=False).flat)
            for ax, key in zip(axes_list, keys):
                finite = np.isfinite(self.scan[key])
                ax.plot(density_axis[finite], self.scan[key][finite], marker="o", ms=3,
                        label=self.LATEX_NAMES[key])
                ax.set_title(f"{self.LATEX_NAMES[key]} [${self.LATEX_UNITS[self.UNITS[key]]}$]")
                ax.set_xlabel(r"$n_e$ [$10^{20}\,\mathrm{m}^{-3}$]")
                ax.grid(alpha=0.3)
                ax.legend()
                if key == "Y_neutron_bb" and not getattr(self.model.plasma, "enable_beam_beam", False):
                    ax.set_ylim(-1.0, 1.0)
                    ax.text(0.5, 0.5, "beam-beam fusion off\n(Losses section)", ha="center", va="center",
                            transform=ax.transAxes, fontsize=9, color="0.45")
            # 4th (otherwise unused) panel of the n_D/n_T/n_b group: the
            # fast-ion / target-ion density ratio -- n_Target is whichever
            # thermal D-T species the DOMINANT NBI beam reacts with (n_T for
            # a D beam, n_D for a T beam), i.e. the species driving the
            # beam-target fusion yield.
            if len(axes_list) > len(keys) and set(keys) == {"n_D", "n_T", "n_b"}:
                major_beam = max(self.model.beams, key=lambda b: b.power_MW, default=None)
                major_species = major_beam.species.upper() if major_beam else "D"
                target_key = "n_T" if major_species == "D" else "n_D"
                target = self.scan[target_key]
                with np.errstate(divide="ignore", invalid="ignore"):
                    ratio = np.where(target > 0.0, self.scan["n_b"] / np.where(target > 0.0, target, 1.0), np.nan)
                ax4 = axes_list[len(keys)]
                finite = np.isfinite(ratio)
                target_sub = target_key[-1]
                ax4.plot(density_axis[finite], ratio[finite], marker="o", ms=3, color="tab:purple",
                         label=fr"$n_{{b0}}/n_{{{target_sub}}}$")
                ax4.set_title(fr"$n_{{b0}}/n_{{{target_sub}}}$  (fast-ion / target-ion ratio)")
                ax4.set_xlabel(r"$n_e$ [$10^{20}\,\mathrm{m}^{-3}$]")
                ax4.grid(alpha=0.3)
                ax4.legend()
        else:
            ax = fig.add_subplot(111)
            for key in keys:
                finite = np.isfinite(self.scan[key])
                ax.plot(density_axis[finite], self.scan[key][finite], marker="o", ms=3,
                        label=f"{self.LATEX_NAMES[key]} [${self.LATEX_UNITS[self.UNITS[key]]}$]")
            ax.set_xlabel(r"$n_e$ [$10^{20}\,\mathrm{m}^{-3}$]")
            ax.grid(alpha=0.3)
            ax.legend()
            if set(keys) == {"Te", "Ti"}:
                all_vals = np.concatenate([self.scan[k][np.isfinite(self.scan[k])] for k in keys])
                top = min(float(all_vals.max()), 100.0) if all_vals.size else 100.0
                ax.set_ylim(0.0, top)
        vmin = self.scan["n_e_valid_min"][0]
        vmax = self.scan["n_e_valid_max"][0]
        plasma = self.model.plasma
        equip_on = plasma.enable_equipartition
        alpha_on = plasma.alpha_heating
        eq_sensitive = selected in self.EQUIP_SENSITIVE
        al_sensitive = selected in self.ALPHA_SENSITIVE
        suptitle = r"valid $n_e=[%.2e, %.2e]$ m$^{-3}$" % (vmin, vmax)
        tags = []
        if eq_sensitive:
            tags.append("equipartition ON (Te, Ti coupled)" if equip_on else "equipartition off")
        if al_sensitive:
            tags.append(f"alpha self-heating ON (f={plasma.f_alpha:.2g})" if alpha_on
                        else "alpha self-heating off")
        if tags:
            suptitle += "\n" + " | ".join(tags)
        fig.suptitle(suptitle, fontsize=10)
        fig.tight_layout()
        self.scan_canvas.draw_idle()

        self.scan_text.delete("1.0", "end")
        eq_tag = "ON (Te,Ti coupled)" if equip_on else "off"
        al_tag = f"ON (f={plasma.f_alpha:.2g})" if alpha_on else "off"
        self.scan_text.insert("end", f"e-i equipartition: {eq_tag}"
                              f"{'  <- affects this group' if eq_sensitive else ''}\n")
        self.scan_text.insert("end", f"alpha self-heating: {al_tag}"
                              f"{'  <- affects this group' if al_sensitive else ''}\n\n")
        self.scan_text.insert("end", "Requested n_e: %.3e .. %.3e m^-3\n"
                              % (self.model.plasma.n_e_min, self.model.plasma.n_e_max))
        self.scan_text.insert("end", "Valid range: %.3e .. %.3e m^-3\n\n" % (vmin, vmax))
        for key in keys:
            vals = self.scan[key][np.isfinite(self.scan[key])]
            if vals.size:
                self.scan_text.insert("end", f"{key} [{self.UNITS[key]}]: {vals.min():.4g} .. {vals.max():.4g}\n")
            else:
                self.scan_text.insert("end", f"{key} [{self.UNITS[key]}]: no feasible values\n")

    def _render_summary(self, result: Result):
        scan = result.scan
        plasma = self.model.plasma
        fig = self.sum_fig
        fig.clear()
        axes = fig.subplots(3, 4, squeeze=False)
        density_axis = scan["n_e"] / 1e20
        theta = np.linspace(0, 2 * np.pi, 400)
        delta = np.clip(plasma.triangularity, -0.999, 0.999)
        panels = [
            ("Geometry [m]", lambda ax: (
                ax.plot(plasma.major_radius + plasma.minor_radius * np.cos(theta + np.arcsin(delta) * np.sin(theta)),
                        plasma.elongation * plasma.minor_radius * np.sin(theta)),
                ax.plot(plasma.major_radius, 0.0, marker="+", color="tab:red", markersize=8))),
            ("Temperature [keV]", lambda ax: (ax.plot(density_axis, scan["Te"], label=r"$T_e$"),
                                              ax.plot(density_axis, scan["Ti"], label=r"$T_i$"))),
            ("Heating [MW]", lambda ax: (ax.plot(density_axis, scan["P_e"], label=r"$P_e$"),
                                         ax.plot(density_axis, scan["P_i"], label=r"$P_i$"))),
            ("Q = P_fus / P_NB", lambda ax: ax.plot(density_axis, scan["Q"])),
            ("Losses [MW]", lambda ax: (ax.plot(density_axis, scan["P_shine-through"], label=r"$P_{shine}$"),
                                        ax.plot(density_axis, scan["P_orbit"], label=r"$P_{orbit}$"),
                                        ax.plot(density_axis, scan["P_cx"], label=r"$P_{cx}$"),
                                        ax.plot(density_axis, scan["P_lost"], label=r"$P_{lost}$", color="k"))),
            (r"Species n [$\mathrm{m}^{-3}$]", lambda ax: (ax.plot(density_axis, scan["n_D"], label=r"$n_D$"),
                                                           ax.plot(density_axis, scan["n_T"], label=r"$n_T$"),
                                                           ax.plot(density_axis, scan["n_b"], label=r"$n_b$"))),
            ("Fusion power [MW]", lambda ax: (ax.plot(density_axis, scan["Pf_tot"], label=r"$P_{f,tot}$"),
                                              ax.plot(density_axis, scan["Pf_th"], label=r"$P_{f,th}$"),
                                              ax.plot(density_axis, scan["Pf_b"], label=r"$P_{f,b}$"),
                                              ax.plot(density_axis, scan["Pf_bb"], label=r"$P_{f,bb}$"))),
            ("Times [s]", lambda ax: (ax.plot(density_axis, scan["tau_S"], label=r"$\tau_S$"),
                                      ax.plot(density_axis, scan["tauE_e"], label=r"$\tau_{E,e}$"),
                                      ax.plot(density_axis, scan["tauE_i"], label=r"$\tau_{E,i}$"),
                                      ax.plot(density_axis, scan["tau_IE"], label=r"$\tau_{IE}$"))),
            ("Neutron rate [1/s]", lambda ax: (ax.plot(density_axis, scan["Y_neutron_b"], label=r"$Y_{beam}$"),
                                               ax.plot(density_axis, scan["Y_neutron_th"], label=r"$Y_{th}$"),
                                               ax.plot(density_axis, scan["Y_neutron_bb"], label=r"$Y_{bb}$"),
                                               ax.plot(density_axis, scan["Y_neutron"], label=r"$Y_{full}$", color="k"))),
            ("R = U_fast / U_th", lambda ax: ax.plot(density_axis, scan["R"])),
            ("Pressure [Pa]", lambda ax: (ax.plot(density_axis, scan["Pr_th"], label=r"$p_{th}$"),
                                          ax.plot(density_axis, scan["Pr_fast"], label=r"$p_{fast}$"))),
            ("Toroidal beta [%]", lambda ax: ax.plot(density_axis, scan["beta_T"])),
        ]
        for ax, (title, draw) in zip(axes.flat, panels):
            draw(ax)
            ax.set_title(title, fontsize=9)
            ax.grid(alpha=0.25)
            ax.tick_params(labelsize=7)
            if title == "Geometry [m]":
                ax.set_aspect("equal")
            else:
                ax.set_xlabel(r"$n_e$ [$10^{20}\,\mathrm{m}^{-3}$]", fontsize=7)
            if title in ("Temperature [keV]", "Heating [MW]", r"Species n [$\mathrm{m}^{-3}$]",
                         "Fusion power [MW]", "Times [s]", "Pressure [Pa]", "Losses [MW]", "Neutron rate [1/s]"):
                ax.legend(fontsize=6, ncol=2)
        fig.text(0.02, 0.005, self._summary_parameters(), fontsize=7, va="bottom", family="monospace")
        fig.subplots_adjust(left=0.06, right=0.98, top=0.95, bottom=0.22, wspace=0.3, hspace=0.42)
        self.sum_canvas.draw_idle()

    # --------------------------------------------------------------- helpers
    def _greenwald(self, n_e0: float | None = None):
        """(<n_e>_line, n_GW, <n_e>/n_GW) in m^-3, for on-axis density n_e0
        (default: central_density).  <n_e> = n_e0 * L, L = line-average of the
        (1-rho^2)^(2 p_n) profile shape.  n_GW = Ip[MA]/(pi a^2) x 1e20.
        """
        p = self.model.plasma
        n_e0 = p.central_density if n_e0 is None else float(n_e0)
        n_gw = (p.plasma_current / 1e6) / (np.pi * max(p.minor_radius, 1e-6) ** 2) * 1e20
        rho = np.linspace(0.0, 1.0, 201)
        shape = float(np.mean(np.maximum(1.0 - rho ** 2, 0.0) ** (2.0 * max(p.density_peaking, 0.0))))
        n_line = n_e0 * shape
        return n_line, n_gw, (n_line / n_gw if n_gw > 0 else float("inf"))

    def _assess(self, model: HotJassModel, op) -> str:
        plasma = model.plasma
        aspect = plasma.major_radius / max(plasma.minor_radius, 1e-6)
        lines = ["ASSUMPTIONS & VALIDITY", "=" * 40, ""]
        lines.append(f"Device: {self.active_device}    aspect ratio A = R0/a = {aspect:.2f}")
        n_line, n_gw, f_gw = self._greenwald()
        lines.append(f"Greenwald (at ne_c={plasma.central_density:.3g} on-axis):"
                     f"  <n_e>_line = n_e0*L = {n_line:.3e}  (L = line-avg of (1-rho^2)^2p_n),"
                     f"  n_GW = Ip/(pi a^2) = {n_gw:.3e} m^-3,  <n_e>/n_GW = {f_gw:.3f}"
                     + ("   ! above the Greenwald limit" if f_gw > 1.0 else ""))
        lines.append(f"Confinement model: {self.confinement_var.get()}")
        lines.append(f"  tau_Ee_mode = {plasma.tau_Ee_mode}, tau_Ei_mode = {plasma.tau_Ei_mode}")
        if plasma.tau_Ee_mode == "fixed":
            lines.append(f"  fixed tauE,e = {plasma.tauE_e:.4g} s")
        if plasma.tau_Ei_mode == "fixed":
            lines.append(f"  fixed tauE,i = {plasma.tauE_i:.4g} s")
        dirs = ", ".join(f"NBI-{i + 1} {'co' if b.co_current else 'counter'}-current"
                         for i, b in enumerate(model.beams))
        orbit_model = getattr(plasma, "orbit_model", "large_aspect")
        st_orbit = orbit_model in ("st_meanshift", "st_pitch")
        orbit_label = {
            "large_aspect": "large-aspect (q* rho_Li, direction-only mean-shift)",
            "st_meanshift": "ST orbits - mean-shift (arbitrary A: q_a, banana + gyro channels)",
            "st_pitch": "ST orbits - pitch-resolved (co/counter/trapped classes)",
        }.get(orbit_model, orbit_model)
        lines.append(f"First-orbit loss: {'ON' if plasma.enable_orbit_loss else 'off'}  ({dirs})")
        lines.append(f"  orbit model: {orbit_label}")
        a = plasma.minor_radius
        wide = False
        for i, b in enumerate(model.beams):
            sp = b.species.upper()
            rho_li = physics.larmor_radius_m(b.beam_energy_keV, plasma.toroidal_field, sp)
            f_orb = op.f_orbit_loss[i] if op.f_orbit_loss else 0.0
            if st_orbit:
                w = physics.st_orbit_widths(
                    b.beam_energy_keV, plasma.toroidal_field, plasma.plasma_current / 1e6,
                    plasma.major_radius, a, plasma.elongation, plasma.triangularity, sp)
                f_c = 1.0 - w["f_trap"]
                lines.append(
                    f"  NBI-{i + 1}: rho_Li = {rho_li / a:.2f} a,  rho_theta = {w['rho_theta'] / a:.2f} a,  "
                    f"w_pass = {w['w_pass'] / a:.2f} a,  w_ban = {w['w_ban'] / a:.2f} a,")
                lines.append(
                    f"          q_a = {w['q_a']:.1f},  passing f_c = {f_c:.2f} / trapped f_t = {w['f_trap']:.2f},  "
                    f"f_orbit = {f_orb:.3f}")
                wide = wide or w["w_ban"] / a > 0.5
            else:
                dr = physics.passing_orbit_width(
                    b.beam_energy_keV, plasma.toroidal_field, plasma.plasma_current / 1e6,
                    plasma.major_radius, a, plasma.elongation, sp)
                lines.append(f"  NBI-{i + 1}: rho_Li = {rho_li * 100:.1f} cm ({rho_li / a:.2f} a), "
                             f"passing orbit dr = {dr * 100:.1f} cm ({dr / a:.2f} a), f_orbit = {f_orb:.3f}")
                wide = wide or rho_li / a > 0.3 or dr / a > 0.3
        if st_orbit:
            lines.append("  ST orbit models: widths from the poloidal gyroradius rho_theta (not q*rho_Li).")
            lines.append("  Drift shift is inward (co) / outward (counter) for both passing and trapped;")
            lines.append("  the gyro channel (born within 1 rho_Li of the LCFS) is direction-independent, so")
            lines.append("  co f_orbit is non-zero but < counter (split ~1.5x, not 0 vs 1).")
            lines.append("  Birth profile: same Beer-Lambert n_e*sigma*exp(-n_e*sigma*x) as the large-aspect")
            lines.append("  model and the Deposition-tab chord plot (stopping_cross_section_m2).")
            lines.append("  Refs: Akers NF (START NBI); Goldston-White-Boozer PRL 47 (1981); Goldston & Rutherford (1995).")
        elif wide:
            lines.append("  ! orbit width >~ 0.3 a: the large-aspect estimate is crude here -- try an ST orbit model.")
        cx_mode = getattr(plasma, "cx_model", "manual_fraction")
        if cx_mode == "manual_fraction":
            lines.append(f"CX loss: manual fraction = {plasma.cx_loss_fraction:.3g} (flat efficiency knob)")
        else:
            cx_label = "manual n0/ne" if cx_mode == "manual_n0" else "penetration (n0_LCFS/ne)"
            ratio = plasma.cx_n0_over_ne if cx_mode == "manual_n0" else plasma.cx_n0_lcfs_over_ne
            ratio_label = "n0/ne" if cx_mode == "manual_n0" else "n0_LCFS/ne"
            lines.append(f"CX loss: {cx_label}, {ratio_label}={ratio:.2g}  "
                         f"(docs/CX_model.tex Sec.2 survival integral; gamma_cx / f_cx,N / f_cx,P per beam on the Dashboard)")
            lines.append("  sigma_cx: Janev & Smith (1993) Sec.2.3.1 form, refit by Swaczyna, Bzowski & Kubiak")
            lines.append("  (arXiv:2411.13174, 2024, Eq. A1) -- 5% accuracy 1-100 keV/amu.")
            if cx_mode == "penetration":
                lines.append("  n0(rho) decay length: INTERIM approximate ionization/CX rate coefficients, not an")
                lines.append("  independently verified ADAS/Voronov fit -- order-of-magnitude only.")
            lines.append("  Escape probability (re-ionization before the wall) is shown but NOT applied: every")
            lines.append("  CX event is conservatively treated as a full loss, per the doc's default.")
        rot_mode = getattr(plasma, "rotation_model", "off")
        if rot_mode == "off":
            lines.append("Bulk toroidal rotation: off (v_phi=0, beam-target reactivity uses a stationary target)")
        else:
            if rot_mode == "manual":
                lines.append(f"Bulk toroidal rotation: manual, v_phi={plasma.manual_v_phi_m_s:.3g} m/s (+=co-current)")
            else:
                lines.append(f"Bulk toroidal rotation: momentum balance, tau_phi={plasma.tau_phi_over_tauEi:.2g}*tau_E,i "
                             f"(v_phi and net NBI torque shown on the Dashboard)")
            lines.append("  Only the beam-target reactivity's relative velocity sees v_phi (per-beam co/counter")
            lines.append("  direction) -- shine-through, first-orbit and CX loss are unaffected (v_phi << v_beam).")
            lines.append("  No validated tau_phi (momentum confinement time) scaling exists; tau_phi~tau_E,i is a")
            lines.append("  crude, commonly-used stand-in, so treat v_phi's absolute scale as order-of-magnitude.")
        if getattr(plasma, "enable_beam_beam", False):
            lines.append("Beam-beam fusion: ON -- reduced, monoenergetic-population estimate (NBI-1 x NBI-2 only,")
            lines.append("  each at its own average fast-ion energy); NOT a full 2-D slowing-down-spectrum integral.")
            lines.append("  Folded into the D-T/D-D fusion power and neutron-rate totals above and into the")
            lines.append("  Scan tab's Neutrons group; T-T pairs (no cross-section fit here) contribute nothing.")
        else:
            lines.append("Beam-beam fusion: off")
        lines.append(f"e-i equipartition: {'ON (coupled Te/Ti solve)' if plasma.enable_equipartition else 'off (decoupled Te, Ti)'}")
        if getattr(plasma, "profile_averaging", False):
            p_ti = plasma.temp_peaking if plasma.temp_peaking_i < 0 else plasma.temp_peaking_i
            pkn = 1.0 + 2.0 * max(plasma.density_peaking, 0.0)
            pkte = 1.0 + 2.0 * max(plasma.temp_peaking, 0.0)
            pkti = 1.0 + 2.0 * max(p_ti, 0.0)
            lines.append(f"Profile-corrected 0-D: ON  -- balance runs on <n_e> = n_e0 / {pkn:.2f} "
                         f"(density peaking {plasma.density_peaking:.2g}).")
            lines.append(f"  Te0 = <Te> x {pkte:.2f} (p_Te={plasma.temp_peaking:.2g}),  "
                         f"Ti0 = <Ti> x {pkti:.2f} (p_Ti={p_ti:.2g}).")
            lines.append(f"  Te0 = {op.Te0_keV:.2f} keV, Ti0 = {op.Ti0_keV:.2f} keV,  "
                         f"neutrons = {op.neutron_rate_s:.2e} n/s.  Fusion / pressure use the")
            lines.append(f"  on-axis values with the separate (1-rho^2) Te and Ti profile shapes.")
        else:
            lines.append("Profile-corrected 0-D: off  -- n_e0 and T treated as uniform "
                         "(reported T is a flat-plasma effective value).")
        if op.pf_dd_w > 0.0:
            lines.append(f"D-D fusion: {op.pf_dd_w * 1e-6:.3g} MW  "
                         f"(thermal {op.pf_dd_thermal_w * 1e-6:.3g} + beam-target {op.pf_dd_beam_w * 1e-6:.3g}); "
                         f"D-T {op.pf_dt_w * 1e-6:.3g} MW.")
        if plasma.alpha_heating:
            lines.append(f"Alpha self-heating: ON  (f_alpha={plasma.f_alpha:.3g}; "
                         f"P_a = f_alpha*(3.5/17.6)*P_fus, fixed-point, e/i split by slowing-down)")
        else:
            lines.append("Alpha self-heating: off")
        aux_e = plasma.p_ecrh_MW * plasma.ecrh_f_e + plasma.p_icrh_MW * plasma.icrh_f_e
        aux_i = plasma.p_ecrh_MW * (1.0 - plasma.ecrh_f_e) + plasma.p_icrh_MW * plasma.icrh_f_i
        if plasma.p_ecrh_MW > 0 or plasma.p_icrh_MW > 0:
            lines.append(f"Aux heating (prescribed source added to the Te/Ti balance):")
            lines.append(f"  ECRH {plasma.p_ecrh_MW:.3g} MW  (f_e={plasma.ecrh_f_e:.3g}),  "
                         f"ICRH {plasma.p_icrh_MW:.3g} MW  (f_e={plasma.icrh_f_e:.3g}, f_i={plasma.icrh_f_i:.3g})")
            lines.append(f"  -> electrons {aux_e:.3g} MW,  ions {aux_i:.3g} MW")
            if plasma.p_icrh_MW > 0 and abs(plasma.icrh_f_e + plasma.icrh_f_i - 1.0) > 1e-6:
                nt = plasma.p_icrh_MW * (1.0 - plasma.icrh_f_e - plasma.icrh_f_i)
                lines.append(f"  ! ICRH f_e+f_i = {plasma.icrh_f_e + plasma.icrh_f_i:.3g} "
                             f"({nt:+.3g} MW not thermalised)")
        else:
            lines.append("Aux heating (ECRH/ICRH): off")
        lines.append("")
        lines.append("Per-beam shine-through model:")
        for i, beam in enumerate(model.beams):
            e_per_amu = beam.beam_energy_keV / physics.beam_mass_number(beam.species.upper())
            lines.append(f"  NBI-{i + 1}: {beam.shine_through_model}   E/A = {e_per_amu:.0f} keV/amu")
            if beam.shine_through_model == "janev_suzuki" and not (100.0 <= e_per_amu <= 1.0e4):
                lines.append("    ! Janev-Suzuki fit valid 100-1e4 keV/amu; input clipped.")
        lines.append("")
        lines.append("Fit ranges / notes:")
        lines.append("  Bosch-Hale DT reactivity valid Ti = 0.2-100 keV.")
        lines.append("  Riviere / Janev-Suzuki stopping: order-of-magnitude fits (see physics.py).")
        if "kaye_nstx_lmode" in (plasma.tau_Ee_mode, plasma.tau_Ei_mode):
            lines.append("  Kaye NSTX L-mode: an ST-appropriate fit (A ~ 1.3-1.5 dataset).")
        if "kaye_nstx_hmode" in (plasma.tau_Ee_mode, plasma.tau_Ei_mode):
            lines.append("  Kaye NSTX H-mode (Kaye 2006, NF 46 848, 'Case 1'): ST fit,")
            lines.append("  tauE ~ Ip^0.57 Bt^1.08 ne^0.44 P^-0.73 (A ~ 1.3-1.5 dataset).")
            lines.append("  P_loss for the scaling = useful beam power + ECRH + ICRH (no alpha).")
        if "iter98y2" in (plasma.tau_Ee_mode, plasma.tau_Ei_mode):
            lines.append("  IPB98(y,2): ELMy H-mode fit to conventional-A devices (A ~ 2.5-4);")
            if aspect < 2.0:
                lines.append(f"  ! A = {aspect:.2f} < 2 -- IPB98(y,2) is being extrapolated well")
                lines.append("    outside its dataset here; a fixed tauE or an ST fit is safer.")
            lines.append("  P_loss for the scaling = useful beam power + ECRH + ICRH (no alpha).")
        if plasma.tau_Ei_mode in ("neoclassical", "neoclassical_arbA"):
            lines.append("  Neoclassical ion transport has no anomalous channel -> Ti can be")
            lines.append("  large / implausible; treat as a lower bound on ion transport.")
        if plasma.tau_Ei_mode == "neoclassical_arbA":
            lines.append("  Arbitrary-A ion channel: chi_i ~ q_a^2 rho_i^2 nu_ii (f_t/f_c) with")
            lines.append("  f_t = Lin-Liu & Miller, q_a = Uckan (Helander PoP 7 (2000); Hinton-Wiley")
            lines.append("  PRL 29 (1972); Satake PoP 9 (2002)). NOT a smooth reduction of the")
            lines.append("  eps^-1.5 form -- can differ by ~1 order of magnitude either way; run both.")
        if aspect < 2.0 and plasma.tau_Ee_mode == "fixed":
            lines.append("  A < 2: conventional-aspect confinement scalings (IPB98) would be")
            lines.append("  extrapolating here; a fixed tauE input sidesteps that.")
        lines.append("")
        lines.append("OPERATING POINT")
        lines.append("-" * 40)
        lines.append(f"feasible: {op.feasible}")
        if not op.feasible:
            lines.append(f"reason: {op.infeasible_reason}")
        else:
            lines.append(f"n_thermal / n_sum = {op.n_thermal_fraction:.3f}")
            if op.n_thermal_fraction < 0.15:
                lines.append("  ! near the feasibility edge: Ti is very sensitive here.")
            if op.Ti_keV and op.Ti_keV > 100.0:
                lines.append(f"  ! Ti = {op.Ti_keV:.0f} keV is implausibly high.")
        return "\n".join(lines)

    def _tau_used_str(self, scan_key: str, mode: str, manual_val: float) -> str:
        """String for an actually-used confinement time: the fixed input value,
        or the min..max range over the current scan for a physics-based mode.
        """
        if mode == "fixed":
            return f"{manual_val:.3g} s (fixed input)"
        sc = getattr(self, "scan", None)
        if sc is not None and scan_key in sc:
            v = sc[scan_key][np.isfinite(sc[scan_key])]
            if v.size:
                return f"[{v.min():.3g} .. {v.max():.3g}] s ({mode})"
        return f"— ({mode})"

    def _summary_parameters(self) -> str:
        p = self.model.plasma
        equip_state = "ON" if p.enable_equipartition else "off"
        alpha_state = f"ON (f_alpha={p.f_alpha:.3g})" if p.alpha_heating else "off"
        orbit_model_short = {"large_aspect": "large-A", "st_meanshift": "ST mean-shift",
                             "st_pitch": "ST pitch-resolved"}.get(
            getattr(p, "orbit_model", "large_aspect"), getattr(p, "orbit_model", ""))
        orbit_state = f"ON [{orbit_model_short}]" if p.enable_orbit_loss else "off"
        prof_state = ("ON (n_e0 -> <n_e>; <T> + T0 reported)"
                      if getattr(p, "profile_averaging", False) else "off (uniform n_e, T)")
        _, n_gw, f_gw_c = self._greenwald()
        f_gw_lo = self._greenwald(p.n_e_min)[2]
        f_gw_hi = self._greenwald(p.n_e_max)[2]
        tau_e_str = self._tau_used_str("tauE_e", p.tau_Ee_mode, p.tauE_e)
        tau_i_str = self._tau_used_str("tauE_i", p.tau_Ei_mode, p.tauE_i)
        lines = [
            f"Device: {self.active_device}    R0={p.major_radius:.3g} m    a={p.minor_radius:.3g} m    "
            f"k={p.elongation:.3g}    delta={p.triangularity:.3g}\n"
            f"B0={p.toroidal_field:.3g} T    Ip={p.plasma_current / 1e6:.3g} MA    "
            f"Zeff={p.effective_charge:.3g}    ne_c(on-axis)={p.central_density:.3g}    "
            f"ne_scan=[{p.n_e_min:.3g}, {p.n_e_max:.3g}] m^-3\n"
            f"n_GW={n_gw:.3g} m^-3    <n_e>/n_GW: {f_gw_c:.2f} at ne_c, "
            f"[{f_gw_lo:.2f} .. {f_gw_hi:.2f}] across the scan\n"
            f"D/T={p.deuterium_fraction:.3g}/{p.tritium_fraction:.3g}    "
            f"conf={self.confinement_var.get()}\n"
            f"tauE,e used = {tau_e_str}    tauE,i used = {tau_i_str}\n"
            f"e-i equipartition={equip_state}    alpha self-heating={alpha_state}    "
            f"first-orbit loss={orbit_state}    CX model={CX_MODELS_INV.get(p.cx_model, p.cx_model)}"
            f"{f' ({p.cx_loss_fraction:.3g})' if p.cx_model == 'manual_fraction' else ''}\n"
            f"profile-corrected 0-D={prof_state}\n"
            f"ECRH={p.p_ecrh_MW:.3g} MW (f_e={p.ecrh_f_e:.3g})    "
            f"ICRH={p.p_icrh_MW:.3g} MW (f_e={p.icrh_f_e:.3g}, f_i={p.icrh_f_i:.3g})",
        ]
        for index, beam in enumerate(self.model.beams, start=1):
            lines.append(f"NBI-{index}: {beam.species.upper()}    P={beam.power_MW:.3g} MW    "
                         f"E={beam.beam_energy_keV:.3g} keV    model={beam.shine_through_model}    "
                         f"{'co' if beam.co_current else 'counter'}-current")
        return "\n".join(lines)

    # ---------------------------------------------------------------- presets
    def _set_entry(self, key: str, value):
        if key in self.entries:
            self.entries[key].delete(0, "end")
            self.entries[key].insert(0, str(value))

    def _apply_preset(self, name: str):
        preset = self.PLASMA_PRESETS[name]
        non_entry_keys = ("profile_averaging", "confinement_mode")
        for field, value in preset.items():
            if field in non_entry_keys:
                continue  # not Plasma text entries -- handled below
            self._set_entry(f"plasma.{field}", value)
        # Full preset overwrite, same as every other field here: devices
        # without their own "profile_averaging" key turn the checkbox off.
        self.profile_var.set(bool(preset.get("profile_averaging", False)))
        # Confinement scaling appropriate to this device's aspect ratio, so
        # e.g. switching DANTE -> ITER doesn't leave an ST-fit Kaye scaling
        # active on a conventional-aspect machine.
        confinement_mode = preset.get("confinement_mode", "Fixed tauE (input)")
        if confinement_mode not in CONFINEMENT_MODES:
            confinement_mode = "Fixed tauE (input)"
        self.confinement_var.set(confinement_mode)
        # Orbit-loss and CX-loss default OFF for every preset (a deliberate
        # simplification for v1, matching Fixed/Kaye/IPB98 launching without
        # tuning a device-specific loss-fraction knob first).
        self.orbit_var.set(False)
        self._set_entry("models.cx_loss_fraction", 0.0)
        major_radius = preset["major_radius"]

        heating = self.MACHINE_HEATING_PRESETS.get(name, {})
        # Tangent radius defaults to on-axis (major_radius) unless a device
        # has a known real, generally off-axis tangency (e.g. TCV's 0.736 m).
        tangent_R_m = heating.get("tangent_R_m", major_radius)
        default_beams = [("D", 0.0, 100.0, True), ("D", 0.0, 100.0, True)]
        for beam_index, (species, power, energy_keV, co_current) in enumerate(
                heating.get("beams", default_beams)):
            self._set_entry(f"beam{beam_index}.species", species)
            self._set_entry(f"beam{beam_index}.power_MW", power)
            self._set_entry(f"beam{beam_index}.beam_energy_keV", energy_keV)
            getattr(self, f"beam_dir_var_{beam_index}").set("co" if co_current else "counter")
            # Reset tangent geometry too, so it doesn't carry over a
            # previously-selected, differently-scaled device's beam aiming.
            self._set_entry(f"beam{beam_index}.tangent_R_m", tangent_R_m)
            self._set_entry(f"beam{beam_index}.tangent_Z_m", 0.0)
        self._set_entry("aux.p_ecrh_MW", heating.get("ecrh_MW", 0.0))
        self._set_entry("aux.p_icrh_MW", heating.get("icrh_MW", 0.0))

        self.active_device = name
        self._highlight_active_preset()
        self._run()

    # --------------------------------------------------------------- export
    def _export_summary(self, kind: str):
        self._ensure_tab_built("Summary")
        extension = ".pdf" if kind == "pdf" else ".png"
        path = filedialog.asksaveasfilename(
            title="Save HI-Jass summary sheet", defaultextension=extension,
            filetypes=[("PDF file", "*.pdf")] if kind == "pdf"
            else [("PNG image", "*.png"), ("JPEG image", "*.jpg")],
        )
        if path:
            self.sum_fig.savefig(path, bbox_inches="tight")

    def _export_json(self):
        if self.last_result is None or not self.last_result.ok:
            self.status.configure(text="Nothing to export yet", text_color="#c0504d")
            return
        path = filedialog.asksaveasfilename(
            title="Save HI-Jass run record", defaultextension=".json",
            filetypes=[("JSON file", "*.json")],
        )
        if not path:
            return
        op = self.last_result.op
        record = {
            "generated": datetime.now().isoformat(timespec="seconds"),
            "git_sha": self._git_sha(),
            "device": self.active_device,
            "mode": self.mode,
            "plasma": dataclasses.asdict(self.model.plasma),
            "beams": [dataclasses.asdict(b) for b in self.model.beams],
            "operating_point": dataclasses.asdict(op),
            "summary": self.last_result.summary,
        }
        with open(path, "w", encoding="utf-8") as handle:
            json.dump(record, handle, indent=2,
                      default=lambda o: float(o) if hasattr(o, "__float__") else str(o))
        self.status.configure(text=f"Wrote {Path(path).name}", text_color="gray")

    # "plasma.<field>" is the entry key for most PlasmaParams fields, but a
    # few live under a different prefix (left over from the Models -> Losses
    # rail reorganisation, or the aux-heating section) -- redirect those.
    _PLASMA_ENTRY_PREFIX_OVERRIDE = {
        "cx_loss_fraction": "models", "cx_n0_over_ne": "models",
        "cx_n0_lcfs_over_ne": "models", "manual_v_phi_m_s": "models",
        "tau_phi_over_tauEi": "models", "f_alpha": "models",
        "p_ecrh_MW": "aux", "ecrh_f_e": "aux", "p_icrh_MW": "aux",
        "icrh_f_e": "aux", "icrh_f_i": "aux",
    }

    def _import_json(self):
        """Load a "Run record (JSON)" previously written by _export_json():
        restores every plasma/beam input field and model toggle, then
        re-solves. Tolerant of a record from an older/newer HI-Jass version
        (a field this version doesn't have is ignored; a field this version
        expects but the record lacks keeps its current value).
        """
        path = filedialog.askopenfilename(
            title="Load HI-Jass run record", filetypes=[("JSON file", "*.json"), ("All files", "*.*")],
        )
        if not path:
            return
        try:
            with open(path, "r", encoding="utf-8") as handle:
                record = json.load(handle)
        except Exception as exc:
            self.status.configure(text=f"Could not read {Path(path).name}: {exc}", text_color="#c0504d")
            return
        plasma_dict = record.get("plasma")
        beam_dicts = record.get("beams")
        if not isinstance(plasma_dict, dict) or not isinstance(beam_dicts, list) or len(beam_dicts) < 2:
            self.status.configure(text="Not a HI-Jass run record (missing plasma/beams)", text_color="#c0504d")
            return

        def set_entry(key, value):
            if key in self.entries:
                self.entries[key].delete(0, "end")
                self.entries[key].insert(0, str(value))

        for field, value in plasma_dict.items():
            if field == "plasma_current":
                set_entry("plasma.plasma_current_MA", value / 1.0e6)
                continue
            prefix = self._PLASMA_ENTRY_PREFIX_OVERRIDE.get(field, "plasma")
            set_entry(f"{prefix}.{field}", value)

        def set_by_label(var, mapping, wanted):
            for label, key in mapping.items():
                if key == wanted:
                    var.set(label)
                    return

        set_by_label(self.confinement_var, CONFINEMENT_MODES,
                     (plasma_dict.get("tau_Ee_mode", "fixed"), plasma_dict.get("tau_Ei_mode", "fixed")))
        self.orbit_var.set(bool(plasma_dict.get("enable_orbit_loss", False)))
        set_by_label(self.orbit_model_var, ORBIT_MODELS, plasma_dict.get("orbit_model", "large_aspect"))
        set_by_label(self.cx_model_var, CX_MODELS, plasma_dict.get("cx_model", "manual_fraction"))
        set_by_label(self.rotation_model_var, ROTATION_MODELS, plasma_dict.get("rotation_model", "off"))
        self.beam_beam_var.set(bool(plasma_dict.get("enable_beam_beam", False)))
        self.equip_var.set(bool(plasma_dict.get("enable_equipartition", False)))
        self.alpha_var.set(bool(plasma_dict.get("alpha_heating", False)))
        self.profile_var.set(bool(plasma_dict.get("profile_averaging", False)))

        for i, beam in enumerate(beam_dicts[:2]):
            prefix = f"beam{i}"
            for field in ("species", "power_MW", "beam_energy_keV", "tangent_R_m",
                          "tangent_Z_m", "manual_shine_through_fraction"):
                if field in beam and beam[field] is not None:
                    set_entry(f"{prefix}.{field}", beam[field])
            set_by_label(getattr(self, f"shine_var_{i}"), SHINE_LABEL_TO_MODEL,
                         beam.get("shine_through_model", "manual"))
            getattr(self, f"beam_dir_var_{i}").set("co" if beam.get("co_current", True) else "counter")

        self.active_device = record.get("device", "custom")
        self._highlight_active_preset()
        self._run()
        self.status.configure(text=f"Loaded {Path(path).name}", text_color="gray")

    @staticmethod
    def _git_sha() -> str:
        try:
            return subprocess.check_output(
                ["git", "rev-parse", "--short", "HEAD"],
                cwd=Path(__file__).parent, stderr=subprocess.DEVNULL,
            ).decode().strip()
        except Exception:
            return "unknown"

    # -------------------------------------------------------------- settings
    def _load_settings(self) -> dict:
        try:
            return json.loads(SETTINGS_PATH.read_text(encoding="utf-8"))
        except Exception:
            return {}

    def _save_settings(self):
        data = {key: entry.get() for key, entry in self.entries.items()}
        data["_mode"] = self.mode_toggle.get()
        data["_ui_scale"] = UI_SCALES[self._ui_scale_idx]
        data["_confinement"] = self.confinement_var.get()
        data["_orbit_loss"] = bool(self.orbit_var.get())
        data["_orbit_model"] = self.orbit_model_var.get()
        data["_equipartition"] = bool(self.equip_var.get())
        data["_alpha_heating"] = bool(self.alpha_var.get())
        data["_profile_averaging"] = bool(self.profile_var.get())
        data["_cx_model"] = self.cx_model_var.get()
        data["_rotation_model"] = self.rotation_model_var.get()
        data["_beam_beam"] = bool(self.beam_beam_var.get())
        data["_active_device"] = self.active_device
        for beam_index in (0, 1):
            data[f"beam{beam_index}._shine"] = getattr(self, f"shine_var_{beam_index}").get()
            data[f"beam{beam_index}._dir"] = getattr(self, f"beam_dir_var_{beam_index}").get()
        try:
            SETTINGS_PATH.parent.mkdir(parents=True, exist_ok=True)
            SETTINGS_PATH.write_text(json.dumps(data, indent=2), encoding="utf-8")
        except Exception:
            pass

    def _on_close(self):
        # CustomTkinter's teardown destroys every widget individually
        # (Python-level, one Tcl call per widget), which for this window's
        # widget count can take tens of seconds and reads as "frozen" or
        # "won't close". Settings are already saved by this point, so there
        # is nothing left to lose by skipping the graceful per-widget
        # destroy: hide the window immediately (instant feedback) and end
        # the process directly -- the OS reclaims every window and the
        # (daemon) solver thread right away.
        self._save_settings()
        self.withdraw()
        os._exit(0)


if __name__ == "__main__":
    try:
        HIJassApp().mainloop()
    except Exception as exc:  # noqa: BLE001
        message = str(exc).lower()
        if "display" in message or "tk" in message or "no $display" in message:
            print("HI-Jass requires a desktop session with a valid DISPLAY.")
        else:
            raise
