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
import queue
import subprocess
import threading
import traceback
from datetime import datetime
from pathlib import Path
from tkinter import filedialog, messagebox

import customtkinter as ctk
import matplotlib
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from matplotlib.sankey import Sankey

from hotjass import physics
from hotjass_core import HotJassModel

matplotlib.use("TkAgg")

SETTINGS_PATH = Path.home() / ".hi_jass" / "settings.json"

CONFINEMENT_MODES = {
    "Fixed tauE (input)": ("fixed", "fixed"),
    "Kaye NSTX L-mode": ("kaye_nstx_lmode", "kaye_nstx_lmode"),
    "Kaye e / neoclassical i": ("kaye_nstx_lmode", "neoclassical"),
    "Fixed e / neoclassical i": ("fixed", "neoclassical"),
}
SHINE_LABEL_TO_MODEL = {"Riviere": "riviere", "Janev": "janev_suzuki", "Manual": "manual"}
SHINE_MODEL_TO_LABEL = {v: k for k, v in SHINE_LABEL_TO_MODEL.items()}


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
        self._expanded = not self._expanded
        self.header.configure(text=self._label())
        if self._expanded:
            self.body.grid(row=1, column=0, sticky="ew", padx=6, pady=6)
        else:
            self.body.grid_remove()


class HIJassApp(ctk.CTk):
    PLASMA_PRESETS = {
        "DANTE": {
            "major_radius": 0.65, "minor_radius": 0.35, "elongation": 2.2,
            "triangularity": -0.35, "effective_charge": 2.0, "toroidal_field": 1.5,
            "plasma_current_MA": 1.5, "central_density": 1.5e20,
            "n_e_min": 1.0e19, "n_e_max": 1.0e20, "density_peaking": 0.1, "temp_peaking": 1.0,
            "deuterium_fraction": 0.5, "tritium_fraction": 0.5, "tauE_e": 0.02, "tauE_i": 0.05,
        },
        "ITER": {
            "major_radius": 6.2, "minor_radius": 2.0, "elongation": 1.85,
            "triangularity": 0.33, "effective_charge": 1.7, "toroidal_field": 5.3,
            "plasma_current_MA": 15.0, "central_density": 1.0e20,
            "n_e_min": 5.0e19, "n_e_max": 1.5e20, "density_peaking": 0.0, "temp_peaking": 0.0,
            "deuterium_fraction": 0.5, "tritium_fraction": 0.5, "tauE_e": 3.7, "tauE_i": 3.7,
        },
        "JET": {
            "major_radius": 2.96, "minor_radius": 1.25, "elongation": 1.7,
            "triangularity": 0.32, "effective_charge": 1.5, "toroidal_field": 3.45,
            "plasma_current_MA": 4.0, "central_density": 6.0e19,
            "n_e_min": 2.0e19, "n_e_max": 1.0e20, "density_peaking": 0.0, "temp_peaking": 0.0,
            "deuterium_fraction": 0.5, "tritium_fraction": 0.5, "tauE_e": 1.5, "tauE_i": 1.5,
        },
        "ST40": {
            "major_radius": 0.45, "minor_radius": 0.30, "elongation": 1.8,
            "triangularity": 0.4, "effective_charge": 1.5, "toroidal_field": 3.0,
            "plasma_current_MA": 2.0, "central_density": 5.0e19,
            "n_e_min": 1.0e19, "n_e_max": 1.0e20, "density_peaking": 0.0, "temp_peaking": 0.0,
            "deuterium_fraction": 0.5, "tritium_fraction": 0.5, "tauE_e": 0.01, "tauE_i": 0.01,
        },
        "T-15MD": {
            "major_radius": 1.5, "minor_radius": 0.67, "elongation": 1.8,
            "triangularity": 0.3, "effective_charge": 1.5, "toroidal_field": 2.0,
            "plasma_current_MA": 2.0, "central_density": 5.0e19,
            "n_e_min": 1.0e19, "n_e_max": 1.0e20, "density_peaking": 0.0, "temp_peaking": 0.0,
            "deuterium_fraction": 0.5, "tritium_fraction": 0.5, "tauE_e": 0.1, "tauE_i": 0.1,
        },
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
        ("Temperature peaking", "temp_peaking", 1.0),
        ("Central n_e [m^-3]", "central_density", 1.5e20),
        ("Scan n_e min [m^-3]", "n_e_min", 1.0e19),
        ("Scan n_e max [m^-3]", "n_e_max", 1.0e20),
        ("D fraction", "deuterium_fraction", 0.5),
        ("T fraction", "tritium_fraction", 0.5),
        ("tauE,e [s]", "tauE_e", 0.02),
        ("tauE,i [s]", "tauE_i", 0.05),
    ]

    OBSERVABLES = {
        "Te, Ti": ["Te", "Ti"],
        "P_e, P_i, Pi_e, P_shine-through": ["P_e", "P_i", "Pi_e", "P_shine-through"],
        "n_D, n_T, n_b": ["n_D", "n_T", "n_b"],
        "Pf_tot, Pf_th, Pf_b": ["Pf_tot", "Pf_th", "Pf_b"],
        "<E_fast>": ["E_fast"],
        "tau_S, tauE_e, tauE_i, tau_IE": ["tau_S", "tauE_e", "tauE_i", "tau_IE"],
        "R = U_fast / U_th": ["R"],
        "Pr_th, Pr_fast (isotropic)": ["Pr_th", "Pr_fast"],
        "beta_T": ["beta_T"],
    }
    EQUIP_SENSITIVE = {"Tₑ, Tᵢ", "Pₑ, Pᵢ, Pᵢₑ, Pshine", "τS, τE,e, τE,i, τIE"}
    ALPHA_SENSITIVE = EQUIP_SENSITIVE | {"Pƒ,tot, Pƒ,th, Pƒ,b", "pₜₕ, pfast", "βt", "R = ufast / Uₜₕ"}
    DISPLAY_GROUPS = {
        "Tₑ, Tᵢ": "Te, Ti",
        "Pₑ, Pᵢ, Pᵢₑ, Pshine": "P_e, P_i, Pi_e, P_shine-through",
        "nᴅ, nₜ, nᵦ": "n_D, n_T, n_b",
        "Pƒ,tot, Pƒ,th, Pƒ,b": "Pf_tot, Pf_th, Pf_b",
        "⟨Efast⟩": "<E_fast>",
        "τS, τE,e, τE,i, τIE": "tau_S, tauE_e, tauE_i, tau_IE",
        "R = ufast / Uₜₕ": "R = U_fast / U_th",
        "pₜₕ, pfast": "Pr_th, Pr_fast (isotropic)",
        "βt": "beta_T",
    }
    UNITS = {
        "Te": "keV", "Ti": "keV", "P_e": "MW", "P_i": "MW", "Pi_e": "MW",
        "P_shine-through": "MW", "n_D": "m^-3", "n_T": "m^-3", "n_b": "m^-3",
        "Pf_tot": "MW", "Pf_th": "MW", "Pf_b": "MW", "E_fast": "keV",
        "tau_S": "s", "tauE_e": "s", "tauE_i": "s", "tau_IE": "s", "R": "1",
        "Pr_th": "Pa", "Pr_fast": "Pa", "beta_T": "%",
    }
    LATEX_UNITS = {
        "keV": r"\mathrm{keV}", "MW": r"\mathrm{MW}", "m^-3": r"\mathrm{m}^{-3}",
        "s": r"\mathrm{s}", "Pa": r"\mathrm{Pa}", "%": r"\%", "1": "1",
    }
    LATEX_NAMES = {
        "Te": r"$T_e$", "Ti": r"$T_i$", "P_e": r"$P_e$", "P_i": r"$P_i$",
        "Pi_e": r"$P_{ie}$", "P_shine-through": r"$P_{shine}$",
        "n_D": r"$n_D$", "n_T": r"$n_T$", "n_b": r"$n_{b0}$",
        "Pf_tot": r"$P_{f,tot}$", "Pf_th": r"$P_{f,th}$", "Pf_b": r"$P_{f,b}$",
        "E_fast": r"$\langle E_{fast}\rangle$", "tau_S": r"$\tau_S$",
        "tauE_e": r"$\tau_{E,e}$", "tauE_i": r"$\tau_{E,i}$", "tau_IE": r"$\tau_{IE}$",
        "R": r"$R = u_{fast}/U_t$", "Pr_th": r"$p_{th}$", "Pr_fast": r"$p_{fast}$",
        "beta_T": r"$\beta_t$",
    }

    DASH_ROWS = [
        ("Feasibility", "feasibility"),
        ("T_e [keV]", "Te"), ("T_i [keV]", "Ti"),
        ("n_e0 [m^-3]", "ne0"), ("n_b0 [m^-3]", "nb0"),
        ("P_NB injected [MW]", "P_NB"), ("P shine-through [MW]", "P_shine"),
        ("P captured [MW]", "P_capt"), ("P first-orbit loss [MW]", "P_orbit"),
        ("  NBI-1  rho_Li / dr / f_orbit", "orb1"),
        ("  NBI-2  rho_Li / dr / f_orbit", "orb2"),
        ("P charge-exchange loss [MW]", "P_cx"), ("P useful (to plasma) [MW]", "P_useful"),
        ("  -> electrons P_e [MW]", "P_e"), ("  -> ions P_i [MW]", "P_i"),
        ("P_ei equipartition (e->i) [MW]", "P_ei"),
        ("P_alpha self-heating [MW]", "P_alpha"),
        ("P_ECRH -> e / i [MW]", "p_ecrh"),
        ("P_ICRH -> e / i [MW]", "p_icrh"),
        ("P_heat total (NBI + aux + alpha) [MW]", "P_heat"),
        ("P_fusion total [MW]", "Pf_tot"), ("  thermal [MW]", "Pf_th"),
        ("  beam-target [MW]", "Pf_b"), ("Q = P_fus / P_NB", "Q"),
        ("<E_fast> [keV]", "E_fast"), ("beta_t [%]", "beta_t"),
        ("<n_e>/n_GW  (Greenwald, at ne_c)", "f_gw"),
        ("Dominant loss", "dominant"),
    ]

    def __init__(self):
        super().__init__()
        self.title("HI-Jass")
        self.geometry("1360x880")
        self.minsize(1100, 720)

        self.model = HotJassModel()
        self.active_device = "DANTE"
        self.mode = "Operating point"
        self.last_result: Result | None = None
        self._result_queue: queue.Queue[Result] = queue.Queue()
        self._running = False
        self._last_infeasible_key: str | None = None
        self._applied_snapshot: dict[str, str] = {}
        self._saved = self._load_settings()

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
        self.run_btn = ctk.CTkButton(run_bar, text="▶  Run", command=self._run)
        self.run_btn.grid(row=1, column=0, sticky="ew", pady=(6, 0))
        self._run_btn_fg = self.run_btn.cget("fg_color")
        self._run_btn_hover = self.run_btn.cget("hover_color")
        self.status = ctk.CTkLabel(run_bar, text="Ready", anchor="w", text_color="gray")
        self.status.grid(row=2, column=0, sticky="ew", pady=(4, 0))

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
        row += 1

        for beam_index in (0, 1):
            beam = self.model.beams[beam_index]
            section = CollapsibleSection(rail, f"NBI-{beam_index + 1}", expanded=beam_index == 0)
            section.grid(row=row, column=0, sticky="ew", pady=4)
            self._add_entries(section.body, f"beam{beam_index}", [
                ("Species", "species", beam.species),
                ("P_NB [MW]", "power_MW", beam.power_MW),
                ("E_b [keV]", "beam_energy_keV", beam.beam_energy_keV),
                ("Tangent R_t [m]", "tangent_R_m", self.model.plasma.major_radius),
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

        models = CollapsibleSection(rail, "Models", expanded=False)
        models.grid(row=row, column=0, sticky="ew", pady=4)
        ctk.CTkLabel(models.body, text="Confinement", anchor="w").grid(
            row=0, column=0, padx=8, pady=4, sticky="w")
        self.confinement_var = ctk.StringVar(
            value=self._saved.get("_confinement", "Fixed tauE (input)"))
        ctk.CTkOptionMenu(
            models.body, values=list(CONFINEMENT_MODES), variable=self.confinement_var,
        ).grid(row=0, column=1, padx=8, pady=4, sticky="ew")
        self.orbit_var = ctk.BooleanVar(value=self._saved.get("_orbit_loss", False))
        ctk.CTkCheckBox(
            models.body, text="First-orbit loss (direction set per NBI)", variable=self.orbit_var,
        ).grid(row=1, column=0, columnspan=2, padx=8, pady=4, sticky="w")
        self._add_entries(models.body, "models", [
            ("CX loss fraction (0-1)", "cx_loss_fraction",
             self._saved.get("models.cx_loss_fraction", 0.0)),
        ], start_row=2)
        self.equip_var = ctk.BooleanVar(value=self._saved.get("_equipartition", False))
        ctk.CTkCheckBox(
            models.body, text="e-i equipartition (couple Te, Ti)", variable=self.equip_var,
        ).grid(row=3, column=0, columnspan=2, padx=8, pady=4, sticky="w")
        self.alpha_var = ctk.BooleanVar(value=self._saved.get("_alpha_heating", False))
        ctk.CTkCheckBox(
            models.body, text="alpha self-heating (P_a into Te, Ti)", variable=self.alpha_var,
        ).grid(row=4, column=0, columnspan=2, padx=8, pady=4, sticky="w")
        self._add_entries(models.body, "models", [
            ("f_alpha (confined fraction 0-1)", "f_alpha", self._saved.get("models.f_alpha", 1.0)),
        ], start_row=5)

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

        for name in ("Dashboard", "Power flow", "Fast ions", "Profiles", "Assumptions"):
            self.tv_op.add(name)
        for name in ("Scan", "Summary"):
            self.tv_scan.add(name)

        self._build_dashboard(self.tv_op.tab("Dashboard"))
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
        self.fi_fig, self.fi_canvas, holder = self._plot_area(self.tv_op.tab("Fast ions"))
        holder.pack(fill="both", expand=True)
        self.prof_fig, self.prof_canvas, holder = self._plot_area(self.tv_op.tab("Profiles"))
        holder.pack(fill="both", expand=True)
        self._build_assumptions(self.tv_op.tab("Assumptions"))
        self._build_scan(self.tv_scan.tab("Scan"))
        self.sum_fig, self.sum_canvas, holder = self._plot_area(self.tv_scan.tab("Summary"), figsize=(11, 8))
        holder.pack(fill="both", expand=True)

        export = ctk.CTkFrame(container, fg_color="transparent")
        export.grid(row=1, column=0, sticky="ew", pady=(6, 0))
        ctk.CTkLabel(export, text="Export:").pack(side="left", padx=(4, 8))
        ctk.CTkButton(export, text="Summary PNG", width=110,
                      command=lambda: self._export_summary("png")).pack(side="left", padx=4)
        ctk.CTkButton(export, text="Summary PDF", width=110,
                      command=lambda: self._export_summary("pdf")).pack(side="left", padx=4)
        ctk.CTkButton(export, text="Run record (JSON)", width=150,
                      command=self._export_json).pack(side="left", padx=4)

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
        self.mode = mode
        if mode == "Scan":
            self.tv_op.grid_remove()
            self.tv_scan.grid()
        else:
            self.tv_scan.grid_remove()
            self.tv_op.grid()

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

        ee_mode, ei_mode = CONFINEMENT_MODES[self.confinement_var.get()]
        plasma.tau_Ee_mode, plasma.tau_Ei_mode = ee_mode, ei_mode
        plasma.enable_orbit_loss = bool(self.orbit_var.get())
        plasma.enable_equipartition = bool(self.equip_var.get())
        if plasma.enable_equipartition and plasma.tau_Ei_mode == "neoclassical":
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
        self.run_btn.configure(state="disabled", text="⏳  Running…",
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
        self.run_btn.configure(state="normal", text="▶  Run",
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
        self._render_powerflow(result)
        self._render_fastions(result)
        self._render_profiles(result)
        self.assump_box.delete("1.0", "end")
        self.assump_box.insert("end", self._assess(self.model, result.op))
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
        self._render_summary(result)

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
        orb = {}
        for i, beam in enumerate(self.model.beams):
            sp = beam.species.upper()
            rho_li = physics.larmor_radius_m(beam.beam_energy_keV, plasma.toroidal_field, sp)
            dr = physics.passing_orbit_width(
                beam.beam_energy_keV, plasma.toroidal_field, plasma.plasma_current / 1e6,
                plasma.major_radius, a, plasma.elongation, sp)
            f_orb = op.f_orbit_loss[i] if op.f_orbit_loss else 0.0
            orb[f"orb{i + 1}"] = (
                f"rho_Li={rho_li * 100:.1f} cm ({rho_li / a:.2f} a),  "
                f"dr={dr * 100:.1f} cm ({dr / a:.2f} a),  f_orbit={f_orb:.3f}")
        n_line, n_gw, f_gw = self._greenwald()

        values = {
            "Te": self._fmt(op.Te_keV), "Ti": self._fmt(op.Ti_keV),
            "ne0": self._fmt(op.ne0_m3, "{:.3e}"), "nb0": self._fmt(op.nb0_m3, "{:.3e}"),
            "P_NB": self._fmt(P_NB), "P_shine": self._fmt(op.P_shine_w * mw),
            "P_capt": self._fmt(op.P_capt_w * mw), "P_orbit": self._fmt(op.P_orbit_loss_w * mw),
            "P_cx": self._fmt(op.P_cx_loss_w * mw), "P_useful": self._fmt(op.P_useful_w * mw),
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
            "Pf_tot": self._fmt(op.pf_total_w * mw), "Pf_th": self._fmt(op.pf_thermal_w * mw),
            "Pf_b": self._fmt(op.pf_beam_w * mw),
            "Q": self._fmt(op.pf_total_w / op.P_NB_total_w if op.P_NB_total_w else None, "{:.3g}"),
            "E_fast": self._fmt(op.avg_fast_energy_keV),
            "beta_t": self._fmt(op.beta_t * 100.0),
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

    def _render_fastions(self, result: Result):
        op = result.op
        model = self.model
        volume = result.volume_m3
        fig = self.fi_fig
        fig.clear()
        ax = fig.add_subplot(111)
        te = op.Te_keV or 1e-3
        any_curve = False
        for i, beam in enumerate(model.beams):
            species = beam.species.upper()
            eb = beam.beam_energy_keV
            f_capt = op.f_capture[i] if op.f_capture else 1.0
            f_orbit = op.f_orbit_loss[i] if op.f_orbit_loss else 0.0
            p_use = beam.power_MW * 1e6 * f_capt * (1.0 - f_orbit) * (1.0 - model.plasma.cx_loss_fraction)
            tau_s = physics.thermalization_time(op.ne0_m3, te, eb, species)
            nb0_i = p_use * tau_s / (eb * 1e3 * physics.E_CHARGE * max(volume, 1e-9))
            grid = np.linspace(1e-3, eb, 400)
            fE = physics.slowing_down_distribution(te, nb0_i, eb, grid, species)
            mean_e = physics.average_fast_energy_keV(te, eb, species)
            line, = ax.plot(grid, fE, label=f"NBI-{i + 1} {species} {eb:.0f} keV   "
                            fr"$\langle E\rangle$={mean_e:.0f} keV,  $\tau_s$={tau_s:.3g} s")
            ax.axvline(eb, color=line.get_color(), linestyle=(0, (4, 3)), linewidth=1.0, alpha=0.8)
            ax.annotate(f"$E_b$ = {eb:.0f}", xy=(eb, 0), xytext=(0, 2),
                        textcoords="offset points", ha="right", va="bottom",
                        fontsize=7, color=line.get_color(), rotation=90)
            any_curve = True
        ax.set_ylim(bottom=0.0)
        ax.set_xlabel("E [keV]")
        ax.set_ylabel(r"$f(E)$ [$\mathrm{m}^{-3}\,\mathrm{keV}^{-1}$]")
        ax.set_title("Per-beam steady-state slowing-down distribution")
        ax.grid(alpha=0.3)
        if any_curve:
            ax.legend(fontsize=8)
        fig.tight_layout()
        self.fi_canvas.draw_idle()

    def _render_profiles(self, result: Result):
        model = self.model
        op = result.op
        plasma = model.plasma
        rho = model.rho_grid()
        density = model.density_profile(rho)
        te_profile = model.temperature_profile(rho, op.Te_keV or 0.0)
        ti_profile = model.temperature_profile(rho, op.Ti_keV or 0.0)

        fig = self.prof_fig
        fig.clear()
        ax_n = fig.add_subplot(221)
        ax_n.plot(rho, density / 1e20, label=r"$n_e$")
        ax_n.set(title=r"$n_e(\rho)$, $p_n=%.2f$" % plasma.density_peaking,
                 xlabel=r"$\rho$", ylabel=r"$n_e$ [$10^{20}\,\mathrm{m}^{-3}$]")
        ax_n.grid(alpha=0.3)
        ax_n.legend()

        ax_t = fig.add_subplot(222)
        ax_t.plot(rho, te_profile, label=r"$T_e$")
        ax_t.plot(rho, ti_profile, label=r"$T_i$")
        ax_t.set(title=r"$T_e(\rho),\ T_i(\rho)$, $p_T=%.2f$" % plasma.temp_peaking,
                 xlabel=r"$\rho$", ylabel=r"$T$ [keV]")
        ax_t.grid(alpha=0.3)
        ax_t.legend()

        ax_s = fig.add_subplot(223)
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

        ax_f = fig.add_subplot(224)
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
            axes = fig.subplots(2, 2, squeeze=False).flat
            for ax, key in zip(axes, keys):
                finite = np.isfinite(self.scan[key])
                ax.plot(density_axis[finite], self.scan[key][finite], marker="o", ms=3,
                        label=self.LATEX_NAMES[key])
                ax.set_title(f"{self.LATEX_NAMES[key]} [${self.LATEX_UNITS[self.UNITS[key]]}$]")
                ax.set_xlabel(r"$n_e$ [$10^{20}\,\mathrm{m}^{-3}$]")
                ax.grid(alpha=0.3)
                ax.legend()
        else:
            ax = fig.add_subplot(111)
            for key in keys:
                finite = np.isfinite(self.scan[key])
                ax.plot(density_axis[finite], self.scan[key][finite], marker="o", ms=3,
                        label=f"{self.LATEX_NAMES[key]} [${self.LATEX_UNITS[self.UNITS[key]]}$]")
            ax.set_xlabel(r"$n_e$ [$10^{20}\,\mathrm{m}^{-3}$]")
            ax.grid(alpha=0.3)
            ax.legend()
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
            (r"$P_{ie}$ [MW]", lambda ax: ax.plot(density_axis, scan["Pi_e"])),
            (r"$P_{shine}$ [MW]", lambda ax: ax.plot(density_axis, scan["P_shine-through"])),
            (r"Species n [$\mathrm{m}^{-3}$]", lambda ax: (ax.plot(density_axis, scan["n_D"], label=r"$n_D$"),
                                                           ax.plot(density_axis, scan["n_T"], label=r"$n_T$"),
                                                           ax.plot(density_axis, scan["n_b"], label=r"$n_b$"))),
            ("Fusion power [MW]", lambda ax: (ax.plot(density_axis, scan["Pf_tot"], label=r"$P_{f,tot}$"),
                                              ax.plot(density_axis, scan["Pf_th"], label=r"$P_{f,th}$"),
                                              ax.plot(density_axis, scan["Pf_b"], label=r"$P_{f,b}$"))),
            ("Times [s]", lambda ax: (ax.plot(density_axis, scan["tau_S"], label=r"$\tau_S$"),
                                      ax.plot(density_axis, scan["tauE_e"], label=r"$\tau_{E,e}$"),
                                      ax.plot(density_axis, scan["tauE_i"], label=r"$\tau_{E,i}$"),
                                      ax.plot(density_axis, scan["tau_IE"], label=r"$\tau_{IE}$"))),
            ("<E_fast> [keV]", lambda ax: ax.plot(density_axis, scan["E_fast"])),
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
                         "Fusion power [MW]", "Times [s]", "Pressure [Pa]"):
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
        lines.append(f"First-orbit loss: {'ON' if plasma.enable_orbit_loss else 'off'}  ({dirs})")
        a = plasma.minor_radius
        wide = False
        for i, b in enumerate(model.beams):
            sp = b.species.upper()
            rho_li = physics.larmor_radius_m(b.beam_energy_keV, plasma.toroidal_field, sp)
            dr = physics.passing_orbit_width(
                b.beam_energy_keV, plasma.toroidal_field, plasma.plasma_current / 1e6,
                plasma.major_radius, a, plasma.elongation, sp)
            f_orb = op.f_orbit_loss[i] if op.f_orbit_loss else 0.0
            lines.append(f"  NBI-{i + 1}: rho_Li = {rho_li * 100:.1f} cm ({rho_li / a:.2f} a), "
                         f"passing orbit dr = {dr * 100:.1f} cm ({dr / a:.2f} a), f_orbit = {f_orb:.3f}")
            wide = wide or rho_li / a > 0.3 or dr / a > 0.3
        if wide:
            lines.append("  ! rho_Li/a or dr/a > 0.3: orbit width ~ machine size; the mean-shift")
            lines.append("    first-orbit-loss estimate is order-unity uncertain here (needs an")
            lines.append("    orbit follower). Co-current f_orbit=0 is a model artefact, not physics.")
        lines.append(f"CX loss fraction: {plasma.cx_loss_fraction:.3g} (flat efficiency knob)")
        lines.append(f"e-i equipartition: {'ON (coupled Te/Ti solve)' if plasma.enable_equipartition else 'off (decoupled Te, Ti)'}")
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
        if plasma.tau_Ee_mode.startswith("kaye"):
            lines.append("  Kaye NSTX L-mode: an ST-appropriate fit (A ~ 1.3-1.5 dataset).")
        if plasma.tau_Ei_mode == "neoclassical":
            lines.append("  Neoclassical ion transport has no anomalous channel -> Ti can be")
            lines.append("  large / implausible; treat as a lower bound on ion transport.")
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
        or the min..max range over the current scan for a physics-based mode
        (the manual input is ignored in that case).
        """
        if mode == "fixed":
            return f"{manual_val:.3g} s (fixed input)"
        sc = getattr(self, "scan", None)
        if sc is not None and scan_key in sc:
            v = sc[scan_key][np.isfinite(sc[scan_key])]
            if v.size:
                return f"[{v.min():.3g} .. {v.max():.3g}] s ({mode}, manual input ignored)"
        return f"— ({mode})"

    def _summary_parameters(self) -> str:
        p = self.model.plasma
        equip_state = "ON" if p.enable_equipartition else "off"
        alpha_state = f"ON (f_alpha={p.f_alpha:.3g})" if p.alpha_heating else "off"
        orbit_state = "ON" if p.enable_orbit_loss else "off"
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
            f"first-orbit loss={orbit_state}    CX loss frac={p.cx_loss_fraction:.3g}\n"
            f"ECRH={p.p_ecrh_MW:.3g} MW (f_e={p.ecrh_f_e:.3g})    "
            f"ICRH={p.p_icrh_MW:.3g} MW (f_e={p.icrh_f_e:.3g}, f_i={p.icrh_f_i:.3g})",
        ]
        for index, beam in enumerate(self.model.beams, start=1):
            lines.append(f"NBI-{index}: {beam.species.upper()}    P={beam.power_MW:.3g} MW    "
                         f"E={beam.beam_energy_keV:.3g} keV    model={beam.shine_through_model}    "
                         f"{'co' if beam.co_current else 'counter'}-current")
        return "\n".join(lines)

    # ---------------------------------------------------------------- presets
    def _apply_preset(self, name: str):
        for field, value in self.PLASMA_PRESETS[name].items():
            key = f"plasma.{field}"
            if key in self.entries:
                self.entries[key].delete(0, "end")
                self.entries[key].insert(0, str(value))
        self.active_device = name
        self._highlight_active_preset()
        self._run()

    # --------------------------------------------------------------- export
    def _export_summary(self, kind: str):
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
        data["_confinement"] = self.confinement_var.get()
        data["_orbit_loss"] = bool(self.orbit_var.get())
        data["_equipartition"] = bool(self.equip_var.get())
        data["_alpha_heating"] = bool(self.alpha_var.get())
        for beam_index in (0, 1):
            data[f"beam{beam_index}._shine"] = getattr(self, f"shine_var_{beam_index}").get()
            data[f"beam{beam_index}._dir"] = getattr(self, f"beam_dir_var_{beam_index}").get()
        try:
            SETTINGS_PATH.parent.mkdir(parents=True, exist_ok=True)
            SETTINGS_PATH.write_text(json.dumps(data, indent=2), encoding="utf-8")
        except Exception:
            pass

    def _on_close(self):
        self._save_settings()
        self.destroy()


if __name__ == "__main__":
    try:
        HIJassApp().mainloop()
    except Exception as exc:  # noqa: BLE001
        message = str(exc).lower()
        if "display" in message or "tk" in message or "no $display" in message:
            print("HI-Jass requires a desktop session with a valid DISPLAY.")
        else:
            raise
