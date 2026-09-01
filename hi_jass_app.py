from __future__ import annotations

from typing import Dict

import customtkinter as ctk
import matplotlib
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

from hotjass_core import HotJassModel

matplotlib.use("TkAgg")


class HIJassApp(ctk.CTk):
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
    DISPLAY_GROUPS = {
        "Tₑ, Tᵢ": "Te, Ti",
        "Pₑ, Pᵢ, Pᵢₑ, Pshine-through": "P_e, P_i, Pi_e, P_shine-through",
        "nᴅ, nₜ, nᵦ": "n_D, n_T, n_b",
        "Pƒ,tot, Pƒ,th, Pƒ,b": "Pf_tot, Pf_th, Pf_b",
        "⟨Efast⟩": "<E_fast>",
        "τS, τE,e, τE,i, τIE": "tau_S, tauE_e, tauE_i, tau_IE",
        "R = ufast / Uₜₕ": "R = U_fast / U_th",
        "pₜₕ, pfast (isotropic)": "Pr_th, Pr_fast (isotropic)",
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

    def __init__(self):
        super().__init__()
        self.title("HI-Jass")
        self.geometry("1280x850")
        self.minsize(1050, 700)
        self.model = HotJassModel()
        self.scan = self.model.density_scan()
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(0, weight=1)
        self.tabview = ctk.CTkTabview(self)
        self.tabview.grid(row=0, column=0, padx=18, pady=18, sticky="nsew")
        for name in ("Plasma", "NBI", "Profiles/Shape", "Results"):
            self.tabview.add(name)
        self._build_plasma_tab()
        self._build_nbi_tab()
        self._build_profiles_tab()
        self._build_results_tab()
        self._run_model()

    def _entry_group(self, frame, fields, start_row=0, inactive=()):
        entries = {}
        for row, (label, attr, default) in enumerate(fields, start_row):
            label_widget = ctk.CTkLabel(frame, text=label, anchor="w")
            label_widget.grid(row=row, column=0, padx=12, pady=5, sticky="ew")
            entry = ctk.CTkEntry(frame)
            entry.insert(0, str(default))
            entry.grid(row=row, column=1, padx=12, pady=5, sticky="ew")
            if attr in inactive:
                entry.configure(state="disabled")
                label_widget.configure(text_color="gray")
            entries[attr] = entry
        return entries

    def _build_plasma_tab(self):
        frame = self.tabview.tab("Plasma")
        frame.grid_columnconfigure(0, weight=1)
        ctk.CTkLabel(frame, text="Low-aspect HotJass case", font=ctk.CTkFont(size=16, weight="bold")).grid(row=0, column=0, columnspan=2, padx=12, pady=(14, 4), sticky="w")
        fields = [
            ("R0 [m]", "major_radius", 0.65), ("a [m]", "minor_radius", 0.35),
            ("Elongation k", "elongation", 2.2), ("Triangularity delta", "triangularity", -0.35),
            ("Zeff", "effective_charge", 2.0), ("B0 [T]", "toroidal_field", 1.5),
            ("Ip [MA]", "plasma_current_MA", 1.5), ("Density peaking", "density_peaking", 1.0),
            ("Temperature peaking", "temp_peaking", 1.0),
            ("n_e_min [m^-3]", "n_e_min", 1.0e19), ("n_e_max [m^-3]", "n_e_max", 2.0e20),
            ("D fraction", "deuterium_fraction", 0.5), ("T fraction", "tritium_fraction", 0.5),
            ("tauE_e [s]", "tauE_e", 0.02), ("tauE_i [s]", "tauE_i", 0.05),
        ]
        self.plasma_entries = self._entry_group(frame, fields, 1, inactive=("density_peaking", "temp_peaking"))
        ctk.CTkLabel(frame, text="Te and Ti are calculated from the 0D balance and shown in Profiles/Shape and Results.", text_color="gray").grid(row=16, column=0, columnspan=2, padx=12, pady=(10, 4), sticky="w")
        self.alpha_var = ctk.BooleanVar(value=False)
        ctk.CTkCheckBox(frame, text="Include alpha heating (not in reference solver)", variable=self.alpha_var, state="disabled", text_color="gray").grid(row=17, column=0, columnspan=2, padx=12, pady=8, sticky="w")
        ctk.CTkButton(frame, text="Apply plasma settings", command=self._apply_plasma_settings).grid(row=18, column=0, columnspan=2, padx=12, pady=16, sticky="ew")

    def _build_nbi_tab(self):
        frame = self.tabview.tab("NBI")
        frame.grid_columnconfigure(0, weight=1)
        self.nbi_entries = {}
        for column, (title, species, power, energy) in enumerate((("NBI-1", "D", 5.0, 120.0), ("NBI-2", "T", 5.0, 180.0))):
            pane = ctk.CTkFrame(frame, fg_color="transparent")
            pane.grid(row=0, column=column, padx=10, pady=10, sticky="nsew")
            frame.grid_columnconfigure(column, weight=1)
            ctk.CTkLabel(pane, text=title, font=ctk.CTkFont(size=16, weight="bold")).grid(row=0, column=0, columnspan=2, pady=8, sticky="w")
            fields = [("Species", "species", species), ("P_NB [MW]", "power_MW", power), ("E_b [keV]", "beam_energy_keV", energy), ("Width", "beam_width", 0.35), ("Shift", "beam_shift", 0.12), ("Angle [deg]", "injection_angle_deg", 25.0), ("Shine-through", "shine_through_fraction", 0.08)]
            self.nbi_entries[column] = self._entry_group(pane, fields, 1, inactive=("beam_width", "beam_shift", "injection_angle_deg", "shine_through_fraction"))
        ctk.CTkButton(frame, text="Apply NBI settings", command=self._apply_nbi_settings).grid(row=1, column=0, columnspan=2, padx=12, pady=18, sticky="ew")

    def _build_profiles_tab(self):
        frame = self.tabview.tab("Profiles/Shape")
        frame.grid_columnconfigure(0, weight=1)
        frame.grid_rowconfigure(0, weight=1)
        self.profile_fig = Figure(figsize=(7, 6), dpi=100)
        self.profile_ax = self.profile_fig.add_subplot(221)
        self.profile_ax2 = self.profile_fig.add_subplot(222)
        self.shape_ax = self.profile_fig.add_subplot(223)
        self.formula_ax = self.profile_fig.add_subplot(224)
        self.profile_canvas = FigureCanvasTkAgg(self.profile_fig, master=frame)
        self.profile_canvas.get_tk_widget().grid(row=0, column=0, padx=10, pady=10, sticky="nsew")

    def _build_results_tab(self):
        frame = self.tabview.tab("Results")
        frame.grid_columnconfigure(1, weight=1)
        frame.grid_rowconfigure(0, weight=1)
        left = ctk.CTkFrame(frame, fg_color="transparent")
        left.grid(row=0, column=0, padx=10, pady=10, sticky="ns")
        ctk.CTkLabel(left, text="n_e scans", font=ctk.CTkFont(size=16, weight="bold")).pack(pady=8)
        self.observable_var = ctk.StringVar(value="Tₑ, Tᵢ")
        ctk.CTkOptionMenu(left, variable=self.observable_var, values=list(self.DISPLAY_GROUPS), command=lambda _: self._plot_results()).pack(pady=8)
        self.result_text = ctk.CTkTextbox(left, width=270, height=420, wrap="word")
        self.result_text.pack(pady=10, fill="y")
        self.results_fig = Figure(figsize=(8, 6), dpi=100)
        self.results_ax = self.results_fig.add_subplot(111)
        self.results_canvas = FigureCanvasTkAgg(self.results_fig, master=frame)
        self.results_canvas.get_tk_widget().grid(row=0, column=1, padx=10, pady=10, sticky="nsew")

    def _apply_plasma_settings(self):
        for attr, entry in self.plasma_entries.items():
            try:
                value = float(entry.get())
                if attr == "plasma_current_MA":
                    value *= 1.0e6
                    attr = "plasma_current"
                setattr(self.model.plasma, attr, value)
            except ValueError:
                pass
        self.model.plasma.alpha_heating = self.alpha_var.get()
        total = self.model.plasma.deuterium_fraction + self.model.plasma.tritium_fraction
        if total > 0:
            self.model.plasma.deuterium_fraction /= total
            self.model.plasma.tritium_fraction /= total
        self._run_model()

    def _apply_nbi_settings(self):
        for index, entries in self.nbi_entries.items():
            beam = self.model.beams[index]
            for attr, entry in entries.items():
                value = entry.get()
                try:
                    setattr(beam, attr, float(value))
                except ValueError:
                    setattr(beam, attr, value.upper())
        self._run_model()

    def _run_model(self):
        self.scan = self.model.density_scan()
        rho = self.model.rho_grid()
        density = self.model.density_profile(rho)
        middle = len(self.scan["Te"]) // 2
        te_profile = self.model.temperature_profile(rho, self.scan["Te"][middle])
        ti_profile = self.model.temperature_profile(rho, self.scan["Ti"][middle])
        p = self.model.plasma
        self.profile_ax.clear()
        self.profile_ax.plot(rho, density / 1e20, label=r"$n_e$")
        self.profile_ax.set(title=r"$n_e(\rho)$, $p_n = %.2f$" % p.density_peaking, xlabel=r"$\rho$ [1]", ylabel=r"$n_e$ [$10^{20}\,\mathrm{m}^{-3}$]")
        self.profile_ax.grid(alpha=0.3); self.profile_ax.legend()
        self.profile_ax2.clear()
        self.profile_ax2.plot(rho, te_profile, label=r"$T_e$"); self.profile_ax2.plot(rho, ti_profile, label=r"$T_i$")
        self.profile_ax2.set(title=r"$T_e(\rho)$, $T_i(\rho)$, $p_T = %.2f$" % p.temp_peaking, xlabel=r"$\rho$ [1]", ylabel=r"$T_e$, $T_i$ [keV]")
        self.profile_ax2.grid(alpha=0.3); self.profile_ax2.legend()
        theta = np.linspace(0, 2 * np.pi, 400)
        self.shape_ax.clear()
        delta = np.clip(p.triangularity, -0.999, 0.999)
        self.shape_ax.plot(p.major_radius + p.minor_radius * np.cos(theta + np.arcsin(delta) * np.sin(theta)), p.elongation * p.minor_radius * np.sin(theta))
        self.shape_ax.axvline(0.0, color="black", linestyle="--", linewidth=0.8, label=r"tokamak center $R=0$")
        self.shape_ax.axvline(p.major_radius, color="tab:red", linestyle=":", linewidth=0.9, label=r"magnetic axis $R_0$")
        self.shape_ax.set(title=r"Shape: $R_0=%.2f\,\mathrm{m}$, $a=%.2f\,\mathrm{m}$, $\kappa=%.2f$, $\delta=%.2f$" % (p.major_radius, p.minor_radius, p.elongation, p.triangularity), xlabel=r"$R$ [m]", ylabel=r"$Z$ [m]"); self.shape_ax.set_aspect("equal"); self.shape_ax.grid(alpha=0.3); self.shape_ax.legend(fontsize=7)
        self.formula_ax.clear(); self.formula_ax.axis("off")
        self.formula_ax.text(0, 0.85, "Profiles / 0D balance", fontsize=11, weight="bold")
        self.formula_ax.text(0, 0.62, r"$n_e(\rho) = n_{e0}(1-\rho^2)^{2p_n}$", fontsize=9)
        self.formula_ax.text(0, 0.42, r"$T_{e,i}(\rho) = T_{e0,i0}(1-\rho^2)^{2p_T}$", fontsize=9)
        self.formula_ax.text(0, 0.22, r"$T_e$, $T_i$ solved at each scanned $n_e$", fontsize=9)
        self.profile_fig.tight_layout(); self.profile_canvas.draw()
        self._plot_results()

    def _plot_results(self):
        selected = self.observable_var.get()
        keys = self.OBSERVABLES[self.DISPLAY_GROUPS[selected]]
        self.results_fig.clear()
        if len(keys) > 2:
            axes = self.results_fig.subplots(2, 2, squeeze=False).flat
            for axis, key in zip(axes, keys):
                finite = np.isfinite(self.scan[key])
                axis.plot(self.scan["n_e"][finite] / 1e20, self.scan[key][finite], marker="o", ms=3, label=self.LATEX_NAMES[key])
                axis.set_title(f"{self.LATEX_NAMES[key]} [${self.LATEX_UNITS[self.UNITS[key]]}$]")
                axis.set_xlabel(r"$n_e$ [$10^{20}\,\mathrm{m}^{-3}$]")
                axis.grid(alpha=0.3)
                axis.legend()
        else:
            self.results_ax = self.results_fig.add_subplot(111)
            for key in keys:
                finite = np.isfinite(self.scan[key])
                self.results_ax.plot(self.scan["n_e"][finite] / 1e20, self.scan[key][finite], marker="o", ms=3, label=f"{self.LATEX_NAMES[key]} [${self.LATEX_UNITS[self.UNITS[key]]}$]")
            self.results_ax.set_xlabel(r"$n_e$ [$10^{20}\,\mathrm{m}^{-3}$]")
            self.results_ax.set_ylabel(", ".join(f"{self.LATEX_NAMES[key]} [${self.LATEX_UNITS[self.UNITS[key]]}$]" for key in keys))
            self.results_ax.grid(alpha=0.3)
            self.results_ax.legend()
        units = ", ".join(f"{self.LATEX_NAMES[key]} [${self.LATEX_UNITS[self.UNITS[key]]}$]" for key in keys)
        valid_min = self.scan["n_e_valid_min"][0]
        valid_max = self.scan["n_e_valid_max"][0]
        self.results_fig.suptitle(r"HotJass scan: %s | valid $n_e=[%.2e, %.2e]$ m$^{-3}$" % (units, valid_min, valid_max))
        self.results_fig.tight_layout()
        self.results_canvas.draw()
        self.result_text.delete("1.0", "end")
        self.result_text.insert("end", "Requested n_e: %.3e to %.3e m^-3\n" % (self.model.plasma.n_e_min, self.model.plasma.n_e_max))
        self.result_text.insert("end", "Valid HotJass range: %.3e to %.3e m^-3\n\n" % (valid_min, valid_max))
        for key in keys:
            values = self.scan[key][np.isfinite(self.scan[key])]
            if values.size:
                self.result_text.insert("end", f"{key} [{self.UNITS[key]}]: {values.min():.4g} .. {values.max():.4g}\n")
            else:
                self.result_text.insert("end", f"{key} [{self.UNITS[key]}]: no feasible values\n")


if __name__ == "__main__":
    try:
        HIJassApp().mainloop()
    except Exception as exc:
        if "DISPLAY" in str(exc) or "tk" in str(exc).lower():
            print("HI-Jass requires a desktop session with a valid DISPLAY variable.")
        else:
            raise
