from __future__ import annotations

from typing import Dict

import customtkinter as ctk
import matplotlib
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

from hotjass_core import HotJassModel, NBIParams, PlasmaParams

matplotlib.use("TkAgg")


class HIJassApp(ctk.CTk):
    def __init__(self):
        super().__init__()
        self.title("HI-Jass")
        self.geometry("1200x800")
        self.minsize(1000, 700)

        self.model = HotJassModel()

        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(0, weight=1)

        self.tabview = ctk.CTkTabview(self)
        self.tabview.grid(row=0, column=0, padx=20, pady=20, sticky="nsew")

        self.tabview.add("Plasma")
        self.tabview.add("NBI")
        self.tabview.add("Profiles")
        self.tabview.add("Geometry")
        self.tabview.add("Summary")

        self._build_plasma_tab()
        self._build_nbi_tab()
        self._build_profiles_tab()
        self._build_geometry_tab()
        self._build_summary_tab()

        self._run_model()

    def _build_plasma_tab(self):
        frame = self.tabview.tab("Plasma")
        frame.grid_columnconfigure(0, weight=1)

        fields = [
            ("Major radius [m]", "major_radius", 6.2),
            ("Minor radius [m]", "minor_radius", 1.8),
            ("Central density [m^-3]", "central_density", 2.5e20),
            ("Density peaking", "density_peaking", 1.2),
            ("Central temperature [keV]", "central_temperature", 8.0),
            ("Temperature peaking", "temp_peaking", 1.4),
            ("Toroidal field [T]", "toroidal_field", 2.5),
            ("Plasma current [A]", "plasma_current", 1.2e6),
        ]

        self.plasma_entries: Dict[str, ctk.CTkEntry] = {}
        for idx, (label, attr, default) in enumerate(fields):
            ctk.CTkLabel(frame, text=label, anchor="w").grid(row=idx, column=0, padx=12, pady=(10, 4), sticky="ew")
            entry = ctk.CTkEntry(frame)
            entry.insert(0, str(default))
            entry.grid(row=idx, column=1, padx=12, pady=(10, 4), sticky="ew")
            self.plasma_entries[attr] = entry

        button = ctk.CTkButton(frame, text="Apply plasma settings", command=self._apply_plasma_settings)
        button.grid(row=len(fields), column=0, columnspan=2, padx=12, pady=18, sticky="ew")

    def _build_nbi_tab(self):
        frame = self.tabview.tab("NBI")
        frame.grid_columnconfigure(0, weight=1)

        fields = [
            ("Injected power [MW]", "injected_power", 5.0),
            ("Beam energy [keV]", "beam_energy", 80.0),
            ("Beam species", "beam_species", "D"),
            ("Beam width", "beam_width", 0.45),
            ("Beam shift", "beam_shift", 0.15),
            ("Injection angle [deg]", "injection_angle_deg", 25.0),
            ("Shine-through fraction", "shine_through_fraction", 0.10),
            ("Slowing-down time [s]", "slowing_down_time", 0.15),
        ]

        self.nbi_entries: Dict[str, ctk.CTkEntry] = {}
        for idx, (label, attr, default) in enumerate(fields):
            ctk.CTkLabel(frame, text=label, anchor="w").grid(row=idx, column=0, padx=12, pady=(10, 4), sticky="ew")
            entry = ctk.CTkEntry(frame)
            entry.insert(0, str(default))
            entry.grid(row=idx, column=1, padx=12, pady=(10, 4), sticky="ew")
            self.nbi_entries[attr] = entry

        button = ctk.CTkButton(frame, text="Apply NBI settings", command=self._apply_nbi_settings)
        button.grid(row=len(fields), column=0, columnspan=2, padx=12, pady=18, sticky="ew")

    def _build_profiles_tab(self):
        frame = self.tabview.tab("Profiles")
        frame.grid_columnconfigure(0, weight=1)
        frame.grid_rowconfigure(0, weight=1)

        self.profile_fig = Figure(figsize=(5.8, 4.0), dpi=100)
        self.profile_ax = self.profile_fig.add_subplot(111)
        self.profile_canvas = FigureCanvasTkAgg(self.profile_fig, master=frame)
        self.profile_canvas.draw()
        self.profile_canvas.get_tk_widget().grid(row=0, column=0, padx=10, pady=10, sticky="nsew")

        self.toolbar = ctk.CTkButton(frame, text="Refresh plots", command=self._run_model)
        self.toolbar.grid(row=1, column=0, padx=10, pady=(0, 10), sticky="ew")

    def _build_geometry_tab(self):
        frame = self.tabview.tab("Geometry")
        frame.grid_columnconfigure(0, weight=1)
        frame.grid_rowconfigure(0, weight=1)

        self.geometry_fig = Figure(figsize=(5.8, 4.0), dpi=100)
        self.geometry_ax = self.geometry_fig.add_subplot(111)
        self.geometry_canvas = FigureCanvasTkAgg(self.geometry_fig, master=frame)
        self.geometry_canvas.draw()
        self.geometry_canvas.get_tk_widget().grid(row=0, column=0, padx=10, pady=10, sticky="nsew")

    def _build_summary_tab(self):
        frame = self.tabview.tab("Summary")
        frame.grid_columnconfigure(0, weight=1)

        self.summary_text = ctk.CTkTextbox(frame, height=22, width=60, wrap="word")
        self.summary_text.grid(row=0, column=0, padx=12, pady=12, sticky="nsew")

    def _apply_plasma_settings(self):
        plasma = self.model.plasma
        for attr, entry in self.plasma_entries.items():
            try:
                value = float(entry.get())
                setattr(plasma, attr, value)
            except ValueError:
                pass
        self._run_model()

    def _apply_nbi_settings(self):
        nbi = self.model.nbi
        for attr, entry in self.nbi_entries.items():
            value = entry.get()
            try:
                setattr(nbi, attr, float(value))
            except ValueError:
                if attr == "beam_species":
                    setattr(nbi, attr, value)
        self._run_model()

    def _run_model(self):
        data = self.model.compute_all()
        rho = data["rho"]
        density = data["density"]
        temperature = data["temperature"]
        fast_ions = data["fast_ions"]

        self.profile_ax.clear()
        self.profile_ax.plot(rho, density / np.max(density), label="Density")
        self.profile_ax.plot(rho, temperature / np.max(temperature), label="Temperature")
        self.profile_ax.plot(rho, fast_ions / np.max(fast_ions), label="Fast ions")
        self.profile_ax.set_title("Normalized plasma and NBI profiles")
        self.profile_ax.set_xlabel("Normalized radius $\rho$")
        self.profile_ax.set_ylabel("Normalized intensity")
        self.profile_ax.legend()
        self.profile_ax.grid(True, alpha=0.3)
        self.profile_canvas.draw()

        self.geometry_ax.clear()
        theta = np.linspace(0, 2 * np.pi, 400)
        r_minor = self.model.plasma.minor_radius
        r_major = self.model.plasma.major_radius
        x = r_major + r_minor * np.cos(theta)
        y = r_minor * np.sin(theta)
        self.geometry_ax.plot(x, y, color="tab:blue", lw=2)
        self.geometry_ax.set_aspect("equal")
        self.geometry_ax.set_title("Plasma geometry")
        self.geometry_ax.grid(True, alpha=0.3)
        self.geometry_canvas.draw()

        summary = self.model.summary()
        lines = [
            "HI-Jass summary",
            "================",
            f"Line-average density: {summary['line_average_density']:.3e} m^-3",
            f"Central temperature: {summary['central_temperature_keV']:.2f} keV",
            f"Beam power to plasma: {summary['beam_power_to_plasma']:.2f} MW",
            f"Injected beam energy: {summary['injected_energy_keV']:.1f} keV",
            f"Toroidal field: {summary['toroidal_field_T']:.2f} T",
            f"Plasma current: {summary['plasma_current_MA']:.2f} MA",
        ]
        self.summary_text.delete("1.0", "end")
        self.summary_text.insert("end", "\n".join(lines))


if __name__ == "__main__":
    app = HIJassApp()
    app.mainloop()
