# HI-Jass User Manual

## Running the application

### Linux

```bash
./run_hi_jass.sh
```

### Windows

Double-click `run_hi_jass.bat`, or run:

```text
.\run_hi_jass.bat
```

The launcher creates the virtual environment, installs the dependencies, and starts the GUI.

## Default NBI configuration

The default configuration contains two neutral beams:

| Beam | Species | Injection energy | Injected power |
|---|---:|---:|---:|
| NBI-1 | D | 120 keV | 5 MW |
| NBI-2 | T | 180 keV | 5 MW |

## Two-beam slowing-down model

The slowing-down quantities are calculated independently for each beam. For beam $j$:

$$
E_{c,j} = \left(\frac{m_j}{m_e}\right)^{1/3} T_e
$$

$$
\tau_{se,j} = 6.27\times10^8\,
\frac{A_j T_e[\mathrm{eV}]^{3/2}}
{n_e[\mathrm{cm}^{-3}]\ln\Lambda}
$$

$$
\tau_{S,j} = \frac{\tau_{se,j}}{3}
\ln\left[1+\left(\frac{E_{b,j}}{E_{c,j}}\right)^{3/2}\right]
$$

The useful beam power includes the HotJass capture and loss cascade. The fast-ion density from each beam is:

$$
 n_{b0,j} =
\frac{P_{\mathrm{useful},j}\tau_{S,j}}
{E_{b,j}(10^3e)V}
$$

The total fast-ion density is:

$$
 n_{b0} = \sum_j n_{b0,j}
$$

## Mean fast-ion energy

Each beam has its own steady-state slowing-down distribution:

$$
 f_j(E) \propto
\frac{\sqrt{E}}
{E^{3/2}+E_{c,j}^{3/2}},
\qquad 0 < E \leq E_{b,j}
$$

Its mean fast-ion energy is:

$$
\langle E_{fast}\rangle_j =
\frac{\int_0^{E_{b,j}} E f_j(E)\,dE}
{\int_0^{E_{b,j}} f_j(E)\,dE}
$$

## Combined values shown in Results

The GUI reports power-weighted effective values for the two beams.

Effective slowing-down time:

$$
\tau_S^{\mathrm{eff}} =
\frac{\sum_j P_{\mathrm{useful},j}\tau_{S,j}}
{\sum_j P_{\mathrm{useful},j}}
$$

Effective mean fast-ion energy:

$$
\langle E_{fast}\rangle =
\frac{\sum_j P_{\mathrm{useful},j}
\langle E_{fast}\rangle_j}
{\sum_j P_{\mathrm{useful},j}}
$$

These effective values are used for the combined fast-ion diagnostics, including fast-ion energy density, pressure, beta, and:

$$
R = \frac{u_{fast}}{U_t}
$$

The underlying HotJass calculation still evaluates each beam separately before forming these combined diagnostics.

## Electron-ion exchange time

The electron-ion exchange time is based on the NRL Plasma Formulary thermal-equilibration rate. For each thermal ion species $s$ (D or T):

$$
\bar{\nu}_{ei,s} =
1.8\times10^{-19}
\frac{\sqrt{m_e m_s}\,Z_e^2 Z_s^2\,n_e\ln\Lambda}
{\left(m_s T_e + m_e T_i\right)^{3/2}}
$$

In this practical-unit expression, masses are in grams, temperatures are in eV, and electron density is in $\mathrm{cm}^{-3}$. The model uses $Z_e=Z_s=1$ for hydrogenic D and T ions.

For a D-T plasma, the effective exchange time is:

$$
τ_{IE} =
\frac{n_D+n_T}
{n_D\bar{\nu}_{ei,D}+n_T\bar{\nu}_{ei,T}}
$$

For a single thermal ion species, this reduces to:

$$
τ_{IE} = \frac{1}{\bar{\nu}_{ei}}
$$

The equivalent energy-over-power form used by the application is:

$$
τ_{IE} =
\frac{\frac{3}{2}(n_D+n_T)|T_i-T_e|(10^3e)V}
{|P_{ie}|}
$$

where $V$ is the plasma volume and $e=1.602176634\times10^{-19}\,\mathrm{J/eV}$. The exchange power is:

$$
P_{ie} = -P_{ei} =
-\frac{3}{2}V(T_e-T_i)
\left(n_D\bar{\nu}_{ei,D}+n_T\bar{\nu}_{ei,T}\right)(10^3e)
$$

The absolute values make $\tau_{IE}$ positive regardless of whether electrons heat ions ($T_e>T_i$) or ions heat electrons ($T_i>T_e$).

## Fusion power

The total fusion power is the sum of the thermal D-T and beam-target contributions:

$$
P_f = P_{f,\mathrm{thermal}} + P_{f,\mathrm{beam}}
$$

### Thermal D-T fusion

The thermal contribution is calculated from the thermal deuterium and tritium densities and the ion-temperature-dependent Bosch-Hale reactivity:

$$
P_{f,\mathrm{thermal}} =
n_{D0}n_{T0}\langle\sigma v\rangle_{DT}(T_i)V E_f
$$

Here $n_{D0}$ and $n_{T0}$ are the thermal species densities, $V$ is the plasma volume, and:

$$
E_f = 17.6\,\mathrm{MeV} = 17.6\times10^6 e
$$

The Bosch-Hale thermal reactivity is evaluated at $T_i$, because the relative velocity is that of the thermal ions, not the electrons.

### Beam-target fusion

Each useful beam is evaluated separately against the stationary thermal ions of the other D-T species:

$$
P_{f,\mathrm{beam}} =
\sum_j V E_f\int_0^{E_{b,j}}
n_{\mathrm{target},j}f_j(E)[\sigma v]_j(E)\,dE
$$

The beam slowing-down distribution is normalized to the fast-ion density of that beam:

$$
f_j(E) \propto
\frac{\sqrt{E}}{E^{3/2}+E_{c,j}^{3/2}},
\qquad
\int_0^{E_{b,j}}f_j(E)\,dE=n_{b0,j}
$$

The beam-target reactivity is:

$$
[\sigma v]_j(E)=\sigma_{DT}(E)v_j(E),
\qquad
v_j(E)=\sqrt{\frac{2E}{m_j}}
$$

For a D beam, $n_{\mathrm{target},j}=n_{T0}$. For a T beam, $n_{\mathrm{target},j}=n_{D0}$. The beam fast-ion density is obtained from the useful beam power:

$$
n_{b0,j} =
\frac{P_{\mathrm{useful},j}\tau_{S,j}}
{E_{b,j}(10^3e)V}
$$

The input shine-through fraction is used to determine captured beam power; the plotted $P_{\mathrm{shine}}$ is then the derived shine-through power from the scan. Charge-exchange and other configured losses similarly reduce the useful power before calculating $n_{b0,j}$.

### Current model limitations

- The beam-target integral uses stationary target ions, so the solved $T_i$ is not included in that reaction integral.
- D-D beam-target reactions are not included.
- Each included D-T reaction is assigned $17.6\,\mathrm{MeV}$ of fusion energy.

## Density scan

The Plasma tab defines the requested scan range with `n_e_min` and `n_e_max`. The Results tab reports the physically valid interval separately. Points below the HotJass charge-neutrality feasibility limit are not assigned physical output values and are omitted from plots.

## Important notation

- $T_e$: electron temperature
- $T_i$: thermal-ion temperature
- $P_e$: beam power deposited into electrons
- $P_i$: beam power deposited into ions
- $P_{ie}$: electron-ion energy exchange power
- $n_D$, $n_T$: thermal deuterium and tritium densities
- $n_{b0}$: total fast-ion density
- $P_{f,tot}$, $P_{f,th}$, $P_{f,b}$: total, thermal, and beam-target fusion power
- $\tau_{E,e}$, $\tau_{E,i}$: electron and ion energy confinement times
- $\tau_{IE}$: effective electron-ion exchange time
- $p_{th}$, $p_{fast}$: thermal and fast-ion pressure
- $\beta_t$: toroidal beta
