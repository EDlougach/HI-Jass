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
