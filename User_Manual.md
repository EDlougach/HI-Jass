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

## Density peaking and shine-through

The density peaking input defines the analytic electron-density profile:

$$
n_e(\rho)=n_{e0}(1-\rho^2)^{2p_n},
\qquad 0\leq\rho\leq1
$$

Shine-through uses the density along the beam chord rather than assuming that
the central density applies everywhere. Each NBI supplies a tangent point
$(R_t,Z_t)$ in the NBI tab. The straight tangential chord at that point is
intersected with the shaped plasma cross-section, and the optical depth is
evaluated as:

$$
	au_b = \sigma_{\mathrm{stop}}(E_b/A_b)
\int_{\mathrm{chord}} n_e(\rho(s))\,ds
$$

and the captured and shine-through fractions are:

$$
f_{\mathrm{capt},b}=1-e^{-\tau_b},
\qquad
f_{\mathrm{shine},b}=e^{-\tau_b}.
$$

The chord is evaluated numerically using the local cylindrical radius along
the tangent path and the plasma elongation. With $p_n=0$, the density is
uniform along the computed chord; with $p_n>0$, the profile is integrated
point by point. A tangent point outside the plasma produces zero captured
power.

### Shine-through model selection

Each NBI has a three-way shine-through selector:

- **Riviere**: the existing Riviere-type analytic stopping cross-section fit.
- **Janev/Suzuki**: the Janev-Boley-Post analytic fit from the supplied paper,
	Nuclear Fusion 29 (1989) 2125, Eqs. (23)-(26), using the published Table III
	coefficients and the carbon-impurity coefficients from Table IV. It models
	total effective stopping, including the paper's density, electron-temperature
	and $Z_{\mathrm{eff}}$ dependence. The fit is published for
	$100\leq E\leq10^4$ keV/u, $10^{12}\leq n_e\leq10^{15}$ cm$^{-3}$ and
	$1\leq T_e\leq50$ keV; inputs outside these ranges are clipped to the
	nearest validity boundary.
- **Manual**: the entered manual shine-through fraction is used directly;
	captured power is $(1-f_{\mathrm{shine}})P_{NB}$.

The selected calculated stopping model is also used by the optional
first-orbit-loss attenuation weighting. The Manual mode bypasses the
cross-section calculation for shine-through.

## First-orbit loss (reduced model)

First-orbit loss is an **opt-in** channel (the *First-orbit loss* checkbox in
the Models section). When enabled, it estimates the fraction of *captured* fast
ions that leave the plasma promptly — on their first drift orbit, before any
slowing-down or charge-exchange — because they are born close enough to the last
closed flux surface (LCFS) that their guiding-centre orbit, or their gyro-orbit,
crosses it.

### Birth profile along the chord

The calculation reuses the shine-through deposition. Along the straight
tangential chord (assumed $R_t=R_0$, length $2a$) the local ionisation rate
density is

$$
\dot n(x) \propto
n_{e0}\,\sigma_{\mathrm{stop}}(E_{b,j}/A_j)\,
e^{-n_{e0}\sigma_{\mathrm{stop}}x},
\qquad x\in[0,2a],
$$

and the path coordinate is mapped to a normalised flux label by

$$
\rho(x) = \frac{|x-a|}{a},
$$

so $\rho=0$ at the chord mid-point and $\rho=1$ at both ends. The stopping
cross-section is the beam's selected shine-through model; **Manual** mode has no
deposition shape, so the Riviere fit is substituted for this weighting only.

### Passing-orbit width

The radial width of a strongly co-passing drift orbit at the full injection
energy is the textbook large-aspect-ratio estimate

$$
\Delta r_j = q_\ast\,\rho_{Li,j},
\qquad
q_\ast = \frac{2\pi a^2 B_t}{\mu_0 R_0 I_p}\,\frac{1+\kappa^2}{2},
\qquad
\rho_{Li,j} = \frac{m_j v_{b,j}}{e B_t}.
$$

$q_\ast$ is the cylindrical, elongation-corrected safety factor and
$\rho_{Li,j}$ the toroidal-field gyroradius at $E_{b,j}$.

### Loss criterion

The orbit is modelled as a **rigid mean shift** of the birth radius: inward for
co-current injection (favourable), outward for counter-current:

$$
\rho_{\mathrm{eff}}(x) = \rho(x) \mp \frac{\Delta r_j}{a}
\qquad(\text{co}\,/\,\text{counter}),
$$

and an ion born at $x$ is counted as lost if $\rho_{\mathrm{eff}}(x) > 1$.
Because $\rho(x)\le 1$, the co-current shift can never trigger this test, so the
mean-shift model gives **exactly zero first-orbit loss for co-current
injection**; for counter-current it loses every ion born with
$\rho(x) > 1 - \Delta r_j/a$.

An optional, direction-independent **Larmor (gyro-orbit) prompt-loss** channel
can be OR-ed in: an ion whose guiding centre is born within one gyro-diameter of
the LCFS,

$$
\rho(x) > 1 - \frac{2\rho_{Li,j}}{a},
$$

is lost regardless of the drift direction (gyration is much faster than the
guiding-centre drift). This channel is disabled by default.

### Loss fraction and power

$$
f_{\mathrm{orbit},j} =
\frac{\int_{\text{lost }x}\dot n(x)\,dx}
     {\int_0^{2a}\dot n(x)\,dx},
\qquad
P_{\mathrm{orbit},j} = f_{\mathrm{orbit},j}\,P_{\mathrm{capt},j}.
$$

First-orbit loss is applied to the captured power *before* the charge-exchange
loss and before the useful power that drives $T_e$, $T_i$ and $n_{b0,j}$.
Injection direction (co-/counter-current) is set per NBI in its input section.

### Reported quantities

The Dashboard and the Assumptions tab report, per beam, $\rho_{Li,j}$ and
$\Delta r_j$ (in cm and as fractions of $a$) and $f_{\mathrm{orbit},j}$.

### Limitations

- **The co-current zero is a model artefact, not physics.** The rigid
	mean-shift ignores the finite radial width of the orbit and the outboard
	excursion that a favourable inward mean shift does not cancel; published
	orbit-following comparisons give of order a few percent co-current
	first-orbit loss, versus tens of percent counter-current.
- The chord label $\rho(x)=|x-a|/a$ folds the inboard and outboard halves
	together, so the model cannot distinguish a high-field-side birth from a
	low-field-side birth. It therefore cannot represent the correct structure —
	co-current losing high-field-side edge births, counter-current losing
	low-field-side edge births, and trapped-ion banana-tip losses on the
	low-field side in either direction.
- Circular cross-section, large-aspect-ratio, $R_t=R_0$ tangential geometry;
	$\Delta r_j$ and $\rho_{Li,j}$ are evaluated at the full injection energy (an
	upper bound — real ions slow down).
- The birth-profile integral in this channel uses the flat central density
	$n_{e0}$, not the peaked profile used elsewhere.
- When $\rho_{Li}/a$ or $\Delta r/a \gtrsim 0.3$ — orbit width comparable to the
	machine size, typical of compact low-aspect-ratio devices — the whole
	estimate is order-unity uncertain and really needs a guiding-centre orbit
	follower. The Assumptions tab flags this case.

## Power balance (electron and ion)

At each central density the model solves a 0-D steady-state power balance for
$T_e$ and $T_i$. The heating that drives each temperature is the sum of the
neutral-beam split, prescribed auxiliary heating (ECRH and ICRH), and — when
enabled — fusion alpha self-heating.

### Beam heating split

Each useful beam $j$ delivers a fraction $L_e(x_j)$ of its power to electrons
and $L_i(x_j)=1-L_e(x_j)$ to ions, with $x_j=E_{b,j}/E_{c,j}(T_e)$ and

$$
L_i(x) = \frac{1}{x}\int_0^{x}\frac{du}{1+u^{3/2}},
\qquad
L_e(x) = 1 - L_i(x).
$$

$L_i\to1$ as $x\to0$ (cold, ion-dominated) and $L_i\to0$ as $x\to\infty$ (fast
beam, electron-dominated). $P_{\mathrm{useful},j}$ is the injected power after
shine-through, first-orbit loss and charge-exchange loss (preceding sections).

### Auxiliary heating (ECRH, ICRH)

ECRH and ICRH are entered as powers with electron/ion fractions. Their
contribution to each channel is

$$
P_{\mathrm{aux},e} = f_e^{\mathrm{EC}}P_{\mathrm{ECRH}}
                   + f_e^{\mathrm{IC}}P_{\mathrm{ICRH}},
\qquad
P_{\mathrm{aux},i} = \left(1-f_e^{\mathrm{EC}}\right)P_{\mathrm{ECRH}}
                   + f_i^{\mathrm{IC}}P_{\mathrm{ICRH}}.
$$

ECRH is electron-dominant, $f_e^{\mathrm{EC}}$ close to 1. For ICRH the electron
and ion fractions are set independently; if
$f_e^{\mathrm{IC}}+f_i^{\mathrm{IC}}<1$ the remainder is treated as power that
does not thermalise (a fast-ion tail or direct loss). These sources are
Te/Ti-independent constants added to the right-hand side of each balance.

### Alpha self-heating (optional)

When enabled, $3.5\,\mathrm{MeV}$ fusion alphas deposit

$$
P_\alpha = f_\alpha\,\frac{E_\alpha}{E_f}\,P_{f,\mathrm{tot}}
         = f_\alpha\,\frac{3.5}{17.6}\,P_{f,\mathrm{tot}},
$$

split between electrons and ions by the alpha slowing-down fraction

$$
L_e^\alpha = L_e\!\left(\frac{E_\alpha}{E_c^\alpha(T_e)}\right),
\qquad
E_c^\alpha = \left(\frac{m_\alpha}{m_e}\right)^{1/3}T_e,
$$

so $P_{\alpha,e}=L_e^\alpha P_\alpha$ and $P_{\alpha,i}=(1-L_e^\alpha)P_\alpha$.
Because $P_\alpha$ depends on $P_{f,\mathrm{tot}}$, which depends on $T_i$, which
depends on $P_{\alpha,i}$, the alpha term is closed by an under-relaxed
fixed-point iteration; with alpha heating off the solve is a single pass.

### Decoupled balance (default)

With electron-ion equipartition off, the two balances are solved in sequence.
The electron balance is a 1-D root-find for $T_e$:

$$
\frac{3}{2}\,n_{e0}\,T_e\,(10^3 e)\,\frac{V}{\tau_{E,e}}
=
\sum_j L_e(x_j)\,P_{\mathrm{useful},j}
+ P_{\mathrm{aux},e} + P_{\alpha,e}.
$$

The ion temperature then follows in closed form (or, for the neoclassical-ion
option, from a second 1-D root-find — see *Confinement time*):

$$
\frac{3}{2}\,n_{\mathrm{thermal}}\,T_i\,(10^3 e)\,\frac{V}{\tau_{E,i}}
=
\sum_j L_i(x_j)\,P_{\mathrm{useful},j}
+ P_{\mathrm{aux},i} + P_{\alpha,i},
$$

where $n_{\mathrm{thermal}}=n_{\mathrm{sum}}(n_{e0},Z_{\mathrm{eff}})-n_{b0}$ is
the thermal-ion density left after the fast ions are subtracted. A non-positive
$n_{\mathrm{thermal}}$ marks the operating point infeasible.

### Coupled balance (equipartition on)

When electron-ion equipartition is enabled, the exchange power
$P_{ei}(T_e,T_i)$ (defined in *Electron-ion exchange time*) is subtracted from
the electron balance and added to the ion balance, and the pair is solved
together:

$$
\frac{3}{2}n_{e0}T_e(10^3 e)\frac{V}{\tau_{E,e}}
= \sum_j L_e(x_j)P_{\mathrm{useful},j} + P_{\mathrm{aux},e} + P_{\alpha,e}
  - P_{ei}(T_e,T_i),
$$

$$
\frac{3}{2}n_{\mathrm{thermal}}T_i(10^3 e)\frac{V}{\tau_{E,i}}
= \sum_j L_i(x_j)P_{\mathrm{useful},j} + P_{\mathrm{aux},i} + P_{\alpha,i}
  + P_{ei}(T_e,T_i).
$$

The reported $P_e$ and $P_i$ remain the pure beam split; $P_{ei}$, $P_\alpha$
and $P_{\mathrm{aux}}$ are reported separately.

## Confinement time

$\tau_{E,e}$ and $\tau_{E,i}$ in the balances above are set per channel by the
*Confinement* selector in the Models section; the electron and ion channels are
chosen independently. When a physics-based scaling is selected the manual
`tauE,e` / `tauE,i` inputs are ignored, and the Summary sheet reports the
computed value (or its range across a scan) instead.

### Fixed (input)

$\tau_{E,e}$ and $\tau_{E,i}$ are the values entered on the Plasma tab,
used unchanged.

### Kaye NSTX L-mode

An NSTX low-aspect-ratio (spherical tokamak) L-mode global fit
(Kaye, *Nucl. Fusion* **46** (2006) 848, Eq. 5 / Fig. 9):

$$
\tau_E = 4.73\times10^{-4}\;
I_p^{1.01}\,B_t^{0.70}\,n_e^{-0.07}\,P_{\mathrm{loss}}^{-0.37}
$$

with $I_p$ in **A**, $B_t$ in **T**, $n_e$ in **m$^{-3}$**, $P_{\mathrm{loss}}$
in **W** and $\tau_E$ in **s**. The loss power is the useful beam power
$P_{\mathrm{loss}}=\sum_j P_{\mathrm{useful},j}$ (auxiliary and alpha heating
are not yet included in $P_{\mathrm{loss}}$). It does not depend on $T_e$ or
$T_i$, so it is evaluated once. Selecting *Kaye NSTX L-mode* uses it for both
channels ($\tau_{E,e}=\tau_{E,i}$).

### Neoclassical ion

A banana-regime order-of-magnitude estimate (Wesson, *Tokamaks*; Chang-Hinton
low-collisionality limit) for the ion channel only:

$$
\chi_{i,\mathrm{neo}} \simeq
\frac{q^2\,\rho_i^2\,\nu_{ii}}{\varepsilon^{3/2}},
\qquad
\tau_{E,i} = \frac{a^2}{\chi_{i,\mathrm{neo}}}
= \frac{a^2\,\varepsilon^{3/2}}{q^2\,\rho_i^2\,\nu_{ii}},
$$

with $\varepsilon=a/R_0$, the elongation-corrected cylindrical safety factor

$$
q = \frac{5\,a^2 B_t}{R_0\,I_p[\mathrm{MA}]}\,\frac{1+\kappa^2}{2},
$$

the thermal-ion Larmor radius

$$
\rho_i = \frac{m_{\mathrm{eff}} v_{th,i}}{e B_t},
\qquad
v_{th,i} = \sqrt{\frac{2 T_i (10^3 e)}{m_{\mathrm{eff}}}},
\qquad
m_{\mathrm{eff}} = x_D m_D + x_T m_T,
$$

and the ion-ion collision frequency summed over D and T,

$$
\nu_{ii} = \sum_{s\in\{D,T\}}
4.80\times10^{-8}\,
\frac{n_s[\mathrm{cm}^{-3}]\,\ln\Lambda}
     {\sqrt{\mu_s}\,T_i[\mathrm{eV}]^{3/2}},
\qquad
\mu_s = m_s/m_p.
$$

Since $\tau_{E,i}\propto T_i^{1/2}$, the ion balance becomes a 1-D fixed-point
root-find for $T_i$ rather than a closed-form substitution. Pure neoclassical
transport has no anomalous channel, so it under-estimates real ion transport
and can predict very hot, sometimes implausible $T_i$; treat it as a lower
bound. It cannot be combined with electron-ion equipartition.

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

When profile peaking is nonzero, the fusion post-processing uses volume
integrals. The thermal contribution becomes:

$$
P_{f,\mathrm{thermal}}=V E_f\left\langle
n_D(\rho)n_T(\rho)\langle\sigma v\rangle_{DT}[T_i(\rho)]\right\rangle_V.
$$

The beam-target contribution similarly uses the local target density and the
local slowing-down distribution evaluated with $T_e(\rho)$. Pressure and
thermal-energy diagnostics use the corresponding volume-averaged $nT$
profiles. The electron and ion power-balance solve remains a 0D central-value
solve; peaking changes those temperatures only indirectly through the
profile-aware shine-through and captured beam power.

### Current model limitations

- The beam-target integral uses stationary target ions, so the solved $T_i$ is not included in that reaction integral.
- D-D beam-target reactions are not included.
- Each included D-T reaction is assigned $17.6\,\mathrm{MeV}$ of fusion energy.
- The radial volume weighting uses the model's elliptical flux-surface
	approximation and does not yet solve a 2D equilibrium.

## Density scan

The Plasma tab defines the requested scan range with `n_e_min` and `n_e_max`. The Results tab reports the physically valid interval separately. Points below the HotJass charge-neutrality feasibility limit are not assigned physical output values and are omitted from plots.

## Important notation

- $T_e$: electron temperature
- $T_i$: thermal-ion temperature
- $P_e$: beam power deposited into electrons
- $P_i$: beam power deposited into ions
- $L_e$, $L_i$: beam-power electron / ion heating fractions ($L_e+L_i=1$)
- $P_{ie}$: electron-ion energy exchange power
- $P_{\mathrm{ECRH}}$, $P_{\mathrm{ICRH}}$: ECRH and ICRH auxiliary heating power
- $P_{\mathrm{aux},e}$, $P_{\mathrm{aux},i}$: auxiliary heating to electrons / ions
- $P_\alpha$, $P_{\alpha,e}$, $P_{\alpha,i}$: fusion alpha self-heating power (total, to electrons, to ions)
- $\tau_{E,e}$, $\tau_{E,i}$: electron / ion energy confinement times used in the balance
- $n_{\mathrm{thermal}}$: thermal-ion density, $n_{\mathrm{sum}}(n_{e0},Z_{\mathrm{eff}})-n_{b0}$
- $n_D$, $n_T$: thermal deuterium and tritium densities
- $n_{b0}$: total fast-ion density
- $\rho_{Li}$: fast-ion toroidal-field Larmor radius at the injection energy
- $\Delta r$: passing drift-orbit radial width, $q_\ast\,\rho_{Li}$
- $f_{\mathrm{orbit}}$: first-orbit loss fraction of captured beam power
- $P_{\mathrm{orbit}}$: first-orbit loss power
- $P_{f,tot}$, $P_{f,th}$, $P_{f,b}$: total, thermal, and beam-target fusion power
- $\tau_{IE}$: effective electron-ion exchange time
- $p_{th}$, $p_{fast}$: thermal and fast-ion pressure
- $\beta_t$: toroidal beta
