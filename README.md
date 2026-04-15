# SIR_-Guias
WepApp for simulation  : Analyzing the Impact of the Parameter Values on the Qualitative Behaviour of the Solutions of a Modified SIR Model with Vaccination and Several Levels of Immunity , Prof. F. Guiaș 
# SIR Immunity Simulator — Prof. F. Guiaș


Interactive browser implementation of the extended SIR epidemiological model developed by **Prof. Dr. Flavius Guiaș** (FH Dortmund, Department of Mechanical Engineering). The model extends the classical SIR framework with **365 discrete immunity levels** and daily immunity decay.

Solves **367 coupled differential equations** in real time — directly in the browser, no server, no build step required. Validated against the original MATLAB implementation with a deviation of **< 0.1%**.

🔗 **Live Demo:** [sir-guias.vercel.app](https://sir-guias.vercel.app)

---

## Quick Start

```bash
git clone https://github.com/Tools00/SIR_-Guias
# Open index.html in any browser — done
```

> **Note:** An internet connection is required. Chart.js 3.9.1 is loaded via CDN.

---

## Features

- **365 immunity levels** — daily immunity decay modeled over a full year
- **Real-time parameters** — R₀, α, γ₀ adjustable via interactive sliders
- **4 visualization modes** — SIR curves, immunity pyramid, phase space plot, stability analysis
- **Two scenarios** — simplified model and full model with vaccination
- **No local dependencies** — only Chart.js 3.9.1 via CDN

---

## The ODE System

The model describes a population split into immunity levels `k = 0, 1, ..., m` with `m = 365`, plus an infected class `I`:

```
// Level k=0 (fully susceptible):
dS/dt   = −β₀·I·S − γ₀·S

// Levels k = 1 ... m−1 (partial immunity):
dSₖ/dt  = βₖ·I·Sₖ₋₁ − βₖ·I·Sₖ − εₖ·Sₖ + εₖ₊₁·Sₖ₊₁ − γₖ·Sₖ

// Level k=m (maximum immunity, vaccination possible):
dSm/dt  = α·I + Σγₖ·Sₖ − εm·Sm − βm·I·Sm

// Infected:
dI/dt   = Σ βₖ·I·Sₖ − α·I
```

The reproduction number decays exponentially with immunity level:

```
Rₖ = R₀ · exp(−α · k)     →     βₖ = a · Rₖ
```

---

## Numerical Solver: 4th-Order Runge-Kutta

The JavaScript implementation replicates MATLAB's `ode45` using a classical RK4 scheme with fixed time step `h = 0.2 days` (5 steps per day):

```
k₁ = f(tₙ,        yₙ)
k₂ = f(tₙ + h/2,  yₙ + h·k₁/2)
k₃ = f(tₙ + h/2,  yₙ + h·k₂/2)
k₄ = f(tₙ + h,    yₙ + h·k₃)

yₙ₊₁ = yₙ + h/6 · (k₁ + 2·k₂ + 2·k₃ + k₄)
```

**MATLAB reference:**
```matlab
options = odeset('RelTol', 2.23e-14, 'AbsTol', 2.23e-14);
[t, y] = ode45(@ode_imm_f_new, tspan, y0, options, m, b, g, e, a);
```

**JavaScript implementation:**
```javascript
function rungeKutta4(f, t, y, h) {
    const k1 = f(t, y);
    const k2 = f(t + h/2, y.map((yi, i) => yi + h * k1[i] / 2));
    const k3 = f(t + h/2, y.map((yi, i) => yi + h * k2[i] / 2));
    const k4 = f(t + h,   y.map((yi, i) => yi + h * k3[i]));
    return y.map((yi, i) => yi + h * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]) / 6);
}
```

---

## Parameters

| Parameter | Symbol | Range | Description |
|---|---|---|---|
| Basic reproduction number | R₀ | 2.0 – 8.0 | Average secondary infections |
| Immunity decay coefficient | α | 0.008 – 0.020 | Exponent of β decay |
| Recovery rate | a = 1/T_inf | 0.05 – 0.5 | Inverse of infection duration |
| Vaccination rate | γ₀ | 0.001 – 0.020 | Daily vaccination quota (levels ≤ k₀) |
| Immunity levels | m | 100 – 365 | Discretization depth |
| Time step | h | 0.2 days | Fixed RK4 step size |

---

## Validation

The JavaScript solution was validated against the original MATLAB implementation by Prof. Dr. Guiaș:

> Maximum deviation **< 0.1%** over 1500 simulation days for R₀ ∈ {2, 4, 6, 8} and α ∈ {0.008, 0.012, 0.016, 0.020}.

---

## Repository Structure

| File | Description |
|---|---|
| `index.html` | Complete application — HTML, CSS and JavaScript (1269 lines) |
| `README.md` | This documentation |
| `.gitignore` | Ignores OS-specific files |

---

## Credits & References

**Theoretical Foundation**

This project implements the model developed by Prof. Dr. Flavius Guiaș. Key publications:

- Guiaș, F. (2023). *Analyzing the Impact of the Parameter Values on the Qualitative Behaviour of the Solutions of a Modified SIR Model with Vaccination and Several Levels of Immunity.* FH Dortmund, Department of Mechanical Engineering.

- Guiaș, F. (2023). *Equilibrium Solutions of a Modified SIR Model with Vaccination and Several Levels of Immunity.* WSEAS Transactions on Systems and Control, Vol. 18, pp. 550–560. [doi:10.37394/23203.2023.18.57](https://doi.org/10.37394/23203.2023.18.57)

- Guiaș, F. (2023). *Epidemic models with several levels of immunity.* In: Skiadas, C.H. (eds.), Quantitative Demography and Health Estimates. Springer Series on Demographic Methods, Vol. 55, pp. 163–174. [doi:10.1007/978-3-031-28697-1](https://doi.org/10.1007/978-3-031-28697-1)

**Classical SIR Model**

- Kermack, W. O., & McKendrick, A. G. (1927). *A contribution to the mathematical theory of epidemics.* Proceedings of the Royal Society A, 115(772), 700–721.

**Visualization**

- Downie, N. (2023). *Chart.js Documentation, Version 3.9.1.* [chartjs.org](https://www.chartjs.org/docs/3.9.1/)

---

**Web implementation:** Engineering Project (Ingenieurmäßige Arbeit), FH Dortmund 2025

---

## License

MIT License — free to use for academic and research purposes.
