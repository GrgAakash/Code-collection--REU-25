# SIHRS Epidemic Model Website

## Overview

This website provides an interactive simulation of the **SIHRS epidemic model** with death and reinfection compartments. The model extends the classic SIR model by adding hospitalization and death, while allowing recovered individuals to become susceptible again.

## Mathematical Model

The SIHRS model defines $s, i, h, r, d : (0, \infty) \to \mathbb{R}_{\geq 0}$, with $s + i + h + r + d = 1$, by:

$$
\begin{aligned}
\dot{s} &= -\beta p_{SI} \, s i + p_{RS} \Lambda r \\
\dot{i} &= \beta p_{SI} \, s i - \gamma (1 - p_{II}) i \\
\dot{h} &= p_{IH} \gamma i - \alpha (1 - p_{HH}) h \\
\dot{r} &= p_{IR} \gamma i + p_{HR} \alpha h - p_{RS} \Lambda r \\
\dot{d} &= p_{ID} \gamma i + p_{HD} \alpha h
\end{aligned}
$$

### Parameters

**Transition Rates:**
- $\beta$: Infection rate
- $\gamma$: I transition rate  
- $\alpha$: H transition rate
- $\Lambda$: R transition rate (reinfection rate)

**Transition Probabilities:**
- $p_{SI}$: S → I probability
- $p_{II}$: I → I probability (stay infected)
- $p_{IH}$: I → H probability
- $p_{IR}$: I → R probability
- $p_{ID}$: I → D probability
- $p_{HH}$: H → H probability (stay hospitalized)
- $p_{HR}$: H → R probability
- $p_{HD}$: H → D probability
- $p_{RR}$: R → R probability (stay recovered)
- $p_{RS}$: R → S probability (reinfection)

**Constraints:**
- $p_{II} + p_{IH} + p_{IR} + p_{ID} = 1$
- $p_{HH} + p_{HR} + p_{HD} = 1$
- $p_{RR} + p_{RS} = 1$

## Key Differences from SIHR Model

| Feature | SIHR Model | SIHRS Model |
|---------|------------|-------------|
| Compartments | 4 (S, I, H, R) | 5 (S, I, H, R, D) |
| Reinfection | No | Yes (R → S) |
| Death | No | Yes (I → D, H → D) |
| Immunity | Permanent | Waning |
| Complexity | Basic | Advanced |

## Technical Implementation

- **Frontend**: HTML5, CSS3, JavaScript
- **Charts**: Chart.js with custom styling
- **Simulation**: Custom Gillespie algorithm implementation
- **ODE Solver**: Euler integration method
- **Validation**: Real-time parameter constraint checking

## Files Structure

```
SIHR Stochastic vs ODE/
├── index_sihrs.html          # Main SIHRS website
├── src/
│   ├── js/
│   │   ├── models/
│   │   │   └── sihrs-model.js    # SIHRS model implementation
│   │   ├── components/
│   │   │   ├── ui-controls.js    # UI control functions
│   │   │   ├── pattern-analysis.js # Pattern analysis
│   │   │   └── download-utils.js # Data export utilities
│   │   └── main-sihrs.js         # Main application logic
│   └── css/
│       └── styles.css            # Styling (shared with SIHR)
└── README_SIHRS.md              # This file
```

Requires JavaScript enabled and Chart.js CDN access.

## Demo
For live demo, please visit https://grgaakash.github.io/SIHR%20Stochastic%20vs%20ODE/index_sihrs.html

## License

This project is for educational and research purposes. 