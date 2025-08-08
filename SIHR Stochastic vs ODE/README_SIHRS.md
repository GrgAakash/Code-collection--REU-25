# SIHRS Epidemic Model Website

## Overview

This website provides an interactive simulation of the **SIHRS epidemic model** with death and reinfection compartments. The model extends the classic SIR model by adding hospitalization and death, while allowing recovered individuals to become susceptible again.

## Mathematical Model

The SIHRS model defines $s, i, h, r, d : (0, \infty) \to \mathbb{R}_{\geq 0}$, with $s + i + h + r + d = 1$, by:

\begin{align}
\dot{s} &= -\beta p_{SI} ~s i +  p_{RS} \Lambda r, \\
\dot{i} &= \beta p_{SI} ~s i - \gamma ~ (1-p_{II}) i, \\
\dot{h} &= p_{IH} \gamma i - \alpha ~(1-p_{HH}) h, \\
\dot{r} &= p_{IR} \gamma i + p_{HR}\alpha h-p_{RS} \Lambda r,\\
\dot{d} &= p_{ID} \gamma i + p_{HD}\alpha h,
\end{align}

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

## Features

### 1. Interactive Parameter Control
- Real-time adjustment of all model parameters
- Validation of probability constraints
- Visual feedback for parameter changes

### 2. Dual Simulation Modes
- **Stochastic**: Multiple runs using Gillespie algorithm
- **Deterministic**: ODE solution using Euler integration
- Side-by-side comparison

### 3. Advanced Visualization
- Multi-compartment plotting (S, I, H, R, D)
- Uncertainty envelopes for stochastic runs
- Interactive charts with zoom and pan

### 4. Threshold Analysis
- R₀ calculation and display
- σ, σ̃, σ̃̃ threshold parameters
- Real-time threshold updates

### 5. Statistics Dashboard
- Peak infected proportion
- Peak time
- Final death proportion
- Total reinfection count

### 6. Data Export
- Download simulation data as CSV
- Export charts as PNG
- Pattern analysis results

## Key Differences from SIHR Model

| Feature | SIHR Model | SIHRS Model |
|---------|------------|-------------|
| Compartments | 4 (S, I, H, R) | 5 (S, I, H, R, D) |
| Reinfection | No | Yes (R → S) |
| Death | No | Yes (I → D, H → D) |
| Immunity | Permanent | Waning |
| Complexity | Basic | Advanced |

## Usage

1. **Open the website**: Navigate to `index_sihrs.html`
2. **Adjust parameters**: Use sliders or number inputs
3. **Run simulation**: Click "Run Simulation"
4. **Analyze results**: View charts and statistics
5. **Export data**: Use download buttons

## Technical Implementation

- **Frontend**: HTML5, CSS3, JavaScript
- **Charts**: Chart.js with custom styling
- **Simulation**: Custom Gillespie algorithm implementation
- **ODE Solver**: Euler integration method
- **Validation**: Real-time parameter constraint checking

## Default Parameters

Based on realistic COVID-19 modeling:
- $\beta = 0.160$ (infection rate)
- $\gamma = 0.126$ (I transition rate)
- $\alpha = 0.100$ (H transition rate)
- $\Lambda = 0.007$ (reinfection rate)
- $p_{IH} = 0.040$ (4% hospitalization rate)
- $p_{ID} = 0.001$ (0.1% death rate from infection)
- $p_{HD} = 0.002$ (0.2% death rate from hospitalization)
- $p_{RS} = 0.980$ (98% reinfection probability)

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

## Browser Compatibility

- Chrome 80+
- Firefox 75+
- Safari 13+
- Edge 80+

Requires JavaScript enabled and Chart.js CDN access.

## License

This project is for educational and research purposes. 