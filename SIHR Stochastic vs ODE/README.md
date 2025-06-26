## 📁 Project Structure

```
StochVsODE/
├── index.html                 # Main HTML file
├── README.md                  # This file
├── src/
│   ├── css/
│   │   └── styles.css         # All CSS styles
│   ├── js/
│   │   ├── main.js            # Main application entry point
│   │   ├── models/
│   │   │   └── sihr-model.js  # SIHR model implementation
│   │   ├── utils/
│   │   │   └── audio.js       # Audio/sound effects
│   │   └── components/
│   │       ├── pattern-analysis.js  # Pattern analysis functionality
│   │       ├── ui-controls.js       # UI controls and animations
│   │       └── download-utils.js    # Download functionality
│   └── assets/
│       ├── images/            # Image assets (future use)
│       └── audio/             # Audio assets (future use)
├── docs/                      # Documentation (future use)
└── public/                    # Public assets (future use)
```


## Mathematical Background

The SIHR model extends the classic SIR model by adding a hospitalized compartment:

```
dS/dt = -p₁βSI
dI/dt = p₁βSI - p₂γI
dH/dt = p₂pₕγI - p₃αH
dR/dt = p₂(1-pₕ)γI + p₃αH
```

Where:
- S: Susceptible individuals
- I: Infected individuals  
- H: Hospitalized individuals
- R: Recovered individuals
