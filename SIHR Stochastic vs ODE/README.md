# Mario's Epidemic Adventure - SIHR Model

A fun and interactive web application that compares stochastic vs deterministic SIHR (Susceptible-Infected-Hospitalized-Recovered) epidemic models with a Mario-themed gaming aesthetic.

## 🎮 Features

- **Interactive Epidemic Simulation** with adjustable parameters (β, γ, α, p₁, p₂, p₃, pₕ)
- **Real-time Visualization** using Chart.js showing both deterministic and stochastic curves
- **Pattern Analysis** that groups similar simulation runs
- **Mario-themed UI** with retro gaming aesthetics and sound effects
- **Matrix-style Background** with falling epidemiological symbols
- **Draggable Mario Coins** for interactive fun
- **Download Functionality** for charts and analysis results

## 📁 Project Structure

```
StochVsODE/
├── index.html                 # Main HTML file
├── old.html                   # Original single-file version
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

## 🚀 Getting Started

1. **Clone or download** the project files
2. **Open `index.html`** in a modern web browser
3. **Interact** with the simulation using the controls

## 🔧 Technical Details

### Models
- **Stochastic SIHR Model**: Uses Gillespie algorithm for discrete-event simulation
- **Deterministic SIHR Model**: Uses RK4 integration for ODE solving
- **Pattern Recognition**: Groups similar simulation runs using correlation analysis

### Technologies Used
- **HTML5** for structure
- **CSS3** for styling and animations
- **JavaScript (ES6+)** for logic
- **Chart.js** for data visualization
- **Web Audio API** for sound effects

### Key Parameters
- **β (beta)**: Infection rate
- **γ (gamma)**: I to H transition rate
- **α (alpha)**: H to R transition rate
- **p₁, p₂, p₃**: Transition probabilities
- **pₕ**: Probability of I to H vs R
- **Population Size**: Number of individuals in simulation
- **Number of Runs**: Stochastic simulation repetitions

## 🎯 Usage

1. **Adjust Parameters**: Use the sliders and number inputs in the sidebar
2. **Apply Changes**: Click "Apply Parameters" to run new simulations
3. **Control Animation**: Use Play/Pause, Reset, and Speed controls
4. **Analyze Patterns**: Click "Analyze" to see pattern analysis
5. **Download Results**: Use download buttons in pattern details

## 🎨 Design Features

- **Mario Theme**: Retro gaming aesthetic with Mario colors and fonts
- **Responsive Design**: Works on desktop and mobile devices
- **Interactive Elements**: Draggable coins, animated backgrounds
- **Sound Effects**: Mario-style audio feedback for interactions
- **Smooth Animations**: CSS transitions and JavaScript animations

## 📊 Analysis Features

- **Real-time Statistics**: Current run, population, R₀, peak hospitalization
- **Pattern Recognition**: Groups similar simulation outcomes
- **Threshold Parameters**: σ, σ̃, σ̃̃, tpi, h(tpi) calculations
- **Download Options**: PNG exports of charts and analysis

## 🔬 Mathematical Background

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

## 🎮 Fun Features

- **Matrix Background**: Falling epidemiological symbols
- **Mario Coins**: Draggable animated coins
- **Sound Effects**: Coin, jump, power-up, and game over sounds
- **Retro Font**: Press Start 2P font for authentic gaming feel
- **Color Scheme**: Mario-inspired color palette

## 📝 License

This project is for educational and research purposes.

## 🤝 Contributing

Feel free to suggest improvements or report issues!

---

**Enjoy exploring epidemic dynamics with Mario! 🍄🏥** 