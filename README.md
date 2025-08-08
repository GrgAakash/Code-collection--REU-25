## Project Overview

# 📁 IPC Modeling Codes

This repository contains MATLAB codes related to the IPC model simulations and case studies presented in our summer 2025 REU project. Each script is accompanied by a brief description of its function and corresponding figures in the report.

---
## 📌 Overview

- **SIR Models** 
- **SIHR Models** 


Please refer to the descriptions below for specific information about each file.

---

## 🚀SIR IPC Simulation Codes

| File | Description | Related Figure |
|------|-------------|----------------|
| `SIHR/siripc.m` | SIR-IPC model code (**agent tracking as of June 4**) — simulates the population-level dynamics with tracking individual agents.Includes comparison with deterministic ODE solution. |Replicates result of Figure 1 of Xia Li's paper|
| `SIHR/sir_combined_model.m` | SIR-IPC model with **agent tracking** - Includes logic to trace individual agent states over time.Includes comparison with deterministic ODE solution. | Replicates result of Figure 1 of Xia Li's paper|
| `SIHR/varpreF.m` | Stochastic SIR model simulation using the Gillespie algorithm with single run across varying population sizes(although the code includes functionality for multiple runs). Includes comparison with deterministic ODE solution and analysis of standard deviation evolution over time. |Replicates result of Figure 2 of Xia Li's paper|


---

## 🚀SIHR IPC Simulation Codes

| File | Description | Related Figure |
|------|-------------|----------------|
| `sihripcv2.m` | SIHR-IPC model code (**agent tracking as of June 4**) — simulates the population-level dynamics with tracking individual agents.Includes comparison with deterministic ODE solution. |Replicates result of Figure 1 of Xia Li's paper for SIHR|
| `sirhipcKeyton.m` | Includes comparison with deterministic ODE solution. |Replicates result of Figure 1 of Xia Li's paper for SIHR|

---

### Code Results

The `Code Results` folder contains images generated from the simulations.

- To view the output from `siripc.m`, navigate to `Code Results/siripc_results`. For example, `R0.95.png` corresponds to a simulation where \( R_0 = 0.95 \).
- For results from `varpreF.m`, filenames like `0.95n3.png` indicate a simulation with \( R_0 = 0.95 \) and 3 runs. For our research, a single run is typically sufficient.

## 🚀Simulations
| File | Description | Related Figure |
|------|-------------|----------------|
| `sihmd_2d.html`|2D animation |The animation can be found at https://cfbpicks.live/reu/|
| `index_sihrs.html`| SIHRS stochastic vs ODE analysis | The simulation can be found at https://grgaakash.github.io/SIHR%20Stochastic%20vs%20ODE/index_sihrs.html|


## 💡 Suggestions Welcome
**Feel free to put your suggestions.**

