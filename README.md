## Project Overview
**CWang**
This repository contains MATLAB codes related to the IPC model simulations and case studies presented in our summer 2025 REU project. For the IPC code, please specify whether the agents are tracked or not.

I remember for example:
- **Aakash's version (`siripc.m`)** models population-level infection dynamics without tracking individual agents.
- **Caden's version** includes agent-tracking features, enabling analysis of individual-level behavior during simulations.

It is good to have both of them here, but with some descriptions to differentiate these two. Now for the ode codes it is may be less writing for now. 

Case study codes should be mentioned with specific geographic regions (e.g., city/county-level scenarios), and descriptions of each case should be included in their respective folders.

As this project is in the early stages, we aim to provide detailed documentation for each file to ensure clarity and collaboration among team members.

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
| `siripc.m` | SIR-IPC model code (**no agent tracking**) — simulates the population-level dynamics without tracking individual agents.Includes comparison with deterministic ODE solution. |Replicates result of Figure 1 of the paper|
| `sir_combined_model.m` | SIR-IPC model with **agent tracking** - Includes logic to trace individual agent states over time.Includes comparison with deterministic ODE solution. | Replicates result of Figure 1 of the paper|
| `varpreF.m` | Stochastic SIR model simulation using the Gillespie algorithm with single run across varying population sizes(although the code includes functionality for multiple runs). Includes comparison with deterministic ODE solution and analysis of standard deviation evolution over time. |Replicates result of Figure 2 of the paper|


---

## 🚀SIHR IPC Simulation Codes

| File | Description | Related Figure |
|------|-------------|----------------|
| 

---

### Code Results

The `Code Results` folder contains images generated from the simulations.

- To view the output from `siripc.m`, navigate to `Code Results/siripc_results`. For example, `R0.95.png` corresponds to a simulation where \( R_0 = 0.95 \).
- For results from `varpreF.m`, filenames like `0.95n3.png` indicate a simulation with \( R_0 = 0.95 \) and 3 runs. For our research, a single run is typically sufficient.




## 💡 Suggestions Welcome
**Feel free to put your suggestions.**

