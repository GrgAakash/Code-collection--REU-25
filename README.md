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


## 💡 Suggestions Welcome

<p>CWang: great idea. I would suggest to write some more detailed descripton of each code. For example, for the IPC code, please specify whether the agents are tracked or not. I remember some version like Caden's can track agents while Aakash's do not. It is good to have both of them here, but with some descriptions to differentiate these two. Now for the ode codes it is may be less writing for now. For the case study codes, I imagine we should at least say things about which city/county this case is about etc. </p>

<p>All in all, it is maybe better to give more details than less, especially at this early stage so that everyone knows what a particular code is about, even if he/she does not write the code. </p>

