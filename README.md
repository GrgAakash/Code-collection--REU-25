# 📁 IPC Modeling Codes

This repository contains MATLAB codes related to the IPC model simulations and case studies presented in our summer 2025 REU project. Each script is accompanied by a brief description of its function and corresponding figures in the report.

---
## 📌 Overview

- **SIR Models** 
- **ODE Simulations**

Please refer to the descriptions below for specific information about each file.

---

## 🧠 IPC Simulation Codes

| File | Description | Related Figure |
|------|-------------|----------------|
| `siripc.m` | SIR-IPC model code (**no agent tracking**) — simulates the population-level dynamics without tracking individual agents. |Replicates result of Figure 1 of the paper|
| `sir_combined_model.m` | SIR-IPC model with **agent tracking** - Includes logic to trace individual agent states over time. | Replicates result of Figure 1 of the paper|

---

## 🔄 ODE-Based Simulations

| File | Description | Related Figure |
|------|-------------|----------------|
| `varpreF.m` | |Replicates result of Figure 2 of the paper|

---


## 💡 Suggestions Welcome

<p>CWang: great idea. I would suggest to write some more detailed descripton of each code. For example, for the IPC code, please specify whether the agents are tracked or not. I remember some version like Caden's can track agents while Aakash's do not. It is good to have both of them here, but with some descriptions to differentiate these two. Now for the ode codes it is may be less writing for now. For the case study codes, I imagine we should at least say things about which city/county this case is about etc. </p>

<p>All in all, it is maybe better to give more details than less, especially at this early stage so that everyone knows what a particular code is about, even if he/she does not write the code. </p>

