# SIHRS Epidemic Model with Death Project

## Overview

This repository contains implementations of the **SIHRS (Susceptible-Infected-Hospitalized-Recovered-Susceptible)** epidemic model with death compartments. The project provides both deterministic (ODE) and stochastic simulation approaches, implemented in multiple programming languages and applied to real-world case studies.

## Project Structure

```
SIHRS/
├── CaseStudies/                    # Real-world applications and data analysis
│   ├── Carson_City,NV/            # Nevada case study
│   │   ├── codes/                 # Data extraction and simulation scripts
│   │   │   ├── Julia/             # Julia implementations
│   │   │   │   ├── SIHRS_hospitalized.jl         # Hospitalization simulation
│   │   │   │   ├── SIHRS_multiple_simulations_infected.jl  # Infected cases simulation
│   │   │   │   ├── calculate_carson_pih.jl       # P(IH) parameter calculation
│   │   │   │   ├── extract_active_cases.jl       # Data extraction script
│   │   │   │   ├── extract_active_deaths.jl      # Active deaths calculation
│   │   │   │   └── extract_daily_deaths.jl       # Daily deaths extraction
│   │   │   ├── Matlab/            # MATLAB implementations
│   │   │   │   ├── SIHRS_hospitalized.m          # Hospitalization simulation
│   │   │   │   └── SIHRS_multiple_simulations_infected.m  # Infected cases simulation
│   │   │   ├── Extras/            # Additional analysis and results
│   │   │   │   ├── carson_city_pih_lagged_14days.csv      # P(IH) analysis results
│   │   │   │   ├── PIH_Calculation_Methodology_Carson.md  # P(IH) documentation
│   │   │   │   ├── carson_city_hospitalization_active_cases_combined.csv  # Combined dataset
│   │   │   │   ├── calculate_pid_average.jl      # PID calculation scripts
│   │   │   │   ├── calculate_pid_lagged.jl       # PID calculation scripts
│   │   │   │   ├── calculate_pid_simple.jl       # PID calculation scripts
│   │   │   │   └── combine_active_cases_hospitalization.jl # Data combination script
│   │   │   ├── carson_city_hospitalization_active_cases_combined.csv  # Combined dataset
│   │   │   ├── hospitalization_Carson_filtered_new.csv  # Cleaned hospitalization data
│   │   │   ├── carson_city_active_cases.csv  # Active cases data
│   │   │   ├── carson_city_daily_deaths.csv  # Daily deaths data
│   │   │   ├── carson_city_combined.csv      # Combined cases and deaths
│   │   │   └── SIHRS_hospitalization_cheating.jl # Additional analysis
│   │   └── Results/               # Generated plots and visualizations
│   ├── Washington, MS/            # Mississippi case study
│   │   ├── codes/                 # Data extraction and simulation scripts
│   │   │   ├── Julia/             # Julia implementations
│   │   │   │   ├── SIHRS_hospitalized.jl         # Hospitalization simulation
│   │   │   │   ├── SIHRS_multiple_simulations_infected.jl  # Infected cases simulation
│   │   │   │   ├── calculate_pih.jl              # P(IH) parameter calculation
│   │   │   │   ├── extract_washington_ms_active_cases.jl   # Data extraction script
│   │   │   │   ├── extract_washington_ms_daily_deaths.jl   # Daily deaths extraction
│   │   │   │   └── test_hospitalization_data.jl            # Data validation script
│   │   │   ├── Matlab/            # MATLAB implementations
│   │   │   │   ├── SIHRS_hospitalized.m          # Hospitalization simulation
│   │   │   │   └── SIHRS_multiple_simulations_infected.m  # Infected cases simulation
│   │   │   ├── Extras/            # Additional analysis and results
│   │   │   │   ├── washington_ms_pih_lagged_14days.csv     # P(IH) analysis results
│   │   │   │   ├── PIH_Calculation_Methodology.md         # P(IH) documentation
│   │   │   │   ├── washington_ms_hospitalization_active_cases_combined.csv  # Combined dataset
│   │   │   │   ├── calculate_pid_average_washington_ms.jl # PID calculation scripts
│   │   │   │   ├── calculate_pid_lagged_washington_ms.jl # PID calculation scripts
│   │   │   │   ├── calculate_pid_simple_washington_ms.jl # PID calculation scripts
│   │   │   │   └── combine_active_cases_hospitalization.jl # Data combination script
│   │   │   ├── hospitalization_MS_filtered.csv        # Cleaned hospitalization data
│   │   │   ├── washington_mississippi_active_cases.csv # Active cases data
│   │   │   ├── washington_mississippi_daily_deaths.csv # Daily deaths data
│   │   │   └── washington_mississippi_combined.csv     # Combined cases and deaths
│   │   └── Results/               # Generated plots and visualizations
│   ├── R_0/                       # R₀ calculation data and analysis
│   │   ├── county_R0_data.csv     # R₀ statistics by county
│   │   ├── filtered_R0_population_fips.csv    # Filtered R₀ data
│   │   └── readme.md              # R₀ analysis documentation
│   ├── case_study.ipynb           # Jupyter notebook for case study analysis
│   ├── readme.md                  # Case studies overview documentation
│   └── hospitalization_readme.md  # Hospitalization analysis documentation
├── SIHRS_main_code/               # Core model implementations
│   ├── Julia/                     # Julia implementations
│   │   ├── Project.toml           # Julia package dependencies
│   │   ├── Manifest.toml          # Julia package manifest
│   │   ├── SIHRS_jl/             # Full SIHRS model (deterministic + stochastic)
│   │   │   ├── SIHRS.jl          # Main simulation script
│   │   │   └── Images/           # Generated plots and results
│   │   ├── SIHRS_I_H_jl/         # Simplified I-H (Infected-Hospitalized) model
│   │   │   ├── SIHRS_I_H_only.jl # I-H simulation script
│   │   │   └── Images/           # Generated plots and results
│   │   ├── sihr.jl               # Basic SIHR model implementation
│   │   └── SIHRS_variance_analysis.jl # Variance analysis for stochastic simulations
│   ├── Matlab/                    # MATLAB implementations
│   │   ├── SIHRS.m               # Full SIHRS model
│   │   ├── SIHRS_I_H_only.m     # I-H only model
│   │   ├── SIHRS_variance_analysis.m # Variance analysis for stochastic simulations
│   │   ├── SIHRS_variance_sqrtroot_analysis.m # Square root variance analysis
│   │   ├── disp_SIHRS_renormalized.m # Renormalized SIHRS display
│   │   ├── justification.md      # Parameter justification documentation
│   │   └── README_SIHRS_Variance_Analysis.md # Variance analysis documentation
│   ├── sihr.m                    # Basic SIHR model implementation
│   └── README.md                  # Detailed mathematical model documentation
└── README.md                      # This file - project overview
```

## Mathematical Model
For a breif overview of SIHRS epidemic model with death compartment, its parameter definitions and constraints, see `SIHRS_main_code/README.md`.

## Implementation Languages

### Julia
- **Full SIHRS Model**: `SIHRS_main_code/Julia/SIHRS_jl/SIHRS.jl`
  - Deterministic ODE solutions
  - Stochastic simulations with configurable population sizes
  - Comprehensive plotting and analysis
- **I-H Only Model**: `SIHRS_main_code/Julia/SIHRS_I_H_jl/SIHRS_I_H_only.jl`
  - Simplified two-compartment model
  - Focused on infection-hospitalization dynamics

### MATLAB
- **Full SIHRS Model**: `SIHRS_main_code/Matlab/SIHRS.m`
- **I-H Only Model**: `SIHRS_main_code/Matlab/SIHRS_I_H_only.m`

## Case Studies

### Carson City, Nevada
- **Location**: `CaseStudies/Carson_City,NV/`
- **Data**: Active cases, daily deaths, hospitalization records
- **Scripts**: Data extraction, multiple simulation runs, hospitalization analysis
- **Results**: Full pandemic trajectories, bandwidth analysis, cumulative deaths

### Washington, Mississippi
- **Location**: `CaseStudies/Washington,MS/`
- **Data**: Active cases, daily deaths, hospitalization records
- **Scripts**: Data extraction, multiple simulation runs, hospitalization analysis
- **Results**: Full pandemic trajectories, bandwidth analysis, cumulative deaths
