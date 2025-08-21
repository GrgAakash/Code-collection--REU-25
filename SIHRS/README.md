# SIHRS Epidemic Model with Death Project

## Overview

This repository contains implementations of the **SIHRS (Susceptible-Infected-Hospitalized-Recovered-Susceptible)** epidemic model with death compartments. The project provides both deterministic (ODE) and stochastic simulation approaches, implemented in multiple programming languages and applied to real-world case studies.

## Project Structure

```
SIHRS/
├── CaseStudies/                    # Real-world applications and data analysis
│   ├── Carson_City,NV/            # Nevada case study
│   │   ├── codes/                 # Data extraction and simulation scripts
│   │   │   ├── Matlab/            # MATLAB implementations
│   │   │   │   ├── SIHRS_hospitalized.m              # Hospitalization simulation
│   │   │   │   ├── SIHRS_hospitalized_carson_aug02.m # Aug 2nd hospitalization simulation
│   │   │   │   ├── SIHRS_hospitalized_carsonparamsBUT_smallN.m # Small N parameter variant
│   │   │   │   ├── SIHRS_multiple_simulations_infected.m      # Infected cases simulation
│   │   │   │   ├── SIHRS_multiple_simulations_infected_aug02.m # Aug 2nd infected simulation
│   │   │   │   ├── carson_city_daily_deaths.csv      # Daily deaths data
│   │   │   │   ├── carson_city_daily_deaths_aug02.csv # Aug 2nd daily deaths data
│   │   │   │   ├── carson_city_combined.csv          # Combined cases and deaths
│   │   │   │   ├── carson_city_combined_aug02.csv    # Aug 2nd combined data
│   │   │   │   ├── carson_city_active_cases.csv      # Active cases data
│   │   │   │   ├── carson_city_active_cases_aug02.csv # Aug 2nd active cases data
│   │   │   │   └── hospitalization_Carson_filtered_new.csv # Cleaned hospitalization data
│   │   │   ├── Extras/            # Additional analysis and results
│   │   │   │   ├── carson_city_pih_lagged_14days.csv      # P(IH) analysis results
│   │   │   │   ├── PIH_Calculation_Methodology_Carson.md  # P(IH) documentation
│   │   │   │   ├── carson_city_hospitalization_active_cases_combined.csv # Combined dataset
│   │   │   │   ├── calculate_carson_pih.jl       # P(IH) parameter calculation
│   │   │   │   ├── calculate_pid_average.jl      # PID average calculation
│   │   │   │   ├── calculate_pid_lagged.jl       # PID lagged calculation
│   │   │   │   ├── calculate_pid_simple.jl       # PID simple calculation
│   │   │   │   ├── combine_active_cases_hospitalization.jl # Data combination script
│   │   │   │   ├── carson_city_pid_lagged_20days.csv      # PID lagged analysis results
│   │   │   │   ├── carson_city_pid_simple.csv             # PID simple analysis results
│   │   │   │   └── SIHRS_hospitalization_cheating.jl     # Additional hospitalization analysis
│   │   │   ├── SIHRS_hospitalized.jl         # Hospitalization simulation (Julia)
│   │   │   ├── SIHRS_multiple_simulations_infected.jl # Infected cases simulation (Julia)
│   │   │   ├── extract_active_cases.jl       # Active cases data extraction
│   │   │   ├── extract_active_deaths.jl      # Active deaths calculation
│   │   │   ├── extract_daily_deaths.jl       # Daily deaths extraction
│   │   │   ├── ratio_analysis_carson.jl      # Ratio analysis for Carson City
│   │   │   ├── carson_ratio_analysis_results.csv # Ratio analysis results
│   │   │   ├── carson_target_ratio_dates.csv     # Target ratio dates
│   │   │   ├── carson_city_daily_deaths.csv  # Daily deaths data
│   │   │   ├── carson_city_combined.csv      # Combined cases and deaths
│   │   │   ├── carson_city_active_cases.csv  # Active cases data
│   │   │   └── hospitalization_Carson_filtered_new.csv # Cleaned hospitalization data
│   │   └── Results/               # Generated plots and visualizations (15 PNG files)
│   ├── Washington, MS/            # Mississippi case study
│   │   ├── codes/                 # Data extraction and simulation scripts
│   │   │   ├── Matlab/            # MATLAB implementations
│   │   │   │   ├── SIHRS_hospitalized.m          # Hospitalization simulation
│   │   │   │   ├── SIHRS_hospitalized_aug30.m    # Aug 30th hospitalization simulation
│   │   │   │   ├── SIHRS_hospitalized_washingtonparamsBUT_smallN.m # Small N variant
│   │   │   │   ├── SIHRS_multiple_simulations_infected.m      # Infected cases simulation
│   │   │   │   ├── SIHRS_multiple_simulations_infected_aug30.m # Aug 30th infected simulation
│   │   │   │   ├── washington_mississippi_daily_deaths.csv    # Daily deaths data
│   │   │   │   ├── washington_mississippi_daily_deaths_aug30.csv # Aug 30th daily deaths
│   │   │   │   ├── washington_mississippi_combined.csv        # Combined cases and deaths
│   │   │   │   ├── washington_mississippi_combined_aug30.csv  # Aug 30th combined data
│   │   │   │   ├── washington_mississippi_active_cases.csv    # Active cases data
│   │   │   │   ├── washington_mississippi_active_cases_aug30.csv # Aug 30th active cases
│   │   │   │   └── hospitalization_MS_filtered.csv           # Cleaned hospitalization data
│   │   │   ├── Extras/            # Additional analysis and results
│   │   │   │   ├── washington_ms_pih_lagged_14days.csv     # P(IH) analysis results
│   │   │   │   ├── PIH_Calculation_Methodology.md         # P(IH) documentation
│   │   │   │   ├── washington_ms_hospitalization_active_cases_combined.csv # Combined dataset
│   │   │   │   ├── calculate_pih.jl                        # P(IH) parameter calculation
│   │   │   │   ├── calculate_pid_average_washington_ms.jl  # PID average calculation
│   │   │   │   ├── calculate_pid_lagged_washington_ms.jl   # PID lagged calculation
│   │   │   │   ├── calculate_pid_simple_washington_ms.jl   # PID simple calculation
│   │   │   │   ├── washington_ms_pid_lagged_20days.csv     # PID lagged analysis results
│   │   │   │   └── washington_ms_pid_simple.csv            # PID simple analysis results
│   │   │   ├── SIHRS_hospitalized.jl         # Hospitalization simulation (Julia)
│   │   │   ├── SIHRS_multiple_simulations_infected.jl # Infected cases simulation (Julia)
│   │   │   ├── extract_washington_ms_active_cases.jl   # Active cases data extraction
│   │   │   ├── extract_washington_ms_daily_deaths.jl   # Daily deaths extraction
│   │   │   ├── ratio_analysis.jl             # Ratio analysis for Washington MS
│   │   │   ├── ratio_analysis_results.csv    # Ratio analysis results
│   │   │   ├── target_ratio_dates.csv        # Target ratio dates
│   │   │   ├── hospitalization_MS_filtered.csv        # Cleaned hospitalization data
│   │   │   ├── washington_mississippi_active_cases.csv # Active cases data
│   │   │   ├── washington_mississippi_daily_deaths.csv # Daily deaths data
│   │   │   └── washington_mississippi_combined.csv     # Combined cases and deaths
│   │   └── Results/               # Generated plots and visualizations (15 PNG files)
│   ├── R_0/                       # R₀ calculation data and analysis
│   │   ├── county_R0_data.csv     # R₀ statistics by county
│   │   ├── filtered_R0_population_fips.csv    # Filtered R₀ data with population and FIPS codes
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
│   │   │   └── Images/           # Generated plots and results (10 PNG files)
│   │   ├── SIHRS_I_H_jl/         # Simplified I-H (Infected-Hospitalized) model
│   │   │   ├── SIHRS_I_H_only.jl # I-H simulation script
│   │   │   └── Images/           # Generated plots and results (4 PNG files)
│   │   └── SIHRS_variance_analysis/ # Variance analysis module
│   │       └── SIHRS_variance_analysis.jl # Variance analysis for stochastic simulations
│   ├── Matlab/                    # MATLAB implementations
│   │   ├── SIHRS.m               # Full SIHRS model
│   │   ├── SIHRS_I_H_only.m     # I-H only model
│   │   ├── SIHRS_variance_analysis.m # Variance analysis for stochastic simulations
│   │   ├── SIHRS_variance_sqrtroot_analysis.m # Square root variance analysis
│   │   ├── disp_SIHRS_renormalized.m # Renormalized SIHRS display
│   │   ├── disp_SIHRS_G_renormalized.m # G-renormalized SIHRS display
│   │   ├── parms_justification.md # Parameter justification documentation
│   │   └── README_SIHRS_Variance_Analysis.md # Variance analysis documentation
│   └── README.md                  # Detailed mathematical model documentation
├── disp_SIHRS_dispersion_CWang.m  # SIHRS dispersion analysis (Prof. Wang)
├── disp_SIHRS_raw_vairance_CWang.m # SIHRS raw variance analysis (Prof. Wang)
├── disp_SIHRS_renormalized_Cwang1_test.m # Renormalized SIHRS test (Prof. Wang)
└── README.md                      # This file - project overview
```

## Mathematical Model
For a brief overview of the SIHRS epidemic model with death compartment, its parameter definitions and constraints, see `SIHRS_main_code/README.md`.

## Implementation Languages

### Julia
- **Full SIHRS Model**: `SIHRS_main_code/Julia/SIHRS_jl/SIHRS.jl`
  - Deterministic ODE solutions
  - Stochastic simulations with configurable population sizes
  - Comprehensive plotting and analysis (10 generated plots)
- **I-H Only Model**: `SIHRS_main_code/Julia/SIHRS_I_H_jl/SIHRS_I_H_only.jl`
  - Simplified two-compartment model focused on infection-hospitalization dynamics
  - Generated visualization plots (4 plots)
- **Variance Analysis**: `SIHRS_main_code/Julia/SIHRS_variance_analysis/SIHRS_variance_analysis.jl`
  - Stochastic variance analysis for model validation

### MATLAB
- **Full SIHRS Model**: `SIHRS_main_code/Matlab/SIHRS.m`
  - Complete SIHRS implementation with all compartments
- **I-H Only Model**: `SIHRS_main_code/Matlab/SIHRS_I_H_only.m`
  - Simplified infection-hospitalization model
- **Variance Analysis Suite**:
  - `SIHRS_variance_analysis.m` - Standard variance analysis
  - `SIHRS_variance_sqrtroot_analysis.m` - Square root variance analysis
  - `disp_SIHRS_renormalized.m` - Renormalized display functions
  - `disp_SIHRS_G_renormalized.m` - G-renormalized display functions

### Additional Analysis Tools
- **Root-level MATLAB scripts** (Prof. Wang's contributions):
  - `disp_SIHRS_dispersion_CWang.m` - Dispersion analysis
  - `disp_SIHRS_raw_vairance_CWang.m` - Raw variance analysis
  - `disp_SIHRS_renormalized_Cwang1_test.m` - Renormalization testing

## Case Studies

### Carson City, Nevada
- **Location**: `CaseStudies/Carson_City,NV/`
- **Data Sources**: 
  - Active cases, daily deaths, and hospitalization records
  - Multiple date snapshots (including August 2nd data)
  - Ratio analysis results and target ratio dates
- **Key Scripts**: 
  - Data extraction (Julia): `extract_active_cases.jl`, `extract_daily_deaths.jl`
  - Simulation runs: `SIHRS_hospitalized.jl`, `SIHRS_multiple_simulations_infected.jl`
  - Ratio analysis: `ratio_analysis_carson.jl`
  - Parameter calculations (P(IH), PID): Located in `Extras/` folder
- **MATLAB Implementations**: Multiple variants including August 2nd and small-N parameter versions
- **Results**: 15 visualization plots including full pandemic trajectories, bandwidth analysis, cumulative deaths, and hospitalization trends

### Washington, Mississippi
- **Location**: `CaseStudies/Washington, MS/`
- **Data Sources**: 
  - Active cases, daily deaths, and hospitalization records
  - Multiple date snapshots (including August 30th data)
  - Ratio analysis results and target ratio dates
- **Key Scripts**: 
  - Data extraction (Julia): `extract_washington_ms_active_cases.jl`, `extract_washington_ms_daily_deaths.jl`
  - Simulation runs: `SIHRS_hospitalized.jl`, `SIHRS_multiple_simulations_infected.jl`
  - Ratio analysis: `ratio_analysis.jl`
  - Parameter calculations (P(IH), PID): Located in `Extras/` folder
- **MATLAB Implementations**: Multiple variants including August 30th and small-N parameter versions
- **Results**: 15 visualization plots including full pandemic trajectories, bandwidth analysis, cumulative deaths, and hospitalization trends

### R₀ Analysis
- **Location**: `CaseStudies/R_0/`
- **Data**: County-level R₀ statistics with population and FIPS code data
- **Purpose**: Basic reproduction number analysis for different geographical regions

## Key Features

### Parameter Estimation
- **P(IH) Calculations**: Probability of infection leading to hospitalization
  - Lagged analysis (14-day lag)
  - Simple and average calculation methods
  - Methodology documentation included
- **PID Analysis**: Parameter estimation using different approaches
  - Simple, average, and lagged (20-day) calculations
  - Results stored as CSV files for both case studies

### Variance and Dispersion Analysis
- **Stochastic Model Validation**: Comprehensive variance analysis tools
- **Dispersion Studies**: Analysis of model behavior under different population sizes
- **Renormalization Techniques**: Advanced mathematical transformations for model analysis

### Data Processing Pipeline
- **Automated Data Extraction**: Julia scripts for processing raw epidemiological data
- **Data Combination**: Scripts to merge active cases with hospitalization data
- **Ratio Analysis**: Tools for analyzing key epidemiological ratios over time

### Visualization and Results
- **Comprehensive Plotting**: Both deterministic and stochastic simulation results
- **Multiple Population Scales**: Analysis across different population sizes (N=316, N=3162, N=10000)
- **Comparative Analysis**: Bandwidth analysis and trajectory comparisons
- **Publication-Ready Figures**: High-quality PNG outputs for research presentations

## Documentation
- **Mathematical Model**: Detailed ODE system documentation in `SIHRS_main_code/README.md`
- **Parameter Justification**: Rationale for parameter choices in `SIHRS_main_code/Matlab/parms_justification.md`
- **Variance Analysis Guide**: Comprehensive documentation in `README_SIHRS_Variance_Analysis.md`
- **Case Study Analysis**: Jupyter notebook for interactive analysis (`case_study.ipynb`)
- **Methodology Documentation**: P(IH) calculation methodologies for both case studies
