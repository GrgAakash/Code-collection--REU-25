# SIHRS Epidemic Model with Death Project

## Overview

This repository contains implementations of the **SIHRS (Susceptible-Infected-Hospitalized-Recovered-Susceptible)** epidemic model with death compartments. The project provides both deterministic (ODE) and stochastic simulation approaches, implemented in multiple programming languages and applied to real-world case studies.

## Project Structure

```
SIHRS/
├── CaseStudies/
│   ├── Carson_City,NV/
│   │   ├── codes/
│   │   │   ├── Matlab/
│   │   │   │   ├── SIHRS_hospitalized.m
│   │   │   │   ├── SIHRS_hospitalized_carson_aug02.m
│   │   │   │   ├── SIHRS_hospitalized_carsonparamsBUT_smallN.m
│   │   │   │   ├── SIHRS_multiple_simulations_infected.m
│   │   │   │   ├── SIHRS_multiple_simulations_infected_aug02.m
│   │   │   │   ├── carson_city_daily_deaths.csv
│   │   │   │   ├── carson_city_daily_deaths_aug02.csv
│   │   │   │   ├── carson_city_combined.csv
│   │   │   │   ├── carson_city_combined_aug02.csv
│   │   │   │   ├── carson_city_active_cases.csv
│   │   │   │   ├── carson_city_active_cases_aug02.csv
│   │   │   │   └── hospitalization_Carson_filtered_new.csv
│   │   │   ├── Extras/
│   │   │   │   ├── carson_city_pih_lagged_14days.csv
│   │   │   │   ├── PIH_Calculation_Methodology_Carson.md
│   │   │   │   ├── carson_city_hospitalization_active_cases_combined.csv
│   │   │   │   ├── calculate_carson_pih.jl
│   │   │   │   ├── calculate_pid_average.jl
│   │   │   │   ├── calculate_pid_lagged.jl
│   │   │   │   ├── calculate_pid_simple.jl
│   │   │   │   ├── combine_active_cases_hospitalization.jl
│   │   │   │   ├── carson_city_pid_lagged_20days.csv
│   │   │   │   ├── carson_city_pid_simple.csv
│   │   │   │   └── SIHRS_hospitalization_cheating.jl
│   │   │   ├── SIHRS_hospitalized.jl
│   │   │   ├── SIHRS_multiple_simulations_infected.jl
│   │   │   ├── extract_active_cases.jl
│   │   │   ├── extract_active_deaths.jl
│   │   │   ├── extract_daily_deaths.jl
│   │   │   ├── ratio_analysis_carson.jl
│   │   │   ├── carson_ratio_analysis_results.csv
│   │   │   ├── carson_target_ratio_dates.csv
│   │   │   ├── carson_city_daily_deaths.csv
│   │   │   ├── carson_city_combined.csv
│   │   │   ├── carson_city_active_cases.csv
│   │   │   └── hospitalization_Carson_filtered_new.csv
│   │   └── Results/
│   ├── Washington, MS/
│   │   ├── codes/
│   │   │   ├── Matlab/
│   │   │   │   ├── SIHRS_hospitalized.m
│   │   │   │   ├── SIHRS_hospitalized_aug30.m
│   │   │   │   ├── SIHRS_hospitalized_washingtonparamsBUT_smallN.m
│   │   │   │   ├── SIHRS_multiple_simulations_infected.m
│   │   │   │   ├── SIHRS_multiple_simulations_infected_aug30.m
│   │   │   │   ├── washington_mississippi_daily_deaths.csv
│   │   │   │   ├── washington_mississippi_daily_deaths_aug30.csv
│   │   │   │   ├── washington_mississippi_combined.csv
│   │   │   │   ├── washington_mississippi_combined_aug30.csv
│   │   │   │   ├── washington_mississippi_active_cases.csv
│   │   │   │   ├── washington_mississippi_active_cases_aug30.csv
│   │   │   │   └── hospitalization_MS_filtered.csv
│   │   │   ├── Extras/
│   │   │   │   ├── washington_ms_pih_lagged_14days.csv
│   │   │   │   ├── PIH_Calculation_Methodology.md
│   │   │   │   ├── washington_ms_hospitalization_active_cases_combined.csv
│   │   │   │   ├── calculate_pih.jl
│   │   │   │   ├── calculate_pid_average_washington_ms.jl
│   │   │   │   ├── calculate_pid_lagged_washington_ms.jl
│   │   │   │   ├── calculate_pid_simple_washington_ms.jl
│   │   │   │   ├── washington_ms_pid_lagged_20days.csv
│   │   │   │   └── washington_ms_pid_simple.csv
│   │   │   ├── SIHRS_hospitalized.jl
│   │   │   ├── SIHRS_multiple_simulations_infected.jl
│   │   │   ├── extract_washington_ms_active_cases.jl
│   │   │   ├── extract_washington_ms_daily_deaths.jl
│   │   │   ├── ratio_analysis.jl
│   │   │   ├── ratio_analysis_results.csv
│   │   │   ├── target_ratio_dates.csv
│   │   │   ├── hospitalization_MS_filtered.csv
│   │   │   ├── washington_mississippi_active_cases.csv
│   │   │   ├── washington_mississippi_daily_deaths.csv
│   │   │   └── washington_mississippi_combined.csv
│   │   └── Results/
│   ├── R_0/
│   │   ├── county_R0_data.csv
│   │   ├── filtered_R0_population_fips.csv
│   │   └── readme.md
│   ├── case_study.ipynb
│   ├── readme.md
│   └── hospitalization_readme.md
├── SIHRS_main_code/               
│   ├── Julia/                     
│   │   ├── Project.toml           
│   │   ├── Manifest.toml         
│   │   ├── SIHRS_jl/             
│   │   │   ├── SIHRS.jl          
│   │   │   └── Images/           
│   │   ├── SIHRS_I_H_jl/         
│   │   │   ├── SIHRS_I_H_only.jl 
│   │   │   └── Images/           
│   │   └── SIHRS_variance_analysis/ 
│   │       └── SIHRS_variance_analysis.jl 
│   ├── Matlab/                    
│   │   ├── SIHRS.m
│   │   ├── SIHRS_I_H_only.m
│   │   ├── SIHRS_variance_analysis.m
│   │   ├── SIHRS_variance_sqrtroot_analysis.m
│   │   ├── disp_SIHRS_DN_renormalized.m
│   │   ├── disp_SIHRS_G_renormalized.m
│   │   ├── SIHRS_cumulative_R_analysis.m
│   │   ├── SIHRS_cumulative_DN_analysis.m
│   │   ├── generate_population_csv.m
│   │   ├── generate_population_csv.py
│   │   ├── intersection_plot_for_hospitalized.m
│   │   ├── intersection_plot_for_infected.m
│   │   ├── intersection_plot_for_recovered.m
│   │   ├── intersection_plot_for_susceptible.m
│   │   ├── parms_justification.md
│   │   ├── explanation of R_N.pdf
│   │   ├── SIHRS_population_and_DN_data_N1600_R0_1.20.csv
│   │   └── extras/
│   │       ├── calculate_G2_V2_zeros.m
│   │       ├── export_R_values_csv.m
│   │       ├── export_sqrt_variance_csv.m
│   │       ├── recreate_csv_export.m
│   │       ├── SIHRS_G_values_N300_T300.csv
│   │       ├── SIHRS_R_values_N300_T300.csv
│   │       └── SIHRS_sqrt_Variance_N300_T300.csv
│   └── README.md
├── disp_SIHRS_dispersion_CWang.m
├── disp_SIHRS_raw_vairance_CWang.m
├── disp_SIHRS_renormalized_Cwang1_test.m
└── README.md
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
  - `disp_SIHRS_renormalized.m` - D_N renormalized display functions
  - `disp_SIHRS_G_renormalized.m` - G-renormalized display functions
  - `SIHRS_cumulative_R_analysis.m` - Cumulative G-renormalized analysis
  - `SIHRS_cumulative_DN_analysis.m` - Cumulative D_N renormalized analysis

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
  - **D_N Renormalization**: Relative fluctuations compared to population size
  - **G-Renormalization**: Relative fluctuations compared to rate of change
- **Cumulative Analysis**: Long-term accumulation of stochastic effects
  - **Cumulative D_N**: `∫₀ᵗ D_N^(l)(τ) dτ` - Accumulated D_N renormalized noise
  - **Cumulative R_N**: `∫₀ᵗ R_N^(l)(τ) dτ` - Accumulated G-renormalized noise

### Data Processing Pipeline
- **Automated Data Extraction**: Julia scripts for processing raw epidemiological data
- **Data Combination**: Scripts to merge active cases with hospitalization data
- **Ratio Analysis**: Tools for analyzing key epidemiological ratios over time

### Visualization and Results
- **Comprehensive Plotting**: Both deterministic and stochastic simulation results
- **Multiple Population Scales**: Analysis across different population sizes (N=1600, N=3000)
- **Comparative Analysis**: Bandwidth analysis and trajectory comparisons
- **Renormalization Visualizations**: 
  - Individual compartment plots (linear and log scale)
  - Combined multi-compartment plots
  - Blow-up behavior detection and visualization
- **Cumulative Analysis Plots**:
  - Linear and logarithmic scale visualizations
  - Individual and combined compartment analysis
  - Long-term noise accumulation patterns
- **Publication-Ready Figures**: High-quality PNG outputs for research presentations

## Documentation
- **Mathematical Model**: Detailed ODE system documentation in `SIHRS_main_code/README.md`
- **Parameter Justification**: Rationale for parameter choices in `SIHRS_main_code/Matlab/parms_justification.md`
- **Variance Analysis Guide**: Comprehensive documentation in `README_SIHRS_Variance_Analysis.md`
- **Case Study Analysis**: Jupyter notebook for interactive analysis (`case_study.ipynb`)
- **Methodology Documentation**: P(IH) calculation methodologies for both case studies
