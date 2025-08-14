# Case Study - SIHRS Model Implementation

The relevant case and death data by county are obtained from the [New York Times COVID-19 Data repository](https://github.com/nytimes/covid-19-data).

## Overview

This directory contains the implementation of the SIHRS (Susceptible-Infected-Hospitalized-Recovered-Susceptible) epidemic model with death compartments for Carson City, Nevada. The models simulate COVID-19 dynamics from March 2020 to December 2021 using both deterministic and stochastic approaches.

## Population & Data Context

- **Location**: Carson City, Nevada & Washington, MS
- **Population**: 56,000 residents & 43,000 residents.
- **Data Sources**: Real COVID-19 case counts, deaths, and hospitalization records

## Code Files Overview

### 1. `SIHRS_multiple_simulations_infected.jl` - Main Pandemic Model
**Purpose**: Simulates the full COVID-19 pandemic trajectory from Patient Zero

**What it does**:
- Starts simulation from **March 25, 2020** (first detected case)
- Runs **59 stochastic simulations** for uncertainty quantification
- Focuses on **infected cases and deaths** dynamics
- Uses real Carson City data for initial conditions and validation
- Generates comprehensive pandemic trajectory analysis

**Key Features**:
- **Typed structs** for optimal performance (`SIHRSParams`)
- **R₀ = 1.23** (targeted reproduction number)
- **Parameter validation** with mathematical constraints
- **Real data integration** for initial conditions and comparison

**Outputs**:
- `SIHRS_Carson_City_Full_Pandemic_bandwidth.png` - 90% prediction intervals
- `SIHRS_Carson_City_Full_Pandemic_trajectories.png` - All stochastic trajectories
- `SIHRS_Carson_City_Full_Pandemic_active_deaths.png` - Active death proportions
- `SIHRS_Carson_City_Full_Pandemic_cumulative_deaths.png` - Cumulative deaths

**Use Case**: Understanding the full pandemic trajectory and uncertainty from the beginning

---

### 2. `SIHRS_hospitalized.jl` - Hospitalization Dynamics Model
**Purpose**: Analyzes hospitalization patterns during the pandemic

**What it does**:
- Simulates from **March 25, 2020** (same start as main model)
- Runs **55 stochastic simulations** 
- Tracks **hospitalization counts** in addition to infections/deaths
- Compares simulation results with real hospitalization data
- Focuses on healthcare resource planning

**Key Features**:
- **Dict-based parameters** (legacy implementation)
- **R₀ = 1.23** (targeted reproduction number)
- **Hospitalization data integration** from real records
- **Active/Daily death calculation** 

**Outputs**:
- `SIHRS_Carson_City_Hospitalization_bandwidth.png` - Hospitalization prediction intervals
- `SIHRS_Carson_City_Hospitalization_trajectories.png` - Hospitalization trajectories

**Use Case**: Healthcare capacity planning and hospitalization forecasting

---

### 3. `SIHRS_hospitalization_cheating.jl` - Mid-Pandemic Forecasting Model
**Purpose**: "Cheating" simulation starting from a known mid-pandemic state. "Cheating" is a tongue in check word as we are using R₀ which is for the early stage exponential growth (Around March - April) but are using active cases, deaths, hospitalized data of August 2. 

**What it does**:
- Starts simulation from **August 2, 2020** (mid-pandemic)
- Uses **pre-calculated initial conditions** based on real data
- Runs **55 stochastic simulations**
- Focuses on **short-term hospitalization forecasting**
- Demonstrates forecasting from known state vs. Patient Zero

**Key Features**:
- **Typed structs** for optimal performance
- **R₀ = 1.23** (targeted reproduction number)
- **Known initial conditions**: S=99.42%, I=0.12%, H=0.01%, R=0.44%, D=0.01%
- **8-month simulation period** (240 days, August 2, 2020 - April 2021)

**Outputs**:
- `SIHRS_Carson_City_Hospitalization_August2_bandwidth.png` - August 2 start bandwidth
- `SIHRS_Carson_City_Hospitalization_August2_trajectories.png` - August 2 start trajectories

**Use Case**: Demonstrating forecasting accuracy when starting from known epidemic state

---

## Data Files

### Input Data
- `carson_city_combined.csv` - Combined cases and deaths data (filtered NYT county-level data for Carson City; columns: `date, county, state, fips, cases, deaths`)

### Output Data
- `carson_city_active_cases.csv` - Processed active case counts
- `carson_city_daily_deaths.csv` - Daily deaths and 7-day moving average

## How the CSVs are generated

- `carson_city_active_cases.csv` (from `extract_active_cases.jl`):
  - Active cases are computed as cumulative cases minus the value 14 days earlier: `active_cases[t] = max(cases[t] - cases[t-14], 0)`.
  - Includes columns: `date, cumulative_cases, cumulative_deaths, active_cases`.

- `carson_city_daily_deaths.csv` (from `extract_daily_deaths.jl`):
  - Daily deaths are computed by differencing cumulative deaths: `daily[t] = deaths[t] - deaths[t-1]` (first day equals cumulative).
  - Adds a 7-day moving average: `moving_avg_7day` for smoothing.
  - Used by plotting code as the "active deaths" reference via the 7-day moving average.
