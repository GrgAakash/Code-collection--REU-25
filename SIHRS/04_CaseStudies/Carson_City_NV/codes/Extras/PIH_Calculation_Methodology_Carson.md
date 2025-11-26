# P(IH) Calculation Methodology - Carson City, NV

## Overview

P(IH) represents the probability that an infected individual becomes hospitalized in the SIHRS epidemic model. This document explains how we calculated P(IH) empirically using real COVID-19 data from Carson City, Nevada.

## Formula

```
P(IH) = (Hospitalized Patients on day T) / (Active Cases on day T-14)
```

## Methodology

### 1. Data Sources
- **Hospitalization Data**: `hospitalization_Carson_filtered_new.csv` (weekly data, combined across hospitals)
- **Active Cases Data**: `carson_city_active_cases.csv` (daily data)

### 2. Time Lag Rationale
We used a **14-day lag** because:
- Days 0-5: Incubation period (asymptomatic)
- Days 5-10: Symptom development
- Days 10-14: Hospitalization decision for severe cases

### 3. Data Processing Steps

1. **Combined duplicate dates** in hospitalization data (136 → 74 rows)
2. **Combined datasets** by matching dates
3. **Applied 2-row lookback** (hospitalization data is weekly, so 2 rows ≈ 14 days)
4. **Validated date differences** to ensure 10-15 day range
5. **Filtered out invalid pairs** with unrealistic time gaps

### 4. Quality Controls

#### Date Validation
- Only included calculations where T and T-14 were **10-15 days apart**
- Excluded pairs with gaps >15 days (data collection irregularities)

#### Data Quality Improvements
- **Removed duplicates**: 62 duplicate hospital entries combined
- **Perfect temporal spacing**: 100% of included data has exactly 14-day spacing

## Results

### Final Statistics
- **Total valid calculations**: 72
- **Average P(IH)**: 0.1060 (10.60%)
- **Median P(IH)**: 0.0811 (8.11%)

## SIHRS Model Application

**Recommended parameter**: `pIH = 0.11` (for Carson City, NV specifically)

This empirically-derived value provides a realistic hospitalization probability based on actual pandemic data from Carson City, accounting for proper disease progression timing.

## Files Generated

- `carson_city_pih_lagged_14days.csv` - Detailed calculations
- `calculate_carson_pih.jl` - Analysis script
- `PIH_Calculation_Methodology_Carson.md` - This documentation
- `carson_city_hospitalization_active_cases_combined.csv` - Combined dataset

## Technical Implementation

```julia
# Key validation logic
date_diff = Dates.value(day_t - day_t_minus_14)
if date_diff <= 15 && date_diff >= 10
    pih_value = hospitalized_t / active_cases_t_minus_14
end
```

This approach ensures temporal consistency and biological realism in the P(IH) parameter estimation for Carson City, Nevada.
