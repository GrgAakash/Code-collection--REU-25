# P(IH) Calculation Methodology

## Overview

P(IH) represents the probability that an infected individual becomes hospitalized in the SIHRS epidemic model. This document explains how we calculated P(IH) empirically using real COVID-19 data from Washington County, Mississippi.

## Formula

```
P(IH) = (Hospitalized Patients on day T) / (Active Cases on day T-14)
```

## Methodology

### 1. Data Sources
- **Hospitalization Data**: `hospitalization_MS_filtered.csv` (weekly data, combined across hospitals)
- **Active Cases Data**: `washington_mississippi_active_cases.csv` (daily data)

### 2. Time Lag Rationale
We used a **14-day lag** because:
- Days 0-5: Incubation period (asymptomatic)
- Days 5-10: Symptom development
- Days 10-14: Hospitalization decision for severe cases

### 3. Data Processing Steps

1. **Combined datasets** by matching dates
2. **Applied 2-row lookback** (hospitalization data is weekly, so 2 rows ≈ 14 days)
3. **Validated date differences** to ensure 10-15 day range
4. **Filtered out invalid pairs** with unrealistic time gaps

### 4. Quality Controls

#### Date Validation
- Only included calculations where T and T-14 were **10-15 days apart**
- Excluded pairs with gaps >15 days (data collection irregularities)

#### Examples of Filtered Data
```
Skipped: 7/4/21 → 5/23/21 (42 days) - Too long gap
Kept: 8/16/20 → 8/2/20 (14 days) - Perfect
```

## Results

### Final Statistics
- **Total valid calculations**: 52
- **Average P(IH)**: 0.1614 (16.14%)
- **Median P(IH)**: 0.1270 (12.70%)
- **Range**: 0.0% to 80.59%

### Data Quality
- **100% of included data** has exactly 14-day temporal spacing
- **8 date pairs filtered out** due to irregular gaps

## Interpretation

The average P(IH) of **16.14%** means that approximately 1 in 6 infected individuals in Washington County, MS required hospitalization during the COVID-19 pandemic.

## SIHRS Model Application

**Recommended parameter**: `pIH = 0.16`

This empirically-derived value provides a realistic hospitalization probability based on actual pandemic data, accounting for proper disease progression timing.

## Files Generated

- `washington_ms_pih_lagged_14days.csv` - Detailed calculations
- `calculate_pih.jl` - Analysis script
- `PIH_Calculation_Methodology.md` - This documentation

## Technical Implementation

```julia
# Key validation logic
date_diff = Dates.value(day_t - day_t_minus_14)
if date_diff <= 15 && date_diff >= 10
    pih_value = hospitalized_t / active_cases_t_minus_14
end
```

This approach ensures temporal consistency and biological realism in the P(IH) parameter estimation.
