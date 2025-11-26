# SIHRS Model: COVID-19 Epidemiological Parameter Analysis

## Project Overview

This repository contains the complete analysis of COVID-19 epidemiological parameters for U.S. counties (March 2020 - December 31, 2021), including:
- **p_IH**: Probability of hospitalization given infection
- **p_ID**: Probability of death given infection
- **CFR**: Case Fatality Rate

The analysis covers **2,400+ counties** with comprehensive data extraction, processing, and statistical analysis.

---

## 📁 Directory Structure

```
SIHRS/
├── 01_RawData/              # Original data sources (779 MB - not in repo)
│   ├── covid_hospital_report.csv    # HealthData.gov (679 MB) - download required
│   └── nyt_covid_data.csv           # NYT COVID-19 data (100 MB) - download required
│
├── 02_Scripts/              # Analysis scripts (Python)
│   ├── extract_hospitalization_data.py    # Extract county-level hospitalization data
│   ├── calculate_pih_all_counties.py      # Calculate p_IH for all counties
│   ├── calculate_pid_all_counties.py      # Calculate p_ID for all counties
│   ├── calculate_cfr_all_counties.py      # Calculate CFR for all counties
│   ├── Script_Methodology.tex/pdf         # Implementation documentation
│   └── requirements.txt                   # Python dependencies
│
├── 03_ProcessedData/        # All analysis results
│   ├── hospitalization_data/   # 2,446 county hospitalization CSVs
│   │   └── by_state/           # Organized by state (50 states + DC)
│   ├── p_IH_analysis/          # Hospitalization probability results
│   │   ├── pih_all_counties.csv           # Individual county results
│   │   ├── pih_by_state.csv               # State-level aggregates
│   │   └── pih_summary_statistics.csv     # National statistics
│   ├── p_ID_analysis/          # Death probability results
│   │   ├── pid_all_counties.csv           # Individual county results
│   │   ├── pid_by_state.csv               # State-level aggregates
│   │   └── pid_summary_statistics.csv     # National statistics
│   ├── CFR_analysis/           # Case fatality rate results
│   │   ├── cfr_all_counties.csv           # Individual county results
│   │   ├── cfr_by_state.csv               # State-level aggregates
│   │   └── cfr_summary_statistics.csv     # National statistics
│   └── R_0/                    # Basic reproduction number data
│
├── 04_CaseStudies/          # Detailed county analyses
│   ├── Carson_City_NV/         # Primary case study (FIPS 32510)
│   │   ├── Parameter_Justification.tex/pdf  # Parameter derivation & justification
│   │   ├── codes/              # Julia/MATLAB scripts & data
│   │   └── Results/            # Simulation results & plots
│   └── Washington_MS/          # Secondary case study
│       ├── codes/              # Julia/MATLAB scripts & data
│       └── Results/            # Simulation results & plots
│
├── 05_Documentation/        # Papers, reports, and documentation
│   ├── methodology_summary.tex/pdf          # Complete methodology
│   └── ALL_COUNTIES_EXTRACTION_REPORT.md   # Extraction report
│
├── 06_ModelCode/            # SIHRS model implementations
│   ├── Cwang/                  # Prof. Wang's dispersion analysis code
│   ├── Julia/                  # Julia implementations
│   └── Matlab/                 # MATLAB implementations
│
└── logs/                    # Execution logs
```

---

## 🚀 Quick Start

### 1. Download Raw Data

The raw data files are not included in this repository due to size (779 MB total). Download them from the original sources:

**Option A: Direct Downloads**
```bash
# Create data directory
mkdir -p 01_RawData

# Download HealthData.gov hospitalization data (679 MB)
# Visit: https://healthdata.gov/Hospital/COVID-19-Reported-Patient-Impact-and-Hospital-Capa/anag-cw7u
# Export as CSV and save to: 01_RawData/covid_hospital_report.csv

# Download NYT COVID-19 data (100 MB)
curl -o 01_RawData/nyt_covid_data.csv https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv
```

**Option B: Use Processed Data**
All processed results are already included in `03_ProcessedData/`. You can skip the data download if you only want to view results.

### 2. Install Dependencies
```bash
cd 02_Scripts
pip install -r requirements.txt
```

### 2. Run Analysis Scripts

**Navigate to scripts directory:**
```bash
cd 02_Scripts
```

**Extract hospitalization data (if needed):**
```bash
python extract_hospitalization_data.py \
  --source ../01_RawData/covid_hospital_report.csv \
  --counties county_list.csv
```

**Calculate epidemiological parameters:**
```bash
python calculate_pih_all_counties.py  # p_IH
python calculate_pid_all_counties.py  # p_ID
python calculate_cfr_all_counties.py  # CFR
```

### 3. View Results

Results are saved in `03_ProcessedData/`:
- `p_IH_analysis/pih_summary_statistics.csv`
- `p_ID_analysis/pid_summary_statistics.csv`
- `CFR_analysis/cfr_summary_statistics.csv`

**Note**: All scripts must be run from the `02_Scripts/` directory as they use relative paths to access data.

---

## 📊 Key Results

### U.S. County Statistics (March 2020 - December 2021)

| Statistic | p_IH | p_ID | CFR |
|-----------|------|------|-----|
| **Mean** | 3.99% | 0.32% | 1.77% |
| **Median** | 1.81% | 0.21% | 1.68% |
| **25th Percentile** | 0.26% | 0.13% | 1.25% |
| **75th Percentile** | 5.67% | 0.36% | 2.17% |
| **95th Percentile** | 13.24% | 0.85% | 3.16% |

### Carson City, NV (Primary Case Study)

| Parameter | Value | National Percentile |
|-----------|-------|---------------------|
| **p_IH** | 10.60% | 91.7th |
| **p_ID** | 0.17% | 39.7th |
| **CFR** | 1.84% | 60.7th |

---

## 📦 Data Sources

1. **HealthData.gov**: COVID-19 Reported Patient Impact and Hospital Capacity
   - Weekly hospitalization data by facility
   - URL: https://healthdata.gov/Hospital/COVID-19-Reported-Patient-Impact-and-Hospital-Capa/anag-cw7u
   - Date range: August 2020 - December 2021
   - File size: 679 MB (2.4 million records)

2. **The New York Times**: COVID-19 Data Repository
   - Daily cases and deaths by county
   - URL: https://github.com/nytimes/covid-19-data
   - Date range: March 2020 - December 2021
   - File size: 100 MB (2.2 million records)

---
