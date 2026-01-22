# SIHRS Model Parameter Estimation

This repository contains Python and Matlab scripts for estimating parameters for the SIHRS (Susceptible-Infected-Hospitalized-Recovered-Susceptible) with death epidemiological model using COVID-19 data from the United States.

## Overview

This repository provides data-driven estimation of key epidemiological parameters for the SIHRS compartmental model, using county-level COVID-19 data from across the United States.

### Parameters Estimated

| Parameter | Description | Method |
|-----------|-------------|--------|
| **p_IH** | Probability an infected person is hospitalized | 7-day lagged ratio with moving average |
| **p_ID** | Probability an infected person dies without hospitalization | 20-day lagged ratio |
| **CFR** | Case fatality rate | Cumulative deaths / cumulative cases |

*Time lags are chosen based on epidemiological literature. See paper for details.*

### Key Features
- Processes **2,300+ US counties** with hospitalization data
- Accounts for **infection-to-hospitalization delay** (7-day lag)
- Applies **7-day moving average** to reduce daily noise
- Includes **Carson City, NV case study** for validation
- **Companion code** for our SIHRS paper — all mathematical formulations and model implementations support our paper.

## Repository Structure

```
SIHRS/
├── 01_RawData/              # Raw COVID-19 data files
├── 02_Scripts/              # Parameter estimation scripts
│   ├── calculate_pih_all_counties_7day_lag_ma.py
│   ├── calculate_pid_all_counties.py
│   ├── calculate_cfr_all_counties.py
│   └── Script_Methodology.tex   # Detailed methodology documentation
├── 03_ProcessedData/        # Processed results and statistics
├── final_code/              # Figure generation code for paper
│   └── Casestudy/           # Carson City validation scripts
├── CITATION.cff             # Citation metadata
├── LICENSE                  # MIT License
└── README.md
```

## Quick Start

See `02_Scripts/Script_Methodology.tex` for detailed methodology and implementation details.

## Citation

**If you use this code, please cite it:**

### BibTeX
```bibtex
@software{sihrs_parameter_estimation_2026,
  author = {Gurung, Aakash and Wagle, Shrinkhal and Carr, Amy and McCann, Caden and Kodatt, Keyton and Song, Yuanyuan and Shao, Yuanzhen and Wang, Chuntian},
  title = {SIHRS Model Parameter Estimation Code},
  year = {2026},
  url = {https://github.com/GrgAakash/Code-collection--REU-25/tree/main/SIHRS},
  version = {1.0.0},
  note = {Paper currently under review}
}
```

### Plain Text
```
Gurung, A., Wagle, S., Carr, A., McCann, C., Kodatt, K., Song, Y., Shao, Y., & Wang, C. (2026). 
SIHRS Model Parameter Estimation Code (Version 1.0.0) [Computer software]. 
https://github.com/GrgAakash/Code-collection--REU-25/tree/main/SIHRS
```

**Note:** The associated paper is currently under review. Once published, please cite the published paper instead.

## License

MIT License - see LICENSE file for details.

## Data Sources

- New York Times COVID-19 Data: https://github.com/nytimes/covid-19-data
- HealthData.gov COVID-19 Hospital Report: https://www.cdc.gov/nchs/covid19/nhcs/hospital-mortality-by-week.htm

