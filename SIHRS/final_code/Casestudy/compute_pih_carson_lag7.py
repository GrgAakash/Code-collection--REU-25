"""
Compute daily ratios H_t / I_{t-7} for Carson City and their 7-day rolling
averages, starting from a user-specified date (default: 2020-08-02).

Uses step function method for daily hospitalization disaggregation (same as
all-counties analysis).

Inputs (expected in the same directory):
  - hospitalization_Carson_daily_interpolated.csv
      columns: date, step_value (smooth_value is ignored)
  - carson_city_active_cases.csv
      columns: date, cumulative_cases, cumulative_deaths, active_cases

Outputs:
  - carson_city_pih_lag7_daily.csv
      columns include:
        date,
        H_t (step function),
        I_t,
        I_t_minus_7,
        H_over_I_lag7,
        H_over_I_lag7_7day_ma
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def load_hospitalizations(path: Path) -> pd.DataFrame:
    """
    Load daily hospitalization data using step function method.
    
    Only uses step_value (smooth_value is ignored for consistency with
    all-counties analysis).
    """
    df = pd.read_csv(path, parse_dates=["date"])
    df = df.sort_values("date").reset_index(drop=True)

    if "step_value" not in df.columns:
        raise ValueError(
            f"Missing required column 'step_value' in {path.name}"
        )

    # Use only step function (same as all-counties analysis)
    df = df[["date", "step_value"]]
    df.rename(columns={"step_value": "H_t"}, inplace=True)
    return df


def load_active_cases(path: Path, recovery_days: int = 7) -> pd.DataFrame:
    """
    Load active cases from CSV using exponential decay method.
    
    If 'cumulative_cases' exists, recalculate using exponential decay:
    I[t] = new_cases[t] + I[t-1] * exp(-gamma) where gamma = 1/recovery_days.
    This matches the continuous-time ODE solution and Poisson clock assumption.
    
    Otherwise, use 'active_cases' column if available.
    """
    import numpy as np
    
    df = pd.read_csv(path, parse_dates=["date"])
    df = df.sort_values("date").reset_index(drop=True)

    if "cumulative_cases" in df.columns:
        # Calculate daily new cases
        df["new_cases"] = df["cumulative_cases"].diff().fillna(df["cumulative_cases"].iloc[0]).clip(lower=0)
        
        # Exponential decay factor: exp(-gamma) where gamma = 1/recovery_days
        gamma = 1.0 / recovery_days
        decay_factor = np.exp(-gamma)
        
        # Recursive calculation: I[t] = new_cases[t] + I[t-1] * exp(-gamma)
        active_cases = []
        for i in range(len(df)):
            if i == 0:
                active_cases.append(df.iloc[i]["new_cases"])
            else:
                current_active = df.iloc[i]["new_cases"] + active_cases[i-1] * decay_factor
                active_cases.append(current_active)
        
        df["I_t"] = active_cases
        df = df[["date", "I_t"]]
    elif "active_cases" in df.columns:
        df = df[["date", "active_cases"]].rename(columns={"active_cases": "I_t"})
    else:
        raise ValueError(
            f"Expected either 'active_cases' or 'cumulative_cases' column in {path.name}."
        )

    return df


def compute_ratios(
    df_hosp: pd.DataFrame,
    df_active: pd.DataFrame,
    start_date: str,
) -> pd.DataFrame:
    """
    Merge hospitalizations and active cases, compute H_t / I_{t-7} and a 7-day
    rolling average of that ratio.
    
    Uses step function method for H_t (consistent with all-counties analysis).
    """
    df = pd.merge(df_hosp, df_active, on="date", how="inner")
    df = df.sort_values("date").reset_index(drop=True)

    # I_{t-7}: 7-day lag of active cases
    df["I_t_minus_7"] = df["I_t"].shift(7)

    # Ratio H_t / I_{t-7} (using step function method)
    df["H_over_I_lag7"] = df["H_t"] / df["I_t_minus_7"]

    # 7-day rolling average of the ratio (require full 7 days)
    df["H_over_I_lag7_7day_ma"] = (
        df["H_over_I_lag7"].rolling(window=7, min_periods=7).mean()
    )

    # Restrict to dates >= start_date and drop rows where I_{t-7} is missing
    start_ts = pd.to_datetime(start_date)
    df = df[df["date"] >= start_ts].copy()
    df = df[df["I_t_minus_7"].notna()].reset_index(drop=True)

    return df


def print_summary_statistics(df: pd.DataFrame, recovery_days: int) -> None:
    """Print summary statistics including 25th-75th percentile range."""
    import numpy as np
    
    # Filter to valid data points (where 7-day MA exists)
    df_valid = df[df["H_over_I_lag7_7day_ma"].notna()].copy()
    
    if len(df_valid) == 0:
        print("No valid data points for summary statistics.")
        return
    
    # Calculate statistics using step function method (consistent with all-counties)
    ratios_ma = df_valid["H_over_I_lag7_7day_ma"]
    
    print("\n" + "=" * 70)
    print(f"CARSON CITY p_IH SUMMARY STATISTICS")
    print(f"Method: Step function (same as all-counties analysis)")
    print(f"Recovery period: {recovery_days} days (gamma = {1.0/recovery_days:.4f})")
    print("=" * 70)
    print(f"Total valid data points: {len(df_valid)}")
    print(f"Date range: {df_valid['date'].min().date()} to {df_valid['date'].max().date()}")
    
    print("\n" + "-" * 70)
    print("STEP FUNCTION METHOD")
    print("-" * 70)
    print(f"  Mean:   {ratios_ma.mean():.6f} ({ratios_ma.mean()*100:.4f}%)")
    print(f"  Median: {ratios_ma.median():.6f} ({ratios_ma.median()*100:.4f}%)")
    p25 = ratios_ma.quantile(0.25)
    p75 = ratios_ma.quantile(0.75)
    print(f"  25th percentile: {p25:.6f} ({p25*100:.4f}%)")
    print(f"  75th percentile: {p75:.6f} ({p75*100:.4f}%)")
    print(f"  IQR range: [{p25:.6f}, {p75:.6f}]")
    print(f"  IQR width: {p75 - p25:.6f} ({(p75 - p25)*100:.4f} percentage points)")
    
    # Check where paper value (0.092) falls
    paper_value = 0.092
    percentile_rank = (ratios_ma <= paper_value).mean() * 100
    print(f"\n  Paper value p_IH = {paper_value:.6f} ({paper_value*100:.4f}%):")
    print(f"    - Percentile rank: {percentile_rank:.2f}th percentile")
    if p25 <= paper_value <= p75:
        print(f"    - Within IQR range: ✓")
    else:
        print(f"    - Within IQR range: ✗")


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Compute daily H_t / I_{t-7} for Carson City and their 7-day "
            "rolling averages."
        )
    )
    parser.add_argument(
        "--hosp-csv",
        type=str,
        default="hospitalization_Carson_daily_interpolated.csv",
        help=(
            "Daily hospitalization CSV with columns 'date', 'step_value'. "
            "Uses step function method (smooth_value is ignored)."
        ),
    )
    parser.add_argument(
        "--active-csv",
        type=str,
        default="carson_city_active_cases.csv",
        help="Daily active cases CSV with columns 'date', 'active_cases' (or 'cumulative_cases' to recalculate).",
    )
    parser.add_argument(
        "--recovery-days",
        type=int,
        default=7,
        help="Recovery window (days) for calculating active cases from cumulative cases (default: 7, consistent with gamma=0.1429).",
    )
    parser.add_argument(
        "--start-date",
        type=str,
        default="2020-08-02",
        help="First date t for which to keep H_t / I_{t-7} (default: 2020-08-02).",
    )
    parser.add_argument(
        "--output-csv",
        type=str,
        default="carson_city_pih_lag7_daily.csv",
        help="Output CSV filename (default: carson_city_pih_lag7_daily.csv).",
    )

    args = parser.parse_args()

    base_dir = Path(".").resolve()
    hosp_path = (base_dir / args.hosp_csv).resolve()
    active_path = (base_dir / args.active_csv).resolve()

    if not hosp_path.exists():
        raise SystemExit(f"Hospitalization CSV not found: {hosp_path}")
    if not active_path.exists():
        raise SystemExit(f"Active cases CSV not found: {active_path}")

    df_hosp = load_hospitalizations(hosp_path)
    df_active = load_active_cases(active_path, recovery_days=args.recovery_days)

    df_out = compute_ratios(df_hosp, df_active, start_date=args.start_date)

    out_path = base_dir / args.output_csv
    df_out.to_csv(out_path, index=False)
    print(f"Wrote H_t / I_(t-7) ratios and 7-day MAs to: {out_path}")
    
    print_summary_statistics(df_out, args.recovery_days)


if __name__ == "__main__":
    main()


