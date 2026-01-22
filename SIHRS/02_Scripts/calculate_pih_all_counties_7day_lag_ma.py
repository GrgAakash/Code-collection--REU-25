#!/usr/bin/env python3
"""
Calculate p_IH for all counties using a daily lagged-ratio + 7-day
moving-average method. 7 days lag is what we used for calculation of gamma so we use it here as well.
A weighted approach with different lag is possible but we didn't implement it here.

For each county:
  1. Build a DAILY hospitalization series H_t by expanding weekly
     7-day-average hospitalization data into 7 identical daily values
     (step function over each week).
  2. Use NYT case data to construct DAILY active cases I_t
  3. For each day t, compute the LAGGED ratio:

         r_t = H_t / I_{t-7}

     whenever I_{t-7} > 0.
  4. Apply a 7-day moving average to r_t.
  5. Take the arithmetic mean of the smoothed ratios over time to obtain
     a single p_IH estimate for that county.

Outputs:
  - ../03_ProcessedData/p_IH_analysis/pih_all_counties_7day_lag_ma.csv
  - ../03_ProcessedData/p_IH_analysis/pih_summary_statistics_7day_lag_ma.
  
Carson city has its own script which you can foundnd in the final_code/Casestudy directory.
"""

from __future__ import annotations

import sys
import logging
from dataclasses import dataclass
from datetime import timedelta
from pathlib import Path
from typing import Optional, Dict, List

import numpy as np
import pandas as pd
from tqdm import tqdm


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


@dataclass
class PIHCalculatorDailyLag7MA:
    """Calculate p_IH using daily H_t / I_{t-7} with a 7-day moving average."""

    recovery_days: int = 7    # for constructing active cases from NYT data (matches gamma = 0.1429)
    lag_days: int = 7         # infection-to-hospitalization lag
    ma_window: int = 7        # 7-day moving average over daily ratios

    def __post_init__(self) -> None:
        self.nyt_data: Optional[pd.DataFrame] = None
        self.active_cases_cache: Dict[str, pd.DataFrame] = {}
        self.results: List[Dict] = []

    # ------------------------------------------------------------------
    # NYT active-case construction
    # ------------------------------------------------------------------
    def load_nyt_data(self, filename: str = "../01_RawData/nyt_covid_data.csv") -> None:
        logger.info(f"Loading NYT data from {filename}...")
        df = pd.read_csv(filename, parse_dates=["date"])
        df["fips"] = df["fips"].fillna(0).astype(int).astype(str).str.zfill(5)
        df = df.sort_values(["fips", "date"])

        end_date = pd.to_datetime("2021-12-31")
        df = df[df["date"] <= end_date]

        self.nyt_data = df
        logger.info(
            f"  Loaded {len(df):,} rows, {df['fips'].nunique()} counties "
            f"(up to {end_date.date()})"
        )

    def calculate_active_cases(self, fips: str) -> Optional[pd.DataFrame]:
        """
        Return daily active cases I_t using exponential decay method.
        
        Uses exact discretization: I[t] = new_cases[t] + I[t-1] * exp(-gamma)
        where gamma = 1/recovery_days. This matches the continuous-time ODE solution
        and the Poisson clock assumption in the SIHRS model.
        """
        if fips in self.active_cases_cache:
            return self.active_cases_cache[fips]

        if self.nyt_data is None:
            raise RuntimeError("NYT data not loaded. Call load_nyt_data() first.")

        county_data = self.nyt_data[self.nyt_data["fips"] == fips].copy()
        if len(county_data) == 0:
            return None

        county_data = county_data.sort_values("date").reset_index(drop=True)
        
        # Calculate daily new cases
        county_data["new_cases"] = county_data["cases"].diff().fillna(county_data["cases"].iloc[0]).clip(lower=0)
        
        # Exponential decay factor: exp(-gamma) where gamma = 1/recovery_days
        # This is the exact probability of remaining active after 1 day (Poisson assumption)
        gamma = 1.0 / self.recovery_days
        decay_factor = np.exp(-gamma)
        
        # Recursive calculation: I[t] = new_cases[t] + I[t-1] * exp(-gamma)
        active_cases = []
        for i in range(len(county_data)):
            if i == 0:
                active_cases.append(county_data.iloc[i]["new_cases"])
            else:
                current_active = county_data.iloc[i]["new_cases"] + active_cases[i-1] * decay_factor
                active_cases.append(current_active)
        
        county_data["active_cases"] = active_cases

        result = county_data[["date", "active_cases"]].copy()
        self.active_cases_cache[fips] = result
        return result

    # ------------------------------------------------------------------
    # Hospitalization data → daily step series
    # ------------------------------------------------------------------
    @staticmethod
    def _build_daily_hosp_from_weekly(
        hosp_df: pd.DataFrame,
        weekly_col: str,
    ) -> pd.DataFrame:
        """
        Expand weekly 7-day-average hospitalization data into a daily step series.

        hosp_df must contain:
          - 'date' (week-end date; originally 'collection_week')
          - weekly_col: the 7-day-average column to use as H_t.
        """
        hosp_df = hosp_df.dropna(subset=["date", weekly_col]).copy()
        hosp_df = hosp_df.sort_values("date").reset_index(drop=True)

        records: List[Dict] = []
        for _, row in hosp_df.iterrows():
            week_end = row["date"]
            value = float(row[weekly_col])
            # Assign the same value to each of the 7 days ending on week_end
            days = pd.date_range(week_end - timedelta(days=6), week_end, freq="D")
            for d in days:
                records.append({"date": d, "H_t": value})

        daily = pd.DataFrame.from_records(records)
        daily = daily.sort_values("date").reset_index(drop=True)
        return daily

    # County-level p_IH calculation

    def calculate_pih_for_county(self, hosp_file: Path) -> Optional[Dict]:
        """
        Calculate p_IH for one county using:
          - daily H_t from weekly hospitalizations,
          - daily I_t from NYT cases,
          - lagged ratio H_t / I_{t-7},
          - 7-day moving average of that ratio,
          - arithmetic mean of the smoothed ratio series.
        """
        try:
            hosp_df = pd.read_csv(hosp_file, parse_dates=["collection_week"])
            if len(hosp_df) < 2:
                return None

            fips = hosp_file.stem.split("_")[-1]
            state = hosp_file.parent.name
            county_name = hosp_file.stem.rsplit("_", 1)[0].replace("_", " ").title()

            hosp_col = (
                "total_adult_patients_hospitalized_confirmed_and_suspected_covid_7_day_avg"
            )
            if hosp_col not in hosp_df.columns:
                return None

            hosp_weekly = (
                hosp_df.groupby("collection_week", as_index=False)
                .agg({hosp_col: "sum"})
                .rename(columns={"collection_week": "date"})
            )

            if len(hosp_weekly) < 2:
                return None

            daily_hosp = self._build_daily_hosp_from_weekly(
                hosp_weekly, weekly_col=hosp_col
            )

            active_df = self.calculate_active_cases(fips)
            if active_df is None or len(active_df) == 0:
                return None

            df = pd.merge(
                daily_hosp,
                active_df.rename(columns={"active_cases": "I_t"}),
                on="date",
                how="inner",
            ).sort_values("date")

            if len(df) == 0:
                return None

            df["I_t_minus_7"] = df["I_t"].shift(self.lag_days)
            df = df[df["I_t_minus_7"] > 0].copy()
            if len(df) < self.ma_window:
                return None

            df["ratio_lag7"] = df["H_t"] / df["I_t_minus_7"]
            df["ratio_lag7_ma"] = df["ratio_lag7"].rolling(
                window=self.ma_window,
                min_periods=self.ma_window,
            ).mean()

            df_valid = df[df["ratio_lag7_ma"].notna()].copy()
            if len(df_valid) == 0:
                return None

            ratios = df_valid["ratio_lag7_ma"].to_numpy()

            p25 = float(np.percentile(ratios, 25))
            p75 = float(np.percentile(ratios, 75))
            
            return {
                "county": county_name,
                "state": state,
                "fips": fips,
                "n_daily_points": len(df),
                "n_ma_points": len(df_valid),
                "mean_pih": float(np.mean(ratios)),
                "median_pih": float(np.median(ratios)),
                "std_pih": float(np.std(ratios)),
                "min_pih": float(np.min(ratios)),
                "max_pih": float(np.max(ratios)),
                "percentile_25": p25,
                "percentile_75": p75,
                "percentile_95": float(np.percentile(ratios, 95)),
                "pih_range_lower": p25,  # 25th percentile (lower bound of IQR)
                "pih_range_upper": p75,  # 75th percentile (upper bound of IQR)
                "pih_iqr": p75 - p25,    # Interquartile range
            }

        except Exception as e:
            logger.debug(f"Error processing {hosp_file}: {e}")
            return None

    # ------------------------------------------------------------------
    # Batch processing and outputs
    # ------------------------------------------------------------------
    def process_all_counties(
        self,
        limit: Optional[int] = None,
        test_county_fips: Optional[str] = None,
    ) -> None:
        """
        Process hospitalization files for all counties (or a subset).

        Args:
            limit: If set, only process the first N counties (for testing).
            test_county_fips: If set, only process this specific county
                              (e.g., '32510' for Carson City).
        """
        base_dir = Path("../03_ProcessedData/hospitalization_data/by_state")
        hosp_files = list(base_dir.rglob("*.csv"))
        logger.info(f"Found {len(hosp_files)} hospitalization files in {base_dir}")

        if test_county_fips:
            hosp_files = [f for f in hosp_files if test_county_fips in f.stem]
            logger.info(
                f"Testing mode: only county FIPS {test_county_fips} "
                f"({len(hosp_files)} file(s))"
            )
        elif limit is not None:
            hosp_files = hosp_files[:limit]
            logger.info(f"Test mode: processing first {limit} counties")

        for hosp_file in tqdm(hosp_files, desc="Calculating P(IH) [7-day lag + MA]"):
            result = self.calculate_pih_for_county(hosp_file)
            if result is not None:
                self.results.append(result)

        logger.info(
            f"Successfully calculated P(IH) for {len(self.results)} counties "
            f"using daily lagged ratio + 7-day MA."
        )

    def save_results(self) -> pd.DataFrame:
        """Save per-county results and a national summary table."""
        output_dir = Path("../03_ProcessedData/p_IH_analysis")
        output_dir.mkdir(parents=True, exist_ok=True)

        df = pd.DataFrame(self.results)
        if df.empty:
            logger.warning("No results to save.")
            return df

        df = df.sort_values("mean_pih", ascending=False)

        out_file = output_dir / "pih_all_counties_7day_lag_ma.csv"
        df.to_csv(out_file, index=False)
        logger.info(f"Saved per-county results: {out_file}")

        # Filter out outliers: counties with mean_pih > 1.0 (clearly erroneous)
        df_filtered = df[df["mean_pih"] <= 1.0].copy()
        n_outliers = len(df) - len(df_filtered)
        if n_outliers > 0:
            logger.info(
                f"Filtered out {n_outliers} counties with mean_pih > 1.0 for summary statistics"
            )

        summary = pd.DataFrame(
            {
                "metric": [
                    "Mean p_IH (7-day lag + MA)",
                    "Median p_IH",
                    "Std Dev",
                    "25th Pctl",
                    "75th Pctl",
                    "95th Pctl",
                    "Counties",
                    "Counties (filtered)",
                ],
                "value": [
                    float(df_filtered["mean_pih"].mean()),
                    float(df_filtered["mean_pih"].median()),
                    float(df_filtered["mean_pih"].std()),
                    float(df_filtered["mean_pih"].quantile(0.25)),
                    float(df_filtered["mean_pih"].quantile(0.75)),
                    float(df_filtered["mean_pih"].quantile(0.95)),
                    len(df),
                    len(df_filtered),
                ],
                "percentage": [
                    f"{df_filtered['mean_pih'].mean() * 100:.4f}%",
                    f"{df_filtered['mean_pih'].median() * 100:.4f}%",
                    f"{df_filtered['mean_pih'].std() * 100:.4f}%",
                    f"{df_filtered['mean_pih'].quantile(0.25) * 100:.4f}%",
                    f"{df_filtered['mean_pih'].quantile(0.75) * 100:.4f}%",
                    f"{df_filtered['mean_pih'].quantile(0.95) * 100:.4f}%",
                    "",
                    "",
                ],
            }
        )

        summary_file = output_dir / "pih_summary_statistics_7day_lag_ma.csv"
        summary.to_csv(summary_file, index=False)
        logger.info(f"Saved summary statistics: {summary_file}")

        return df

    @staticmethod
    def print_statistics(df: pd.DataFrame) -> None:
        """Print a short summary of the national distribution."""
        if df.empty:
            print("No results to summarize.")
            return

        # Filter out outliers: counties with mean_pih > 1.0
        df_filtered = df[df["mean_pih"] <= 1.0].copy()
        n_outliers = len(df) - len(df_filtered)

        print("\n" + "=" * 70)
        print("P(IH) STATISTICS - DAILY LAGGED RATIO + 7-DAY MOVING AVERAGE")
        print("=" * 70)
        print(f"Total counties: {len(df)}")
        if n_outliers > 0:
            print(f"Counties filtered (mean_pih > 1.0): {n_outliers}")
            print(f"Counties used for statistics: {len(df_filtered)}")
        print(
            f"\nMean p_IH:   {df_filtered['mean_pih'].mean():.6f} "
            f"({df_filtered['mean_pih'].mean() * 100:.4f}%)"
        )
        print(
            f"Median p_IH: {df_filtered['mean_pih'].median():.6f} "
            f"({df_filtered['mean_pih'].median() * 100:.4f}%)"
        )
        print(
            f"25th pctl:   {df_filtered['mean_pih'].quantile(0.25):.6f} "
            f"({df_filtered['mean_pih'].quantile(0.25) * 100:.4f}%)"
        )
        print(
            f"75th pctl:   {df_filtered['mean_pih'].quantile(0.75):.6f} "
            f"({df_filtered['mean_pih'].quantile(0.75) * 100:.4f}%)"
        )


def main() -> None:
    # Simple CLI parsing (to keep dependencies minimal)
    limit: Optional[int] = None
    test_county: Optional[str] = None

    if "--limit" in sys.argv:
        idx = sys.argv.index("--limit")
        if idx + 1 < len(sys.argv):
            try:
                limit = int(sys.argv[idx + 1])
            except ValueError:
                logger.warning("Ignoring invalid --limit value")

    if "--test-county" in sys.argv:
        idx = sys.argv.index("--test-county")
        if idx + 1 < len(sys.argv):
            test_county = sys.argv[idx + 1]

    calc = PIHCalculatorDailyLag7MA(recovery_days=7, lag_days=7, ma_window=7)
    calc.load_nyt_data("../01_RawData/nyt_covid_data.csv")
    calc.process_all_counties(limit=limit, test_county_fips=test_county)
    df = calc.save_results()
    calc.print_statistics(df)

    print("\n✓ Complete!")
    print("\nUsage examples:")
    print("  python calculate_pih_all_counties_7day_lag_ma.py")
    print("  python calculate_pih_all_counties_7day_lag_ma.py --limit 10")
    print("  python calculate_pih_all_counties_7day_lag_ma.py --test-county 32510")


if __name__ == "__main__":
    main()


