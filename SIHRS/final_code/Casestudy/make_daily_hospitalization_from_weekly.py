"""
Construct daily hospitalization time series for Carson City from weekly 7-day
averages using step function method (same as all-counties analysis).

Method:
    - Expand each weekly average into 7 identical daily values (step function).
    - Assign the weekly 7-day average value to all 7 days of that week.

Output:
    - hospitalization_Carson_daily_interpolated.csv
        Columns: date, step_value (smooth_value is not generated)
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def load_weekly_data(csv_path: Path) -> pd.DataFrame:
    """
    Load the weekly hospitalization data.

    Expected columns (as in hospitalization_Carson_filtered_new.csv):
        - collection_week: week-end date (string, e.g. '11/22/20')
        - total_adult_patients_hospitalized_confirmed_and_suspected_covid_7_day_avg
    """
    df = pd.read_csv(csv_path)

    # Normalize column names just in case
    df.columns = [c.strip() for c in df.columns]

    if "collection_week" not in df.columns:
        raise ValueError("Expected a 'collection_week' column in the input CSV.")

    # Try to detect the 7-day average column
    avg_cols = [
        c
        for c in df.columns
        if "7_day_avg" in c.lower() or "7_day" in c.lower()
        or "7 day" in c.lower()
    ]
    if not avg_cols:
        raise ValueError(
            "Could not find a 7-day average hospitalization column. "
            "Expected something like '..._7_day_avg' in the header."
        )
    if len(avg_cols) > 1:
        # Use the first match but warn in case of ambiguity
        print(
            f"Warning: multiple candidate 7-day average columns found {avg_cols}. "
            f"Using '{avg_cols[0]}'."
        )
    avg_col = avg_cols[0]

    df["collection_week"] = pd.to_datetime(df["collection_week"])
    df = df.sort_values("collection_week").reset_index(drop=True)
    df = df[["collection_week", avg_col]].rename(
        columns={avg_col: "weekly_7day_avg"}
    )

    return df


def make_step_function_daily(df_weekly: pd.DataFrame) -> pd.DataFrame:
    """
    Expand each weekly observation into 7 daily observations, each equal to the
    weekly 7-day average (a "step function" over each week).

    We treat 'collection_week' as the week-end date and assign the same value
    to each of the 7 days ending on that date.
    """
    records = []
    for _, row in df_weekly.iterrows():
        week_end = row["collection_week"]
        value = row["weekly_7day_avg"]

        # 7 days ending on week_end: week_end-6, ..., week_end
        days = pd.date_range(week_end - pd.Timedelta(days=6), week_end, freq="D")
        for d in days:
            records.append({"date": d, "step_value": value})

    df_daily_step = pd.DataFrame.from_records(records)
    df_daily_step = df_daily_step.sort_values("date").reset_index(drop=True)
    return df_daily_step


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Disaggregate weekly 7-day average hospitalization data into "
            "daily step-function series (same as all-counties analysis)."
        )
    )
    parser.add_argument(
        "--input",
        type=str,
        default="hospitalization_Carson_filtered_new.csv",
        help="Path to the weekly hospitalization CSV (default: hospitalization_Carson_filtered_new.csv).",
    )
    parser.add_argument(
        "--output-csv",
        type=str,
        default="hospitalization_Carson_daily_interpolated.csv",
        help="Output CSV file for daily step series.",
    )

    args = parser.parse_args()

    input_path = Path(args.input).resolve()
    if not input_path.exists():
        raise SystemExit(f"Input file not found: {input_path}")

    df_weekly = load_weekly_data(input_path)
    df_daily = make_step_function_daily(df_weekly)

    # Add smooth_value column (same as step_value) for backward compatibility
    # with existing CSV files, though only step_value is used
    df_daily["smooth_value"] = df_daily["step_value"]

    out_csv_path = input_path.parent / args.output_csv
    df_daily.to_csv(out_csv_path, index=False)
    print(f"Daily step-function series written to: {out_csv_path}")


if __name__ == "__main__":
    main()


