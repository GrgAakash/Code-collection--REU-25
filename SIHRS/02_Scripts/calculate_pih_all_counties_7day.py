#!/usr/bin/env python3

import pandas as pd
import numpy as np
from pathlib import Path
from datetime import timedelta
from tqdm import tqdm
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('pih_calculation_7day.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


class PIHCalculator7Day:
    def __init__(self, recovery_days=14, lag_rows=1):
        self.recovery_days = recovery_days
        self.lag_rows = lag_rows
        self.nyt_data = None
        self.active_cases_cache = {}
        self.results = []

    def load_nyt_data(self, filename='../01_RawData/nyt_covid_data.csv'):
        logger.info(f"Loading NYT data from {filename}...")
        self.nyt_data = pd.read_csv(filename, parse_dates=['date'])
        self.nyt_data['fips'] = self.nyt_data['fips'].fillna(0).astype(int).astype(str).str.zfill(5)
        self.nyt_data = self.nyt_data.sort_values(['fips', 'date'])
        end_date = pd.to_datetime('2021-12-31')
        self.nyt_data = self.nyt_data[self.nyt_data['date'] <= end_date]
        logger.info(f"  Loaded {len(self.nyt_data):,} rows, {self.nyt_data['fips'].nunique()} counties")

    def calculate_active_cases(self, fips):
        if fips in self.active_cases_cache:
            return self.active_cases_cache[fips]
        county_data = self.nyt_data[self.nyt_data['fips'] == fips].copy()
        if len(county_data) == 0:
            return None
        county_data = county_data.sort_values('date')
        county_data['active_cases'] = county_data['cases'] - county_data['cases'].shift(self.recovery_days).fillna(0)
        county_data['active_cases'] = county_data['active_cases'].clip(lower=0)
        result = county_data[['date', 'active_cases']].copy()
        self.active_cases_cache[fips] = result
        return result

    def calculate_pih_for_county(self, hosp_file):
        try:
            hosp_df = pd.read_csv(hosp_file, parse_dates=['collection_week'])
            if len(hosp_df) < 2:
                return None
            fips = hosp_file.stem.split('_')[-1]
            state = hosp_file.parent.name
            county_name = hosp_file.stem.rsplit('_', 1)[0].replace('_', ' ').title()
            hosp_col = 'total_adult_patients_hospitalized_confirmed_and_suspected_covid_7_day_avg'
            if hosp_col not in hosp_df.columns:
                return None
            hosp_df = hosp_df.groupby('collection_week', as_index=False).agg({hosp_col: 'sum'})
            hosp_df = hosp_df.sort_values('collection_week')
            hosp_df.rename(columns={'collection_week': 'date', hosp_col: 'hospitalized'}, inplace=True)
            active_df = self.calculate_active_cases(fips)
            if active_df is None or len(active_df) == 0:
                return None
            pih_values = []
            for i in range(self.lag_rows, len(hosp_df)):
                date_t = hosp_df.iloc[i]['date']
                hospitalized_t = hosp_df.iloc[i]['hospitalized']
                date_t_minus_7 = hosp_df.iloc[i - self.lag_rows]['date']
                date_diff = (date_t - date_t_minus_7).days
                if 5 <= date_diff <= 10:
                    matching = active_df[
                        (active_df['date'] >= date_t_minus_7 - timedelta(days=1)) &
                        (active_df['date'] <= date_t_minus_7 + timedelta(days=1))
                    ]
                    if len(matching) > 0:
                        active_cases_t_minus_7 = matching.iloc[0]['active_cases']
                        if active_cases_t_minus_7 > 0:
                            pih = hospitalized_t / active_cases_t_minus_7
                            pih_values.append(pih)
            if len(pih_values) == 0:
                return None
            pih_array = np.array(pih_values)
            return {
                'county': county_name, 'state': state, 'fips': fips,
                'n_calculations': len(pih_values),
                'mean_pih': np.mean(pih_array), 'median_pih': np.median(pih_array),
                'std_pih': np.std(pih_array), 'min_pih': np.min(pih_array), 'max_pih': np.max(pih_array),
                'percentile_25': np.percentile(pih_array, 25), 'percentile_75': np.percentile(pih_array, 75),
                'percentile_95': np.percentile(pih_array, 95),
            }
        except Exception as e:
            return None

    def process_all_counties(self):
        logger.info("Processing all counties (7-DAY LAG)...")
        hosp_files = list(Path('../03_ProcessedData/hospitalization_data/by_state').rglob('*.csv'))
        logger.info(f"Found {len(hosp_files)} hospitalization files")
        for hosp_file in tqdm(hosp_files, desc="Calculating P(IH) [7-day]"):
            result = self.calculate_pih_for_county(hosp_file)
            if result is not None:
                self.results.append(result)
        logger.info(f"Successfully calculated P(IH) for {len(self.results)} counties")

    def save_results(self):
        output_dir = Path('../03_ProcessedData/p_IH_analysis')
        output_dir.mkdir(parents=True, exist_ok=True)
        df = pd.DataFrame(self.results)
        df = df.sort_values('mean_pih', ascending=False)
        df.to_csv(output_dir / 'pih_all_counties_7day.csv', index=False)
        logger.info(f"Saved: {output_dir / 'pih_all_counties_7day.csv'}")
        
        summary = pd.DataFrame({
            'metric': ['Mean p_IH (7-day)', 'Median p_IH', 'Std Dev', '25th Pctl', '75th Pctl', '95th Pctl', 'Counties'],
            'value': [df['mean_pih'].mean(), df['median_pih'].median(), df['mean_pih'].std(),
                     df['mean_pih'].quantile(0.25), df['mean_pih'].quantile(0.75), 
                     df['mean_pih'].quantile(0.95), len(df)],
            'percentage': [f"{df['mean_pih'].mean()*100:.2f}%", f"{df['median_pih'].median()*100:.2f}%",
                          f"{df['mean_pih'].std()*100:.2f}%", f"{df['mean_pih'].quantile(0.25)*100:.2f}%",
                          f"{df['mean_pih'].quantile(0.75)*100:.2f}%", f"{df['mean_pih'].quantile(0.95)*100:.2f}%", '']
        })
        summary.to_csv(output_dir / 'pih_summary_statistics_7day.csv', index=False)
        return df

    def print_statistics(self, df):
        print("\n" + "="*70)
        print("P(IH) STATISTICS - 7-DAY LAG")
        print("="*70)
        print(f"Total counties: {len(df)}")
        print(f"\nMean p_IH:   {df['mean_pih'].mean():.4f} ({df['mean_pih'].mean()*100:.2f}%)")
        print(f"Median p_IH: {df['median_pih'].median():.4f} ({df['median_pih'].median()*100:.2f}%)")
        print(f"25th pctl:   {df['mean_pih'].quantile(0.25):.4f} ({df['mean_pih'].quantile(0.25)*100:.2f}%)")
        print(f"75th pctl:   {df['mean_pih'].quantile(0.75):.4f} ({df['mean_pih'].quantile(0.75)*100:.2f}%)")
        print("\n" + "="*70)
        print("COMPARISON: 7-DAY vs 14-DAY LAG")
        print("="*70)
        print(f"14-day mean p_IH: 0.0399 (3.99%)")
        print(f"7-day mean p_IH:  {df['mean_pih'].mean():.4f} ({df['mean_pih'].mean()*100:.2f}%)")
        print(f"Difference:       {(df['mean_pih'].mean() - 0.0399)*100:.2f} percentage points")
        carson = df[df['fips'] == '32510']
        if len(carson) > 0:
            print(f"\nCarson City: 7-day={carson.iloc[0]['mean_pih']:.4f} vs 14-day=0.1060")


def main():
    calc = PIHCalculator7Day(recovery_days=14, lag_rows=1)
    calc.load_nyt_data('../01_RawData/nyt_covid_data.csv')
    calc.process_all_counties()
    df = calc.save_results()
    calc.print_statistics(df)
    print("\n✓ Complete!")


if __name__ == '__main__':
    main()
