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
        logging.FileHandler('pih_calculation_final.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


class PIHCalculator:
    def __init__(self, recovery_days=14):
        self.recovery_days = recovery_days
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
        
        logger.info(f"  Loaded {len(self.nyt_data):,} rows")
        logger.info(f"  Date range: {self.nyt_data['date'].min().date()} to {self.nyt_data['date'].max().date()}")
        logger.info(f"  Unique counties: {self.nyt_data['fips'].nunique()}")

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
            if len(hosp_df) < 3:
                return None

            fips = hosp_file.stem.split('_')[-1]
            state = hosp_file.parent.name
            county_name = hosp_file.stem.rsplit('_', 1)[0].replace('_', ' ').title()

            hosp_col = 'total_adult_patients_hospitalized_confirmed_and_suspected_covid_7_day_avg'
            if hosp_col not in hosp_df.columns:
                return None

            hosp_df = hosp_df.groupby('collection_week', as_index=False).agg({
                hosp_col: 'sum'
            })
            hosp_df = hosp_df.sort_values('collection_week')
            hosp_df.rename(columns={'collection_week': 'date', hosp_col: 'hospitalized'}, inplace=True)

            active_df = self.calculate_active_cases(fips)
            if active_df is None or len(active_df) == 0:
                return None

            pih_values = []
            for i in range(2, len(hosp_df)):
                date_t = hosp_df.iloc[i]['date']
                hospitalized_t = hosp_df.iloc[i]['hospitalized']
                date_t_minus_14 = hosp_df.iloc[i-2]['date']
                date_diff = (date_t - date_t_minus_14).days

                if 10 <= date_diff <= 15:
                    matching = active_df[
                        (active_df['date'] >= date_t_minus_14 - timedelta(days=1)) &
                        (active_df['date'] <= date_t_minus_14 + timedelta(days=1))
                    ]

                    if len(matching) > 0:
                        active_cases_t_minus_14 = matching.iloc[0]['active_cases']
                        if active_cases_t_minus_14 > 0:
                            pih = hospitalized_t / active_cases_t_minus_14
                            pih_values.append(pih)

            if len(pih_values) == 0:
                return None

            pih_array = np.array(pih_values)
            return {
                'county': county_name,
                'state': state,
                'fips': fips,
                'n_calculations': len(pih_values),
                'mean_pih': np.mean(pih_array),
                'median_pih': np.median(pih_array),
                'std_pih': np.std(pih_array),
                'min_pih': np.min(pih_array),
                'max_pih': np.max(pih_array),
                'percentile_25': np.percentile(pih_array, 25),
                'percentile_50': np.percentile(pih_array, 50),
                'percentile_75': np.percentile(pih_array, 75),
                'percentile_90': np.percentile(pih_array, 90),
                'percentile_95': np.percentile(pih_array, 95),
            }

        except Exception as e:
            logger.warning(f"Error processing {hosp_file.name}: {e}")
            return None

    def process_all_counties(self):
        logger.info("="*70)
        logger.info("Processing all counties...")
        logger.info("="*70)

        hosp_files = list(Path('../03_ProcessedData/hospitalization_data/by_state').rglob('*.csv'))
        logger.info(f"Found {len(hosp_files)} hospitalization files")

        for hosp_file in tqdm(hosp_files, desc="Calculating P(IH)"):
            result = self.calculate_pih_for_county(hosp_file)
            if result is not None:
                self.results.append(result)

        logger.info(f"\nSuccessfully calculated P(IH) for {len(self.results)} counties")

    def save_results(self):
        output_dir = Path('../03_ProcessedData/p_IH_analysis')
        output_dir.mkdir(parents=True, exist_ok=True)

        df = pd.DataFrame(self.results)
        df = df.sort_values('mean_pih', ascending=False)

        df.to_csv(output_dir / 'pih_all_counties_14day.csv', index=False)
        logger.info(f"Saved: {output_dir / 'pih_all_counties_14day.csv'}")

        summary = pd.DataFrame({
            'metric': [
                'Overall Mean p_IH',
                'Overall Median p_IH',
                'Standard Deviation',
                '5th Percentile',
                '25th Percentile',
                '50th Percentile',
                '75th Percentile',
                '95th Percentile',
                'Minimum',
                'Maximum',
                'Total Counties',
                'Total Calculations'
            ],
            'value': [
                df['mean_pih'].mean(),
                df['median_pih'].median(),
                df['mean_pih'].std(),
                df['mean_pih'].quantile(0.05),
                df['mean_pih'].quantile(0.25),
                df['mean_pih'].quantile(0.50),
                df['mean_pih'].quantile(0.75),
                df['mean_pih'].quantile(0.95),
                df['mean_pih'].min(),
                df['mean_pih'].max(),
                len(df),
                df['n_calculations'].sum()
            ],
            'percentage': [
                f"{df['mean_pih'].mean()*100:.4f}%",
                f"{df['median_pih'].median()*100:.4f}%",
                f"{df['mean_pih'].std()*100:.4f}%",
                f"{df['mean_pih'].quantile(0.05)*100:.4f}%",
                f"{df['mean_pih'].quantile(0.25)*100:.4f}%",
                f"{df['mean_pih'].quantile(0.50)*100:.4f}%",
                f"{df['mean_pih'].quantile(0.75)*100:.4f}%",
                f"{df['mean_pih'].quantile(0.95)*100:.4f}%",
                f"{df['mean_pih'].min()*100:.4f}%",
                f"{df['mean_pih'].max()*100:.4f}%",
                '',
                ''
            ]
        })
        summary.to_csv(output_dir / 'pih_summary_statistics_14day.csv', index=False)
        logger.info(f"Saved: {output_dir / 'pih_summary_statistics_14day.csv'}")

        return df

    def print_statistics(self, df):
        logger.info("\n" + "="*70)
        logger.info("P(IH) STATISTICS - ALL U.S. COUNTIES")
        logger.info("="*70)

        logger.info(f"\nTotal counties: {len(df)}")
        logger.info(f"Total calculations: {df['n_calculations'].sum():,}")

        logger.info(f"\nMEAN P(IH) ACROSS ALL COUNTIES:")
        logger.info(f"  Mean:             {df['mean_pih'].mean():.4f} ({df['mean_pih'].mean()*100:.2f}%)")
        logger.info(f"  Median:           {df['mean_pih'].median():.4f} ({df['mean_pih'].median()*100:.2f}%)")
        logger.info(f"  Std Dev:          {df['mean_pih'].std():.4f}")
        logger.info(f"  25th percentile:  {df['mean_pih'].quantile(0.25):.4f} ({df['mean_pih'].quantile(0.25)*100:.2f}%)")
        logger.info(f"  75th percentile:  {df['mean_pih'].quantile(0.75):.4f} ({df['mean_pih'].quantile(0.75)*100:.2f}%)")

        logger.info(f"\nCARSON CITY VALIDATION:")
        carson = df[df['fips'] == '32510']
        if len(carson) > 0:
            carson_pih = carson.iloc[0]['mean_pih']
            logger.info(f"  Calculated: {carson_pih:.4f} ({carson_pih*100:.2f}%)")
            logger.info(f"  Expected:   0.1060 (10.60%)")
            logger.info(f"  Difference: {abs(carson_pih - 0.1060):.4f}")
            
            percentile = (df['mean_pih'] < carson_pih).mean() * 100
            logger.info(f"  Percentile: {percentile:.1f}th")
        else:
            logger.warning("  Carson City not found in results!")

        logger.info("="*70)


def main():
    logger.info("="*70)
    logger.info("P(IH) CALCULATION - CARSON CITY METHODOLOGY (14-DAY LAG)")
    logger.info("="*70)

    calc = PIHCalculator(recovery_days=14)
    calc.load_nyt_data('../01_RawData/nyt_covid_data.csv')
    calc.process_all_counties()
    df = calc.save_results()
    calc.print_statistics(df)

    logger.info("\n✓ Complete!")


if __name__ == '__main__':
    main()
