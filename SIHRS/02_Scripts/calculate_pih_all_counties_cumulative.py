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
        logging.FileHandler('pih_cumulative_calculation.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


class PIHCalculatorCumulative:
    def __init__(self, lag_days=14):
        self.lag_days = lag_days
        self.nyt_data = None
        self.new_infections_cache = {}
        self.results = []

    def load_nyt_data(self, filename='../01_RawData/nyt_covid_data.csv'):
        logger.info(f"Loading NYT data from {filename}...")
        self.nyt_data = pd.read_csv(filename, parse_dates=['date'])
        self.nyt_data['fips'] = self.nyt_data['fips'].fillna(0).astype(int).astype(str).str.zfill(5)
        self.nyt_data = self.nyt_data.sort_values(['fips', 'date'])
        
        start_date = pd.to_datetime('2020-03-01')
        end_date = pd.to_datetime('2021-12-31')
        self.nyt_data = self.nyt_data[
            (self.nyt_data['date'] >= start_date) & 
            (self.nyt_data['date'] <= end_date)
        ]
        
        logger.info(f"  Loaded {len(self.nyt_data):,} rows")
        logger.info(f"  Date range: {self.nyt_data['date'].min().date()} to {self.nyt_data['date'].max().date()}")
        logger.info(f"  Unique counties: {self.nyt_data['fips'].nunique()}")

    def calculate_new_infections(self, fips):
        if fips in self.new_infections_cache:
            return self.new_infections_cache[fips]

        county_data = self.nyt_data[self.nyt_data['fips'] == fips].copy()
        if len(county_data) == 0:
            return None

        county_data = county_data.sort_values('date')
        county_data['new_infections'] = county_data['cases'].diff().fillna(0)
        county_data['new_infections'] = county_data['new_infections'].clip(lower=0)

        result = county_data[['date', 'new_infections', 'cases']].copy()
        self.new_infections_cache[fips] = result
        return result

    def calculate_pih_cumulative_for_county(self, hosp_file):
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

            hosp_df = hosp_df.groupby('collection_week', as_index=False).agg({
                hosp_col: 'sum'
            })
            hosp_df = hosp_df.sort_values('collection_week')
            hosp_df.rename(columns={'collection_week': 'date', hosp_col: 'hospitalized'}, inplace=True)

            infections_df = self.calculate_new_infections(fips)
            if infections_df is None or len(infections_df) == 0:
                return None

            start_date = pd.to_datetime('2020-03-01')
            end_date = pd.to_datetime('2021-12-31')
            
            hosp_df = hosp_df[
                (hosp_df['date'] >= start_date) & 
                (hosp_df['date'] <= end_date)
            ]
            
            infections_df = infections_df[
                (infections_df['date'] >= start_date) & 
                (infections_df['date'] <= end_date)
            ]

            if len(hosp_df) == 0 or len(infections_df) == 0:
                return None

            total_hospitalizations = hosp_df['hospitalized'].sum()
            min_hosp_date = hosp_df['date'].min()
            max_hosp_date = hosp_df['date'].max()
            
            infections_start = min_hosp_date - timedelta(days=self.lag_days + 7)
            infections_end = max_hosp_date - timedelta(days=self.lag_days - 7)
            
            infections_for_hosp = infections_df[
                (infections_df['date'] >= infections_start) &
                (infections_df['date'] <= infections_end)
            ]
            
            total_new_infections = infections_for_hosp['new_infections'].sum()

            if total_new_infections > 0:
                p_IH = total_hospitalizations / total_new_infections
            else:
                return None

            daily_ratios = []
            for i in range(2, len(hosp_df)):
                date_t = hosp_df.iloc[i]['date']
                hospitalized_t = hosp_df.iloc[i]['hospitalized']
                date_t_minus_lag = hosp_df.iloc[i-2]['date']
                
                matching = infections_df[
                    (infections_df['date'] >= date_t_minus_lag - timedelta(days=2)) &
                    (infections_df['date'] <= date_t_minus_lag + timedelta(days=2))
                ]
                
                if len(matching) > 0:
                    infections_window = matching['new_infections'].sum()
                    if infections_window > 0:
                        ratio = hospitalized_t / infections_window
                        daily_ratios.append(ratio)

            if len(daily_ratios) > 0:
                daily_array = np.array(daily_ratios)
                return {
                    'county': county_name,
                    'state': state,
                    'fips': fips,
                    'p_IH_cumulative': p_IH,
                    'total_hospitalizations': total_hospitalizations,
                    'total_new_infections': total_new_infections,
                    'n_daily_calculations': len(daily_ratios),
                    'mean_pih_daily': np.mean(daily_array) if len(daily_array) > 0 else None,
                    'median_pih_daily': np.median(daily_array) if len(daily_array) > 0 else None,
                    'std_pih_daily': np.std(daily_array) if len(daily_array) > 0 else None,
                }
            else:
                return {
                    'county': county_name,
                    'state': state,
                    'fips': fips,
                    'p_IH_cumulative': p_IH,
                    'total_hospitalizations': total_hospitalizations,
                    'total_new_infections': total_new_infections,
                    'n_daily_calculations': 0,
                    'mean_pih_daily': None,
                    'median_pih_daily': None,
                    'std_pih_daily': None,
                }

        except Exception as e:
            logger.warning(f"Error processing {hosp_file.name}: {e}")
            import traceback
            logger.debug(traceback.format_exc())
            return None

    def process_all_counties(self):
        logger.info("="*70)
        logger.info("Processing all counties using CUMULATIVE RATIO method...")
        logger.info(f"Lag period: {self.lag_days} days")
        logger.info("="*70)

        hosp_files = list(Path('../03_ProcessedData/hospitalization_data/by_state').rglob('*.csv'))
        logger.info(f"Found {len(hosp_files)} hospitalization files")

        for hosp_file in tqdm(hosp_files, desc="Calculating P(IH) cumulative"):
            result = self.calculate_pih_cumulative_for_county(hosp_file)
            if result is not None:
                self.results.append(result)

        logger.info(f"\nSuccessfully calculated P(IH) for {len(self.results)} counties")

    def save_results(self):
        output_dir = Path('../03_ProcessedData/p_IH_analysis')
        output_dir.mkdir(parents=True, exist_ok=True)

        df = pd.DataFrame(self.results)
        df = df.sort_values('p_IH_cumulative', ascending=False)

        output_file = output_dir / f'pih_all_counties_cumulative_{self.lag_days}day.csv'
        df.to_csv(output_file, index=False)
        logger.info(f"Saved: {output_file}")

        summary_data = []
        valid_pih = df[df['p_IH_cumulative'].notna()]['p_IH_cumulative']
        valid_daily = df[df['mean_pih_daily'].notna()]['mean_pih_daily']
        
        summary_data.append({
            'metric': 'Method',
            'cumulative': 'Cumulative Ratio',
            'daily_mean': 'Daily Mean (for comparison)'
        })
        
        summary_data.append({
            'metric': f'p_IH ({self.lag_days}-day lag)',
            'cumulative': valid_pih.mean() if len(valid_pih) > 0 else None,
            'daily_mean': valid_daily.mean() if len(valid_daily) > 0 else None
        })
        
        summary_data.append({
            'metric': 'Median',
            'cumulative': valid_pih.median() if len(valid_pih) > 0 else None,
            'daily_mean': valid_daily.median() if len(valid_daily) > 0 else None
        })
        
        summary_data.append({
            'metric': 'Standard Deviation',
            'cumulative': valid_pih.std() if len(valid_pih) > 0 else None,
            'daily_mean': valid_daily.std() if len(valid_daily) > 0 else None
        })
        
        summary_data.append({
            'metric': '5th Percentile',
            'cumulative': valid_pih.quantile(0.05) if len(valid_pih) > 0 else None,
            'daily_mean': valid_daily.quantile(0.05) if len(valid_daily) > 0 else None
        })
        
        summary_data.append({
            'metric': '25th Percentile',
            'cumulative': valid_pih.quantile(0.25) if len(valid_pih) > 0 else None,
            'daily_mean': valid_daily.quantile(0.25) if len(valid_daily) > 0 else None
        })
        
        summary_data.append({
            'metric': '75th Percentile',
            'cumulative': valid_pih.quantile(0.75) if len(valid_pih) > 0 else None,
            'daily_mean': valid_daily.quantile(0.75) if len(valid_daily) > 0 else None
        })
        
        summary_data.append({
            'metric': '95th Percentile',
            'cumulative': valid_pih.quantile(0.95) if len(valid_pih) > 0 else None,
            'daily_mean': valid_daily.quantile(0.95) if len(valid_daily) > 0 else None
        })
        
        summary_data.append({
            'metric': 'Min',
            'cumulative': valid_pih.min() if len(valid_pih) > 0 else None,
            'daily_mean': valid_daily.min() if len(valid_daily) > 0 else None
        })
        
        summary_data.append({
            'metric': 'Max',
            'cumulative': valid_pih.max() if len(valid_pih) > 0 else None,
            'daily_mean': valid_daily.max() if len(valid_daily) > 0 else None
        })
        
        summary_data.append({
            'metric': 'Total Counties',
            'cumulative': len(valid_pih),
            'daily_mean': len(valid_daily)
        })

        summary_df = pd.DataFrame(summary_data)
        summary_file = output_dir / f'pih_summary_cumulative_{self.lag_days}day.csv'
        summary_df.to_csv(summary_file, index=False)
        logger.info(f"Saved: {summary_file}")

        logger.info("\n" + "="*70)
        logger.info("SUMMARY COMPARISON")
        logger.info("="*70)
        if len(valid_pih) > 0:
            logger.info(f"Cumulative Ratio Method:")
            logger.info(f"  Mean p_IH: {valid_pih.mean():.6f} ({valid_pih.mean()*100:.2f}%)")
            logger.info(f"  Median p_IH: {valid_pih.median():.6f} ({valid_pih.median()*100:.2f}%)")
            logger.info(f"  Std Dev: {valid_pih.std():.6f} ({valid_pih.std()*100:.2f}%)")
        if len(valid_daily) > 0:
            logger.info(f"\nDaily Mean Method (for comparison):")
            logger.info(f"  Mean p_IH: {valid_daily.mean():.6f} ({valid_daily.mean()*100:.2f}%)")
            logger.info(f"  Median p_IH: {valid_daily.median():.6f} ({valid_daily.median()*100:.2f}%)")
            logger.info(f"  Std Dev: {valid_daily.std():.6f} ({valid_daily.std()*100:.2f}%)")


def main():
    logger.info("="*70)
    logger.info("p_IH Calculation: Cumulative Ratio Method")
    logger.info("="*70)
    
    calculator_14 = PIHCalculatorCumulative(lag_days=14)
    calculator_14.load_nyt_data()
    calculator_14.process_all_counties()
    calculator_14.save_results()
    
    logger.info("\n" + "="*70)
    logger.info("Also calculating with 7-day lag for comparison...")
    logger.info("="*70)
    
    calculator_7 = PIHCalculatorCumulative(lag_days=7)
    calculator_7.nyt_data = calculator_14.nyt_data
    calculator_7.new_infections_cache = calculator_14.new_infections_cache
    calculator_7.process_all_counties()
    calculator_7.save_results()
    
    logger.info("\n" + "="*70)
    logger.info("Calculation complete!")
    logger.info("="*70)


if __name__ == '__main__':
    main()

