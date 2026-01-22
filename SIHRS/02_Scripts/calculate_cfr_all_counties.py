#!/usr/bin/env python3

import pandas as pd
import numpy as np
from pathlib import Path
from tqdm import tqdm
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('cfr_calculation_final.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


class CFRCalculator:
    def __init__(self):
        self.nyt_data = None
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

    def calculate_cfr_for_county(self, fips, county_name, state):
        try:
            county_data = self.nyt_data[self.nyt_data['fips'] == fips].copy()
            if len(county_data) == 0:
                return None
            
            county_data = county_data.sort_values('date')
            final_row = county_data.iloc[-1]
            total_cases = final_row['cases']
            total_deaths = final_row['deaths']
            
            if total_cases <= 0:
                return None
            
            cfr = total_deaths / total_cases
            
            county_data['daily_cfr'] = county_data['deaths'] / county_data['cases'].replace(0, np.nan)
            valid_cfr = county_data['daily_cfr'].dropna()
            
            if len(valid_cfr) == 0:
                return None
            
            return {
                'county': county_name,
                'state': state,
                'fips': fips,
                'total_cases': int(total_cases),
                'total_deaths': int(total_deaths),
                'cfr': cfr,
                'cfr_percent': cfr * 100,
                'mean_daily_cfr': valid_cfr.mean(),
                'median_daily_cfr': valid_cfr.median(),
                'std_daily_cfr': valid_cfr.std(),
                'min_daily_cfr': valid_cfr.min(),
                'max_daily_cfr': valid_cfr.max(),
            }

        except Exception as e:
            logger.warning(f"Error processing {county_name}, {state}: {e}")
            return None

    def process_all_counties(self):
        logger.info("="*70)
        logger.info("Processing all counties for CFR...")
        logger.info("="*70)
        
        counties = self.nyt_data.groupby(['fips', 'county', 'state']).size().reset_index()
        counties = counties[counties['fips'] != '00000']
        
        logger.info(f"Found {len(counties)} unique counties")
        
        for _, row in tqdm(counties.iterrows(), total=len(counties), desc="Calculating CFR"):
            fips = row['fips']
            county_name = row['county']
            state = row['state']
            
            result = self.calculate_cfr_for_county(fips, county_name, state)
            if result is not None:
                self.results.append(result)
        
        logger.info(f"\nSuccessfully calculated CFR for {len(self.results)} counties")

    def save_results(self):
        output_dir = Path('../03_ProcessedData/CFR_analysis')
        output_dir.mkdir(parents=True, exist_ok=True)

        df = pd.DataFrame(self.results)
        df = df.sort_values('cfr', ascending=False)

        df.to_csv(output_dir / 'cfr_all_counties.csv', index=False)
        logger.info(f"Saved: {output_dir / 'cfr_all_counties.csv'}")

        state_summary = df.groupby('state').agg({
            'cfr': ['mean', 'median', 'std', 'count'],
            'total_cases': 'sum',
            'total_deaths': 'sum'
        }).round(6)
        state_summary.columns = ['_'.join(col).strip() for col in state_summary.columns]
        state_summary = state_summary.reset_index()
        state_summary['state_cfr'] = state_summary['total_deaths_sum'] / state_summary['total_cases_sum']
        state_summary = state_summary.sort_values('state_cfr', ascending=False)
        state_summary.to_csv(output_dir / 'cfr_by_state.csv', index=False)
        logger.info(f"Saved: {output_dir / 'cfr_by_state.csv'}")

        total_cases = df['total_cases'].sum()
        total_deaths = df['total_deaths'].sum()
        overall_cfr = total_deaths / total_cases
        
        summary = pd.DataFrame({
            'metric': [
                'Overall CFR (U.S.)',
                'Mean County CFR',
                'Median County CFR',
                'Standard Deviation',
                '5th Percentile',
                '25th Percentile',
                '50th Percentile',
                '75th Percentile',
                '95th Percentile',
                'Minimum',
                'Maximum',
                'Total Counties',
                'Total U.S. Cases',
                'Total U.S. Deaths'
            ],
            'value': [
                overall_cfr,
                df['cfr'].mean(),
                df['cfr'].median(),
                df['cfr'].std(),
                df['cfr'].quantile(0.05),
                df['cfr'].quantile(0.25),
                df['cfr'].quantile(0.50),
                df['cfr'].quantile(0.75),
                df['cfr'].quantile(0.95),
                df['cfr'].min(),
                df['cfr'].max(),
                len(df),
                total_cases,
                total_deaths
            ],
            'percentage': [
                f"{overall_cfr*100:.4f}%",
                f"{df['cfr'].mean()*100:.4f}%",
                f"{df['cfr'].median()*100:.4f}%",
                f"{df['cfr'].std()*100:.4f}%",
                f"{df['cfr'].quantile(0.05)*100:.4f}%",
                f"{df['cfr'].quantile(0.25)*100:.4f}%",
                f"{df['cfr'].quantile(0.50)*100:.4f}%",
                f"{df['cfr'].quantile(0.75)*100:.4f}%",
                f"{df['cfr'].quantile(0.95)*100:.4f}%",
                f"{df['cfr'].min()*100:.4f}%",
                f"{df['cfr'].max()*100:.4f}%",
                '',
                '',
                ''
            ]
        })
        summary.to_csv(output_dir / 'cfr_summary_statistics.csv', index=False)
        logger.info(f"Saved: {output_dir / 'cfr_summary_statistics.csv'}")

        return df

    def print_statistics(self, df):
        total_cases = df['total_cases'].sum()
        total_deaths = df['total_deaths'].sum()
        overall_cfr = total_deaths / total_cases
        
        logger.info("\n" + "="*70)
        logger.info("CASE FATALITY RATE (CFR) - ALL U.S. COUNTIES")
        logger.info("Date Range: March 2020 - December 31, 2021")
        logger.info("="*70)

        logger.info(f"\nTOTAL U.S. STATISTICS:")
        logger.info(f"  Total Cases:  {total_cases:,}")
        logger.info(f"  Total Deaths: {total_deaths:,}")
        logger.info(f"  Overall CFR:  {overall_cfr:.6f} ({overall_cfr*100:.4f}%)")

        logger.info(f"\nCOUNTY-LEVEL CFR STATISTICS:")
        logger.info(f"  Total counties: {len(df)}")
        logger.info(f"  Mean CFR:       {df['cfr'].mean():.6f} ({df['cfr'].mean()*100:.4f}%)")
        logger.info(f"  Median CFR:     {df['cfr'].median():.6f} ({df['cfr'].median()*100:.4f}%)")
        logger.info(f"  Std Dev:        {df['cfr'].std():.6f}")
        logger.info(f"  25th %ile:      {df['cfr'].quantile(0.25):.6f} ({df['cfr'].quantile(0.25)*100:.4f}%)")
        logger.info(f"  75th %ile:      {df['cfr'].quantile(0.75):.6f} ({df['cfr'].quantile(0.75)*100:.4f}%)")

        logger.info(f"\nCARSON CITY VALIDATION:")
        carson = df[df['fips'] == '32510']
        if len(carson) > 0:
            carson_cfr = carson.iloc[0]['cfr']
            carson_cases = carson.iloc[0]['total_cases']
            carson_deaths = carson.iloc[0]['total_deaths']
            logger.info(f"  Cases:  {carson_cases:,}")
            logger.info(f"  Deaths: {carson_deaths:,}")
            logger.info(f"  CFR:    {carson_cfr:.6f} ({carson_cfr*100:.4f}%)")
            
            percentile = (df['cfr'] < carson_cfr).mean() * 100
            logger.info(f"  Percentile: {percentile:.1f}th")
        else:
            logger.warning("  Carson City not found!")

        logger.info("="*70)


def main():
    logger.info("="*70)
    logger.info("CASE FATALITY RATE (CFR) CALCULATION")
    logger.info("CFR = Total Deaths / Total Cases")
    logger.info("Date Range: March 2020 - December 31, 2021")
    logger.info("="*70)

    calc = CFRCalculator()
    calc.load_nyt_data('../01_RawData/nyt_covid_data.csv')
    calc.process_all_counties()
    df = calc.save_results()
    calc.print_statistics(df)

    logger.info("\n✓ Complete!")


if __name__ == '__main__':
    main()
