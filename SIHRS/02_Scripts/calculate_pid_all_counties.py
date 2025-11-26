#!/usr/bin/env python3
"""
Complete P(ID) Calculation for All U.S. Counties

CRITICAL METHODOLOGY (Matching Carson City):
- Formula: P(ID) = Daily_Deaths[T] / Active_Cases[T-20]
- 20-day lag accounts for longer time from infection to death
- Uses daily data (not weekly like p_IH)
- Includes zeros in calculation

Data Sources:
- NYT COVID-19 data: cases and deaths by county (daily)
- Active cases: cases[t] - cases[t-14]
"""

import pandas as pd
import numpy as np
from pathlib import Path
from datetime import timedelta
from tqdm import tqdm
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('pid_calculation_final.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


class PIDCalculator:
    """Calculate P(ID) for all U.S. counties using Carson City methodology."""

    def __init__(self, lag_days=20, recovery_days=14):
        """
        Initialize calculator.
        
        Args:
            lag_days: Days between infection and death (default 20)
            recovery_days: Days for active case calculation (default 14)
        """
        self.lag_days = lag_days
        self.recovery_days = recovery_days
        self.nyt_data = None
        self.results = []

    def load_nyt_data(self, filename='../01_RawData/nyt_covid_data.csv'):
        """Load NYT COVID-19 data."""
        logger.info(f"Loading NYT data from {filename}...")
        self.nyt_data = pd.read_csv(filename, parse_dates=['date'])
        
        # Clean FIPS codes
        self.nyt_data['fips'] = self.nyt_data['fips'].fillna(0).astype(int).astype(str).str.zfill(5)
        self.nyt_data = self.nyt_data.sort_values(['fips', 'date'])
        
        # Filter to end at 2021-12-31 to match paper timeline
        end_date = pd.to_datetime('2021-12-31')
        self.nyt_data = self.nyt_data[self.nyt_data['date'] <= end_date]
        
        logger.info(f"  Loaded {len(self.nyt_data):,} rows")
        logger.info(f"  Date range: {self.nyt_data['date'].min().date()} to {self.nyt_data['date'].max().date()}")
        logger.info(f"  Unique counties: {self.nyt_data['fips'].nunique()}")

    def calculate_pid_for_county(self, fips, county_name, state):
        """
        Calculate P(ID) for a single county.
        
        Args:
            fips: 5-digit FIPS code
            county_name: County name
            state: State abbreviation
            
        Returns:
            dict with P(ID) statistics or None if calculation fails
        """
        try:
            # Get county data
            county_data = self.nyt_data[self.nyt_data['fips'] == fips].copy()
            
            if len(county_data) < self.lag_days + 10:
                return None
            
            county_data = county_data.sort_values('date')
            
            # Calculate daily deaths
            county_data['daily_deaths'] = county_data['deaths'].diff().fillna(0)
            
            # Calculate active cases (cases[t] - cases[t-14])
            county_data['active_cases'] = county_data['cases'] - county_data['cases'].shift(self.recovery_days).fillna(0)
            county_data['active_cases'] = county_data['active_cases'].clip(lower=0)
            
            # Calculate P(ID) with 20-day lag
            pid_values = []
            
            for i in range(self.lag_days, len(county_data)):
                date_t = county_data.iloc[i]['date']
                daily_deaths_t = county_data.iloc[i]['daily_deaths']
                
                # Get active cases from 20 days earlier
                active_cases_t_minus_20 = county_data.iloc[i - self.lag_days]['active_cases']
                
                # Include all calculations where active_cases > 0 (including when deaths = 0)
                if active_cases_t_minus_20 > 0:
                    pid = daily_deaths_t / active_cases_t_minus_20
                    pid_values.append(pid)
            
            if len(pid_values) == 0:
                return None
            
            pid_array = np.array(pid_values)
            
            return {
                'county': county_name,
                'state': state,
                'fips': fips,
                'n_calculations': len(pid_values),
                'mean_pid': np.mean(pid_array),
                'median_pid': np.median(pid_array),
                'std_pid': np.std(pid_array),
                'min_pid': np.min(pid_array),
                'max_pid': np.max(pid_array),
                'percentile_25': np.percentile(pid_array, 25),
                'percentile_50': np.percentile(pid_array, 50),
                'percentile_75': np.percentile(pid_array, 75),
                'percentile_90': np.percentile(pid_array, 90),
                'percentile_95': np.percentile(pid_array, 95),
            }

        except Exception as e:
            logger.warning(f"Error processing {county_name}, {state}: {e}")
            return None

    def process_all_counties(self):
        """Process all counties in the NYT dataset."""
        logger.info("="*70)
        logger.info("Processing all counties for P(ID)...")
        logger.info("="*70)
        
        # Get unique county-state-fips combinations
        counties = self.nyt_data.groupby(['fips', 'county', 'state']).size().reset_index()
        counties = counties[counties['fips'] != '00000']  # Remove invalid FIPS
        
        logger.info(f"Found {len(counties)} unique counties")
        
        for _, row in tqdm(counties.iterrows(), total=len(counties), desc="Calculating P(ID)"):
            fips = row['fips']
            county_name = row['county']
            state = row['state']
            
            result = self.calculate_pid_for_county(fips, county_name, state)
            if result is not None:
                self.results.append(result)
        
        logger.info(f"\nSuccessfully calculated P(ID) for {len(self.results)} counties")

    def save_results(self):
        """Save all results."""
        output_dir = Path('../03_ProcessedData/p_ID_analysis')
        output_dir.mkdir(parents=True, exist_ok=True)

        df = pd.DataFrame(self.results)
        df = df.sort_values('mean_pid', ascending=False)

        # Save county-level results
        df.to_csv(output_dir / 'pid_all_counties.csv', index=False)
        logger.info(f"Saved: {output_dir / 'pid_all_counties.csv'}")

        # Save state-level summary
        state_summary = df.groupby('state').agg({
            'mean_pid': ['mean', 'median', 'std', 'count'],
            'n_calculations': 'sum'
        }).round(6)
        state_summary.columns = ['_'.join(col).strip() for col in state_summary.columns]
        state_summary = state_summary.reset_index().sort_values('mean_pid_mean', ascending=False)
        state_summary.to_csv(output_dir / 'pid_by_state.csv', index=False)
        logger.info(f"Saved: {output_dir / 'pid_by_state.csv'}")

        # Save overall summary
        summary = pd.DataFrame({
            'metric': [
                'Overall Mean p_ID',
                'Overall Median p_ID',
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
                df['mean_pid'].mean(),
                df['median_pid'].median(),
                df['mean_pid'].std(),
                df['mean_pid'].quantile(0.05),
                df['mean_pid'].quantile(0.25),
                df['mean_pid'].quantile(0.50),
                df['mean_pid'].quantile(0.75),
                df['mean_pid'].quantile(0.95),
                df['mean_pid'].min(),
                df['mean_pid'].max(),
                len(df),
                df['n_calculations'].sum()
            ],
            'percentage': [
                f"{df['mean_pid'].mean()*100:.4f}%",
                f"{df['median_pid'].median()*100:.4f}%",
                f"{df['mean_pid'].std()*100:.4f}%",
                f"{df['mean_pid'].quantile(0.05)*100:.4f}%",
                f"{df['mean_pid'].quantile(0.25)*100:.4f}%",
                f"{df['mean_pid'].quantile(0.50)*100:.4f}%",
                f"{df['mean_pid'].quantile(0.75)*100:.4f}%",
                f"{df['mean_pid'].quantile(0.95)*100:.4f}%",
                f"{df['mean_pid'].min()*100:.4f}%",
                f"{df['mean_pid'].max()*100:.4f}%",
                '',
                ''
            ]
        })
        summary.to_csv(output_dir / 'pid_summary_statistics.csv', index=False)
        logger.info(f"Saved: {output_dir / 'pid_summary_statistics.csv'}")

        return df

    def print_statistics(self, df):
        """Print summary statistics."""
        logger.info("\n" + "="*70)
        logger.info("P(ID) STATISTICS - ALL U.S. COUNTIES")
        logger.info("="*70)

        logger.info(f"\nTotal counties: {len(df)}")
        logger.info(f"Total calculations: {df['n_calculations'].sum():,}")

        logger.info(f"\nMEAN P(ID) ACROSS ALL COUNTIES:")
        logger.info(f"  Mean:             {df['mean_pid'].mean():.6f} ({df['mean_pid'].mean()*100:.4f}%)")
        logger.info(f"  Median:           {df['mean_pid'].median():.6f} ({df['mean_pid'].median()*100:.4f}%)")
        logger.info(f"  Std Dev:          {df['mean_pid'].std():.6f}")
        logger.info(f"  25th percentile:  {df['mean_pid'].quantile(0.25):.6f} ({df['mean_pid'].quantile(0.25)*100:.4f}%)")
        logger.info(f"  75th percentile:  {df['mean_pid'].quantile(0.75):.6f} ({df['mean_pid'].quantile(0.75)*100:.4f}%)")

        # Carson City validation
        logger.info(f"\nCARSON CITY VALIDATION:")
        carson = df[df['fips'] == '32510']
        if len(carson) > 0:
            carson_pid = carson.iloc[0]['mean_pid']
            logger.info(f"  Calculated: {carson_pid:.6f} ({carson_pid*100:.4f}%)")
            
            # Carson's percentile
            percentile = (df['mean_pid'] < carson_pid).mean() * 100
            logger.info(f"  Percentile: {percentile:.1f}th")
        else:
            logger.warning("  Carson City not found in results!")

        logger.info("="*70)


def main():
    logger.info("="*70)
    logger.info("P(ID) CALCULATION - CARSON CITY METHODOLOGY")
    logger.info("20-day lag: infection → death")
    logger.info("="*70)

    calc = PIDCalculator(lag_days=20, recovery_days=14)
    calc.load_nyt_data('../01_RawData/nyt_covid_data.csv')
    calc.process_all_counties()
    df = calc.save_results()
    calc.print_statistics(df)

    logger.info("\n✓ Complete!")


if __name__ == '__main__':
    main()

