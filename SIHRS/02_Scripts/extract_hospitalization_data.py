#!/usr/bin/env python3
import argparse
import logging
import sys
from pathlib import Path
from typing import Optional, Dict, List
from datetime import datetime
from urllib.parse import urlparse

import pandas as pd
from tqdm import tqdm


# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


class HospitalizationDataExtractor:
    """Extract and filter COVID-19 hospitalization data for multiple counties."""

    def __init__(
        self,
        source_path: str,
        counties_path: str,
        start_date: str,
        end_date: str,
        output_dir: str
    ):
        """
        Initialize the extractor with configuration parameters.

        Args:
            source_path: Path or URL to the main hospitalization CSV
            counties_path: Path to the county list CSV
            start_date: Start date for filtering (YYYY-MM-DD)
            end_date: End date for filtering (YYYY-MM-DD)
            output_dir: Base directory for output files
        """
        self.source_path = source_path
        self.counties_path = counties_path
        self.start_date = pd.to_datetime(start_date)
        self.end_date = pd.to_datetime(end_date)
        self.output_dir = Path(output_dir)

        self.hospital_data: Optional[pd.DataFrame] = None
        self.counties_data: Optional[pd.DataFrame] = None
        self.stats = {
            'total_counties': 0,
            'counties_with_data': 0,
            'counties_without_data': 0,
            'total_rows_processed': 0
        }

    def validate_inputs(self) -> bool:
        """
        Validate all input parameters and files.

        Returns:
            True if all validations pass, False otherwise
        """
        logger.info("Validating inputs...")

        # Check if counties file exists
        if not Path(self.counties_path).exists():
            logger.error(f"Counties file not found: {self.counties_path}")
            return False

        # Check if source is URL or local file
        if self._is_url(self.source_path):
            logger.info(f"Source is URL: {self.source_path}")
        else:
            if not Path(self.source_path).exists():
                logger.error(f"Source file not found: {self.source_path}")
                return False

        # Validate date range
        if self.start_date > self.end_date:
            logger.error(f"Start date ({self.start_date}) is after end date ({self.end_date})")
            return False

        logger.info("All inputs validated successfully")
        return True

    @staticmethod
    def _is_url(path: str) -> bool:
        """Check if the given path is a URL."""
        try:
            result = urlparse(path)
            return all([result.scheme, result.netloc])
        except Exception:
            return False

    def load_hospital_data(self) -> bool:
        """
        Load and filter the main hospitalization dataset.

        Returns:
            True if data loaded successfully, False otherwise
        """
        try:
            logger.info(f"Loading hospital data from: {self.source_path}")

            # Load data from URL or local file
            if self._is_url(self.source_path):
                logger.info("Downloading data from URL (this may take a while)...")
                self.hospital_data = pd.read_csv(self.source_path, low_memory=False)
            else:
                self.hospital_data = pd.read_csv(self.source_path, low_memory=False)

            initial_rows = len(self.hospital_data)
            logger.info(f"Loaded {initial_rows:,} rows")

            # Convert collection_week to datetime
            # Handle various date formats
            if 'collection_week' not in self.hospital_data.columns:
                logger.error("Column 'collection_week' not found in hospital data")
                return False

            self.hospital_data['collection_week'] = pd.to_datetime(
                self.hospital_data['collection_week'],
                errors='coerce'
            )

            # Check for required columns
            required_col = 'total_adult_patients_hospitalized_confirmed_and_suspected_covid_7_day_avg'
            if required_col not in self.hospital_data.columns:
                logger.warning(f"Column '{required_col}' not found. Available columns:")
                logger.warning(f"{list(self.hospital_data.columns)}")

            # Filter by date range
            logger.info(f"Filtering data between {self.start_date.date()} and {self.end_date.date()}")
            self.hospital_data = self.hospital_data[
                (self.hospital_data['collection_week'] >= self.start_date) &
                (self.hospital_data['collection_week'] <= self.end_date)
            ]

            # Filter out null and -999999 values in the main hospitalization column
            if required_col in self.hospital_data.columns:
                before_filter = len(self.hospital_data)
                self.hospital_data = self.hospital_data[
                    self.hospital_data[required_col].notna() &
                    (self.hospital_data[required_col] != -999999)
                ]
                after_filter = len(self.hospital_data)
                logger.info(f"Removed {before_filter - after_filter:,} rows with null/-999999 values")

            # Convert fips_code to string format for matching (handle float FIPS codes)
            if 'fips_code' in self.hospital_data.columns:
                if self.hospital_data['fips_code'].dtype == 'float64':
                    # Convert float to int, then to string, handling NaN
                    self.hospital_data['fips_code'] = self.hospital_data['fips_code'].fillna(0).astype(int).astype(str).str.zfill(5)
                    # Replace '00000' (from NaN) back to empty string for filtering
                    self.hospital_data['fips_code'] = self.hospital_data['fips_code'].replace('00000', '')
                else:
                    self.hospital_data['fips_code'] = self.hospital_data['fips_code'].fillna('').astype(str).str.zfill(5)

            final_rows = len(self.hospital_data)
            logger.info(f"After filtering: {final_rows:,} rows ({(final_rows/initial_rows)*100:.1f}% of original)")

            if final_rows == 0:
                logger.error("No data remaining after filtering")
                return False

            return True

        except Exception as e:
            logger.error(f"Error loading hospital data: {e}")
            return False

    def load_counties_data(self) -> bool:
        """
        Load the county list CSV.

        Returns:
            True if data loaded successfully, False otherwise
        """
        try:
            logger.info(f"Loading counties data from: {self.counties_path}")
            self.counties_data = pd.read_csv(self.counties_path)

            # Check for required columns
            required_cols = ['fips']
            missing_cols = [col for col in required_cols if col not in self.counties_data.columns]

            if missing_cols:
                logger.error(f"Missing required columns in counties file: {missing_cols}")
                logger.error(f"Available columns: {list(self.counties_data.columns)}")
                return False

            # Ensure fips is string with leading zeros preserved
            self.counties_data['fips'] = self.counties_data['fips'].astype(str).str.zfill(5)

            self.stats['total_counties'] = len(self.counties_data)
            logger.info(f"Loaded {self.stats['total_counties']} counties")

            return True

        except Exception as e:
            logger.error(f"Error loading counties data: {e}")
            return False

    def create_county_slug(self, county_name: str, state_abbr: str) -> str:
        """
        Create a URL-friendly slug for the county.

        Args:
            county_name: Name of the county
            state_abbr: State abbreviation

        Returns:
            Slugified county name
        """
        # Remove common suffixes and clean the name
        name = county_name.lower()
        name = name.replace(' county', '').replace(' city', '').replace(' parish', '')
        name = name.replace(' ', '_').replace(',', '').replace('.', '')
        return f"{name}_{state_abbr.lower()}"

    def extract_county_data(self, county_row: pd.Series) -> Optional[pd.DataFrame]:
        """
        Extract hospitalization data for a single county.

        Args:
            county_row: Row from the counties DataFrame

        Returns:
            Filtered DataFrame for the county, or None if no data
        """
        fips = str(county_row['fips']).zfill(5)

        # Ensure fips_code in hospital data is also string with leading zeros
        if 'fips_code' not in self.hospital_data.columns:
            logger.error("Column 'fips_code' not found in hospital data")
            return None

        # FIPS codes should already be converted in load_hospital_data()

        # Filter by FIPS code
        county_data = self.hospital_data[self.hospital_data['fips_code'] == fips].copy()

        if len(county_data) == 0:
            return None

        # Aggregate by date (sum across multiple hospitals in the same county)
        # This is crucial for counties with multiple hospitals reporting
        numeric_columns = [
            'total_adult_patients_hospitalized_confirmed_and_suspected_covid_7_day_avg',
            'total_adult_patients_hospitalized_confirmed_covid_7_day_avg',
            'total_adult_patients_hospitalized_suspected_covid_7_day_avg',
            'staffed_adult_icu_bed_occupancy_7_day_avg',
            'total_pediatric_patients_hospitalized_confirmed_and_suspected_covid_7_day_avg',
            'inpatient_beds_used_covid_7_day_avg',
            'percentage_of_inpatient_beds_used_covid_7_day_avg'
        ]
        
        # Only aggregate columns that exist
        agg_columns = {col: 'sum' for col in numeric_columns if col in county_data.columns}
        
        # Keep state and fips_code as well
        groupby_cols = ['collection_week']
        if 'state' in county_data.columns:
            groupby_cols.append('state')
        if 'fips_code' in county_data.columns:
            groupby_cols.append('fips_code')
        
        # Aggregate by date, summing across hospitals
        county_data = county_data.groupby(groupby_cols).agg(agg_columns).reset_index()
        
        # Calculate total_adult_and_pediatric_covid_patients
        if 'total_adult_patients_hospitalized_confirmed_and_suspected_covid_7_day_avg' in county_data.columns:
            county_data['total_adult_and_pediatric_covid_patients'] = (
                county_data['total_adult_patients_hospitalized_confirmed_and_suspected_covid_7_day_avg'].fillna(0) +
                county_data.get('total_pediatric_patients_hospitalized_confirmed_and_suspected_covid_7_day_avg', 0).fillna(0)
            )
        
        # Sort by date
        county_data = county_data.sort_values('collection_week')

        return county_data

    def process_all_counties(self) -> bool:
        """
        Process all counties and extract their hospitalization data.

        Returns:
            True if processing completed (even if some counties had no data)
        """
        logger.info(f"Processing {self.stats['total_counties']} counties...")

        # Create output directory if it doesn't exist
        self.output_dir.mkdir(parents=True, exist_ok=True)

        for idx, county_row in tqdm(
            self.counties_data.iterrows(),
            total=len(self.counties_data),
            desc="Processing counties"
        ):
            county_name = county_row.get('name', 'Unknown')
            state_abbr = county_row.get('state', 'XX')
            fips = str(county_row['fips']).zfill(5)

            # Extract data for this county
            county_data = self.extract_county_data(county_row)

            if county_data is None or len(county_data) == 0:
                logger.warning(f"No data for {county_name}, {state_abbr} (FIPS: {fips})")
                self.stats['counties_without_data'] += 1
                continue

            # Create output directory structure to match existing pattern:
            # CaseStudies/<County>,<State>/codes/hospitalization_<county>_filtered_new.csv
            county_dir_name = f"{county_name},{state_abbr}"
            county_dir = self.output_dir / county_dir_name / "codes"
            county_dir.mkdir(parents=True, exist_ok=True)

            # Create output filename matching existing pattern
            # e.g., "hospitalization_Carson_filtered_new.csv" or "hospitalization_MS_filtered.csv"
            county_name_clean = county_name.replace(' County', '').replace(' City', '').replace(' ', '_')
            output_file = county_dir / f"hospitalization_{county_name_clean}_filtered_new.csv"

            # Save to CSV
            county_data.to_csv(output_file, index=False)

            self.stats['counties_with_data'] += 1
            self.stats['total_rows_processed'] += len(county_data)

            logger.info(
                f"✓ {county_name}, {state_abbr}: {len(county_data)} rows → {output_file}"
            )

        return True

    def print_summary(self):
        """Print a summary of the extraction process."""
        logger.info("\n" + "="*70)
        logger.info("EXTRACTION SUMMARY")
        logger.info("="*70)
        logger.info(f"Total counties processed:     {self.stats['total_counties']}")
        logger.info(f"Counties with data:           {self.stats['counties_with_data']}")
        logger.info(f"Counties without data:        {self.stats['counties_without_data']}")
        logger.info(f"Total rows extracted:         {self.stats['total_rows_processed']:,}")
        logger.info(f"Date range:                   {self.start_date.date()} to {self.end_date.date()}")
        logger.info(f"Output directory:             {self.output_dir.absolute()}")
        logger.info("="*70)

    def run(self) -> bool:
        """
        Run the complete extraction pipeline.

        Returns:
            True if successful, False otherwise
        """
        if not self.validate_inputs():
            return False

        if not self.load_counties_data():
            return False

        if not self.load_hospital_data():
            return False

        if not self.process_all_counties():
            return False

        self.print_summary()
        return True


def main():
    """Main entry point for the script."""
    parser = argparse.ArgumentParser(
        description='Extract COVID-19 hospitalization data for multiple counties',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Using a local file
  python extract_hospitalization_data.py \\
    --source ../01_RawData/covid_hospital_report.csv \\
    --counties county_list.csv

  # Using a URL
  python extract_hospitalization_data.py \\
    --source "https://healthdata.gov/resource/data.csv" \\
    --counties county_list.csv \\
    --start-date 2020-03-01 \\
    --end-date 2021-12-31

  # Custom output directory
  python extract_hospitalization_data.py \\
    --source ../01_RawData/covid_hospital_report.csv \\
    --counties county_list.csv \\
    --output-dir ../03_ProcessedData
        """
    )

    parser.add_argument(
        '--source',
        required=True,
        help='Path or URL to the COVID-19 Hospital Report CSV file'
    )

    parser.add_argument(
        '--counties',
        required=True,
        help='Path to the CSV file containing county list (must have columns: name, state, fips, population, R0, rating)'
    )

    parser.add_argument(
        '--start-date',
        default='2020-03-12',
        help='Start date for filtering (YYYY-MM-DD). Default: 2020-03-12'
    )

    parser.add_argument(
        '--end-date',
        default='2021-12-31',
        help='End date for filtering (YYYY-MM-DD). Default: 2021-12-31'
    )

    parser.add_argument(
        '--output-dir',
        default='../03_ProcessedData',
        help='Base directory for output files. Default: ../03_ProcessedData'
    )

    parser.add_argument(
        '--verbose',
        '-v',
        action='store_true',
        help='Enable verbose logging'
    )

    args = parser.parse_args()

    # Set logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Create extractor and run
    extractor = HospitalizationDataExtractor(
        source_path=args.source,
        counties_path=args.counties,
        start_date=args.start_date,
        end_date=args.end_date,
        output_dir=args.output_dir
    )

    success = extractor.run()

    if success:
        logger.info("\n✓ Extraction completed successfully!")
        sys.exit(0)
    else:
        logger.error("\n✗ Extraction failed!")
        sys.exit(1)


if __name__ == '__main__':
    main()
