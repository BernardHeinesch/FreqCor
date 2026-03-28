# -*- coding: utf-8 -*-
"""
FREQCOR_Start.py: Main entry point for FREQCOR spectral/cospectral correction analysis
=====================================================================================

Overview:
---------
FREQCOR is a tool for applying spectral and cospectral correction procedures to eddy covariance
datasets, supporting both stable and unstable atmospheric conditions. This script orchestrates
the full workflow, from reading and filtering raw data to applying corrections and generating
summary outputs.

Key Features:
-------------
- Batch processing of multiple sites and configurations
- Automatic reading of (co)spectra, meteorological, and quality control files
- Flexible selection and filtering of data (e.g., by wind sector, quality)
- Calculation of correction factors for both stable and unstable regimes
- Application of correction factors to flux data and computation of uncertainties
- Generation of summary reports for each site and for all processed sites
- Command-line interface for automation and integration

Workflow:
---------
1. **Parameter Setup:**
   - Reads user-specified parameters from .ini configuration files.
2. **Data Reading:**
   - Loads (co)spectra, meteorological, and data quality files for each site.
3. **Data Selection:**
   - Filters out files based on quality control and user-defined criteria.
4. **Correction Factor Computation:**
   - Computes correction factors for unstable and stable conditions.
5. **Flux Correction:**
   - Applies correction factors to the full dataset and calculates uncertainties.
6. **Reporting:**
   - Outputs results and summary tables for each site and for all processed sites.

Usage:
------
- **Manual mode:** Set `manual = True` and specify an `.ini` file to process a single
  configuration interactively.

- **Batch mode:** Run from the command line to process all sites found in the input directory,
  each with its own subfolder of `.ini` files under the ini directory:

    python FREQCOR_Start.py [--sites SITE1 SITE2 ...] [--list-sites]
                            [--input-dir DIR] [--ini-dir DIR] [--results-dir DIR]

- **Flat mode:** Pass a single directory containing any number of `.ini` files directly.
  No site-name pattern or subfolder structure is required. Mutually exclusive with
  --sites, --list-sites, --input-dir and --ini-dir:

    python FREQCOR_Start.py --flat-ini-dir DIR [--results-dir DIR]

  If --results-dir is omitted in flat mode, the summary CSV is written to flat-ini-dir itself.

- See the README or documentation for more details on configuration and options.

Authorship & Acknowledgments:
-----------------------------
- Original MATLAB version by Marc Aubinet
- Python translation and extension: Ariane Faures
- Contributions: [see repository contributors]

References & Documentation:
---------------------------
- For full documentation, variable explanations, and workflow schematics, see the project README
  and supplementary docs.
- For issues or contributions, visit the project repository.

Created: June 28, 2022
Last updated: August 28, 2025
"""

# %% Package and functions import
import os
import glob
import time
import pandas as pd
import argparse
import sys
import re
import traceback
from FREQCOR_Main import FREQCOR_Main, FREQCORMainError
from FREQCOR_Read_EP import FREQCORReadError
from FREQCOR_Sel_cof import FREQCORSelCofError
from FREQCOR_Compute import FREQCORComputeError
from FREQCOR_LUT_cof import FREQCORLUTCofError
from FREQCOR_LUT_CF import FREQCORLUTCFError
import configparser


## MANUAL MODE ##
# Runs the code directly with the .ini file provided.
# Insert ini file to be processed
manual = True
if manual == True:
    _src_dir = os.path.dirname(os.path.abspath(__file__))
    ini_file = os.path.normpath(os.path.join(_src_dir, '..', 'examples', 'metadata', 'FREQCOR_config_example.ini'))
    os.chdir(_src_dir)
    config = configparser.ConfigParser()
    config.read(ini_file)
    LUTcof, LUTCF_u, LUTCF_s = FREQCOR_Main(config)
    sys.exit()

# %% Site Pattern
# looks for this pattern in sites_dir path to detect the list of sites to process
# when changing it, use the format ^() as the code will then look for a group
# in the pattern
# site_pattern = re.compile(r'^([A-Z]{2}-[A-Za-z0-9]{3})_output_EP.*') #ETC
site_pattern = re.compile(r'^([A-Z]{2}-[A-Za-z0-9]{3})')

# %% Define default paths
def get_default_paths():
    # sites_dir=r'D:\Ariane\spectral_analysis\a_inputs'
    sites_dir=r'C:\ICOS\Multi-site\spectral_corrections_paper\data\icos_sites\b_runs_ini\b\2024'
    output_base_dir=r'C:\ICOS\Multi-site\spectral_corrections_paper\data\icos_sites'
    # sites_dir= r'/home/g.nicolini/dev/ETCPS'
    # output_base_dir =  r'/home/a.faures/dev/multisite_outputs'

    paths = {
        # Path to site folders
        'sites_dir': sites_dir,
        # Path to working and output directory: change run letter (eg. 'a') if
        # needed, according to the ini used
        'ini_base_dir': os.path.join(output_base_dir, 'b_runs_ini', 'b', '2024'),  #,'2024'
        'results_base_dir': os.path.join(output_base_dir, 'd_results', '2024')     #,'2024'
    }
    return paths


# %% Core processing function
def process_ini_files(ini_files, results_dir, label):
    """Run FREQCOR_Main on a list of .ini files and write a summary CSV.

    This is the single processing core shared by both batch and flat modes.
    It iterates over the provided list of .ini files, runs FREQCOR_Main on
    each, handles all known exceptions, and saves a summary CSV to results_dir.

    Parameters
    ----------
    ini_files : list of str
        Full paths to the .ini files to process, in the order they will be run.
    results_dir : str
        Directory where the summary CSV will be written. Created if absent.
    label : str
        Human-readable label used in console output and in the summary CSV
        (e.g. a site name such as 'BE-Dor', or 'flat' for flat mode).

    Returns
    -------
    list of dict
        One entry per .ini file with keys: label, ini_file, gas, method,
        status, duration_seconds.
    """
    print(f"\n{'='*50}")
    print(f"Processing: {label}")
    print(f"{'='*50}")

    os.makedirs(results_dir, exist_ok=True)

    results = []
    for ini_file in ini_files:
        ini_filename = os.path.basename(ini_file)
        print(f"\nProcessing configuration: {ini_filename}")

        # Read the configuration
        config = configparser.ConfigParser()
        config.read(ini_file)

        # Extract configuration details for results tracking
        try:
            gas = 'co2' if config.get('PROCEDURE_OPTIONS', 'gss') == '1' else 'h2o'
            method = 'cosp' if config.get('PROCEDURE_OPTIONS', 'sps') == '1' else 'sp'
        except (configparser.NoSectionError, configparser.NoOptionError):
            gas = 'unknown'
            method = 'unknown'

        # Run FREQCOR_Main and handle all known error types
        start_time = time.time()
        try:
            LUTcof, LUTCF_u, LUTCF_s = FREQCOR_Main(config)
            status = "Success"
        except FREQCORReadError as e:
            print(f"Data reading error for {ini_filename}: {str(e)}")
            status = f"Read Error: {str(e)}"
        except FREQCORMainError as e:
            print(f"Processing error for {ini_filename}: {str(e)}")
            status = f"Processing Error: {str(e)}"
        except FREQCORSelCofError as e:
            print(f"Selection error for {ini_filename}: {str(e)}")
            status = f"Selection Error: {str(e)}"
        except FREQCORComputeError as e:
            print(f"Computation error for {ini_filename}: {str(e)}")
            status = f"Computation Error: {str(e)}"
        except FREQCORLUTCofError as e:
            print(f"LUT computation error for {ini_filename}: {str(e)}")
            status = f"LUT Error: {str(e)}"
        except FREQCORLUTCFError as e:
            print(f"LUT CF computation error for {ini_filename}: {str(e)}")
            status = f"LUT CF Error: {str(e)}"
        except Exception as e:
            print(f"Unexpected error processing {ini_filename}: {str(e)}")
            traceback.print_exc()
            status = f"Error: {str(e)}"
        duration = time.time() - start_time

        results.append({
            'label': label,
            'ini_file': ini_filename,
            'gas': gas,
            'method': method,
            'status': status,
            'duration_seconds': duration
        })

    # Write per-label summary CSV
    if results:
        results_df = pd.DataFrame(results)
        report_path = os.path.join(results_dir, f"{label}_processing_summary.csv")
        results_df.to_csv(report_path, index=False)
        print(f"\nSummary report saved to: {report_path}")

    return results


# %% Main execution
def main(mode, paths=None, specific_sites=None):
    """Orchestrate processing for batch or flat mode.

    Parameters
    ----------
    mode : str
        'batch' to process sites discovered via site-pattern matching, or
        'flat' to process all .ini files directly from paths['flat_ini_dir'].
    paths : dict or None
        Path configuration. When None, get_default_paths() is used for batch
        mode. Expected keys depend on mode:
        - batch: 'sites_dir', 'ini_base_dir', 'results_base_dir'
        - flat:  'flat_ini_dir', 'results_base_dir'
    specific_sites : list of str or None
        Batch mode only. When provided, restricts processing to the listed
        site codes. Ignored in flat mode.
    """
    if paths is None:
        paths = get_default_paths()

    # ------------------------------------------------------------------ #
    # FLAT MODE: process every .ini file in a single directory directly   #
    # ------------------------------------------------------------------ #
    if mode == 'flat':
        flat_ini_dir = paths['flat_ini_dir']
        results_base_dir = paths['results_base_dir']

        ini_files = sorted(glob.glob(os.path.join(flat_ini_dir, "*.ini")))
        if not ini_files:
            print(f"Warning: No .ini files found in {flat_ini_dir}")
            return

        print(f"Found {len(ini_files)} .ini file(s) to process.")
        process_ini_files(ini_files, results_base_dir, label='flat')

    # ------------------------------------------------------------------ #
    # BATCH MODE: discover sites via pattern, process each one in turn    #
    # ------------------------------------------------------------------ #
    else:
        # Look for available sites in the input folder
        all_site_dirs = []
        site_folder_mapping = {}

        for folder in os.listdir(paths['sites_dir']):
            folder_path = os.path.join(paths['sites_dir'], folder)
            if os.path.isdir(folder_path):
                match = site_pattern.match(folder)
                if match:
                    site_name = match.group(1)  # Extract just the XX-Ccc part
                    all_site_dirs.append(site_name)
                    site_folder_mapping[site_name] = folder

        if not all_site_dirs:
            print("No site directories found!")
            return

        # Filter to requested sites if provided
        if specific_sites:
            site_dirs = [d for d in all_site_dirs if d in specific_sites]
            if not site_dirs:
                print(f"None of the specified sites {specific_sites} were found. "
                      f"Available sites: {all_site_dirs}")
                return
        else:
            site_dirs = all_site_dirs

        print(f"Found {len(site_dirs)} sites to process: {', '.join(site_dirs)}")

        # Process each site
        all_results = []
        for site in site_dirs:
            site_ini_dir = os.path.join(paths['ini_base_dir'], site)

            if not os.path.exists(site_ini_dir):
                print(f"Warning: No ini directory found for site {site} at {site_ini_dir}")
                continue

            ini_files = glob.glob(os.path.join(site_ini_dir, "*.ini"))
            if not ini_files:
                print(f"Warning: No ini files found for site {site}")
                continue

            site_results_dir = os.path.join(paths['results_base_dir'], site)
            site_results = process_ini_files(ini_files, site_results_dir, label=site)
            if site_results:
                all_results.extend(site_results)

        # Write overall summary across all sites
        if all_results:
            all_results_df = pd.DataFrame(all_results)
            overall_report_path = os.path.join(paths['results_base_dir'],
                                               "overall_processing_summary.csv")
            all_results_df.to_csv(overall_report_path, index=False)
            print(f"\nOverall summary report saved to: {overall_report_path}")

    print("\nAll processing complete!")


# %% Run the main function
if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Process multiple sites with FREQCOR')
    parser.add_argument('--sites', nargs='+',
                        help='Batch mode: specific sites to process (e.g. BE-Lon BE-Dor)')
    parser.add_argument('--list-sites', action='store_true',
                        help='Batch mode: list available sites and exit')
    parser.add_argument('--input-dir',
                        help='Batch mode: base directory for site inputs '
                             '(default: from get_default_paths)')
    parser.add_argument('--ini-dir',
                        help='Batch mode: base directory for ini files '
                             '(default: from get_default_paths)')
    parser.add_argument('--results-dir',
                        help='Base directory for results (default: from get_default_paths). '
                             'In flat mode, defaults to flat-ini-dir if not specified.')
    parser.add_argument('--flat-ini-dir',
                        help='Flat mode: directory containing .ini files to process '
                             'sequentially, with no site-name pattern or subfolder structure. '
                             'Mutually exclusive with --sites, --list-sites, --input-dir '
                             'and --ini-dir.')
    args = parser.parse_args()

    # Validate flat mode exclusivity
    if args.flat_ini_dir and (args.sites or args.list_sites or args.input_dir or args.ini_dir):
        parser.error("--flat-ini-dir cannot be combined with --sites, --list-sites, "
                     "--input-dir or --ini-dir.")

    # ------------------------------------------------------------------ #
    # FLAT MODE                                                           #
    # ------------------------------------------------------------------ #
    if args.flat_ini_dir:
        paths = {
            'flat_ini_dir': args.flat_ini_dir,
            # Fall back to the ini directory itself if no results dir is given
            'results_base_dir': args.results_dir if args.results_dir else args.flat_ini_dir
        }
        main(mode='flat', paths=paths)

    # ------------------------------------------------------------------ #
    # BATCH MODE                                                          #
    # ------------------------------------------------------------------ #
    else:
        paths = get_default_paths()
        if args.input_dir:
            paths['sites_dir'] = args.input_dir
        if args.ini_dir:
            paths['ini_base_dir'] = args.ini_dir
        if args.results_dir:
            paths['results_base_dir'] = args.results_dir

        # Just list sites if requested, then exit without processing
        if args.list_sites:
            available_sites = []
            for folder in os.listdir(paths['sites_dir']):
                if os.path.isdir(os.path.join(paths['sites_dir'], folder)):
                    match = site_pattern.match(folder)
                    if match:
                        available_sites.append(match.group(1))
            if available_sites:
                print(f"Available sites: {', '.join(available_sites)}")
            else:
                print("No sites found.")
            sys.exit(0)

        main(mode='batch', paths=paths, specific_sites=args.sites)


