
# -*- coding: utf-8 -*-

"""Input reading utilities for the FREQCOR processing chain."""

import datetime
import glob
import json
import numpy as np
import os
import pandas as pd
import re
import sys
from .FREQCOR_VM_flag import FREQCOR_VM_flag

# Custom exception class for FREQCOR_Read errors
class FREQCORReadError(Exception):
    """Exception raised for errors in the FREQCOR_Read function.
    
    Attributes:
        message -- explanation of the error
    """
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)
        
def FREQCOR_Read_EP(config, jsite):
    """
    Read data files from EddyPro and return ideal and real (co)spectra dataframes,
    relevant meteo and fluxes data as well as correspondent qc flags.
    Compatible input files are EddyPro binned cospectra csv files.

    Parameters
    ----------
    config : configparser.ConfigParser
        Parsed configuration.
    jsite : str
        Station identifier, ICOS format (CO-Unt).

    Returns
    -------
    list
     [nspec, Icon_raw, Rcon_raw, Icosp_raw, nfreq, freq, freqn, freqcount, meteo_df,
     WS, WD, Ustar, Zeta, FlagF_H, FlagF_g, FlagVM_w, FlagVM_T, FlagVM_g,
     flag_wd, massman_coef].
     nspec : int
         Number of valid (co)spectra files in the directory.
     Icon_raw : DataFrame
         Containing all valid ideal content (raw) (co)spectra (either cospectra or spectra; depends on `sps`).
         Column names are the half-hours names.
     Rcon_raw : DataFrame
         Containing all valid real content (raw) (co)spectra (either cospectra or spectra; depends on `sps`).
         Column names are the half-hours names.
     Icosp_raw : DataFrame
         Containing cospectra used for CF computation (raw).
         Column names are the half-hours names.
     nfreq : int
         Actual number of frequencies with spectral content.
     freq : Series
         Natural frequencies.
     freqn : DataFrame
         Normalised frequencies. Column names are the half-hours names.
     freqcount : Series
         Number of frequencies per bin for weighting calculations.
     meteo_df : DataFrame
         Contains relevant fluxes and meteo data.
     WS : Series
         Wind speed.
     WD : Series
         Wind direction.
     Ustar : Series
         Friction velocity (u*).
     Zeta : Series
         Stability parameter (z-d)/L.
     FlagF_H : Series
         Mauder & Foken (2004) quality flag on sensible heat flux.
     FlagF_g : Series
         Mauder & Foken (2004) quality flag on gas flux.
     FlagVM_w : Array of int
         Vickers & Mahrt (1997) flag on vertical wind speed.
     FlagVM_T : Array of int
         Vickers & Mahrt (1997) flag on sonic temperature.
     FlagVM_g : Array of int
         Vickers & Mahrt (1997) flag on gas concentration.
     massman_coef : array
         Massman coefficients for the site [A0, kf0, mu].

     Raises
     ------
     FREQCORReadError
        If inputs cannot be read or validated (e.g., no files found 
        or file counts don't match between different data sources).
    """
    
    # Extract variables
    dirname = config['IO']['input_path']
    dirname_sp = config['IO']['input_path_sp']
    massmanpath = config['IO']['massman_path']
    fcname = config['IO']['binned_cosp']
    foname = config['IO']['flx_meteo']
    sps = int(config['PROCEDURE_OPTIONS']['sps'])
    gss = int(config['PROCEDURE_OPTIONS']['gss'])
    output_path = config['IO']['output_path']
    output_pattern = config['IO']['output_file']
    vitale_qc_flags = int(config['PROCEDURE_OPTIONS']['vitale_qc_flags'])
    if vitale_qc_flags == 1:
        vitale_path = config['IO']['vitale_path']
    #%% Load data
    # Existing data from previous run
    if 'PROCESSED_DATA' in config and int(config['PROCESSED_DATA']['enable_loading']) == 1:
        try:
            data_dir = os.path.join(output_path,f'{output_pattern}_intermediate_data')
            file_prefix = config['PROCESSED_DATA']['file_prefix']
            
            print(f"Attempting to load previously processed data from {data_dir}...")
            
            # Check if metadata file exists to verify all required files are present
            metadata_path = os.path.join(data_dir, f"{file_prefix}_metadata.json")
            if not os.path.exists(metadata_path):
                print(f"Metadata file not found at {metadata_path}. Will process data from scratch.")
            else:
                # Load metadata to check what files should be available
                with open(metadata_path, 'r') as f:
                    metadata = json.load(f)
                
                # Check if all required files exist
                all_files_exist = True
                for filename in metadata['files_saved']:
                    if not os.path.exists(os.path.join(data_dir, filename)):
                        print(f"Missing file: {filename}. Will process data from scratch.")
                        all_files_exist = False
                        break
                
                if all_files_exist:
                    print("Loading all data from saved files...")
                    
                    # Load scalar values
                    with open(os.path.join(data_dir, f"{file_prefix}_scalar_values.json"), 'r') as f:
                        scalar_values = json.load(f)
                    nspec = scalar_values['nspec']
                    nfreq = scalar_values['nfreq']
                    
                    # Load DataFrames
                    Icon_raw = pd.read_pickle(os.path.join(data_dir, f"{file_prefix}_Icon_raw.pkl"))
                    Rcon_raw = pd.read_pickle(os.path.join(data_dir, f"{file_prefix}_Rcon_raw.pkl"))
                    Icosp_raw = pd.read_pickle(os.path.join(data_dir, f"{file_prefix}_Icosp_raw.pkl"))
                    freqn = pd.read_pickle(os.path.join(data_dir, f"{file_prefix}_freqn.pkl"))
                    meteo_df = pd.read_pickle(os.path.join(data_dir, f"{file_prefix}_meteo_df.pkl"))
                    
                    # Load Series
                    freq = pd.read_pickle(os.path.join(data_dir, f"{file_prefix}_freq.pkl"))
                    WS = pd.read_pickle(os.path.join(data_dir, f"{file_prefix}_WS.pkl"))
                    WD = pd.read_pickle(os.path.join(data_dir, f"{file_prefix}_WD.pkl"))
                    Ustar = pd.read_pickle(os.path.join(data_dir, f"{file_prefix}_Ustar.pkl"))
                    Zeta = pd.read_pickle(os.path.join(data_dir, f"{file_prefix}_Zeta.pkl"))
                    massman_coef = np.load(os.path.join(data_dir, f"{file_prefix}_massman_coef.npy"), allow_pickle=True)
                    if vitale_qc_flags == 1:
                        FlagF_H = pd.read_pickle(os.path.join(data_dir, f"{file_prefix}_FlagF_H_vitale.pkl"))
                        FlagF_g = pd.read_pickle(os.path.join(data_dir, f"{file_prefix}_FlagF_g_vitale.pkl"))
                        flag_wd = pd.read_pickle(os.path.join(data_dir, f"{file_prefix}_flag_wd_vitale.pkl"))
                        
                        # Load arrays
                        FlagVM_w = np.load(os.path.join(data_dir, f"{file_prefix}_FlagVM_w_vitale.npy"))
                        FlagVM_T = np.load(os.path.join(data_dir, f"{file_prefix}_FlagVM_T_vitale.npy"))
                        FlagVM_g = np.load(os.path.join(data_dir, f"{file_prefix}_FlagVM_g_vitale.npy"))
                    else: 
                        FlagF_H = pd.read_pickle(os.path.join(data_dir, f"{file_prefix}_FlagF_H.pkl"))
                        FlagF_g = pd.read_pickle(os.path.join(data_dir, f"{file_prefix}_FlagF_g.pkl"))
                        
                        # Load arrays
                        FlagVM_w = np.load(os.path.join(data_dir, f"{file_prefix}_FlagVM_w.npy"))
                        FlagVM_T = np.load(os.path.join(data_dir, f"{file_prefix}_FlagVM_T.npy"))
                        FlagVM_g = np.load(os.path.join(data_dir, f"{file_prefix}_FlagVM_g.npy"))
                        try : # to make this version compatible with older intermediate data outputs
                        # where the flag_wd variable was not considered
                            flag_wd = pd.read_pickle(os.path.join(data_dir, f"{file_prefix}_flag_wd.pkl"))
                        except:
                            flag_wd = pd.Series([np.nan] * len(FlagF_H), index=FlagF_H.index) 
                    print("Successfully loaded all processed data!")
                    print(f"Loaded {nspec} spectra with {nfreq} frequency bins.")
                    
                    # Create the array with count of frequencies per bin if absent
                    try:
                        freqcount = pd.read_pickle(os.path.join(data_dir, f"{file_prefix}_freqcount.pkl"))
                    except FileNotFoundError:
                        freqcount = pd.Series(320*freq)  # Default weighting factor taken from example site
                    # Return the loaded data
                    return [nspec, Icon_raw, Rcon_raw, Icosp_raw, nfreq, freq, freqn, freqcount, meteo_df, WS, WD, Ustar, 
                            Zeta, FlagF_H, FlagF_g, FlagVM_w, FlagVM_T, FlagVM_g, flag_wd, massman_coef]
                    
        except Exception as e:
            print(f"Error loading processed data: {e}")
    # Read time window parameters if enabled
    enable_time_window = 0
    exclusion_windows = []          # list of (start, end) datetime pairs to exclude
    if 'TIME_WINDOW' in config:
        if 'enable_time_window' in config['TIME_WINDOW']:
            enable_time_window = int(config['TIME_WINDOW']['enable_time_window'])
        
        if enable_time_window == 1:
            start_datetime_str = config['TIME_WINDOW']['start_datetime']
            end_datetime_str = config['TIME_WINDOW']['end_datetime']
            
            # Parse datetime strings to datetime objects
            start_datetime = datetime.datetime.strptime(start_datetime_str, '%Y-%m-%d %H:%M')
            end_datetime = datetime.datetime.strptime(end_datetime_str, '%Y-%m-%d %H:%M')
            
            print(f"  Time window enabled: {start_datetime} to {end_datetime}")
        
        # ------------------------------------------------------------------
        # Parse optional exclusion windows
        # ------------------------------------------------------------------
        enable_exclusion = 0
        if 'enable_exclusion_windows' in config['TIME_WINDOW']:
            enable_exclusion = int(config['TIME_WINDOW']['enable_exclusion_windows'])
        
        if enable_exclusion == 1:
            if 'date_exclusion_start' not in config['TIME_WINDOW'] or \
               'date_exclusion_end'   not in config['TIME_WINDOW']:
                raise FREQCORReadError(
                    "enable_exclusion_windows = 1 but date_exclusion_start and/or "
                    "date_exclusion_end are missing from the [TIME_WINDOW] section."
                )
            
            # Split comma-separated lists and strip whitespace
            excl_starts_raw = [s.strip() for s in
                               config['TIME_WINDOW']['date_exclusion_start'].split(',')
                               if s.strip()]
            excl_ends_raw   = [s.strip() for s in
                               config['TIME_WINDOW']['date_exclusion_end'].split(',')
                               if s.strip()]
            
            # Paired count check
            if len(excl_starts_raw) != len(excl_ends_raw):
                raise FREQCORReadError(
                    f"date_exclusion_start has {len(excl_starts_raw)} entries but "
                    f"date_exclusion_end has {len(excl_ends_raw)}. "
                    "Each exclusion window must have both a start and an end date."
                )
            
            for i, (es_raw, ee_raw) in enumerate(zip(excl_starts_raw, excl_ends_raw), start=1):
                try:
                    es = datetime.datetime.strptime(es_raw, '%Y-%m-%d %H:%M')
                except ValueError:
                    raise FREQCORReadError(
                        f"Exclusion window {i}: cannot parse start date '{es_raw}'. "
                        "Expected format: YYYY-MM-DD HH:MM"
                    )
                try:
                    ee = datetime.datetime.strptime(ee_raw, '%Y-%m-%d %H:%M')
                except ValueError:
                    raise FREQCORReadError(
                        f"Exclusion window {i}: cannot parse end date '{ee_raw}'. "
                        "Expected format: YYYY-MM-DD HH:MM"
                    )
                
                # Warn (not error) if exclusion window is fully outside the main window
                if enable_time_window == 1:
                    if ee <= start_datetime or es >= end_datetime:
                        print(
                            f"Warning: Exclusion window {i} ({es} â€“ {ee}) lies entirely "
                            f"outside the main processing window ({start_datetime} â€“ "
                            f"{end_datetime}) and will have no effect."
                        )
                    elif es < start_datetime or ee > end_datetime:
                        print(
                            f"Warning: Exclusion window {i} ({es} â€“ {ee}) partially "
                            f"extends beyond the main processing window "
                            f"({start_datetime} â€“ {end_datetime}). "
                            "Only the overlapping portion will be excluded."
                        )
                
                exclusion_windows.append((es, ee))
            
            if exclusion_windows:
                print(f"Exclusion windows enabled: {len(exclusion_windows)} window(s) defined.")
                for i, (es, ee) in enumerate(exclusion_windows, start=1):
                    print(f"  Exclusion {i}: {es} to {ee}")
    
    # Read fluxes, weather and qc data.
    # fluxes and weather data
    fileb = glob.glob(os.path.join(dirname,'*'+foname+'*_adv*.csv'))[0]
    with open(fileb, mode='r') as file:
        meteo_data = pd.read_csv(file, header=1, na_values="-9999", skiprows=[2],
                                low_memory=False, memory_map=True)
    # Convert timestamp strings to datetime objects for filtering
    # If timestamp column doesn't exist, try to create it from date and time columns
    if 'timestamp' in meteo_data.columns:
        meteo_data['datetime'] = pd.to_datetime(meteo_data['timestamp'])
    else:
        meteo_data['datetime'] = pd.to_datetime(meteo_data['date'] + ' ' + meteo_data['time'])
    
    if vitale_qc_flags == 1:
        try:
            qc_files = glob.glob(os.path.join(vitale_path, jsite + '_vitale_qc_flags*.csv'))
            filec = qc_files[0]
            with open(filec, mode='r') as file:
                flag_data = pd.read_csv(
                    file,
                    header=0,
                    parse_dates=["TIMESTAMP_END"],
                    low_memory=False,
                    memory_map=True,
                )
                flag_data['datetime'] = flag_data['TIMESTAMP_END']
        except Exception as e:
            print(f"Error reading Vitale QC details file: {e}. Will switch to EddyPro output flags.")
            vitale_qc_flags = 0
            
    # Remove missing data rows
    meteo_data = meteo_data[meteo_data.filename != 'not_enough_data']
    # flag_data = flag_data[flag_data.filename[:] != 'not_enough_data']
    # meteo_data = meteo_data.reset_index(drop=True)
    # flag_data = flag_data.reset_index(drop=True)
    
    # Filter by time window if enabled
    all_cospectra_files = glob.glob(os.path.join(dirname_sp,'eddypro_binned_cospectra','*' + fcname + '*_adv.csv'))
    filtered_cospectra_files = all_cospectra_files.copy() 
    if enable_time_window == 1:
        meteo_data = meteo_data[(meteo_data['datetime'] >= start_datetime) & 
                              (meteo_data['datetime'] <= end_datetime)]
        # flag_data = flag_data[(flag_data['datetime'] >= start_datetime) & 
        #                     (flag_data['datetime'] <= end_datetime)]
        # meteo_data = meteo_data.reset_index(drop=True)
        # flag_data = flag_data.reset_index(drop=True)
        
        
        
    # ------------------------------------------------------------------
    # Apply exclusion windows to meteo_data (within or without main window)
    # Exclusion is applied as a mask: timestamps inside ANY exclusion interval
    # are removed.  The cospectra file list is derived from the surviving
    # meteo_data timestamps afterwards, keeping the existing logic intact.
    # ------------------------------------------------------------------
    if exclusion_windows:
        before_excl = len(meteo_data)
        # Build a boolean mask: True = keep, False = falls inside an exclusion window
        keep_mask = pd.Series(True, index=meteo_data.index)
        for i, (es, ee) in enumerate(exclusion_windows, start=1):
            in_excl = (meteo_data['datetime'] >= es) & (meteo_data['datetime'] <= ee)
            n_excl  = in_excl.sum()
            keep_mask = keep_mask & ~in_excl
            print(f"  Exclusion window {i} ({es} â€“ {ee}): {n_excl} record(s) removed.")
        meteo_data = meteo_data.loc[keep_mask]
        print(
            f"Exclusion windows applied. Records before: {before_excl}, "
            f"after: {len(meteo_data)} "
            f"({before_excl - len(meteo_data)} excluded in total)."
        )

    if enable_time_window == 1 or exclusion_windows:
        # Rebuild filtered_timestamps from meteo_data surviving all filters
        filtered_timestamps = set()
        date_parts = meteo_data['datetime'].dt.strftime('%Y%m%d').values
        time_parts = meteo_data['datetime'].dt.strftime('%H%M').values
        for date_str, time_str in zip(date_parts, time_parts):
            filtered_timestamps.add(date_str + '-' + time_str)
        
        # Match the set of filtered timestamps to the timestamps in cospectra
        # file names
        filtered_cospectra_files = []
        for filen in all_cospectra_files:
            filename = filen.split('\\')[-1]
            timestamp_match = re.search(r'(\d{8}-\d{4})_', filename)
            if timestamp_match:
                file_timestamp = timestamp_match.group(1)
                if file_timestamp in filtered_timestamps:
                    filtered_cospectra_files.append(filen)
        
        if len(filtered_cospectra_files) == 0:
            raise FREQCORReadError('No cospectra files found for the specified time window')
                
        print(f"  Found {len(filtered_cospectra_files)} cospectra files matching the time window")
    
    if vitale_qc_flags == 1:
        # Apply filtering to flag_data as well:
        # to do: do not reset index and keep datetime indexing when possible
        flag_data = (
            flag_data
            .set_index('datetime')
            .loc[meteo_data['datetime']]
            .reset_index()
        )
    # reset meteo_data index for later consistency and correct indexing
    # to do : remove integer-indexing and keep datetime indexing when possible
    meteo_data = meteo_data.reset_index(drop=True)
    
    # Use the filtered cospectra files sorted in alphabetical order
    total_to_process = filtered_cospectra_files
    total_to_process.sort()
    nspec = len(total_to_process)
    
    # Verify the number of selected (co)spectra
    if vitale_qc_flags == 1:
        if not (meteo_data['datetime'] == flag_data['datetime']).all():
            raise FREQCORReadError(
                "Timestamp misalignment detected between meteo_data and Vitale QC flag data."
            )
        if len(meteo_data) != len(flag_data) or len(meteo_data) != nspec:
            raise FREQCORReadError(
                'The number of binned_cospectra files must match the number of rows in the full_output and Vitale qc flags files. '
                'If using a time window, ensure there are cospectra files for the specified period.'
            )
    else:
        if len(meteo_data) != nspec:
            raise FREQCORReadError(
                'The number of binned_cospectra files must match the number of rows in the full_output file. '
                'If using a time window, ensure there are cospectra files for the specified period.'
            )

    # Infer number of bins directly from the first cospectra file
    if nspec == 0:
        raise FREQCORReadError('No cospectra files found to process')
    try:
        with open(total_to_process[0], mode='r') as file:
            first_file_data = pd.read_csv(file, header=11, na_values="-9999", low_memory=False)
    except Exception as e:
        raise FREQCORReadError(f"Could not read first cospectra file to infer bin count: {total_to_process[0]} ({e})")

    jfbin_inferred = int(first_file_data.shape[0])
    if jfbin_inferred <= 0:
        raise FREQCORReadError(f"Invalid inferred bin count from first cospectra file: {jfbin_inferred}")

    # Number of frequencies per bin
    numfreq = pd.DataFrame(np.zeros(shape=(jfbin_inferred, nspec)))
    # Natural frequency of each bin
    freqm = pd.DataFrame(np.zeros(shape=(jfbin_inferred, nspec)))
    # Normalised frequency of each bin
    freqn = pd.DataFrame(np.zeros(shape=(jfbin_inferred, nspec)))
    # Ideal (co)spectral data from file (sensible heat)
    Icosbrut = pd.DataFrame(np.zeros(shape=(jfbin_inferred, nspec)))
    # Real (co)spectral data from file (gas-specific)
    Rcosbrut = pd.DataFrame(np.zeros(shape=(jfbin_inferred, nspec)))
    # Ideal (co)spectral data from file (sensible heat) - for correction factor
    # computation
    Ccosbrut = pd.DataFrame(np.zeros(shape=(jfbin_inferred, nspec)))
    
    # Process files sequentially
    print(f"  Reading {len(total_to_process)} cospectra files...")
    for k, filen in enumerate(total_to_process):
        try:
            with open(filen, mode='r') as file:
                file_data = pd.read_csv(file, header=11, na_values="-9999", low_memory=False)

            if int(file_data.shape[0]) != jfbin_inferred:
                raise FREQCORReadError(
                    f"Inconsistent number of bins across cospectra files. "
                    f"First file has {jfbin_inferred} rows, but file '{os.path.basename(filen)}' has {int(file_data.shape[0])} rows. "
                    f"Please re-export cospectra with consistent binning."
                )
            if sps == 1:
                Icosbrut.iloc[:,k]=file_data.loc[:,'f_nat*cospec(w_ts)/cov(w_ts)'] 
                if gss == 1:
                    Rcosbrut.iloc[:,k]=file_data.loc[:,'f_nat*cospec(w_co2)/cov(w_co2)']   
                elif gss == 2:
                    Rcosbrut.iloc[:,k]=file_data.loc[:,'f_nat*cospec(w_h2o)/cov(w_h2o)'] 
                elif gss == 3:
                    try:
                        Rcosbrut.iloc[:,k] = file_data.loc[:,'f_nat*cospec(w_o3)/cov(w_o3)']
                    except KeyError:
                        print(f"Warning: Column for O3 not found in file {filen}")
                elif gss == 4:
                    Rcosbrut.iloc[:,k]=file_data.loc[:,'f_nat*cospec(w_ch4)/cov(w_ch4)']
                elif gss == 5:
                    Rcosbrut.iloc[:,k]=file_data.loc[:,'f_nat*cospec(w_n2o)/cov(w_n2o)']
            else:
                Icosbrut.iloc[:,k]=file_data.loc[:,'f_nat*spec(ts)/var(ts)']    
                if gss == 1:
                    Rcosbrut.iloc[:,k]=file_data.loc[:,'f_nat*spec(co2)/var(co2)']     
                elif gss == 2:
                    Rcosbrut.iloc[:,k]=file_data.loc[:,'f_nat*spec(h2o)/var(h2o)']    
                elif gss == 3:
                    try:
                        Rcosbrut.iloc[:,k] = file_data.loc[:,'f_nat*spec(o3)/var(o3)']
                    except KeyError:
                        print(f"Warning: Column for O3 not found in file {filen}")
                elif gss == 4:
                    Rcosbrut.iloc[:,k]=file_data.loc[:,'f_nat*spec(ch4)/var(ch4)']
                elif gss == 5:
                    Rcosbrut.iloc[:,k]=file_data.loc[:,'f_nat*spec(n2o)/var(n2o)']
            
            
            Ccosbrut.iloc[:,k]=file_data.loc[:,'f_nat*cospec(w_ts)/cov(w_ts)']    
            
            # Name the columns with datetime information extracted from file name
            filename = os.path.basename(filen)
            timestamp_match = re.search(r'(\d{8}-\d{4})', filename)
            if timestamp_match:
                file_timestamp = timestamp_match.group(1)
            else:
                raise FREQCORReadError(f"Could not extract timestamp from filename: {filename}. Expected format YYYYMMDD-HHMM (e.g., 20220628-1430). Check that your cospectra files follow the standard EddyPro naming convention.")
            numfreq.rename(columns={k: file_timestamp}, inplace=True)
            freqm.rename(columns={k: file_timestamp}, inplace=True)
            freqn.rename(columns={k: file_timestamp}, inplace=True)
            Icosbrut.rename(columns={k: file_timestamp}, inplace=True)
            Rcosbrut.rename(columns={k: file_timestamp}, inplace=True)
            Ccosbrut.rename(columns={k: file_timestamp}, inplace=True)
            
            # Fill numfreq with the number of frequencies per bin
            numfreq.iloc[:, k] = file_data.iloc[:, 0]
            
            # Fill freqm and freqn with the natural and normal frequency values
            freqm.iloc[:, k] = file_data.iloc[:, 1]
            freqn.iloc[:, k] = file_data.iloc[:, 2]
            
            # Check if all files have the same number of bins. If not, it could
            # lead to mismatch issues in the rest of the procedure. 
            # Warning printed.
            if k > 0:
                if not freqm.iloc[:, k].equals(freqm.iloc[:, 0]):
                    print(f'Warning: file {file_timestamp} has a different number of bins')
                if not numfreq.iloc[:, k].equals(numfreq.iloc[:, 0]):
                    print(f'File {file_timestamp} has a different number of frequencies per bin')

            # Show progress
            progress_interval = max(1, min(10, int(len(total_to_process) * 0.05)))
            if (k+1) % progress_interval == 0 or k+1 == len(total_to_process):
                progress_percent = int((k+1) / len(total_to_process) * 100)
                width = 30
                left = width * progress_percent // 100
                right = width - left
                sys.stdout.write(
                    f'\rReading cospectra files: [{("#" * left) + (" " * right)}] {k+1}/{len(total_to_process)} ({progress_percent}%)'
                )
                sys.stdout.flush()
        except Exception as e:
            print(f"\nError reading file {filen}: {e}")
            continue

    print()

    # Prepare variables
    # Filter out columns where WS, Zeta, or WD is missing
    valid_ts = ~meteo_data[['wind_speed', '(z-d)/L', 'wind_dir']].isna().any(axis=1)
    meteo_data = meteo_data.loc[valid_ts,:]
    if vitale_qc_flags == 1:
        flag_data = flag_data.loc[valid_ts,:]
    
    # match boolean index to (co)spectra dataframes column names
    valid_ts.index = Icosbrut.columns
    Icosbrut = Icosbrut.loc[:, valid_ts]
    Rcosbrut = Rcosbrut.loc[:, valid_ts]
    Ccosbrut = Ccosbrut.loc[:, valid_ts]
    freqn = freqn.loc[:, valid_ts]
    freqm = freqm.loc[:, valid_ts]
    numfreq = numfreq.loc[:, valid_ts]

    nspec = Icosbrut.shape[1]
    
    
    # Suppression of missing frequencies lines.
    # Assumption that all files have the same frequency bins (checked above).
    # Therefore, the natural frequency dataframe contains identical columns and
    # the first column can be used to represent all the others.
    freq = freqm.copy()
    freq = freq.iloc[:,0]
    
    valid_rows = (numfreq.iloc[:,0] != 0)
    freqm = freqm.loc[valid_rows,:]
    freq = freq.loc[valid_rows]
    Icosbrut = Icosbrut.loc[valid_rows,:]
    Rcosbrut = Rcosbrut.loc[valid_rows,:]
    Ccosbrut = Ccosbrut.loc[valid_rows,:]
    freqn = freqn.loc[valid_rows]
    if (freqn <= 0).any().any():
        raise FREQCORReadError('Error: freqn contains non-positive values. Normalized frequency and (z-d) cannot be < 0. Please check your data.')
    numfreq = numfreq.loc[valid_rows,:]
    
    # Reset the indexes so to have a continuous integer series of indexes
    freqm = freqm.reset_index(drop=True)    
    freq = freq.reset_index(drop=True)
    Icosbrut = Icosbrut.reset_index(drop=True)    
    Rcosbrut = Rcosbrut.reset_index(drop=True)    
    Ccosbrut = Ccosbrut.reset_index(drop=True)    
    
    freqn = freqn.reset_index(drop=True)
    numfreq = numfreq.reset_index(drop=True)
    # Define one variable to represent the overall weighting of bins
    freqcount = numfreq.iloc[:,0]
    # Number of available bins after cleaning
    nfreq = freqm.shape[0]
    
    # Division of (co)spectral density by frequency (because Eddypro files are 
    # in format f_nat*(co)spec/(co)var) 
    try:
        if freqm.empty:
            print("Warning: freqm is empty. Cannot perform division.")
        elif Icosbrut.empty or Rcosbrut.empty or Ccosbrut.empty:
            print("Warning: One or more spectral dataframes are empty.")
        else:
            Icon_raw = Icosbrut.div(freqm, fill_value=0)
            Rcon_raw = Rcosbrut.div(freqm, fill_value=0)
            Icosp_raw = Ccosbrut.div(freqm, fill_value=0)
    except Exception as e:
        print(f"Error during division operation: {e}")
        Icon_raw = pd.DataFrame(index=Icosbrut.index, columns=Icosbrut.columns)
        Rcon_raw = pd.DataFrame(index=Rcosbrut.index, columns=Rcosbrut.columns)
        Icosp_raw = pd.DataFrame(index=Ccosbrut.index, columns=Ccosbrut.columns)
    
    meteo_data = meteo_data.reset_index(drop=True)
    if vitale_qc_flags == 1:
        flag_data = flag_data.reset_index(drop=True)
    
    # Rename variables
    H=meteo_data.loc[:,'H'].copy()                                   # Sensible heat
    if gss == 1 or gss == 2 :
        Fc=(meteo_data.loc[:,'un_co2_flux'].copy()).rename('Fc')     # CO2 Flux (uncorrected)
        FcEP=(meteo_data.loc[:,'co2_flux'].copy()).rename('FcEP')    # CO2 Flux (EddyPro corrected)        
        LE=(meteo_data.loc[:,'un_LE'].copy()).rename('LE')           # Latent heat (uncorrected)
        LEEP=(meteo_data.loc[:,'LE'].copy()).rename('LEEP')          # Latent heat (EP corrected)
        RH=meteo_data.loc[:,'RH'].copy()                             # Relative humidity
        meteo_df = pd.concat([H,Fc,FcEP,LE,LEEP,RH],axis=1)
    elif gss==3 : 
        Fc=(meteo_data.loc[:,'un_o3_flux'].copy()).rename('Fc')      # O3 flux (uncorrected)
        FcEP=(meteo_data.loc[:,'o3_flux'].copy()).rename('FcEP')     # O3 flux (EP corrected) 
        meteo_df = pd.concat([H,Fc,FcEP],axis=1)
    elif gss==4 : 
        Fc=(meteo_data.loc[:,'un_ch4_flux'].copy()).rename('Fc')      # CH4 flux (uncorrected)
        FcEP=(meteo_data.loc[:,'ch4_flux'].copy()).rename('FcEP')     # CH4 flux (EP corrected) 
        meteo_df = pd.concat([H,Fc,FcEP],axis=1)
    elif gss==5 : 
        Fc=(meteo_data.loc[:,'un_n2o_flux'].copy()).rename('Fc')      # N2O flux (uncorrected)
        FcEP=(meteo_data.loc[:,'n2o_flux'].copy()).rename('FcEP')     # N2O flux (EP corrected) 
        meteo_df = pd.concat([H,Fc,FcEP],axis=1)
    
    WS=meteo_data.loc[:,'wind_speed'].copy()                   # Wind speed
    WD=meteo_data.loc[:,'wind_dir'].copy()                     # Wind direction
    Ustar=meteo_data.loc[:,'u*'].copy()                        # Friction velocity
    Zeta=meteo_data.loc[:,'(z-d)/L'].copy()                    # Stability parameter
    FlagVM_sp=meteo_data.loc[:,'spikes_hf'].copy()             # Vickers and Mahrt flag spikes
    FlagVM_do=meteo_data.loc[:,'drop_out_hf'].copy()           # Vickers and Mahrt flag drop out
    FlagVM_skh=meteo_data.loc[:,'skewness_kurtosis_hf'].copy() # Vickers and Mahrt flag skewness hard
    FlagVM_dh=meteo_data.loc[:,'discontinuities_hf'].copy()    # Vickers and Mahrt flag discontinuitites hard
    
    # Get global flags for later filtering:
        # When Vitale outputs are available, the global flag corresponds to the
        # _DATA_FLAG variable for each flux
        # When the output from EddyPro is used, the global flag is computed with
        # a specific function to put together results from different VM tests
    if vitale_qc_flags == 1 :
        FlagF_H = flag_data.loc[:,'H_SSITC_TEST'].copy()
        # Both "w" and "t_sonic" quality flags are included in the "H_DATA_FLAG"
        # to be consistent with the current code structure both variables are 
        # defined in the same way.
        FlagVM_w = flag_data.loc[:,'H_DATA_FLAG'].copy()
        FlagVM_T = flag_data.loc[:,'H_DATA_FLAG'].copy()
        flag_wd = flag_data.loc[:,'WSECT_FLAG'].copy()
        if gss==1:
            FlagF_g = flag_data.loc[:,'FC_SSITC_TEST'].copy()
            FlagVM_g = flag_data.loc[:,'NEE_DATA_FLAG'].copy()
        elif gss==2:
            FlagF_g = flag_data.loc[:,'LE_SSITC_TEST'].copy()
            FlagVM_g = flag_data.loc[:,'LE_DATA_FLAG'].copy()                 
        # do they exist for other gas?
        # "VM" flags are expected to be arrays in the rest of the code: convert
        FlagVM_w = np.array(FlagVM_w)
        FlagVM_T = np.array(FlagVM_T)
        FlagVM_g = np.array(FlagVM_g)
    else:
        FlagF_H = meteo_data.loc[:,'qc_H'].copy()
        # flag_wd absent if vitale outputs are unavailable: create empty Series
        # for compatiblity with the rest of the routine
        flag_wd = FlagF_H.copy()
        flag_wd[:] = np.nan
            
        # Foken flags
        if gss==1:
            FlagF_g=meteo_data.loc[:,'qc_co2_flux']          
        elif gss==2:
            FlagF_g=meteo_data.loc[:,'qc_LE']                  
        elif gss==3:
            FlagF_g=meteo_data.loc[:,'qc_o3_flux']
        elif gss==4:
            FlagF_g=meteo_data.loc[:,'qc_ch4_flux']
        elif gss==5:
            FlagF_g=meteo_data.loc[:,'qc_n2o_flux']
        
        
        # Create a global quality flag 
        # for VM test
        NumberDataFLAGF = len(meteo_data)
        FlagVM_w, FlagVM_T, FlagVM_g = FREQCOR_VM_flag(NumberDataFLAGF, FlagVM_sp, FlagVM_do, FlagVM_skh, FlagVM_dh, gss)
        
    #%% Read Massman coefficients 
    # from file when available
    try:
        massman_coef = pd.read_csv(massmanpath, header=0, comment='#')
        massman_coef = massman_coef.set_index('var')
        massman_coef = massman_coef.loc[:,jsite]
        
        A0 = [float(x) for x in massman_coef['A0'].split(',')]
        kf0 = [float(x) for x in massman_coef['kf0'].split(',')]
        mu = [float(x) for x in massman_coef['mu'].split(',')]
    
        massman_coef = [A0, kf0, mu]
    except:
        massman_coef=np.array([np.nan]*3)
    
    #%% Save the processed data 
    # to files for faster reloading if requested
    if 'PROCESSED_DATA' in config and int(config['PROCESSED_DATA']['enable_saving']) == 1:
        try:
            output_dir = os.path.join(output_path, f'{output_pattern}_intermediate_data')
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            if 'use_timestamp' in config['PROCESSED_DATA'] and int(config['PROCESSED_DATA']['use_timestamp']) == 1:
                timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
                file_prefix = f"{jsite}_{timestamp}"
            else:
                file_prefix = f"{jsite}"
                
            print(f"Saving processed data to {output_dir} for faster reloading...")
            
            # Save DataFrames
            Icon_raw.to_pickle(os.path.join(output_dir, f"{file_prefix}_Icon_raw.pkl"))
            Rcon_raw.to_pickle(os.path.join(output_dir, f"{file_prefix}_Rcon_raw.pkl"))
            Icosp_raw.to_pickle(os.path.join(output_dir, f"{file_prefix}_Icosp_raw.pkl"))
            freqn.to_pickle(os.path.join(output_dir, f"{file_prefix}_freqn.pkl"))
            meteo_df.to_pickle(os.path.join(output_dir, f"{file_prefix}_meteo_df.pkl"))
            
            # Save Series
            freq.to_pickle(os.path.join(output_dir, f"{file_prefix}_freq.pkl"))
            WS.to_pickle(os.path.join(output_dir, f"{file_prefix}_WS.pkl"))
            WD.to_pickle(os.path.join(output_dir, f"{file_prefix}_WD.pkl"))
            Ustar.to_pickle(os.path.join(output_dir, f"{file_prefix}_Ustar.pkl"))
            Zeta.to_pickle(os.path.join(output_dir, f"{file_prefix}_Zeta.pkl"))
            freqcount.to_pickle(os.path.join(output_dir, f"{file_prefix}_freqcount.pkl"))
            # Save arrays
            np.save(os.path.join(output_dir, f"{file_prefix}_massman_coef.npy"), massman_coef)
            
            if vitale_qc_flags == 1:
                FlagF_H.to_pickle(os.path.join(output_dir, f"{file_prefix}_FlagF_H_vitale.pkl"))
                FlagF_g.to_pickle(os.path.join(output_dir, f"{file_prefix}_FlagF_g_vitale.pkl"))
                flag_wd.to_pickle(os.path.join(output_dir, f"{file_prefix}_flag_wd_vitale.pkl"))
                
                # Save arrays 
                np.save(os.path.join(output_dir, f"{file_prefix}_FlagVM_w_vitale.npy"), FlagVM_w)
                np.save(os.path.join(output_dir, f"{file_prefix}_FlagVM_T_vitale.npy"), FlagVM_T)
                np.save(os.path.join(output_dir, f"{file_prefix}_FlagVM_g_vitale.npy"), FlagVM_g)
            else:
                
                FlagF_H.to_pickle(os.path.join(output_dir, f"{file_prefix}_FlagF_H.pkl"))
                FlagF_g.to_pickle(os.path.join(output_dir, f"{file_prefix}_FlagF_g.pkl"))
                flag_wd.to_pickle(os.path.join(output_dir, f"{file_prefix}_flag_wd.pkl"))
                
                # Save arrays
                np.save(os.path.join(output_dir, f"{file_prefix}_FlagVM_w.npy"), FlagVM_w)
                np.save(os.path.join(output_dir, f"{file_prefix}_FlagVM_T.npy"), FlagVM_T)
                np.save(os.path.join(output_dir, f"{file_prefix}_FlagVM_g.npy"), FlagVM_g)
                
            # Save scalar values in a JSON file
            scalar_values = {
                'nspec': nspec,
                'nfreq': nfreq
            }
            with open(os.path.join(output_dir, f"{file_prefix}_scalar_values.json"), 'w') as f:
                json.dump(scalar_values, f)
                
            print(f"Successfully saved processed data to {output_dir}")
            
            # Save a metadata file with information about what was saved
            # vitale files but only if the option is activated: 
                
                    # f"{file_prefix}_FlagF_H_vitale.pkl",
                    # f"{file_prefix}_FlagF_g_vitale.pkl",
                    # f"{file_prefix}_FlagVM_w_vitale.npy",
                    # f"{file_prefix}_FlagVM_T_vitale.npy",
                    # f"{file_prefix}_FlagVM_g_vitale.npy",
                    # f"{file_prefix}_flag_wd_vitale.npy",
            metadata = {
                'timestamp': datetime.datetime.now().isoformat(),
                'site': jsite,
                'gas': gss,
                 'files_saved': [
                     f"{file_prefix}_Icon_raw.pkl",
                     f"{file_prefix}_Rcon_raw.pkl",
                     f"{file_prefix}_Icosp_raw.pkl",
                     f"{file_prefix}_freqn.pkl",
                     f"{file_prefix}_meteo_df.pkl",
                     f"{file_prefix}_freq.pkl",
                    f"{file_prefix}_freqcount.pkl",
                    f"{file_prefix}_WS.pkl",
                    f"{file_prefix}_WD.pkl",
                    f"{file_prefix}_Ustar.pkl",
                    f"{file_prefix}_Zeta.pkl",
                    f"{file_prefix}_FlagF_H.pkl",
                    f"{file_prefix}_FlagF_g.pkl",
                    f"{file_prefix}_FlagVM_w.npy",
                    f"{file_prefix}_FlagVM_T.npy",
                    f"{file_prefix}_FlagVM_g.npy",
                    f"{file_prefix}_massman_coef.npy",
                    f"{file_prefix}_scalar_values.json"
                ]
            }
            with open(os.path.join(output_dir, f"{file_prefix}_metadata.json"), 'w') as f:
                json.dump(metadata, f, indent=4)
                
        except Exception as e:
            print(f"Error saving processed data: {e}")
    
    del freqm
    del FlagVM_do
    del FlagVM_skh
    del FlagVM_dh
    del Icosbrut
    del Rcosbrut
    return [nspec, Icon_raw, Rcon_raw, Icosp_raw, nfreq, freq, freqn, freqcount, meteo_df, WS, WD, Ustar, 
             Zeta, FlagF_H, FlagF_g, FlagVM_w, FlagVM_T, FlagVM_g, flag_wd, massman_coef]

