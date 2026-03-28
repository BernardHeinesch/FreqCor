
# -*- coding: utf-8 -*-
"""
Created on Dec 28 2022

@author: Bernard Heinesch

This subroutine is the FREQCOR_Read.py routine, adapted for data coming from GEddySoft 

Inputs :
      - GEddySoft daily output files
Outputs : 
      -(co)spectra matrices: Icon_raw(nfreq, nspec), Rcon_raw(nfreq, nspec) (either cospectra or spectra; depends on `sps`);
      - meteo and flux: x vectors(nspec);
      - General Flags: x vectors(nspec);

Comments :
    - filtering is applied :
        - for H, using Vitale 2020 tests and user limit on absolute value (from the config)
        - for N2O flux, using Vitale 2020 tests and user limit on absolute value (from the config)
        - for VOC flux, using only abs(flux) > percentile 90 of absolute value and removal of plume events
    - ! the mz choice is hard-coded, at the beginning of the routine !
    - some variables are missing in the GEddySoft output but are needed by FREQCOR for filtering steps. 
      They are therefore created, filled with values that will end in no-filtering, and returned, to mimic
      FREQCOR_Read
TO DO:
    - target_mz should not be hard-coded but aim was to not modify the FREQCOR structure
"""

import glob
import numpy as np
from .FREQCOR_VM_flag import FREQCOR_VM_flag
import pandas as pd
import h5py
import math
import datetime
import re
import os
import matplotlib.pyplot as plt

def calculate_freqcount(freq_array, T=1800.0):
    """
    Calculate the frequency count (population in each bin) from a frequency array.
    
    Parameters:
        freq_array (numpy array): Array of frequencies (e.g., from an FFT analysis).
        T (float): Time window in seconds. Default is 1800.0 seconds (30 minutes).
    
    Returns:
        pd.DataFrame: A DataFrame with the frequency count for each bin.
    """
    # Round the last frequency to the upper integer for Nyquist
    nyquist_rounded = math.ceil(freq_array[-1])
    fs = 2 * nyquist_rounded  # Sampling frequency (2x Nyquist frequency)
    
    # Number of samples
    N = int(T * fs)
    
    # FFT frequency axis (positive frequencies only)
    fft_freqs = np.fft.rfftfreq(N, d=1/fs)
    
    # Number of bins = number of lines in freq_array
    NUM_BINS = len(freq_array)
    
    # Define bin edges (geometric mean between consecutive freq_array elements)
    bin_edges = np.zeros(NUM_BINS + 1)
    bin_edges[0] = 0  # First bin starts at 0
    for i in range(1, NUM_BINS):
        bin_edges[i] = np.sqrt(freq_array[i] * freq_array[i - 1])
    
    # Last bin: include all remaining FFT frequencies
    bin_edges[-1] = fft_freqs[-1] + 1e-12  # Tiny epsilon to ensure the last bin is inclusive
    
    # Count FFT frequencies per bin, ensuring at least 1 per bin
    n_freq = np.zeros(NUM_BINS, dtype=int)
    for i in range(NUM_BINS):
        if i < NUM_BINS - 1:
            idx = np.where((fft_freqs > bin_edges[i]) & (fft_freqs <= bin_edges[i + 1]))[0]
        else:
            idx = np.where(fft_freqs > bin_edges[i])[0]  # Last bin: all remaining frequencies
        n_freq[i] = len(idx) if len(idx) > 0 else 1
    
    # Build DataFrame with frequency counts
    freqcount = pd.DataFrame({'n_freq': n_freq})
    
    return freqcount

def filtering(HF5, HD5, KID, DIP, DDI, SST_M98, ITC_w, OOR_w):

    base_thresholds = {
        'HF5': 4,
        'HD5': 4,
        'KID': 50,
        'DIP': 0.05,  # way to ground this threshold because not used by the carbon portal (if used, should be 0.05) !
        'DDI': 3000,
        'ITC_w': 0.6,
        'SST_M98': 3,
    }

    # Initialize the flags as pandas Series with zeros (matching ITC_w's shape)
    FlagVM_w = pd.Series(np.zeros_like(ITC_w, dtype=int))  # Flag for ITC_w
    FlagVM_T = pd.Series(np.zeros_like(ITC_w, dtype=int))  # Not used in the task
    FlagVM_g = pd.Series(np.zeros_like(ITC_w, dtype=int))  # Flag for DIP and DDI
    FlagF_H = pd.Series(np.zeros_like(ITC_w, dtype=int))   # Not directly used
    FlagF_g = pd.Series(np.zeros_like(ITC_w, dtype=int))   # Flag for SST_M98
    
    # Apply filtering logic based on thresholds
    FlagVM_w = pd.Series(np.where(ITC_w > base_thresholds['ITC_w'], 1, FlagVM_w))
    FlagVM_w = pd.Series(np.where(OOR_w > 0, 1, FlagVM_w))
    
    FlagVM_g = pd.Series(np.where(HF5 > base_thresholds['HF5'], 1, 0))
    FlagVM_g = pd.Series(np.where(HD5 > base_thresholds['HD5'], 1, FlagVM_g))
    FlagVM_g = pd.Series(np.where(KID > base_thresholds['KID'], 1, FlagVM_g))
    FlagVM_g = pd.Series(np.where(DIP < base_thresholds['DIP'], 1, FlagVM_g))
    FlagVM_g = pd.Series(np.where(DDI > base_thresholds['DDI'], 1, FlagVM_g))
    
    FlagF_g = pd.Series(np.where(SST_M98 > base_thresholds['SST_M98'], 2, 0))
    
    # Return the flags as pandas Series so you can access them with .values
    return FlagVM_w, FlagVM_T, FlagVM_g, FlagF_H, FlagF_g

def FREQCOR_Read_GEddySoft(config, jsite):

    # %%  Initialisation

    # USER CHOICE !
    # mz to be selected (will choose the closest one)
    # 33.033 : proposed flux treshold of 0.1
    # 45.033 : proposed flux treshold of 0.01
    # 69.0699 : proposed flux treshold of 0.5
    # 137.132 : proposed flux treshold of 0.2
    target_mz = 44.013
    apply_filtering = True
    plot = True

    # Collect all candidate files
    total_to_process = glob.glob(config['IO']['input_path'] + config['IO']['binned_cosp'])

    # Filter files by date in filename if time window enabled
    if config['TIME_WINDOW']['enable_time_window'] == '1':
        # Read time window parameters
        start_datetime_str = config['TIME_WINDOW']['start_datetime']
        end_datetime_str = config['TIME_WINDOW']['end_datetime']
        # Parse datetime strings to datetime objects
        start_datetime = datetime.datetime.strptime(start_datetime_str, '%Y-%m-%d %H:%M')
        end_datetime = datetime.datetime.strptime(end_datetime_str, '%Y-%m-%d %H:%M')

        # Filter files by date in filename
        # inclusive date-range based on filename day (since filename has no time)
        start_day = start_datetime.date()
        end_day = end_datetime.date()
        date_pattern = re.compile(r'QCL_BE-Lon_fluxes_(\d{8})\.hdf5$', re.IGNORECASE)
        filtered = []
        skipped_badname = 0
        for filen in total_to_process:
            base = os.path.basename(filen)
            m = date_pattern.search(base)
            if not m:
                skipped_badname += 1
                continue
            yyyymmdd = m.group(1)
            try:
                file_day = datetime.datetime.strptime(yyyymmdd, '%Y%m%d').date()
            except ValueError:
                skipped_badname += 1
                continue
            if start_day <= file_day <= end_day:
                filtered.append(filen)
        total_to_process = sorted(filtered)
        print(f"Selected {len(total_to_process)} file(s) within window days [{start_day} .. {end_day}]")
        if skipped_badname:
            print(f"Skipped {skipped_badname} file(s) due to unexpected name/date format")

    # get total number of available (co)spectra
    nspec = 0
    for filen in total_to_process:
        with h5py.File(filen, 'r') as hdf5_f:
            if 'TRACER' in hdf5_f:
                name_key = hdf5_f['TRACER'].keys()
                if len(name_key):  # file with TRACER values
                    ds = hdf5_f['time']
                    nspec = nspec + ds.shape[0]
    print(str(nspec) + ' nspec available')

    # initialize vectors and matrices
    timestamp = pd.Series(dtype='float64')
    freqm = pd.DataFrame(dtype='float64')
    Icon_raw = pd.DataFrame(dtype='float64')
    Rcon_raw = pd.DataFrame(dtype='float64')
    H = pd.DataFrame(dtype='float64')
    Fc = pd.DataFrame(dtype='float64')
    FcEP = pd.DataFrame(dtype='float64')
    LE = pd.DataFrame(dtype='float64')
    RH = pd.DataFrame(dtype='float64')
    meteo_df = pd.DataFrame(dtype='float64')
    Rcon_raw = pd.DataFrame(dtype='float64')
    WS = pd.Series(dtype='float64', name='wind_speed')
    WD = pd.Series(dtype='float64', name='wind_dir')
    Ustar = pd.Series(dtype='float64', name='u*')
    Zeta = pd.Series(dtype='float64', name='(z-d)/L')
    OOR_w = pd.Series(dtype='float64', name='OOR_w')
    num_spikes_w = pd.Series(dtype='float64', name='spikes_w')
    HF5 = pd.Series(dtype='float64', name='HF5')
    HD5 = pd.Series(dtype='float64', name='HD5')
    KID = pd.Series(dtype='float64', name='KID')
    DIP = pd.Series(dtype='float64', name='DIP')
    DDI = pd.Series(dtype='float64', name='DDI')
    SST_M98 = pd.Series(dtype='float64', name='SST_M98')
    ITC_w = pd.Series(dtype='float64', name='ITC_w')

    # %% loop over GEddySoft daily output files

    # loop over the daily files in the given directory matching the name pattern
    k = 0
    for filen in total_to_process:
    
        with h5py.File(filen, 'r') as hdf5_f:

            if 'TRACER' in hdf5_f:
                name_key = hdf5_f['TRACER'].keys()
                if len(name_key):  # file with TRACER values
    
                    # find the mz closest to target_mz 
                    for x in list(name_key)[:-1]:
                        if abs(float(hdf5_f['TRACER'][x]['mz'][()])-target_mz) < 0.01:
                            target_channel = x
                            break
                        
                    # find the number of 1/2h 
                    ds = hdf5_f['time']
                    n_hh = ds.shape[0]
        
                    # Create all the tables that will contain the 1/2 (co)spectra or frequency and the FLUX, METEO and QUALITY FLAGS
                    # Each colum will correspond to one specific 1/2 hour
                    # Dimensions: max nb of bins x nb of (co)spectra
            
                    # timestamp = timestamp.append(pd.Series(np.transpose(np.array(hdf5_f['time']).astype('U17'))),ignore_index=True)
                    timestamp = pd.concat([timestamp, pd.Series(np.array(hdf5_f['time']).astype('U17'))], ignore_index=True)                    

                    # Natural frequency of each bin
                    freqm_d = pd.DataFrame(np.transpose(np.array(hdf5_f['freq'])))
                    freqm_d = freqm_d[np.repeat(freqm_d.columns.values,n_hh)]
                    freqm = pd.concat([freqm,freqm_d],axis=1,ignore_index=True)
            
                    # Ideal cospectral data from file (sensible heat)
                    Icon_raw = pd.concat([Icon_raw,pd.DataFrame(np.transpose(np.array(hdf5_f['MET']['cospec_wT_scaled'])))],axis=1,ignore_index=True)
                    if config['PROCEDURE_OPTIONS']['sps'] == '1':
                        # Real cospectral data from file (VOC)
                        Rcon_raw = pd.concat([Rcon_raw,pd.DataFrame(np.transpose(np.array(hdf5_f['TRACER'][target_channel]['cospec_scaled'])))],axis=1,ignore_index=True)
                    elif config['PROCEDURE_OPTIONS']['sps'] == '2':
                        # Real spectral data from file (VOC)
                        Rcon_raw = pd.concat([Rcon_raw,pd.DataFrame(np.transpose(np.array(hdf5_f['TRACER'][target_channel]['spec_scaled'])))],axis=1,ignore_index=True)

                    # Extract timestamps (assuming 1D array, adjust if necessary)
                    timestamps = np.array(hdf5_f['time'])
                    # Make sure timestamps length matches n_hh (number of rows in meteo_df_d)
                    assert len(timestamps) == n_hh, "Timestamps length does not match data rows!"
                    
                    # Create meteo_df_d with timestamps as index
                    meteo_df_d = pd.concat([
                        pd.DataFrame(np.zeros(shape=(n_hh,))),
                        pd.DataFrame(np.nan_to_num(np.transpose(np.array(hdf5_f['MET']['wT']) * 1206))),
                        pd.DataFrame(np.nan_to_num(np.transpose(np.array(hdf5_f['TRACER'][target_channel]['flux'])))),
                        pd.DataFrame(np.zeros(shape=(n_hh,))),
                        pd.DataFrame(np.zeros(shape=(n_hh,))),
                        pd.DataFrame(np.zeros(shape=(n_hh,)))
                    ], axis=1)
                    
                    meteo_df_d.iloc[:, 3:7] = 999
                    
                    meteo_df_d.iloc[:,0] = timestamps
                    
                    # Convert byte strings to strings, then to datetime
                    meteo_df_d.iloc[:,0] = pd.to_datetime([x.decode() if isinstance(x, bytes) else str(x) for x in meteo_df_d.iloc[:,0]],
                                                      format='%Y-%m-%d %H-%M-%S')
                    
                    # When concatenating, preserve index (do not use ignore_index=True)
                    meteo_df = pd.concat([meteo_df, meteo_df_d], axis=0,ignore_index=True)

                    WS = pd.concat([WS, pd.DataFrame(np.transpose(np.array(hdf5_f['MET']['wsh'])), columns=['wind_speed']).squeeze()], ignore_index=True)
                    WD = pd.concat([WD, pd.DataFrame(np.transpose(np.array(hdf5_f['MET']['wdir'])), columns=['wind_dir']).squeeze()], ignore_index=True)
                    Ustar = pd.concat([Ustar, pd.DataFrame(np.transpose(np.array(hdf5_f['MET']['ust'])), columns=['u*']).squeeze()], ignore_index=True)
                    Zeta = pd.concat([Zeta, pd.DataFrame(np.transpose(np.array(hdf5_f['MET']['zoL'])), columns=['(z-d)/L']).squeeze()], ignore_index=True)
                    OOR_w = pd.concat([OOR_w, pd.DataFrame(np.transpose(np.array(hdf5_f['MET']['qaqc']['IPT_w'][:,11])), columns=['OOR_w']).squeeze()], ignore_index=True)
                    num_spikes_w = pd.concat([num_spikes_w, pd.DataFrame(np.transpose(np.array(hdf5_f['MET']['qaqc']['num_spikes_w'])), columns=['num_spikes_w']).squeeze()], ignore_index=True)                    
                    HF5 = pd.concat([HF5, pd.DataFrame(np.transpose(np.array(hdf5_f['TRACER'][target_channel]['qaqc']['IPT'][:,4])), columns=['HF5']).squeeze()], ignore_index=True)                    
                    HD5 = pd.concat([HD5, pd.DataFrame(np.transpose(np.array(hdf5_f['TRACER'][target_channel]['qaqc']['IPT'][:,6])), columns=['HD5']).squeeze()], ignore_index=True)                    
                    KID = pd.concat([KID, pd.DataFrame(np.transpose(np.array(hdf5_f['TRACER'][target_channel]['qaqc']['IPT'][:,3])), columns=['KID']).squeeze()], ignore_index=True)                    
                    DIP = pd.concat([DIP, pd.DataFrame(np.transpose(np.array(hdf5_f['TRACER'][target_channel]['qaqc']['IPT'][:,10])), columns=['DIP']).squeeze()], ignore_index=True)                    
                    DDI = pd.concat([DDI, pd.DataFrame(np.transpose(np.array(hdf5_f['TRACER'][target_channel]['qaqc']['IPT'][:,9])), columns=['DDI']).squeeze()], ignore_index=True)                    
                    SST_M98 = pd.concat([SST_M98, pd.DataFrame(np.transpose(np.array(hdf5_f['TRACER'][target_channel]['qaqc']['SST_M98'])), columns=['SST_M98']).squeeze()], ignore_index=True)                    
                    ITC_w = pd.concat([ITC_w, pd.DataFrame(np.transpose(np.array(hdf5_f['MET']['qaqc']['ITC_w'])), columns=['ITC_w']).squeeze()], ignore_index=True)

                    k = k + n_hh
                
    # Rename the columns in the DataFrames with the info on the date and
    # time of the 1/2 hour to which it corresponds.
    timestamp = timestamp.replace({'-':''}, regex=True)
    timestamp = timestamp.replace({' ':'-'}, regex=True)
    freqm.columns = timestamp.iloc[:]
    Icon_raw.columns = timestamp.iloc[:]
    Rcon_raw.columns = timestamp.iloc[:]

    # %% Preparation of variables

    meteo_df.columns = ['timestamp', 'H','Fc', 'FcEP', 'LE', 'RH']

    # Natural frequency of each bin
    freqn = np.divide(freqm * (float(config['EC_SETUP']['zmeas']) - float(config['EC_SETUP']['disph'])), WS.values.reshape(1, -1))

    # freqm: your log-spaced bin frequencies (1D, from your DataFrame)
    freq_array = freqm.iloc[:, 0].to_numpy()
    
    # create freqcount, the population in each bin. Present in the EP outputs but not in the GEddySoft outputs
    freqcount = calculate_freqcount(freq_array, 1800.0)

    # Division of (co)spectral density by frequency (because GEddysoft output gives 
    # f_nat*cospec/cov or f_nat*spec/var) 
    Icon_raw=Icon_raw.div(freqm)
    Rcon_raw=Rcon_raw.div(freqm)

    # %% Filtering and cleaning

    # based on raw data flagging
    if apply_filtering:
        FlagVM_w, FlagVM_T, FlagVM_g, FlagF_H, FlagF_g = filtering(HF5, HD5, KID, DIP, DDI, SST_M98, ITC_w, OOR_w)

    # Remove freq bins fully filled with Nans
    idx = Icon_raw.notna().any(axis=1) & Rcon_raw.notna().any(axis=1)
    freqm = freqm[idx]
    freqn = freqn[idx]
    Icon_raw = Icon_raw[idx]
    Rcon_raw = Rcon_raw[idx]
    
    # filtering based on full NaNs or null data
    idx = (Icon_raw.notna().any(axis=0) & 
           Rcon_raw.notna().any(axis=0) & 
           ~WS.isnull().values)

    freqm = freqm.loc[:,idx]
    freqn = freqn.loc[:,idx]
    Icon_raw = Icon_raw.loc[:,idx]
    Rcon_raw = Rcon_raw.loc[:,idx]
    idx = idx.reset_index(drop=True)
    meteo_df = meteo_df.loc[idx,:]
    WS = WS[idx]
    WD = WD[idx]
    Ustar = Ustar[idx]
    Zeta = Zeta[idx]
    OOR_w = OOR_w[idx]
    num_spikes_w = num_spikes_w[idx]
    FlagVM_w = FlagVM_w[idx]
    FlagVM_T = FlagVM_T[idx]
    FlagVM_g = FlagVM_g[idx]
    FlagF_H = FlagF_H[idx]
    FlagF_g = FlagF_g[idx]

    nspec = sum(idx)
    
    # Reset the indexes so to have a continuous integer series of indexes
    freqm = freqm.reset_index(drop=True)
    freqn = freqn.reset_index(drop=True)
    Icon_raw = Icon_raw.reset_index(drop=True)    
    Rcon_raw = Rcon_raw.reset_index(drop=True)    
    meteo_df = meteo_df.reset_index(drop=True)
    WS = WS.reset_index(drop=True)
    WD = WD.reset_index(drop=True)
    Ustar = Ustar.reset_index(drop=True)
    Zeta = Zeta.reset_index(drop=True)

    #%% 
    
    # Series of natural frequencies, taken as the first column of the freqm
    # DataFrame where all the columns should contain the same values (see here above)
    freq = freqm.copy()
    freq = freq.iloc[:,0]
    # Same for normalised frequencies
    # freqn = freqn.iloc[:,0]
   
    # Cleaned (co)spectra vector length (number of actually present frequency bins)
    nfreq = freq.shape[0]

    #%% Read Massman coefficients from file
    massman_coef = pd.read_csv(config['IO']['massman_path'], header=0, comment='#')
    massman_coef = massman_coef.set_index('var')
    massman_coef = massman_coef.loc[:,jsite]
    
    A0 = [float(x) for x in massman_coef['A0'].split(',')]
    kf0 = [float(x) for x in massman_coef['kf0'].split(',')]
    mu = [float(x) for x in massman_coef['mu'].split(',')]

    massman_coef = [A0, kf0, mu]

    #%%
    # Suppression of useless intermediate variables
    del freqm
    del freqm_d
    del meteo_df_d
    # del FlagVM_sp
    # del FlagVM_do
    # del FlagVM_skh
    # del FlagVM_sks
    # del FlagVM_ds
    # del FlagVM_dh

    # create a copy of Icon_raw for further use
    Icosp_raw = Icon_raw.copy()
    
    if plot:    
        # Use the first column as the time axis
        time = meteo_df.iloc[:, 0]
    
        fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(16, 8), sharex=True, sharey=False)
    
        # First subplot: H flux
        axes[0].scatter(time, meteo_df.loc[:, 'H'], s=10, color='black', label='H flux')
        axes[0].scatter(time[FlagVM_w > 0], meteo_df.loc[FlagVM_w > int(config['USER_LIMITS']['fvm']), 'H'], s=10, color='green', label=f'FlagVM_w > {config["USER_LIMITS"]["fvm"]}', alpha=0.7)
        axes[0].scatter(time[FlagVM_T > 0], meteo_df.loc[FlagVM_T > int(config['USER_LIMITS']['fvm']), 'H'], s=10, color='purple', label=f'FlagVM_T > {config["USER_LIMITS"]["fvm"]}', alpha=0.7)
        axes[0].scatter(time[FlagF_H > 0], meteo_df.loc[FlagF_H > int(config['USER_LIMITS']['ff']), 'H'], s=10, color='red', label=f'FlagF_H > {config["USER_LIMITS"]["ff"]}', alpha=0.7)

        percentile_90 = np.percentile(abs(meteo_df.loc[:, 'H']), 90)
        axes[0].axhline(y=percentile_90, color='red', linestyle='--', label='90th percentile')
        axes[0].axhline(y=-percentile_90, color='red', linestyle='--', label='-90th percentile')
        axes[0].axhline(y=float(config['USER_LIMITS']['hlim']), color='green', linestyle='--', label='user limit')
        axes[0].set_title('H flux')
        axes[0].legend()
        axes[0].set_xlabel('Time')
        axes[0].set_ylabel('H')
    
        # Second subplot: gas flux
        axes[1].scatter(time, meteo_df.loc[:, 'Fc'], s=10, color='black', label='Fgas flux')

        axes[1].scatter(time[FlagF_H > 0], meteo_df.loc[FlagF_H > int(config['USER_LIMITS']['ff']), 'Fc'], s=10, color='red', label=f'FlagF_H > {config["USER_LIMITS"]["ff"]}', alpha=0.7)
        axes[1].scatter(time[FlagF_g > 0], meteo_df.loc[FlagF_g > int(config['USER_LIMITS']['ff']), 'Fc'], s=10, color='blue', label=f'FlagF_g > {config["USER_LIMITS"]["ff"]}', alpha=0.7)
        axes[1].scatter(time[FlagVM_w > 0], meteo_df.loc[FlagVM_w > int(config['USER_LIMITS']['fvm']), 'Fc'], s=10, color='green', label=f'FlagVM_w > {config["USER_LIMITS"]["fvm"]}', alpha=0.7)
        axes[1].scatter(time[FlagVM_T > 0], meteo_df.loc[FlagVM_T > int(config['USER_LIMITS']['fvm']), 'Fc'], s=10, color='purple', label=f'FlagVM_T > {config["USER_LIMITS"]["fvm"]}', alpha=0.7)
        axes[1].scatter(time[FlagVM_g > 0], meteo_df.loc[FlagVM_g > int(config['USER_LIMITS']['fvm']), 'Fc'], s=10, color='orange', label=f'FlagVM_g > {config["USER_LIMITS"]["fvm"]}', alpha=0.7)

        percentile_90 = np.percentile(abs(meteo_df.loc[:, 'Fc']), 90)
        axes[1].axhline(y=percentile_90, color='red', linestyle='--', label='90th percentile')
        axes[1].axhline(y=-percentile_90, color='red', linestyle='--', label='-90th percentile')
        axes[1].axhline(y=float(config['USER_LIMITS']['fclim']), color='green', linestyle='--', label='user limit')
        axes[1].axhline(y=-float(config['USER_LIMITS']['fclim']), color='green', linestyle='--', label='user limit')
        axes[1].set_title('Fgas flux')
        axes[1].legend()
        axes[1].set_xlabel('Time')
        axes[1].set_ylabel('Fgas')
    
        plt.tight_layout()

        if int(config['PROCEDURE_OPTIONS']['plot_save']) == 1:
            plt.savefig(config['IO']['output_path'] + '/flux filtering.png')
            plt.close()
        else:
            plt.show()
        
    # remove the timestamp column
    meteo_df = meteo_df.iloc[:, 1:]
    return [nspec, Icon_raw, Rcon_raw, Icosp_raw, nfreq, freq, freqn, freqcount, meteo_df, WS, WD,
            Ustar, Zeta, FlagF_H, FlagF_g, FlagVM_w, FlagVM_T, FlagVM_g, massman_coef]

