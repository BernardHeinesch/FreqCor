# -*- coding: utf-8 -*-
"""
Created on Dec 28 2022 - merged version Nov 10 2026

@author: Bernard Heinesch

Unified reader for GEddySoft daily HDF5 output files. It merges the former
``FREQCOR_Read_TOF`` and ``FREQCOR_Read_IRGA`` routines:

The source is chosen explicitly in the INI via
[IO].source_type:

- ``source_type = tof``: reads the ``TRACER`` group (PTR-TOF-MS data); the
  channel is selected by its ``mz`` value, closest to the required
  ``target_mz`` setting (tolerance 0.01).
- ``source_type = irga``: reads the ``IRGA`` group; the channel is selected
  by matching the channel ``name`` dataset (e.g. 'CO2 dry', 'H2O dry')
  against the gas key derived from ``gss``.

Inputs :
      - GEddySoft daily HDF5 output files (TRACER and/or IRGA groups)
Outputs :
      -(co)spectra matrices: Icon_raw(nfreq, nspec), Rcon_raw(nfreq, nspec) (either cospectra or spectra; depends on `sps`);
      - meteo and flux: x vectors(nspec);
      - General Flags: x vectors(nspec);

Comments :
    - meteo_df columns follow the EP convention:
        H    : sensible heat flux [Wm-2]          <- MET/wT * 1206
        Fc   : CO2 flux [umolm-2s-1] if a CO2 IRGA channel exists, else flux of
               the selected channel
        FcEP : placeholder (999)
        LE   : latent heat flux [Wm-2]            <- IRGA H2O flux [mmolm-2s-1] * 44.1
               (999 when no H2O channel is available)
        LEEP : copy of LE (GEddySoft provides a single flux; the column is
               required by FREQCOR_Sel_general for gss=2)
        RH   : relative humidity [%]              <- MET/RH_air (999 if absent)
    - filtering is applied using GEddySoft qaqc tests (HF5, HD5, KID, DIP, DDI,
      SST_M98 on the selected channel when its qaqc datasets exist; ITC_w,
      OOR_w on the sonic) through the `filtering` helper. When the channel has
      no qaqc group, no-filtering flags are created (former TOF behaviour).
    - TIME_WINDOW filtering is supported: files are kept/discarded according to
      the YYYYMMDD date in their filename.
"""

import glob
import numpy as np
import pandas as pd
import h5py
import math
import datetime
import re
import os
import matplotlib.pyplot as plt
from FREQCOR_Read_EP import FREQCORReadError

# Conversion mmol H2O m-2 s-1 -> W m-2 (latent heat flux)
# lambda_v ~= 2.45e6 J/kg * 18.02e-3 kg/mol = 44.1 J/mmol
LE_PER_MMOL_H2O = 44.1

# Gas code (PROCEDURE_OPTIONS.gss) -> substring looked for in IRGA channel names
GAS_CHANNEL_NAMES = {1: 'co2', 2: 'h2o', 3: 'o3', 4: 'ch4', 5: 'n2o'}

# Tolerance on mz matching for TRACER channels
MZ_TOLERANCE = 0.01


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


def _find_named_channel(group, gas_key):
    """Return the key of the channel whose ``name`` dataset contains ``gas_key``.

    Only h5py groups are considered (scalar datasets such as ``cf_lpf`` are
    skipped). Returns None if no channel matches.
    """
    for x in group.keys():
        grp = group[x]
        if not isinstance(grp, h5py.Group) or 'name' not in grp:
            continue
        nm = grp['name'][()]
        nm = nm.decode() if isinstance(nm, bytes) else str(nm)
        if gas_key in nm.lower():
            return x
    return None


def _find_mz_channel(group, target_mz, tol=MZ_TOLERANCE):
    """Return the key of the channel whose ``mz`` is closest to ``target_mz``.

    Only channels within ``tol`` of the target are considered. Returns None if
    no channel matches.
    """
    target_channel = None
    best = tol
    for x in group.keys():
        grp = group[x]
        if not isinstance(grp, h5py.Group) or 'mz' not in grp:
            continue
        diff = abs(float(grp['mz'][()]) - target_mz)
        if diff < best:
            target_channel = x
            best = diff
    return target_channel


def FREQCOR_Read_GEddySoft(config, jsite):

    # %%  Initialisation

    apply_filtering = True
    plot = int(config['PROCEDURE_OPTIONS'].get('plot_aux', '0')) == 1

    sps = int(config['PROCEDURE_OPTIONS']['sps'])
    gss = int(config['PROCEDURE_OPTIONS']['gss'])
    gas_key = GAS_CHANNEL_NAMES.get(gss, 'co2')

    # Source selection (INI [IO]): 'tof' -> TRACER group, 'irga' -> IRGA group
    src_type = config['IO'].get('source_type', '').strip().lower()
    if src_type in ('tof', 'tracer'):
        src_group = 'TRACER'
    elif src_type == 'irga':
        src_group = 'IRGA'
    else:
        raise FREQCORReadError(
            "[IO].source_type must be 'irga' or 'tof' "
            f"(got '{src_type}')"
        )

    # Target mz (INI [IO]): required in 'tof' mode, ignored in 'irga' mode.
    target_mz_str = config['IO'].get('target_mz', '').strip()
    if src_group == 'TRACER' and not target_mz_str:
        raise FREQCORReadError(
            "[IO].target_mz is required when source_type = 'tof'"
        )
    target_mz = float(target_mz_str) if target_mz_str else None

    # Map config names to variables used in the script
    massmanpath = config['IO']['massman_path']
    zmeas = float(config['EC_SETUP']['zmeas'])
    disph = float(config['EC_SETUP']['disph'])

    # Collect all candidate files
    total_to_process = glob.glob(config['IO']['binned_cosp'])

    # Filter files by date in filename if time window enabled
    if 'TIME_WINDOW' in config and config['TIME_WINDOW'].get('enable_time_window', '0') == '1':
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
        date_pattern = re.compile(r'(\d{8})\.hdf5$', re.IGNORECASE)
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
            if src_group in hdf5_f:
                ds = hdf5_f['time']
                nspec = nspec + ds.shape[0]
    print(str(nspec) + ' nspec available')

    # initialize vectors and matrices
    timestamp = pd.Series(dtype='float64')
    freqm = pd.DataFrame(dtype='float64')
    Icon_raw = pd.DataFrame(dtype='float64')
    Rcon_raw = pd.DataFrame(dtype='float64')
    meteo_df = pd.DataFrame(dtype='float64')
    WS = pd.Series(dtype='float64', name='wind_speed')
    WD = pd.Series(dtype='float64', name='wind_dir')
    Ustar = pd.Series(dtype='float64', name='u*')
    Zeta = pd.Series(dtype='float64', name='(z-d)/L')
    OOR_w = pd.Series(dtype='float64', name='OOR_w')
    num_spikes_w = pd.Series(dtype='float64', name='num_spikes_w')
    HF5 = pd.Series(dtype='float64', name='HF5')
    HD5 = pd.Series(dtype='float64', name='HD5')
    KID = pd.Series(dtype='float64', name='KID')
    DIP = pd.Series(dtype='float64', name='DIP')
    DDI = pd.Series(dtype='float64', name='DDI')
    SST_M98 = pd.Series(dtype='float64', name='SST_M98')
    ITC_w = pd.Series(dtype='float64', name='ITC_w')
    has_channel_qaqc = False

    # %% loop over GEddySoft daily output files

    # loop over the daily files in the given directory matching the name pattern
    for filen in total_to_process:

        with h5py.File(filen, 'r') as hdf5_f:

            if src_group not in hdf5_f:
                print(f"Warning: no '{src_group}' group in "
                      f"{os.path.basename(filen)} - file skipped")
                continue

            if src_group == 'TRACER':
                # find the mz closest to target_mz (required in tof mode)
                target_channel = _find_mz_channel(hdf5_f['TRACER'], target_mz)
                if target_channel is None:
                    print(f"Warning: no TRACER channel with mz within "
                          f"{MZ_TOLERANCE} of {target_mz} in "
                          f"{os.path.basename(filen)} - file skipped")
                    continue
                fc_source = target_channel
                h2o_channel = None
            else:  # 'IRGA'
                # find the IRGA channel matching the gas selected by gss
                target_channel = _find_named_channel(hdf5_f['IRGA'], gas_key)
                if target_channel is None:
                    print(f"Warning: no IRGA channel matching '{gas_key}' "
                          f"(gss={gss}) in {os.path.basename(filen)} - file skipped")
                    continue
                # CO2 channel feeds Fc when present (EP convention);
                # H2O channel feeds LE (mmol m-2 s-1 -> W m-2)
                co2_channel = _find_named_channel(hdf5_f['IRGA'], 'co2')
                h2o_channel = _find_named_channel(hdf5_f['IRGA'], 'h2o')
                fc_source = co2_channel if co2_channel is not None else target_channel

            chan = hdf5_f[src_group][target_channel]

            # find the number of 1/2h
            ds = hdf5_f['time']
            n_hh = ds.shape[0]

            # Create all the tables that will contain the 1/2 (co)spectra or frequency and the FLUX, METEO and QUALITY FLAGS
            # Each colum will correspond to one specific 1/2 hour
            # Dimensions: max nb of bins x nb of (co)spectra

            timestamp = pd.concat([timestamp, pd.Series(np.array(hdf5_f['time']).astype('U17'))], ignore_index=True)

            # Natural frequency of each bin
            freqm_d = pd.DataFrame(np.transpose(np.array(hdf5_f['freq'])))
            freqm_d = freqm_d[np.repeat(freqm_d.columns.values, n_hh)]
            freqm = pd.concat([freqm, freqm_d], axis=1, ignore_index=True)

            # Ideal cospectral data from file (sensible heat)
            Icon_raw = pd.concat([Icon_raw, pd.DataFrame(np.transpose(np.array(hdf5_f['MET']['cospec_wT_scaled'])))], axis=1, ignore_index=True)
            if sps == 2 and 'spec_scaled' in chan:
                # Real spectral data from file (gas)
                Rcon_raw = pd.concat([Rcon_raw, pd.DataFrame(np.transpose(np.array(chan['spec_scaled'])))], axis=1, ignore_index=True)
            else:
                # Real cospectral data from file (gas)
                Rcon_raw = pd.concat([Rcon_raw, pd.DataFrame(np.transpose(np.array(chan['cospec_scaled'])))], axis=1, ignore_index=True)

            # Extract timestamps
            timestamps = np.array(hdf5_f['time'])
            # Make sure timestamps length matches n_hh (number of rows in meteo_df_d)
            assert len(timestamps) == n_hh, "Timestamps length does not match data rows!"

            # Flux/meteo columns (EP convention)
            flux_fc = np.nan_to_num(np.transpose(np.array(hdf5_f[src_group][fc_source]['flux'])))
            if h2o_channel is not None:
                le_vals = np.nan_to_num(np.transpose(np.array(hdf5_f['IRGA'][h2o_channel]['flux']))) * LE_PER_MMOL_H2O
            else:
                le_vals = np.full(n_hh, 999.0)
            if 'RH_air' in hdf5_f['MET']:
                rh_vals = np.nan_to_num(np.transpose(np.array(hdf5_f['MET']['RH_air'])))
            else:
                rh_vals = np.full(n_hh, 999.0)

            # Create meteo_df_d with timestamps as first column
            meteo_df_d = pd.concat([
                pd.DataFrame(np.zeros(shape=(n_hh,))),                                   # timestamp placeholder
                pd.DataFrame(np.nan_to_num(np.transpose(np.array(hdf5_f['MET']['wT']) * 1206))),  # H
                pd.DataFrame(flux_fc),                                                    # Fc
                pd.DataFrame(np.full(n_hh, 999.0)),                                       # FcEP placeholder
                pd.DataFrame(le_vals),                                                    # LE
                pd.DataFrame(le_vals),                                                    # LEEP: copy of LE
                pd.DataFrame(rh_vals),                                                    # RH
            ], axis=1)

            meteo_df_d.iloc[:, 0] = timestamps

            # Convert byte strings to strings, then to datetime
            meteo_df_d.iloc[:, 0] = pd.to_datetime([x.decode() if isinstance(x, bytes) else str(x) for x in meteo_df_d.iloc[:, 0]],
                                                   format='%Y-%m-%d %H-%M-%S')

            # When concatenating, preserve index (do not use ignore_index=True)
            meteo_df = pd.concat([meteo_df, meteo_df_d], axis=0, ignore_index=True)

            WS = pd.concat([WS, pd.DataFrame(np.transpose(np.array(hdf5_f['MET']['wsh'])), columns=['wind_speed']).squeeze()], ignore_index=True)
            WD = pd.concat([WD, pd.DataFrame(np.transpose(np.array(hdf5_f['MET']['wdir'])), columns=['wind_dir']).squeeze()], ignore_index=True)
            Ustar = pd.concat([Ustar, pd.DataFrame(np.transpose(np.array(hdf5_f['MET']['ust'])), columns=['u*']).squeeze()], ignore_index=True)
            Zeta = pd.concat([Zeta, pd.DataFrame(np.transpose(np.array(hdf5_f['MET']['zoL'])), columns=['(z-d)/L']).squeeze()], ignore_index=True)
            OOR_w = pd.concat([OOR_w, pd.DataFrame(np.transpose(np.array(hdf5_f['MET']['qaqc']['IPT_w'][:, 11])), columns=['OOR_w']).squeeze()], ignore_index=True)
            num_spikes_w = pd.concat([num_spikes_w, pd.DataFrame(np.transpose(np.array(hdf5_f['MET']['qaqc']['num_spikes_w'])), columns=['num_spikes_w']).squeeze()], ignore_index=True)
            ITC_w = pd.concat([ITC_w, pd.DataFrame(np.transpose(np.array(hdf5_f['MET']['qaqc']['ITC_w'])), columns=['ITC_w']).squeeze()], ignore_index=True)

            # Channel-level qaqc (IRGA and possibly TRACER channels)
            if 'qaqc' in chan and 'IPT' in chan['qaqc'] and 'SST_M98' in chan['qaqc']:
                has_channel_qaqc = True
                HF5 = pd.concat([HF5, pd.DataFrame(np.transpose(np.array(chan['qaqc']['IPT'][:, 4])), columns=['HF5']).squeeze()], ignore_index=True)
                HD5 = pd.concat([HD5, pd.DataFrame(np.transpose(np.array(chan['qaqc']['IPT'][:, 6])), columns=['HD5']).squeeze()], ignore_index=True)
                KID = pd.concat([KID, pd.DataFrame(np.transpose(np.array(chan['qaqc']['IPT'][:, 3])), columns=['KID']).squeeze()], ignore_index=True)
                DIP = pd.concat([DIP, pd.DataFrame(np.transpose(np.array(chan['qaqc']['IPT'][:, 10])), columns=['DIP']).squeeze()], ignore_index=True)
                DDI = pd.concat([DDI, pd.DataFrame(np.transpose(np.array(chan['qaqc']['IPT'][:, 9])), columns=['DDI']).squeeze()], ignore_index=True)
                SST_M98 = pd.concat([SST_M98, pd.DataFrame(np.transpose(np.array(chan['qaqc']['SST_M98'])), columns=['SST_M98']).squeeze()], ignore_index=True)
            else:
                # No channel qaqc: no-filtering values
                HF5 = pd.concat([HF5, pd.Series(np.zeros(n_hh))], ignore_index=True)
                HD5 = pd.concat([HD5, pd.Series(np.zeros(n_hh))], ignore_index=True)
                KID = pd.concat([KID, pd.Series(np.zeros(n_hh))], ignore_index=True)
                DIP = pd.concat([DIP, pd.Series(np.full(n_hh, 1.0))], ignore_index=True)  # DIP < 0.05 flags; keep at 1 to not flag
                DDI = pd.concat([DDI, pd.Series(np.zeros(n_hh))], ignore_index=True)
                SST_M98 = pd.concat([SST_M98, pd.Series(np.zeros(n_hh))], ignore_index=True)

    if meteo_df.empty:
        raise FREQCORReadError(
            "No usable data found: no file matching binned_cosp contains a "
            f"'{src_group}' group with a channel matching "
            f"{'mz=' + str(target_mz) if target_mz is not None else repr(gas_key) + ' (gss=' + str(gss) + ')'}"
        )

    # Rename the columns in the DataFrames with the info on the date and
    # time of the 1/2 hour to which it corresponds.
    timestamp = timestamp.replace({'-': ''}, regex=True)
    timestamp = timestamp.replace({' ': '-'}, regex=True)
    freqm.columns = timestamp.iloc[:]
    Icon_raw.columns = timestamp.iloc[:]
    Rcon_raw.columns = timestamp.iloc[:]

    # %% Preparation of variables

    meteo_df.columns = ['timestamp', 'H', 'Fc', 'FcEP', 'LE', 'LEEP', 'RH']

    # Natural frequency of each bin
    freqn = np.divide(freqm * (zmeas - disph), WS.values.reshape(1, -1))

    # freqm: your log-spaced bin frequencies (1D, from your DataFrame)
    freq_array = freqm.iloc[:, 0].to_numpy()

    # create freqcount, the population in each bin. Present in the EP outputs but not in the GEddySoft outputs
    freqcount = calculate_freqcount(freq_array, 1800.0)

    # Division of (co)spectral density by frequency (because GEddysoft output gives
    # f_nat*cospec/cov or f_nat*spec/var)
    Icon_raw = Icon_raw.div(freqm)
    Rcon_raw = Rcon_raw.div(freqm)

    # %% Filtering and cleaning

    # based on raw data flagging
    if apply_filtering:
        FlagVM_w, FlagVM_T, FlagVM_g, FlagF_H, FlagF_g = filtering(HF5, HD5, KID, DIP, DDI, SST_M98, ITC_w, OOR_w)
    if not has_channel_qaqc:
        print("Note: no channel-level qaqc found - using no-filtering flags "
              "for HF5/HD5/KID/DIP/DDI/SST_M98")

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

    freqm = freqm.loc[:, idx]
    freqn = freqn.loc[:, idx]
    Icon_raw = Icon_raw.loc[:, idx]
    Rcon_raw = Rcon_raw.loc[:, idx]
    idx = idx.reset_index(drop=True)
    meteo_df = meteo_df.loc[idx, :]
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

    # %%

    # Series of natural frequencies, taken as the first column of the freqm
    # DataFrame where all the columns should contain the same values (see here above)
    freq = freqm.copy()
    freq = freq.iloc[:, 0]
    # Same for normalised frequencies
    # freqn = freqn.iloc[:,0]

    # Cleaned (co)spectra vector length (number of actually present frequency bins)
    nfreq = freq.shape[0]

    # %% Read Massman coefficients
    # from file when available (optional: used for reference plots only)
    try:
        massman_coef = pd.read_csv(massmanpath, header=0, comment='#')
        massman_coef = massman_coef.set_index('var')
        massman_coef = massman_coef.loc[:, jsite]

        A0 = [float(x) for x in massman_coef['A0'].split(',')]
        kf0 = [float(x) for x in massman_coef['kf0'].split(',')]
        mu = [float(x) for x in massman_coef['mu'].split(',')]

        massman_coef = [A0, kf0, mu]
    except:
        massman_coef = np.array([np.nan]*3)

    # %%
    # Suppression of useless intermediate variables
    del freqm
    del freqm_d
    del meteo_df_d

    # create a copy of Icon_raw for further use
    Icosp_raw = Icon_raw.copy()

    if plot:
        # Use the first column as the time axis
        time = meteo_df.iloc[:, 0]
        # Column plotted as "gas flux": LE for H2O (gss=2), Fc otherwise
        gas_col = 'LE' if gss == 2 else 'Fc'
        # User limits used for the reference lines (from [USER_LIMITS])
        hlim = float(config['USER_LIMITS']['hl_cof_min'])
        fclim = float(config['USER_LIMITS']['lel_cof_min'] if gss == 2 else config['USER_LIMITS']['fcl_cof_min'])

        fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(16, 8), sharex=True, sharey=False)

        # First subplot: H flux
        axes[0].scatter(time, meteo_df.loc[:, 'H'], s=10, color='black', label='H flux')
        axes[0].scatter(time[FlagVM_w > 0], meteo_df.loc[FlagVM_w > int(config['USER_LIMITS']['fvm']), 'H'], s=10, color='green', label=f'FlagVM_w > {config["USER_LIMITS"]["fvm"]}', alpha=0.7)
        axes[0].scatter(time[FlagVM_T > 0], meteo_df.loc[FlagVM_T > int(config['USER_LIMITS']['fvm']), 'H'], s=10, color='purple', label=f'FlagVM_T > {config["USER_LIMITS"]["fvm"]}', alpha=0.7)
        axes[0].scatter(time[FlagF_H > 0], meteo_df.loc[FlagF_H > int(config['USER_LIMITS']['ff']), 'H'], s=10, color='red', label=f'FlagF_H > {config["USER_LIMITS"]["ff"]}', alpha=0.7)

        percentile_90 = np.percentile(abs(meteo_df.loc[:, 'H']), 90)
        axes[0].axhline(y=percentile_90, color='red', linestyle='--', label='90th percentile')
        axes[0].axhline(y=-percentile_90, color='red', linestyle='--', label='-90th percentile')
        axes[0].axhline(y=hlim, color='green', linestyle='--', label='user limit')
        axes[0].set_title('H flux')
        axes[0].legend()
        axes[0].set_xlabel('Time')
        axes[0].set_ylabel('H')

        # Second subplot: gas flux
        axes[1].scatter(time, meteo_df.loc[:, gas_col], s=10, color='black', label=f'{gas_col} flux')

        axes[1].scatter(time[FlagF_H > 0], meteo_df.loc[FlagF_H > int(config['USER_LIMITS']['ff']), gas_col], s=10, color='red', label=f'FlagF_H > {config["USER_LIMITS"]["ff"]}', alpha=0.7)
        axes[1].scatter(time[FlagF_g > 0], meteo_df.loc[FlagF_g > int(config['USER_LIMITS']['ff']), gas_col], s=10, color='blue', label=f'FlagF_g > {config["USER_LIMITS"]["ff"]}', alpha=0.7)
        axes[1].scatter(time[FlagVM_w > 0], meteo_df.loc[FlagVM_w > int(config['USER_LIMITS']['fvm']), gas_col], s=10, color='green', label=f'FlagVM_w > {config["USER_LIMITS"]["fvm"]}', alpha=0.7)
        axes[1].scatter(time[FlagVM_T > 0], meteo_df.loc[FlagVM_T > int(config['USER_LIMITS']['fvm']), gas_col], s=10, color='purple', label=f'FlagVM_T > {config["USER_LIMITS"]["fvm"]}', alpha=0.7)
        axes[1].scatter(time[FlagVM_g > 0], meteo_df.loc[FlagVM_g > int(config['USER_LIMITS']['fvm']), gas_col], s=10, color='orange', label=f'FlagVM_g > {config["USER_LIMITS"]["fvm"]}', alpha=0.7)

        percentile_90 = np.percentile(abs(meteo_df.loc[:, gas_col]), 90)
        axes[1].axhline(y=percentile_90, color='red', linestyle='--', label='90th percentile')
        axes[1].axhline(y=-percentile_90, color='red', linestyle='--', label='-90th percentile')
        axes[1].axhline(y=fclim, color='green', linestyle='--', label='user limit')
        axes[1].axhline(y=-fclim, color='green', linestyle='--', label='user limit')
        axes[1].set_title(f'{gas_col} flux')
        axes[1].legend()
        axes[1].set_xlabel('Time')
        axes[1].set_ylabel(gas_col)

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
