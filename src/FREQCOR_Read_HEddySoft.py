# -*- coding: utf-8 -*-
"""
Created on Oct 28 2026

@author: Matthias Cuntz

This subroutine is the FREQCOR_Read.py routine,
adapted for data coming from HEddySoft

Inputs :
      - HEddySoft (half-)hourly output files
Outputs :
      - (co)spectra matrices:
        icon(nfreq, nspec),
        rcon(nfreq, nspec) (either cospectra or spectra; depends on `sps`);
      - meteo and flux: x vectors(nspec);
      - General Flags: x vectors(nspec);

Comments :
    - filtering is applied :
        - for H, using Vitale 2020 tests and user limit on absolute value
          (from the config)
        - for N2O flux, using Vitale 2020 tests and user limit on absolute
          value (from the config)
        - for VOC flux, using only abs(flux) > percentile 90 of absolute
          value and removal of plume events
    - ! the mz choice is hard-coded, at the beginning of the routine !
    - some variables are missing in the HEddySoft output but are needed
      by FREQCOR for filtering steps. They are therefore created, filled with
      values that will end in no-filtering, and returned, to mimic FREQCOR_Read
TO DO:
    - target_mz should not be hard-coded but aim was to not modify the FREQCOR
      structure
"""
import os
import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from FREQCOR_Read_GEddySoft import calculate_freqcount


__all__ = ['FREQCOR_Read_HEddySoft']


def _filtering(HF5, HD5, KID, DIP, DDI, SST_M98, ITC_w, OOR_w):

    base_thresholds = {
        'HF5': 4,
        'HD5': 4,
        # ignore for the moment
        # 'KID': 50,
        'KID': 1000,
        # way to ground this threshold because not used by the carbon portal;
        # if used, should be 0.05
        # 'DIP': 0.05,
        'DIP': 0.0,
        'DDI': 3000,
        'ITC_w': 0.6,
        'SST_M98': 3,
    }

    # Initialize the flags as pandas Series with zeros (matching ITC_w's shape)
    FlagVM_w = np.zeros_like(ITC_w, dtype=int)  # Flag for ITC_w
    FlagVM_T = np.zeros_like(ITC_w, dtype=int)  # Not used in the task
    FlagVM_g = np.zeros_like(ITC_w, dtype=int)  # Flag for DIP and DDI
    FlagF_H = np.zeros_like(ITC_w, dtype=int)   # Not directly used
    FlagF_g = np.zeros_like(ITC_w, dtype=int)   # Flag for SST_M98

    # Apply filtering logic based on thresholds
    FlagVM_w = np.where(ITC_w > base_thresholds['ITC_w'], 1, FlagVM_w)
    FlagVM_w = np.where(OOR_w > 0, 1, FlagVM_w)

    FlagVM_g = np.where(HF5 > base_thresholds['HF5'], 1, 0)
    FlagVM_g = np.where(HD5 > base_thresholds['HD5'], 1, FlagVM_g)
    FlagVM_g = np.where(KID > base_thresholds['KID'], 1, FlagVM_g)
    FlagVM_g = np.where(DIP < base_thresholds['DIP'], 1, FlagVM_g)
    FlagVM_g = np.where(DDI > base_thresholds['DDI'], 1, FlagVM_g)

    FlagF_g = np.where(SST_M98 > base_thresholds['SST_M98'], 2, 0)

    # Return the flags as pandas Series so you can access them with .values
    return (pd.Series(FlagVM_w), pd.Series(FlagVM_T), pd.Series(FlagVM_g),
            pd.Series(FlagF_H), pd.Series(FlagF_g))


def FREQCOR_Read_HEddySoft(config, jsite):
    # %%  Initialisation
    apply_filtering = True
    plot = True

    # get spectral file

    input_path = config['IO']['input_path']
    input_file = config['IO']['binned_cosp']
    if not (input_file.startswith('/') or (input_file[1] == ':')):
        input_file = f'{input_path}/{input_file}'

    if not os.path.exists(input_file):
        raise ValueError(f'Input file does not exists: {input_file}')

    # %% loop over HEddySoft daily output files

    ispecies = int(config['PROCEDURE_OPTIONS']['gss'])
    if ispecies == 1:
        species = 'co2'
    elif ispecies == 2:
        species = 'h2o'
    else:
        raise ValueError('Species not implemented in HEddySoft')

    # read input file
    sdf = pd.DataFrame()
    idf = pd.DataFrame()
    freqm = pd.DataFrame()
    icon = pd.DataFrame()
    rcon = pd.DataFrame()

    ldf = pickle.load(open(input_file, 'rb'))

    # ['TIMESTAMP_END', 'wT', 'wind_speed', 'wind_direction', ...]
    sdf = ldf['scalar']

    # ipt_w, ipt_co2, ipt_h2o
    ipt_w = ldf['ipt_w']
    ipt_gas = ldf[f'ipt_{species}']
    idf = pd.DataFrame(
        {'OOR_w': ipt_w[11],
         'HF5': ipt_gas[4],
         'HD5': ipt_gas[6],
         'KID': ipt_gas[3],
         'DIP': ipt_gas[10],
         'DDI': ipt_gas[9]},
        index=sdf.index)
    idf.index.name = sdf.index.name

    # freq, cospec_wT_scaled, spec_co2_scaled, cospec_wc_co2_scaled,
    # spec_h2o_scaled, cospec_wc_h2o_scaled
    freqm = ldf['freq']
    icon = ldf['cospec_wT']

    if config['PROCEDURE_OPTIONS']['sps'] == '1':
        # Real cospectral data
        rcon = ldf[f'cospec_wc_{species}']
    elif config['PROCEDURE_OPTIONS']['sps'] == '2':
        # Real spectral data
        rcon = ldf[f'spec_{species}']

    sdf = pd.concat([sdf, idf], axis=1)
    sdf['Fc'] = sdf[f'F_{species}']
    sdf['FcEP'] = sdf['Fc']
    sdf['LEEP'] = sdf['LE']
    sdf['wind_dir'] = sdf['wind_direction']

    sdf = sdf.set_index(sdf.index.strftime('%Y%m%d-%H%M'))
    freqm = freqm.set_index(freqm.index.strftime('%Y%m%d-%H%M'))
    icon = icon.set_index(icon.index.strftime('%Y%m%d-%H%M'))
    rcon = rcon.set_index(rcon.index.strftime('%Y%m%d-%H%M'))

    # sdf = sdf.T
    freqm = freqm.T
    icon = icon.T
    rcon = rcon.T

    # %% Preparation of variables

    # Natural frequency of each bin
    freqn = np.divide(
        freqm * (float(config['EC_SETUP']['zmeas']) -
                 float(config['EC_SETUP']['disph'])),
        sdf['wind_speed'].values.reshape(1, -1))

    # create freqcount, the population in each bin.
    # Present in the EP outputs but not in the GEddySoft outputs
    freq_array = freqm.iloc[:, 0].to_numpy()
    freqcount = calculate_freqcount(freq_array, 1800.)

    # Division of (co)spectral density by frequency
    # (because GEddysoft output gives  f_nat*cospec/cov  or f_nat*spec/var)
    icon = icon.div(freqm)
    rcon = rcon.div(freqm)

    # %% Filtering and cleaning

    # based on raw data flagging
    if apply_filtering:
        FlagVM_w, FlagVM_T, FlagVM_g, FlagF_H, FlagF_g = _filtering(
            sdf['HF5'].to_numpy(),
            sdf['HD5'].to_numpy(),
            sdf['KID'].to_numpy(),
            sdf['DIP'].to_numpy(),
            sdf['DDI'].to_numpy(),
            sdf[f'steady_state_wc_{species}_Mahrt1998'].to_numpy(),
            sdf['itc_w'].to_numpy(),
            sdf['OOR_w'].to_numpy())
    else:
        FlagVM_w = pd.Series(np.zeros_like(sdf['itc_w'], dtype=int))
        FlagVM_T = pd.Series(np.zeros_like(sdf['itc_w'], dtype=int))
        FlagVM_g = pd.Series(np.zeros_like(sdf['itc_w'], dtype=int))
        FlagF_H = pd.Series(np.zeros_like(sdf['itc_w'], dtype=int))
        FlagF_g = pd.Series(np.zeros_like(sdf['itc_w'], dtype=int))

    # Remove freq bins fully filled with Nans
    idx = icon.notna().any(axis=1) & rcon.notna().any(axis=1)
    freqm = freqm[idx]
    freqn = freqn[idx]
    icon = icon[idx]
    rcon = rcon[idx]

    # filtering based on full NaNs or null data
    idx = (icon.notna().any(axis=0) &
           rcon.notna().any(axis=0) &
           (~sdf['wind_speed'].isnull().values))
    freqm = freqm.loc[:, idx]
    freqn = freqn.loc[:, idx]
    icon = icon.loc[:, idx]
    rcon = rcon.loc[:, idx]
    idx = idx.to_numpy()
    sdf = sdf.loc[idx, :]
    FlagVM_w = FlagVM_w[idx]
    FlagVM_T = FlagVM_T[idx]
    FlagVM_g = FlagVM_g[idx]
    FlagF_H = FlagF_H[idx]
    FlagF_g = FlagF_g[idx]
    nspec = sum(idx)

    # Reset the indexes so to have a continuous integer series of indexes
    freqm = freqm.reset_index(drop=True)
    freqn = freqn.reset_index(drop=True)
    icon = icon.reset_index(drop=True)
    rcon = rcon.reset_index(drop=True)
    # sdf = sdf.reset_index(drop=True)

    # %%

    # Series of natural frequencies, taken as the first column of the freqm
    # DataFrame where all the columns should contain the same values (see here above)
    freq = freqm.copy()
    freq = freq.iloc[:, 0]
    # Same for normalised frequencies
    # freqn = freqn.iloc[:,0]

    # Cleaned (co)spectra vector length (number of actually present frequency bins)
    nfreq = freq.shape[0]

    #%% Read Massman coefficients from file
    if config['IO']['massman_path']:
        massman_coef = pd.read_csv(config['IO']['massman_path'], header=0, comment='#')
        massman_coef = massman_coef.set_index('var')
        massman_coef = massman_coef.loc[:,jsite]

        A0 = [float(x) for x in massman_coef['A0'].split(',')]
        kf0 = [float(x) for x in massman_coef['kf0'].split(',')]
        mu = [float(x) for x in massman_coef['mu'].split(',')]

        massman_coef = [A0, kf0, mu]
    else:
        massman_coef = [np.nan, np.nan, np.nan]

    # create a copy of icon for further use
    icosp = icon.copy()

    if plot:
        # Use the first column as the time axis
        time = sdf.index

        fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(16, 8),
                                 sharex=True, sharey=False)

        # First subplot: H flux
        axes[0].scatter(time, sdf['H'], s=10, color='black', label='H flux')
        axes[0].scatter(
            time[FlagVM_w > 0],
            sdf.loc[(FlagVM_w > int(config['USER_LIMITS']['fvm'])).values, 'H'],
            s=10, color='green',
            label=f'FlagVM_w > {config["USER_LIMITS"]["fvm"]}',
            alpha=0.7)
        axes[0].scatter(
            time[FlagVM_T > 0],
            sdf.loc[(FlagVM_T > int(config['USER_LIMITS']['fvm'])).values, 'H'],
            s=10, color='purple',
            label=f'FlagVM_T > {config["USER_LIMITS"]["fvm"]}',
            alpha=0.7)
        axes[0].scatter(
            time[FlagF_H > 0],
            sdf.loc[(FlagF_H > int(config['USER_LIMITS']['ff'])).values, 'H'],
            s=10, color='red',
            label=f'FlagF_H > {config["USER_LIMITS"]["ff"]}',
            alpha=0.7)

        percentile_90 = np.percentile(abs(sdf.loc[:, 'H']), 90)
        axes[0].axhline(y=percentile_90, color='red', linestyle='--',
                        label='90th percentile')
        axes[0].axhline(y=-percentile_90, color='red', linestyle='--',
                        label='-90th percentile')
        axes[0].axhline(y=float(config['USER_LIMITS']['hlim']), color='green',
                        linestyle='--', label='user limit')
        axes[0].set_title('H flux')
        axes[0].legend()
        axes[0].set_xlabel('Time')
        axes[0].set_ylabel('H')

        # Second subplot: gas flux
        axes[1].scatter(time, sdf.loc[:, 'Fc'], s=10, color='black',
                        label='Fgas flux')
        axes[1].scatter(
            time[FlagF_H > 0],
            sdf.loc[(FlagF_H > int(config['USER_LIMITS']['ff'])).values, 'Fc'],
            s=10, color='red',
            label=f'FlagF_H > {config["USER_LIMITS"]["ff"]}',
            alpha=0.7)
        axes[1].scatter(
            time[FlagF_g > 0],
            sdf.loc[(FlagF_g > int(config['USER_LIMITS']['ff'])).values, 'Fc'],
            s=10, color='blue',
            label=f'FlagF_g > {config["USER_LIMITS"]["ff"]}',
            alpha=0.7)
        axes[1].scatter(
            time[FlagVM_w > 0],
            sdf.loc[(FlagVM_w > int(config['USER_LIMITS']['fvm'])).values, 'Fc'],
            s=10, color='green',
            label=f'FlagVM_w > {config["USER_LIMITS"]["fvm"]}',
            alpha=0.7)
        axes[1].scatter(
            time[FlagVM_T > 0],
            sdf.loc[(FlagVM_T > int(config['USER_LIMITS']['fvm'])).values, 'Fc'],
            s=10, color='purple',
            label=f'FlagVM_T > {config["USER_LIMITS"]["fvm"]}',
            alpha=0.7)
        axes[1].scatter(
            time[FlagVM_g > 0],
            sdf.loc[(FlagVM_g > int(config['USER_LIMITS']['fvm'])).values, 'Fc'],
            s=10, color='orange',
            label=f'FlagVM_g > {config["USER_LIMITS"]["fvm"]}',
            alpha=0.7)

        percentile_90 = np.percentile(abs(sdf.loc[:, 'Fc']), 90)
        axes[1].axhline(y=percentile_90, color='red', linestyle='--',
                        label='90th percentile')
        axes[1].axhline(y=-percentile_90, color='red', linestyle='--',
                        label='-90th percentile')
        axes[1].axhline(y=float(config['USER_LIMITS']['fclim']),
                        color='green', linestyle='--', label='user limit')
        axes[1].axhline(y=-float(config['USER_LIMITS']['fclim']),
                        color='green', linestyle='--', label='user limit')
        axes[1].set_title('Fgas flux')
        axes[1].legend()
        axes[1].set_xlabel('Time')
        axes[1].set_ylabel('Fgas')

        plt.tight_layout()

        if int(config['PROCEDURE_OPTIONS']['plot_save']) == 1:
            plt.savefig(config['IO']['output_path'] + '/flux_filtering.png')
            plt.close()
        else:
            plt.show()

    icon = icon.reset_index(drop=True)
    rcon = rcon.reset_index(drop=True)
    icosp = icosp.reset_index(drop=True)
    freq = freq.reset_index(drop=True)
    freqn = freqn.reset_index(drop=True)
    freqcount = freqcount.reset_index(drop=True)
    sdf = sdf.reset_index(drop=True)

    return [nspec, icon, rcon, icosp, nfreq, freq, freqn, freqcount,
            sdf, sdf['wind_speed'], sdf['wind_dir'], sdf['u_star'],
            sdf['z/L'], FlagF_H, FlagF_g, FlagVM_w, FlagVM_T, FlagVM_g,
            massman_coef]
