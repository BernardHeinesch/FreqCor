# -*- coding: utf-8 -*-

"""Output writer utilities for saving FREQCOR results to disk."""

import os
import pandas as pd
from datetime import datetime
from . import  as fqfct

def FREQCOR_write_outputs(gss, LUT_cof, LUTCF_u, LUTCF_s, outputpath, meteo_df, FccorL, FccorG,
                          freqn, WS, WD, Zeta, sps, CF_L, 
                          CF_G, CFHL, config, run_tag=None):
    """
    Write FREQCOR outputs to disk.

    This function writes:
    - 6_flux_out__{run_tag}.csv
    - 4_LUT_cof__{run_tag}.csv
    - 5_LUT_CF__{run_tag}.csv
    - 0_run_input__{run_tag}.ini

    Parameters
    ----------
    gss : int
        Gas species selector (e.g., 1 for CO2, 2 for H2O).
    LUT_cof : dict
        Cut-off frequency lookup tables (structure depends on gas species).
    LUTCF_u : dict
        Correction factor look-up table (unstable conditions).
    LUTCF_s : dict
        Correction factor look-up table (stable conditions).
    outputpath : str
        Output directory.
    meteo_df : pandas.DataFrame
        Half-hourly meteorological and flux data.
    FccorL : pandas.Series
        Corrected flux (Lorentz model).
    FccorG : pandas.Series
        Corrected flux (Gauss model).
    freqn : pandas.DataFrame
        Frequency index or timestamps (used to create output index).
    WS : pandas.Series
        Wind speed time series.
    WD : pandas.Series
        Wind direction time series.
    Zeta : pandas.Series
        Stability parameter time series.
    sps : int
        Approach selector (cospectral or spectral).
    CF_L : numpy.ndarray or pandas.Series
        Correction factor (Lorentz model).
    CF_G : numpy.ndarray or pandas.Series
        Correction factor (Gauss model).
    CFHL : numpy.ndarray or pandas.Series
        Correction factor for sensor separation according to Horst and Lenschow 2009.
    config : configparser.ConfigParser
        Run configuration to write back to disk.
    run_tag : str
        Run identifier used to generate output filenames.

    Returns
    -------
    None
        Outputs are written to disk.

    Raises
    ------
    ValueError
        If `run_tag` is missing.
    """
    if not run_tag:
        raise ValueError("run_tag is required for output naming")
    run_suffix = run_tag
    
    # --- Flux and meteorological data ---
    meteo_df.Fc.name = 'FCuncor'
    if gss==2:
        meteo_df.LE.name = 'LEuncor'
    dates = freqn.columns.to_series(name='Timestamp') # get a date column to save it to the excel file
    dates = pd.to_datetime(dates, format='%Y%m%d-%H%M')
    
    if gss == 2:
        FccorL = FccorL.rename('LEcorL')
        FccorG = FccorG.rename('LEcorG')
        meteo_cols = ['LE', 'LEEP', 'RH']
        out_cols = [WS, WD, Zeta, meteo_df.loc[:, meteo_cols], FccorL, FccorG, pd.Series(CF_L, name='CF_L'), pd.Series(CF_G, name='CF_G')]
    else:
        FccorL = FccorL.rename('FCcorL')
        FccorG = FccorG.rename('FCcorG')
        meteo_cols = ['Fc', 'FcEP']
        out_cols = [WS, WD, Zeta, meteo_df.loc[:, meteo_cols], FccorL, FccorG, pd.Series(CF_L, name='CF_L'), pd.Series(CF_G, name='CF_G')]
    
    df = pd.concat(out_cols, axis=1)
    if sps == 2:
        df = pd.concat([df, pd.Series(CFHL, name='CF_HL')], axis=1)
    df.index = dates

    if gss == 2:
        df = df.rename(
            columns={
                'wind_speed': 'wind_speed (m s-1)',
                'wind_dir': 'wind_dir (deg)',
                '(z-d)/L': '(z-d)/L (-)',
                'LE': 'LE (W m-2)',
                'LEEP': 'LEEP (W m-2)',
                'RH': 'RH (%)',
                'LEcorL': 'LEcorL (W m-2)',
                'LEcorG': 'LEcorG (W m-2)',
                'CF_L': 'CF_L (-)',
                'CF_G': 'CF_G (-)',
                'CF_HL': 'CF_HL (-)',
            }
        )
    else:
        df = df.rename(
            columns={
                'wind_speed': 'wind_speed (m s-1)',
                'wind_dir': 'wind_dir (deg)',
                '(z-d)/L': '(z-d)/L (-)',
                'Fc': 'Fc (Âµmol m-2 s-1)',
                'FcEP': 'FcEP (Âµmol m-2 s-1)',
                'FCcorL': 'FCcorL (Âµmol m-2 s-1)',
                'FCcorG': 'FCcorG (Âµmol m-2 s-1)',
                'CF_L': 'CF_L (-)',
                'CF_G': 'CF_G (-)',
                'CF_HL': 'CF_HL (-)',
            }
        )
    fileflx = os.path.join(outputpath, f"6_flux_out__{run_suffix}.csv")
    df.to_csv(fileflx, encoding='utf-8-sig')
    
    # --- cof LUT writing to file ---
    fqfct.write_coflut_to_csv(LUT_cof, gss, outputpath, run_tag=run_tag)
    

    # --- CF LUT writing to file ---
    cf_file = os.path.join(outputpath, f"5_LUT_CF__{run_suffix}.csv")
    with open(cf_file, "w", newline="") as f:
        f.write("unstable\n")
        fqfct.write_cflut_to_csv(f, LUTCF_u, gss)
        f.write("\n\n")
        f.write("stable\n")
        fqfct.write_cflut_to_csv(f, LUTCF_s, gss)
    

    # Config (ini) file: to have a record of the options used for the run
    config_path = os.path.join(outputpath, f"0_run_input__{run_suffix}.ini")
    with open(config_path, 'w') as configfile:
        config.write(configfile)
