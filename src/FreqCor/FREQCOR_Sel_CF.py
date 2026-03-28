# -*- coding: utf-8 -*-

"""Selection utilities for filtering (co)spectra used in correction factor computation."""

import numpy as np
import pandas as pd
from . import  as fqplt

# Custom exception class for FREQCOR_Sel_CF errors
class FREQCORSelCFError(Exception):
    """Exception raised for errors in `FREQCOR_Sel_CF`."""
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)



def FREQCOR_Sel_CF(Icos_CF, freq, sts, meteo_df, Hlcfu, Hlcfs, 
                   Hlim_max, FlagVM_T, FlagF_H, FlagVM_w, FVM, FF, outputpath, plot, run_tag=None):
    """
    Select cospectra used for correction factor (CF) computation.

    Tests for all (co)spectra
        Stationarity test (tau, H)
        VM tests (T)
        Min value (H, u*)
    Additional test for cospectra
        VM test (w) 

    Parameters
    ----------
    Icos_CF : pandas.DataFrame
        Ideal cospectra for CF (stability-filtered).
    freq : pandas.Series
        Natural frequencies (Hz). Used for plotting only.
    sts : int
        Stability class selector (1: unstable, 2: stable).
    meteo_df : pandas.DataFrame
        Half-hourly meteorological and flux data (used for H thresholding).
    Hlcfu : float
        Minimum sensible heat flux threshold for unstable conditions.
    Hlcfs : float
        Maximum sensible heat flux threshold for stable conditions.
    FlagVM_T : pandas.Series
        Vickers & Mahrt (1997) global flag for temperature.
    FlagF_H : pandas.Series
        Mauder & Foken (2004) quality flag for sensible heat flux.
    FlagVM_w : pandas.Series
        Vickers & Mahrt (1997) global flag for vertical wind speed.
    FVM : float
        Threshold applied to Vickers & Mahrt flags.
    FF : float
        Threshold applied to Mauder & Foken flags.
    outputpath : str
        Output directory (used for plots).
    plot : list[int]
        Plotting options selector.
    run_tag : str, optional
        Run identifier used to generate plot filenames.

    Returns
    -------
    Icos_CF : pandas.DataFrame
        Quality-filtered cospectra used for CF computation.

    Raises
    ------
    ValueError
        If `run_tag` is missing while plotting is enabled.
    FREQCORSelCFError
        If all cospectra are filtered out.
    """
    
    stability = ['unst','st']        
    Icos_CF_plot = Icos_CF.copy()
    Icos_CF_filt = Icos_CF.copy()

    fluxes_df = meteo_df.copy()
    fluxes_df.index = Icos_CF_filt.columns 
    if sts == 1:
    # Unstable: select H >= Hlcfu and H <= Hlim_max
        flux_mask = (
                fluxes_df['H'].between(Hlcfu, Hlim_max)
                ).fillna(False)
        
    else:
        # Stable: select H <= Hlcfs and >= -Hlim_max
        flux_mask = (
                fluxes_df['H'].between(-Hlim_max, Hlcfs)
                ).fillna(False)
        
    flagged = (
        (FlagVM_T > FVM) |
        (FlagVM_w > FVM) |
        (FlagF_H > FF) 
    )
    
    flag_mask = pd.Series(~flagged,index=Icos_CF_filt.columns)
    
    valid_mask = flux_mask & flag_mask
    Icos_CF_filt.iloc[:, ~valid_mask] = np.nan
    
    mask_history = [flux_mask,
             flag_mask]
    step_names = ['Flux thresholds',
                  'VM and FF flags']
    
    if plot[2] == 1:
        if not run_tag:
            raise ValueError("run_tag is required for output naming")
        run_suffix = run_tag
        step_plot_info = {
            "Flux thresholds": {
                "flux_limits": [Hlcfu, Hlim_max]
            }
        }
        
        fqplt.plot_stepwise_filtering(
            mask_history=mask_history,
            step_names=step_names,
            fluxes_df=fluxes_df,
            flux_col='H',
            cospectra_df=Icos_CF_plot,
            freq=freq,
            outputpath=outputpath,
            plot=plot,
            init_valid_mask=None,
            step_plot_info=step_plot_info,
            filename=f"3_filtering_CF_H_{stability[sts-1]}__{run_suffix}.png",
        )
       
    if Icos_CF_filt.isnull().sum().sum() == Icos_CF_filt.size:
        raise FREQCORSelCFError(f'No (co)spectra after CF filtering ({stability[sts-1]}able conditions)')
    
    return Icos_CF_filt
