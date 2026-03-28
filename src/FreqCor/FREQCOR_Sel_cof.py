# -*- coding: utf-8 -*-

"""Selection utilities for filtering (co)spectra used in cut-off frequency computation."""

import numpy as np
import pandas as pd
from . import  as fqplt
from .FREQCOR_functions import spectral_outlier_mask_iqr

# Custom exception class for FREQCOR_Sel_cof errors
class FREQCORSelCofError(Exception):
    """Exception raised for errors in `FREQCOR_Sel_cof`."""
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)

def FREQCOR_Sel_cof(Icon_sel, Rcon_sel, freq, gss, sps, meteo_df, Hlim, 
                    Fclim, LElim, Hlim_max, Fclim_max, LElim_max, FlagVM_T,
                    FlagVM_g, FlagVM_w, FlagF_H, FlagF_g, FVM, FF, outputpath,
                    plot, cospec_abs_limit, run_tag=None):
    """
    Select (co)spectra used for cut-off frequency computation.

    Tests for spectra all gases

    - Stationarity test (tau, H)
    - VM tests (T)
    - Min value (H, u*)

    Tests for spectra specific to gas

    - Stationarity test (Fc or LE)
    - VM tests (c or h)
    - Min value (Fc or LE)

    Additional test for CO2 and heat cospectra

    - VM test (w)

    Absolute values (depend on stability conditions and gas type):

    Minimum value for H (FH), Fc or LE (FG)

    Parameters
    ----------
    Icon_sel, Rcon_sel : pandas.DataFrame
        Ideal content (sel) and real content (sel) (co)spectral data (either cospectra or spectra; depends on `sps`).
    freq : pandas.Series
        Natural frequencies (Hz). Used for outlier filtering and plotting.
    gss : int
        Gas species selector.
    sps : int
        Approach selector (cospectral or spectral).
    meteo_df : pandas.DataFrame
        Half-hourly meteorological and flux data.
    Hlim, Fclim, LElim : float
        Flux thresholds used for selection (H + gas flux or LE depending on `gss`).
    FlagVM_T, FlagVM_g, FlagVM_w : pandas.Series
        Vickers & Mahrt (1997) global flags for temperature, gas concentration, and wind.
    FlagF_H, FlagF_g : pandas.Series
        Mauder & Foken (2004) quality flags for sensible heat flux and gas flux.
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
    Icon_cof : pandas.DataFrame
        Quality-filtered ideal (co)spectral data [nfreq, nspec] (either cospectra or spectra; depends on `sps`).
    Rcon_cof : pandas.DataFrame
        Quality-filtered real (co)spectral data [nfreq, nspec] (either cospectra or spectra; depends on `sps`).

    Raises
    ------
    ValueError
        If `run_tag` is missing while plotting is enabled.
    FREQCORSelCofError
        If all (co)spectra are filtered out.
    """

    # (co)spectral absolute limit: if a single point of a hh (co)spectra if outside + or - this treshold, 
    # the whole hh is set to nan
    Icon_cof = Icon_sel.copy()
    Rcon_cof = Rcon_sel.copy()
   
    Icos_plot = Icon_sel.copy()
    Rcos_plot = Rcon_sel.copy()
    
    # if not indexed with 1/2 h names yet:
    fluxes_df = meteo_df.copy()
    fluxes_df.index = Icon_cof.columns 
    
    # Vectorized selection based on flux limits
    if gss != 2:
        flux_mask = (
                    (fluxes_df['H'].abs().between(abs(Hlim), abs(Hlim_max)))
                    & (fluxes_df['Fc'].abs().between(abs(Fclim), abs(Fclim_max)))
                ).fillna(False)
        flux_col_gas = 'Fc'
    else:
        flux_mask = (
                    (fluxes_df['H'].abs().between(abs(Hlim), abs(Hlim_max)))
                    & (fluxes_df['LE'].abs().between(abs(LElim), abs(LElim_max)))
                ).fillna(False)
        flux_col_gas = 'LE'
    
    # Vectorized flag filtering (Vickers & Mahrt; Foken)
    flagged = (
        (FlagVM_T > FVM) |
        (FlagVM_g > FVM) |
        (FlagF_H > FF) |
        (FlagF_g > FF)
    )
    if sps == 1:
        flagged = flagged | (FlagVM_w > FVM)
    
    flag_mask = pd.Series(~flagged,index=Icon_cof.columns)
    
    # Absolue values outliers mask
    mask_i = Icon_cof.ge(-cospec_abs_limit).all() & Icon_cof.le(cospec_abs_limit).all()
    mask_r = Rcon_cof.ge(-cospec_abs_limit).all() & Rcon_cof.le(cospec_abs_limit).all()
    
    # Overall validity mask
    valid_mask = flux_mask & flag_mask & mask_i & mask_r
    
    Icon_cof.loc[:, ~valid_mask] = np.nan
    Rcon_cof.loc[:, ~valid_mask] = np.nan

    keep_mask_i, stats_i = spectral_outlier_mask_iqr(Icon_cof.copy(), 
                                          freq,
                                       fmin=0.005,
                                       fmax=0.1,
                                       k=2.5,
                                       min_fraction=0.15,
                                       )
    keep_mask_r, stats_r = spectral_outlier_mask_iqr(Rcon_cof.copy(), 
                                          freq,
                                       fmin=0.005,
                                       fmax=0.1,
                                       k=2.5,
                                       min_fraction=0.15,
                                       )
    iqr_mask = keep_mask_i & keep_mask_r
    
    Icon_cof.loc[:,~iqr_mask] = np.nan 
    Rcon_cof.loc[:,~iqr_mask] = np.nan
    
    mask_history = [flux_mask,
             flag_mask,
             mask_i & mask_r,
             iqr_mask]
    step_names = ['Flux thresholds',
                  'VM and FF flags',
                  'Absolute (co)spectral value',
                  'IQR over frequency band']
    
    if plot[2] == 1:
        if not run_tag:
            raise ValueError("run_tag is required for output naming")
        run_suffix = run_tag
        
        step_plot_info = {
            "Flux thresholds": {
                "flux_limits": [Hlim, Hlim_max]
            },
            "Absolute (co)spectral value": {
                "cospec_abs_limit": cospec_abs_limit
            },
            "IQR over frequency band": {
                "iqr_stats": stats_i
            }
        }
    
        
        fqplt.plot_stepwise_filtering(
            mask_history=mask_history,
            step_names=step_names,
            fluxes_df=fluxes_df,
            flux_col='H',
            cospectra_df=Icos_plot,
            freq=freq,
            outputpath=outputpath,
            plot=plot,
            init_valid_mask=None,
            step_plot_info=step_plot_info,
            filename=f"3_filtering_cof_H__{run_suffix}.png",
        )
        
        step_plot_info = {
            "Flux thresholds": {
                "flux_limits": [Fclim, Fclim_max]
            },
            "Absolute (co)spectral value": {
                "cospec_abs_limit": cospec_abs_limit
            },
            "IQR over frequency band": {
                "iqr_stats": stats_r
            }
        }
    
        fqplt.plot_stepwise_filtering(
            mask_history=mask_history,
            step_names=step_names,
            fluxes_df=fluxes_df,
            flux_col=flux_col_gas,
            cospectra_df=Rcos_plot,
            freq=freq,
            outputpath=outputpath,
            plot=plot,
            init_valid_mask=None,
            step_plot_info=step_plot_info,
            filename=f"3_filtering_cof_gas__{run_suffix}.png",
        )
    
    if Icon_cof.isnull().sum().sum() == Icon_cof.size:
        raise FREQCORSelCofError('No (co)spectra after cof filtering')
        
        
    return Icon_cof, Rcon_cof
