# -*- coding: utf-8 -*-

"""Core computation orchestration for cut-off frequencies and correction factors."""

from .FREQCOR_LUT_CF import FREQCOR_LUT_CF
from .FREQCOR_LUT_cof import FREQCOR_LUT_cof
from .FREQCOR_Sel_CF import FREQCOR_Sel_CF
from .FREQCOR_Sel_cof import FREQCOR_Sel_cof
from .FREQCOR_Sel_stunst import FREQCOR_Sel_stunst

# Custom exception class for FREQCOR_Compute errors
class FREQCORComputeError(Exception):
    """Exception raised for errors in `FREQCOR_Compute`."""
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)

def FREQCOR_Compute(Icon_sel, Rcon_sel, Icosp_sel, nspec, nfreq, gss, sps, 
                    classnum_first, classnum, classnumCF_u, classnumCF_s, WS, WD, Zeta, 
                    meteo_df, Hlim, Fclim, LElim, Hlcfu, Hlcfs,
                    Hlim_max, Fclim_max, LElim_max, Hlcf_abs_max,
                    FlagVM_T, FlagVM_g, FlagVM_w, FlagF_H, FlagF_g, FVM, FF, 
                    jtmin_hz, jtmax_hz, varclim, freq, freqn, freqcount, 
                    tfmin_hz, tfmax_hz, plot, outputpath, 
                    massman_coef, tf_sonic, tf_peltola, cospec_abs_limit, run_tag=None):
    """
    Compute cut-off frequency and correction factor lookup tables.

    The routine orchestrates the full computation pipeline:
    - selection of very good (co)spectra for transfer function (TF) estimation,
    - computation of TF cut-off frequencies and their look-up table,
    - selection of fair quality (co)spectra for correction factor (CF) estimation,
    - computation of CF look-up tables for unstable and stable conditions.

    Parameters
    ----------
    Icon_sel : pandas.DataFrame
        Ideal content (sel) (co)spectral data (either cospectra or spectra; depends on `sps`).
    Rcon_sel : pandas.DataFrame
        Real content (sel) (co)spectral data (either cospectra or spectra; depends on `sps`).
    Icosp_sel : pandas.DataFrame
        Ideal cospectra for CF (sel).
    nspec : int
        Number of total initial valid (co)spectra files.
    nfreq : int
        Number of frequency bins.
    gss : int
        Gas species selector.
    sps : int
        Approach selector (cospectral or spectral).
    classnum_first : int
        Number of classes for the first sorting variable (e.g., wind direction or relative humidity).
    classnum : int
        Number of wind speed classes for cut-off frequency computation.
    classnumCF_u : int
        Number of wind speed classes for correction factors (unstable conditions).
    classnumCF_s : int
        Number of wind speed classes for correction factors (stable conditions).
    WS, WD, Zeta : pandas.Series
        Wind speed, wind direction, and stability parameter time series.
    meteo_df : pandas.DataFrame
        Half-hourly meteorological and flux data.
    Hlim, Fclim, LElim : float
        Flux thresholds used during (co)spectra selection for cut-off frequency computation.
    Hlcfu, Hlcfs : float
        Sensible heat flux thresholds used during (co)spectra selection for CF computation
        (unstable and stable conditions).
    FlagVM_T, FlagVM_g, FlagVM_w : array-like
        Vickers & Mahrt (1997) flags for sonic temperature, gas concentration, and vertical wind speed.
    FlagF_H, FlagF_g : array-like
        Mauder & Foken (2004) quality flags for sensible heat flux and gas flux.
    FVM : int
        Threshold applied to Vickers & Mahrt (1997) flags.
    FF : int
        Threshold applied to Mauder & Foken flags.
    jtmin_hz, jtmax_hz : float
        Frequency range (Hz) for transfer function normalization.
    varclim : float
        Threshold for coefficient of variation of the normalization parameter.
    freq : pandas.Series
        Natural frequencies (Hz).
    freqn : pandas.DataFrame
        Normalized frequencies per time step.
    freqcount : pandas.DataFrame
        Frequency-bin population counts per time step.
    tfmin_hz, tfmax_hz : float
        Frequency range (Hz) for transfer function fitting.
    plot : list[int]
        Plotting options selector.
    outputpath : str
        Output directory.
    massman_coef : object
        Massman coefficients or configuration required by downstream routines.
    tf_sonic, tf_peltola : object
        Transfer function configuration for the sonic and Peltola response.
    run_tag : str, optional
        Run identifier propagated to writers/plots.

    Returns
    -------
    LUT_cof : dict
        Cut-off frequency look-up table.
    classize : int
        Number of (co)spectra used per class during cut-off frequency computation.
    LUTCF_u : dict
        Correction factor look-up table (unstable conditions).
    matsortCF_u : pandas.DataFrame
        Sorting matrix used to build the unstable CF look-up table.
    classizeCF_u : int
        Number of cospectra used per class for the unstable CF look-up table.
    LUTCF_s : dict
        Correction factor look-up table (stable conditions).
    matsortCF_s : pandas.DataFrame
        Sorting matrix used to build the stable CF look-up table.
    classizeCF_s : int
        Number of cospectra used per class for the stable CF look-up table.

    """
    
    # %% (co)spectra selection
    Icon_cof, Rcon_cof = FREQCOR_Sel_cof(Icon_sel, Rcon_sel, freq, gss, sps, 
                                   meteo_df, Hlim, Fclim, LElim, Hlim_max, Fclim_max, LElim_max, FlagVM_T,
                                   FlagVM_g, FlagVM_w, FlagF_H, FlagF_g, FVM, 
                                   FF, outputpath, plot, cospec_abs_limit, run_tag=run_tag)
     
          
    # %% Cut off frequency computation 
    
    LUT_cof, classize = \
        FREQCOR_LUT_cof(sps, gss, classnum_first, classnum, Icon_cof, Rcon_cof, WS, WD, Zeta, 
                        meteo_df, nfreq, jtmin_hz,jtmax_hz,
                        varclim,freq,freqn,tfmin_hz,tfmax_hz,plot,
                        outputpath, tf_peltola, run_tag)
    
    # %% Correction factor computation
    # Unstable conditions
    sts=1
    Icos_CF = FREQCOR_Sel_stunst(Icosp_sel, sts, Zeta)
    # Data selection 
    Icos_CF = FREQCOR_Sel_CF(Icos_CF, freq, sts, meteo_df, Hlcfu, 
                           Hlcfs, Hlcf_abs_max, FlagVM_T, FlagF_H, FlagVM_w, FVM, FF,
                           outputpath, plot, run_tag=run_tag)
    
    # Correction factor computation
    LUTCF_u, matsortCF_u, classizeCF_u = FREQCOR_LUT_CF(Icos_CF, nspec, sts, gss, WS, 
                                                   Zeta, WD, meteo_df, LUT_cof, 
                                                   classnum, classnumCF_u, nfreq,  
                                                   freq, freqn, sps, 
                                                   tf_sonic, tf_peltola,
                                                   outputpath, plot, massman_coef, run_tag)
    
    
    
    # Stable conditions
    sts=2
    Icos_CF = FREQCOR_Sel_stunst(Icosp_sel, sts, Zeta)
    # Data selection 
    Icos_CF = FREQCOR_Sel_CF(Icos_CF, freq, sts, meteo_df, Hlcfu, 
                           Hlcfs, Hlcf_abs_max, FlagVM_T, FlagF_H, FlagVM_w, FVM, FF,
                           outputpath, plot, run_tag=run_tag)
    
    # Correction factor computation
    LUTCF_s, matsortCF_s, classizeCF_s = FREQCOR_LUT_CF(Icos_CF, nspec, sts, gss, WS, 
                                                   Zeta, WD, meteo_df, LUT_cof, 
                                                   classnum, classnumCF_s, nfreq,  
                                                   freq, freqn, sps, 
                                                   tf_sonic, tf_peltola,
                                                   outputpath, plot, massman_coef, run_tag)
    
    return LUT_cof, classize, LUTCF_u, matsortCF_u, classizeCF_u, LUTCF_s, matsortCF_s, classizeCF_s
