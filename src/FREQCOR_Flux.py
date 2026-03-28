# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 14:49:58 2022

@author: Ariane Faurès
"""

import numpy as np
import pandas as pd
import FREQCOR_functions as fqfct
from FREQCOR_Sensor_Separation import FREQCOR_Sensor_Separation

def FREQCOR_Flux(nspec, sps, gss, meteo_df, WS, WD, Zeta, LUTCF_u,
             classnumCF_u, LUTCF_s, classnumCF_s, zmeas, disph, rsz, rseast,
             rsnorth, eq):
    """
    Compute corrected fluxes using spectral correction factors and LUTs.

    Parameters
    ----------
    nspec : int
        Number of half-hourly periods (number of spectra).
    sps : int
        Correction method selection (1: standard, 2: sensor separation correction).
    gss : int
        Gas species selection (1: CO2, 2: H2O, etc.).
    meteo_df : pandas.DataFrame
        Meteorological and flux data for each period.
    WS : array-like
        Wind speed series.
    WD : array-like
        Wind direction series.
    Zeta : array-like
        Stability parameter series.
    LUTCF_u : dict
        LUT for unstable conditions.
    classnumCF_u : int
        Number of LUT classes for unstable conditions.
    LUTCF_s : dict
        LUT for stable conditions.
    classnumCF_s : int
        Number of LUT classes for stable conditions.
    zmeas : float
        Measurement height (m).
    disph : float
        Displacement height (m).
    rsz : float
        Sensor separation (z direction, m).
    rseast : float
        Sensor separation (east direction, m).
    rsnorth : float
        Sensor separation (north direction, m).
    eq : int
        Equation selection (13 or 16).

    Returns
    -------
    list
        [FccorL, FccorG, CFHL, CF_L, CF_G]
        FccorL: pandas.Series of corrected fluxes (Lorentzian).
        FccorG: pandas.Series of corrected fluxes (Gaussian).
        CFHL: Sensor separation correction factor.
        CF_L: numpy.ndarray of correction factors (Lorentzian).
        CF_G: numpy.ndarray of correction factors (Gaussian).
    """
    

    CFHL = np.full(nspec, np.nan)
    
    # Initialisation of variables
    FccorL = pd.Series(np.full(nspec, np.nan))
    FcuncL = pd.Series(np.full(nspec, np.nan))
    FccorG = pd.Series(np.full(nspec, np.nan))
    FcuncG = pd.Series(np.full(nspec, np.nan))
    
    CF_L = np.full(nspec, np.nan)
    CF_G = np.full(nspec, np.nan)
    CFunc_L = np.full(nspec, np.nan)
    CFunc_G = np.full(nspec, np.nan)
    
    XVR=[np.nan,np.nan] # when 2-levels LUT
    if gss==2:
        RH=meteo_df.loc[:,'RH']
        
    for js in range(nspec):
        # Application of the correction
        # Selection of LUT_CF class and application of factor 
        # Selection of sorting variable (user entry)
        XV=WS[js]
            
        if Zeta[js]<=0: 
            LUTCF = LUTCF_u
            classnumCF = classnumCF_u
        else: 
            LUTCF = LUTCF_s  
            classnumCF = classnumCF_s
        
        if not(gss==2):
            if np.isfinite(WD[js]):
                XVR=fqfct.ReadLUT(XV,LUTCF,classnumCF,gss,WD[js])
                CF_L[js]=LUTCF[int(XVR[1])]['ws'].iloc[int(XVR[0]),2]
                CF_G[js]=LUTCF[int(XVR[1])]['ws'].iloc[int(XVR[0]),4]
                CFunc_L[js]=LUTCF[int(XVR[1])]['ws'].iloc[int(XVR[0]),3]
                CFunc_G[js]=LUTCF[int(XVR[1])]['ws'].iloc[int(XVR[0]),5]
            else:
                CF_L[js]=np.nan
                CF_G[js]=np.nan
                CFunc_L[js]=np.nan
                CFunc_G[js]=np.nan
            
            FccorL[js]=meteo_df.loc[js,'Fc']*CF_L[js]
            FcuncL[js]=meteo_df.loc[js,'Fc']*CFunc_L[js]
            FccorG[js]=meteo_df.loc[js,'Fc']*CF_G[js]
            FcuncG[js]=meteo_df.loc[js,'Fc']*CFunc_G[js]
            
        elif gss==2:
            XVR=fqfct.ReadLUT(XV,LUTCF,classnumCF,gss,RH[js])
            if np.isfinite(RH[js]):
                CF_L[js]=LUTCF[int(XVR[1])]['ws'].iloc[int(XVR[0]),2]
                CF_G[js]=LUTCF[int(XVR[1])]['ws'].iloc[int(XVR[0]),4]
                CFunc_L[js]=LUTCF[int(XVR[1])]['ws'].iloc[int(XVR[0]),3]
                CFunc_G[js]=LUTCF[int(XVR[1])]['ws'].iloc[int(XVR[0]),5]
            else:
                CF_L[js]=np.nan
                CF_G[js]=np.nan
                CFunc_L[js]=np.nan
                CFunc_G[js]=np.nan
                
            FccorL[js]=meteo_df.loc[js,'LE']*CF_L[js]
            FcuncL[js]=meteo_df.loc[js,'LE']*CFunc_L[js]
            FccorG[js]=meteo_df.loc[js,'LE']*CF_G[js]
            FcuncG[js]=meteo_df.loc[js,'LE']*CFunc_G[js]
        
    # Computation of sensor separation correction factor (spectral approach)
    if sps==2:
        CFHL = FREQCOR_Sensor_Separation(zmeas,disph,rsz,WD,rseast,rsnorth,Zeta,eq)
            
    # Application of sensor separation correction
    if sps==2:
        FccorL=FccorL*CFHL 
        FccorG=FccorG*CFHL
        CF_L = CF_L*CFHL
        CF_G = CF_G*CFHL

    return [FccorL, FccorG, CFHL, CF_L, CF_G]



