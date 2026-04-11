# -*- coding: utf-8 -*-

"""General filtering utilities for (co)spectra and associated meteorological time series."""

import numpy as np

def FREQCOR_Sel_general(Icon_raw, Rcon_raw, Icosp_raw, meteo_df, WS, WD, 
                        Ustar, Zeta,  FlagF_H, FlagF_g, FlagVM_T, 
                        FlagVM_w, FlagVM_g, flag_wd, config):
    """
    Apply general filtering to all (co)spectra and associated meteorological series.
    
    The routine removes time steps with missing meteorology/fluxes and excludes wind
    direction sectors configured in `config['USER_LIMITS']`. Invalid columns in
    (co)spectral matrices are set to NaN.
    
    Parameters
    ----------
    Icon_raw, Rcon_raw, Icosp_raw : pandas.DataFrame
        Ideal content (raw) and real content (raw) (either cospectra or spectra; depends on `sps`), and ideal cospectra for CF (raw)
        with shape ``[nfreq, nspec]``.
    meteo_df : pandas.DataFrame
        Meteorological and flux time series table with length ``nspec``.
    WS, WD, Ustar, Zeta : pandas.Series
        Wind speed, wind direction, friction velocity, and stability parameter.
    FlagF_H, FlagF_g : pandas.Series
        Mauder & Foken (2004) quality flags for sensible heat flux and gas flux.
    FlagVM_T, FlagVM_w, FlagVM_g : array-like
        Vickers & Mahrt (1997) quality flags for temperature, vertical wind speed,
        and gas concentration.
    flag_wd : pandas.Series or numpy.ndarray
        Additional wind-direction flag (ETC/ICOS Vitale). Values equal to 2 are
        excluded.
    config : configparser.ConfigParser
        Run configuration. Uses:
        - `USER_LIMITS.wdmin`, `USER_LIMITS.wdmax`
        - `PROCEDURE_OPTIONS.gss`

    Returns
    -------
    Icon_sel, Rcon_sel, Icosp_sel : pandas.DataFrame
        Filtered ideal content (sel) and real content (sel) (either cospectra or spectra; depends on `sps`), and ideal cospectra for CF (sel)
        with shape ``[nfreq,nspec]``.
    meteo_df : pandas.DataFrame
        Filtered meteorological/flux table.
    WS, WD, Ustar, Zeta : pandas.Series
        Filtered time series.
    FlagF_H, FlagF_g, FlagVM_T, FlagVM_w, FlagVM_g : numpy.ndarray
        Filtered quality flags (converted to float arrays).

    Notes
    -----
    Wind direction exclusion is applied as follows:
    - if ``wdmin <= wdmax``: exclude ``(wdmin, wdmax)``
    - else: exclude ``(wdmin, 360)`` and ``(0, wdmax)`` (wrap-around)

    Data with missing sensible heat flux (FlagF_H == -9999) is set to NaN.
    """
    wdmax = float(config['USER_LIMITS']['wdmax'])
    wdmin = float(config['USER_LIMITS']['wdmin'])
    gss = int(config['PROCEDURE_OPTIONS']['gss'])
    
    # Convert negative wind directions to 0-360° range
    if wdmin < 0:
        wdmin = wdmin + 360
    if wdmax < 0:
        wdmax = wdmax + 360

    # Copy all input data
    Icon_sel = Icon_raw.copy()
    Rcon_sel = Rcon_raw.copy()
    Icosp_sel = Icosp_raw.copy()
    
    # Series to array of float
    FlagF_g = FlagF_g.values 
    FlagF_H = FlagF_H.values
    
    # Remove data when H is absent
    idx = FlagF_H==-9999
    meteo_df.loc[idx,'H'] = np.nan
    meteo_df.loc[idx,'Fc']=np.nan
    meteo_df.loc[idx,'FcEP']=np.nan
    if gss == 1 or gss == 2:
        meteo_df.loc[idx,'LE']=np.nan
        meteo_df.loc[idx,'RH']=np.nan

    WS[idx]=np.nan
    WD[idx]=np.nan
    Ustar[idx]=np.nan
    Zeta[idx]=np.nan

    if FlagF_g.dtype == 'int64' or FlagF_g.dtype == 'int32':
        FlagF_g = FlagF_g.astype('float')
    if FlagF_H.dtype == 'int64' or FlagF_H.dtype == 'int32':
        FlagF_H = FlagF_H.astype('float')
    
    FlagF_H[idx]=np.nan
    FlagF_g[idx]=np.nan
    FlagVM_T = FlagVM_T.astype('float')
    FlagVM_T[idx]=np.nan
    FlagVM_w=FlagVM_w.astype('float')
    FlagVM_w[idx]=np.nan
    FlagVM_g= FlagVM_g.astype('float')
    FlagVM_g[idx]=np.nan

    # Remove data with missing meteorological parameters
    if gss != 2:    
        meteo_valid = meteo_df[['H','Fc','FcEP']].notna().all(axis=1) & WS.notna() & WD.notna() & Ustar.notna() & Zeta.notna()
        # meteo_valid = meteo_df.drop(columns="RH").notna().all(axis=1) & WS.notna() & WD.notna() & Ustar.notna() & Zeta.notna()
    else: 
        meteo_valid = meteo_df[['H','LE','LEEP','RH']].notna().all(axis=1) & WS.notna() & WD.notna() & Ustar.notna() & Zeta.notna()
    # Wind direction filter
    if wdmin <= wdmax:
        wd_exclude = (WD > wdmin) & (WD < wdmax)
    else:
        wd_exclude = (WD > wdmin) | (WD < wdmax)
    
    # add etc-icos vitale wind direction flag. If the flag is not available,
    # it is a Series full of nan: the boolean will be full of "False" and 
    # therefore the filtering will still be valid
    wd_exclude_vitale = flag_wd == 2
    # Combined valid mask
    valid_mask = np.array(meteo_valid & ~wd_exclude & ~wd_exclude_vitale)
    
    # Set invalid data to nan
    Icon_sel.iloc[:, ~valid_mask] = np.nan
    Rcon_sel.iloc[:, ~valid_mask] = np.nan
    Icosp_sel.iloc[:, ~valid_mask] = np.nan
    
    del Icon_raw
    del Rcon_raw
    del Icosp_raw
    
    return Icon_sel,Rcon_sel,Icosp_sel,meteo_df,WS,WD,Ustar,Zeta,FlagF_H,FlagF_g,FlagVM_T,FlagVM_w,FlagVM_g

