# -*- coding: utf-8 -*-

"""Look-up table construction for spectral correction factors (CF)."""

from FREQCOR_functions import ReadLUT, Simpson, TFsonic, av_Kaimal, remove_outliers
import FREQCOR_plot as fqplt
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

class FREQCORLUTCFError(Exception):
    """Exception raised for errors in `FREQCOR_LUT_CF`."""
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)

_AV_IDEALS_CACHE = {}

def FREQCOR_LUT_CF(Icos_CF, nspec, sts, gss, WS, Zeta, WD, meteo_df, LUT_cof, 
                   classnum, classnumCF, nfreq, freq, freqn, sps, 
                   tf_sonic, tf_peltola, outputpath, plot, massman_coef, run_tag):
    """
    Build the correction factor (CF) look-up table.

    For each half-hourly cospectrum, the routine retrieves cut-off frequencies from
    `LUT_cof`, builds transfer functions (Lorentz/Gauss, optionally Peltola form),
    computes correction factors via numerical integration, then aggregates them into
    a LUT sorted by a first variable (wind direction for CO2 or relative humidity for
    H2O) and by wind speed classes.

    CF LUT uncertainties:
    - ``unc_*_tf``: transfer-function-related uncertainty component, computed as the
      maximum deviation from the central CF value ("M") within each WS class.
    - ``unc_*_sd``: within-class spread of CF values, computed as the standard
      deviation of CF values after outlier removal.
    Total uncertainties are not produced; users can recompute them from the component
    uncertainties if needed.

    Parameters
    ----------
    Icos_CF : pandas.DataFrame
        Cospectra used for CF computation.
    nspec : int
        Number of total initial valid (co)spectra files.
    sts : int
        Stability class selector (1: unstable, 2: stable).
    gss : int
        Gas species selector.
    WS, Zeta, WD : pandas.Series
        Wind speed, stability parameter, and wind direction time series.
    meteo_df : pandas.DataFrame
        Half-hourly meteorological and flux data (used for RH when `gss == 2`).
    LUT_cof : dict
        Cut-off frequency look-up table.
    classnum : int
        Number of wind speed classes used for cut-off frequency LUT construction.
    classnumCF : int
        Number of wind speed classes used for CF LUT construction.
    nfreq : int
        Number of frequency bins.
    freq : pandas.Series
        Natural frequencies (Hz).
    freqn : pandas.DataFrame
        Normalized frequencies per time step.
    sps : int
        Approach selector (cospectral or spectral).
    tf_sonic : int
        If 1, apply theoretical correction to sonic cospectra.
    tf_peltola : int
        If 1, use the Peltola et al. (2021) transfer function form where applicable.
    outputpath : str
        Output directory.
    plot : list[int]
        Plotting options selector.
    massman_coef : array-like
        Massman coefficients for the site (e.g., ``[A0, kf0, mu]``).
    run_tag : str
        Run identifier used to generate output filenames.

    Returns
    -------
    LUT_CF : dict
        Correction factor look-up table.
    matsortCF : pandas.DataFrame
        Sorting matrix with meteorology and correction factors (after filtering missing values).
    classize_ws : int
        Number of (co)spectra used per wind speed class in the LUT.

    Raises
    ------
    ValueError
        If `run_tag` is missing.
    FREQCORLUTCFError
        If there are not enough cospectra to populate the wind speed classes.
    """

    if not run_tag:
        raise ValueError("run_tag is required for output naming")
    stability=['unst','st']
    print('  Stability class: ' + stability[sts-1] + 'able. Computing spectral correction factors')
    cols = ['CFL_M','CFL_L','CFL_H','CFG_M','CFG_L','CFG_H']
    matCF = pd.DataFrame(np.full(shape=(nspec,len(cols)),fill_value=np.nan), index = None, columns=cols)
    WSR = pd.DataFrame(np.zeros(shape=[nspec,2]))
    WSR[:] = np.nan
    if gss != 2:
        RH = pd.Series(np.zeros(shape=meteo_df.shape[0]))
    elif gss==2:
        RH = meteo_df.RH
        
    # %% Correction factor computation from individual valid H cospectra
    for js in range(nspec):
        Icoscf=Icos_CF.copy()
        Icoscf = Icoscf.iloc[:,js]
        if np.isfinite(sum(Icoscf)) and np.isfinite(RH[js]):
        
            # Selecting cut off frequency in LUT_cof: according to WD/RH and WS
            if not(gss==2): 
                WSR.iloc[js,0],WSR.iloc[js,1]=ReadLUT(WS[js],LUT_cof,classnum,gss,WD[js])
                ws_row = LUT_cof[int(WSR.iloc[js,1])]['ws'].iloc[int(WSR.iloc[js,0])]
                cof_L = ws_row['cof_L']
                unc_L_tf = ws_row['unc_L_tf']
                cof_G = ws_row['cof_G']
                unc_G_tf = ws_row['unc_G_tf']
            elif gss == 2: 
                WSR.iloc[js,0],WSR.iloc[js,1]=ReadLUT(WS[js],LUT_cof,classnum,gss,RH[js])
                ws_row = LUT_cof[int(WSR.iloc[js,1])]['ws'].iloc[int(WSR.iloc[js,0])]
                cof_L = ws_row['cof_L']
                unc_L_tf = ws_row['unc_L_tf']
                cof_G = ws_row['cof_G']
                unc_G_tf = ws_row['unc_G_tf']
                
            if not(np.isfinite(cof_L) and np.isfinite(cof_G)):
                matCF.iloc[js,0]=np.nan
                matCF.iloc[js,1]=np.nan
                matCF.iloc[js,2]=np.nan
                matCF.iloc[js,3]=np.nan
                matCF.iloc[js,4]=np.nan
                matCF.iloc[js,5]=np.nan
                continue
            
            
            # Correct H cospectra for high spectral losses (path averaging and time response)
            # if requested by the user            
            if tf_sonic == 1 : 
                # Sensible heat cospectra to be corrected and transfer function computation
                Rcoscfs = Icoscf.copy() 
                TFs = TFsonic(freq,WS[js]) 
                # Building of ideal cospectra: divide the degraded by TF
                Icoscfs=Rcoscfs/TFs
                Icoscf = Icoscfs.copy()
                
                
            # Building of transfer functions and their confidence limits
            # Cospectra: 
            #       Lorentzian/Gaussian : TF
            #       Peltola Lorentzian : sqrt(TF)
            # Spectra : 
            #       Lorentzian/Gaussian : sqrt(TF)    
            if sps==1:
                if tf_peltola ==1: 
                    TraFunL_M=(1/(1+(freq/cof_L)**2))**0.5                                 
                    TraFunL_L=(1/(1+(freq/(cof_L+unc_L_tf))**2))**0.5  
                    TraFunL_H=(1/(1+(freq/(cof_L-unc_L_tf))**2))**0.5             
                else:
                    TraFunL_M=1/(1+(freq/cof_L)**2)                                             
                    TraFunL_L=1/(1+(freq/(cof_L+unc_L_tf))**2)
                    if unc_L_tf < cof_L:
                        TraFunL_H=1/(1+(freq/(cof_L-unc_L_tf))**2)              
                    else:
                        TraFunL_H=1/(1+(freq/0.005))**2                
            else:   
                TraFunL_M=(1/(1+(freq/cof_L)**2))**0.5                                 
                TraFunL_L=(1/(1+(freq/(cof_L+unc_L_tf))**2))**0.5  
                TraFunL_H=(1/(1+(freq/(cof_L-unc_L_tf))**2))**0.5            
            if sps==1:
                TraFunG_M=np.exp(-np.log(2)*(freq/cof_G)**2)            
                TraFunG_L=np.exp(-np.log(2)*(freq/(cof_G+unc_G_tf))**2)  
                TraFunG_H=np.exp(-np.log(2)*(freq/(cof_G-unc_G_tf))**2)  
            else:
                TraFunG_M=(np.exp(-np.log(2)*(freq/cof_G)**2))**0.5
                TraFunG_L=(np.exp(-np.log(2)*(freq/(cof_G+unc_G_tf))**2))**0.5  
                TraFunG_H=(np.exp(-np.log(2)*(freq/(cof_G-unc_G_tf))**2))**0.5  
    
            # Building of degraded cospectra          
            RcoscfL_M=Icoscf*TraFunL_M
            RcoscfL_L=Icoscf*TraFunL_L
            RcoscfL_H=Icoscf*TraFunL_H
            RcoscfG_M=Icoscf*TraFunG_M
            RcoscfG_L=Icoscf*TraFunG_L
            RcoscfG_H=Icoscf*TraFunG_H
            # Preparation of vectors for numerical integration (Simpson method)   
            IcoscfS=Simpson(Icoscf)
            RcoscfL_MS=Simpson(RcoscfL_M)
            RcoscfL_LS=Simpson(RcoscfL_L)
            RcoscfL_HS=Simpson(RcoscfL_H)
            RcoscfG_MS=Simpson(RcoscfG_M)
            RcoscfG_LS=Simpson(RcoscfG_L)
            RcoscfG_HS=Simpson(RcoscfG_H)
    
            # Elimination of vector extremities        
            IcoscfS = IcoscfS[1:-1]
            RcoscfL_MS = RcoscfL_MS[1:-1]
            RcoscfL_LS = RcoscfL_LS[1:-1]
            RcoscfL_HS = RcoscfL_HS[1:-1]
            RcoscfG_MS = RcoscfG_MS[1:-1]
            RcoscfG_LS = RcoscfG_LS[1:-1]
            RcoscfG_HS = RcoscfG_HS[1:-1]
            
            # Definition of the frequency increment     
            df = pd.Series(np.zeros(shape = nfreq-1)) 
            for jf in range(1,nfreq-1):
                df[jf]=(freq[jf+1]-freq[jf-1])/2
            df = df.drop(0)
            
            # Correction factor computation    
            CFLL_M = np.nansum(IcoscfS*df.T)/np.nansum(RcoscfL_MS*df.T)
            CFLL_L = np.nansum(IcoscfS*df.T)/np.nansum(RcoscfL_LS*df.T)
            CFLL_H = np.nansum(IcoscfS*df.T)/np.nansum(RcoscfL_HS*df.T)
            CFLG_M = np.nansum(IcoscfS*df.T)/np.nansum(RcoscfG_MS*df.T)
            CFLG_L = np.nansum(IcoscfS*df.T)/np.nansum(RcoscfG_LS*df.T)
            CFLG_H = np.nansum(IcoscfS*df.T)/np.nansum(RcoscfG_HS*df.T)
            
            # Putting correction factors and their confidence limits in a matrix
            matCF.loc[js, cols] = [CFLL_M, CFLL_L, CFLL_H, CFLG_M, CFLG_L, CFLG_H]
    
    
    # %% Building the CF look-up table
    outlier_info_L = {}
    outlier_info_G = {}
    valid_mask = np.isfinite(matCF).all(axis=1)
    mat_meteo_cf =  pd.concat([WS,WD,RH,Zeta,matCF], axis = 1)
  
    # Suppression of missing CF rows
    matCF = matCF[valid_mask].reset_index(drop=True)
    mat_meteo_cf = mat_meteo_cf[valid_mask].reset_index(drop=True)
    Zeta = Zeta[valid_mask].reset_index(drop=True)
    Cos_i=Icos_CF.iloc[:,valid_mask.to_numpy()]
    freqn = freqn.loc[:,valid_mask.to_numpy()]
    
    if not(gss==2):
        classnum_wd = len(LUT_cof)
        li = list(np.linspace(1, classnum_wd, classnum_wd, dtype=int))
        LUT_CF = dict.fromkeys(li)
        all_data_cf=dict.fromkeys(li)
        
        # WD, apply symmetry: map [180deg-360deg] to [0deg-180deg])
        mat_meteo_cf['wd_reduced'] = mat_meteo_cf['wind_dir'] % 360
        mat_meteo_cf['wd_reduced'] = mat_meteo_cf['wd_reduced'].apply(lambda wd: wd - 180 if wd > 180 else wd)
        bin_edges = np.linspace(0, 180, classnum_wd + 1)
        wd_bin_labels = list(range(1, classnum_wd + 1))
        mat_meteo_cf['wd_bin'] = pd.cut(mat_meteo_cf['wd_reduced'], bins=bin_edges, labels=wd_bin_labels, include_lowest=True)
        
        for jc in range(1, classnum_wd + 1):
            wd_class_df = mat_meteo_cf[mat_meteo_cf['wd_bin'] == jc].copy().reset_index(drop=True)         
            cf_cols = [
                'ws_mean',
                'ws_max',
                'CF_L',
                'unc_L_tf',
                'unc_L_sd',
                'CF_G',
                'unc_G_tf',
                'unc_G_sd',
                'unc_wsCFclass',
            ]
            lut_cf_temp = {
                'wd': pd.DataFrame(),
                'ws': pd.DataFrame(np.zeros(shape=(classnumCF, len(cf_cols))), columns=cf_cols),
                'ws_n': pd.Series(np.full(shape=(classnumCF,), fill_value=np.nan))
            }
        
            lut_cf_temp['wd'].loc[0, 'wdmean'] = wd_class_df['wd_reduced'].mean()
            lut_cf_temp['wd'].loc[0, 'wdmax'] = wd_class_df['wd_reduced'].max()
            lut_cf_temp['wd'].loc[0, 'uncclass'] = wd_class_df['wd_reduced'].std()
        
            matsortws = wd_class_df.sort_values(by='wind_speed').reset_index(drop=True)
            
            classize_wd = len(matsortws)
            classize_ws = int(np.floor(classize_wd / classnumCF))
            if classize_ws < 2:
                raise FREQCORLUTCFError('Not enough cospectra in each LUT wind speed class (within WD class) for LUT construction.')
            
            unc_CFL_tf = pd.Series(np.full(shape=[classnumCF], fill_value=np.nan))
            unc_CFG_tf = pd.Series(np.full(shape=[classnumCF], fill_value=np.nan))
            unc_CFL_sd = pd.Series(np.full(shape=[classnumCF], fill_value=np.nan))
            unc_CFG_sd = pd.Series(np.full(shape=[classnumCF], fill_value=np.nan))
            
            all_cf_ws = {
                'ws_l': [],
                'ws_g': [],
                'cf_l': [],
                'cf_g': [],
                'mean_cf_l': [],
                'mean_cf_g': [],
                'wd': lut_cf_temp['wd'].loc[0,'wdmean'],
                'mean_ws':[],
            }
            
            bins = np.array_split(matsortws.index, classnumCF)
            for jws, bin_idx in enumerate(bins, 1):
                bin_df = matsortws.loc[bin_idx]
                ws_values = bin_df['wind_speed']
                cf_lorentz = bin_df['CFL_M']
                cf_gauss = bin_df['CFG_M']
            
                class_id = f"wd{jc}_ws{jws}"
                cf_lorentz_clean = remove_outliers(cf_lorentz, class_id, ws_values, outlier_info_L)
                cf_gauss_clean = remove_outliers(cf_gauss, class_id, ws_values, outlier_info_G)

                ws_clean = ws_values.loc[cf_lorentz_clean.index.intersection(cf_gauss_clean.index)]
                if len(ws_clean) == 0:
                    ws_clean = ws_values.loc[cf_lorentz_clean.index]
                if len(ws_clean) == 0:
                    ws_clean = ws_values.loc[cf_gauss_clean.index]
                if len(ws_clean) == 0:
                    ws_clean = ws_values

                lut_cf_temp['ws'].loc[jws-1,'ws_mean']=np.mean(ws_clean)
                lut_cf_temp['ws'].loc[jws-1,'ws_max']=np.max(ws_clean)            
                lut_cf_temp['ws'].loc[jws-1,'unc_wsCFclass']=1.96*np.std(ws_clean)
                lut_cf_temp['ws_n'].iloc[jws-1]=int(len(ws_clean))
                lut_cf_temp['ws'].loc[jws-1,'CF_L']=np.median(cf_lorentz_clean)
                lut_cf_temp['ws'].loc[jws-1,'CF_G']=np.median(cf_gauss_clean)
                
                # uncertainty due to transfer function
                mean_CFL_M = np.mean(bin_df['CFL_M'])
                mean_CFL_L = np.mean(bin_df['CFL_L'])
                mean_CFL_H = np.mean(bin_df['CFL_H'])
                unc_CFL_tf[jws-1] = max(abs(mean_CFL_M - mean_CFL_L), abs(mean_CFL_H - mean_CFL_M))

                mean_CFG_M = np.mean(bin_df['CFG_M'])
                mean_CFG_L = np.mean(bin_df['CFG_L'])
                mean_CFG_H = np.mean(bin_df['CFG_H'])
                unc_CFG_tf[jws-1] = max(abs(mean_CFG_M - mean_CFG_L), abs(mean_CFG_H - mean_CFG_M))
                
                unc_CFL_sd[jws-1] = np.std(cf_lorentz_clean)
                unc_CFG_sd[jws-1] = np.std(cf_gauss_clean)

                lut_cf_temp['ws'].loc[jws-1,'unc_L_tf']=unc_CFL_tf[jws-1]
                lut_cf_temp['ws'].loc[jws-1,'unc_L_sd']=unc_CFL_sd[jws-1]
                lut_cf_temp['ws'].loc[jws-1,'unc_G_tf']=unc_CFG_tf[jws-1]
                lut_cf_temp['ws'].loc[jws-1,'unc_G_sd']=unc_CFG_sd[jws-1]

                all_cf_ws['ws_l'].append(ws_values.loc[cf_lorentz_clean.index])
                all_cf_ws['cf_l'].append(cf_lorentz_clean)
                all_cf_ws['ws_g'].append(ws_values.loc[cf_gauss_clean.index])
                all_cf_ws['cf_g'].append(cf_gauss_clean)
            
            all_cf_ws['ws_l'] = pd.concat(all_cf_ws['ws_l']).reset_index(drop=True) if len(all_cf_ws['ws_l']) else pd.Series(dtype=float)
            all_cf_ws['cf_l'] = pd.concat(all_cf_ws['cf_l']).reset_index(drop=True) if len(all_cf_ws['cf_l']) else pd.Series(dtype=float)
            all_cf_ws['ws_g'] = pd.concat(all_cf_ws['ws_g']).reset_index(drop=True) if len(all_cf_ws['ws_g']) else pd.Series(dtype=float)
            all_cf_ws['cf_g'] = pd.concat(all_cf_ws['cf_g']).reset_index(drop=True) if len(all_cf_ws['cf_g']) else pd.Series(dtype=float)

            all_cf_ws['mean_cf_l'].append(lut_cf_temp['ws'].loc[:,'CF_L'])
            all_cf_ws['mean_cf_g'].append(lut_cf_temp['ws'].loc[:,'CF_G'])
            all_cf_ws['mean_ws'].append(lut_cf_temp['ws'].loc[:,'ws_mean'])
            
            all_data_cf[jc]=all_cf_ws
            LUT_CF[jc]=lut_cf_temp
        
    elif gss==2:
        classnum_rh=len(LUT_cof)
        li = list(np.linspace(1,classnum_rh,classnum_rh,dtype=int))
        LUT_CF= dict.fromkeys(li)
        all_data_cf=dict.fromkeys(li)
        
        mat_meteo_cf=mat_meteo_cf.sort_values(by = ['RH'])
        bins_rh = np.array_split(mat_meteo_cf.index, classnum_rh)
        for jc, bin_idx in enumerate(bins_rh, 1):
            bin_rh = mat_meteo_cf.loc[bin_idx]
            cf_cols = [
                'ws_mean',
                'ws_max',
                'CF_L',
                'unc_L_tf',
                'unc_L_sd',
                'CF_G',
                'unc_G_tf',
                'unc_G_sd',
                'unc_wsCFclass',
            ]
            lut_cf_temp={
                'rh':pd.DataFrame(),
                'ws':pd.DataFrame(np.zeros(shape = (classnumCF,len(cf_cols))), columns=cf_cols),
                'ws_n': pd.Series(np.full(shape=(classnumCF,), fill_value=np.nan))
                }
            
            lut_cf_temp['rh'].loc[0,'rhmean']=np.mean(bin_rh['RH'])           
            lut_cf_temp['rh'].loc[0,'rhmax']=np.max(bin_rh['RH'])            
            lut_cf_temp['rh'].loc[0,'uncclass']=1.96*np.std(bin_rh['RH']) 
            
            matsortws=bin_rh.sort_values(by = ['wind_speed']).reset_index(drop=True)
            
            classize_rh = len(matsortws)
            classize_ws = int(np.floor(classize_rh / classnumCF))
            if classize_ws < 2:
                raise FREQCORLUTCFError('Not enough cospectra in each LUT wind speed class (within RH class) for LUT construction.')
            
            unc_CFL_tf = pd.Series(np.full(shape=[classnumCF], fill_value=np.nan))
            unc_CFG_tf = pd.Series(np.full(shape=[classnumCF], fill_value=np.nan))
            unc_CFL_sd = pd.Series(np.full(shape=[classnumCF], fill_value=np.nan))
            unc_CFG_sd = pd.Series(np.full(shape=[classnumCF], fill_value=np.nan))
            
            all_cf_ws = {
                'ws_l': [],
                'ws_g': [],
                'cf_l': [],
                'cf_g': [],
                'mean_cf_l': [],
                'mean_cf_g': [],
                'rh': lut_cf_temp['rh'].loc[0,'rhmean'],
                'mean_ws':[],
            }
            
            bins = np.array_split(matsortws, classnumCF)
            for jws, bin_df in enumerate(bins, 1):
                ws_values = bin_df['wind_speed']
                cf_lorentz = bin_df['CFL_M']
                cf_gauss = bin_df['CFG_M']
           
                class_id = f"rh{jc}_ws{jws}"
                cf_lorentz_clean = remove_outliers(cf_lorentz, class_id, ws_values, outlier_info_L)
                cf_gauss_clean = remove_outliers(cf_gauss, class_id, ws_values, outlier_info_G)

                ws_clean = ws_values.loc[cf_lorentz_clean.index.intersection(cf_gauss_clean.index)]
                if len(ws_clean) == 0:
                    ws_clean = ws_values.loc[cf_lorentz_clean.index]
                if len(ws_clean) == 0:
                    ws_clean = ws_values.loc[cf_gauss_clean.index]
                if len(ws_clean) == 0:
                    ws_clean = ws_values

                
                lut_cf_temp['ws'].loc[jws-1,'ws_mean']=np.mean(ws_clean)            
                lut_cf_temp['ws'].loc[jws-1,'ws_max']=np.max(ws_clean)             
                lut_cf_temp['ws'].loc[jws-1,'unc_wsCFclass']=1.96*np.std(ws_clean)
                lut_cf_temp['ws_n'].iloc[jws-1]=int(len(ws_clean))
                lut_cf_temp['ws'].loc[jws-1,'CF_L']=np.median(cf_lorentz_clean)
                lut_cf_temp['ws'].loc[jws-1,'CF_G']=np.median(cf_gauss_clean) 
               
                # uncertainty due to transfer function
                mean_CFL_M = np.mean(bin_df['CFL_M'])
                mean_CFL_L = np.mean(bin_df['CFL_L'])
                mean_CFL_H = np.mean(bin_df['CFL_H'])
                unc_CFL_tf[jws-1] = max(abs(mean_CFL_M - mean_CFL_L), abs(mean_CFL_H - mean_CFL_M))

                mean_CFG_M = np.mean(bin_df['CFG_M'])
                mean_CFG_L = np.mean(bin_df['CFG_L'])
                mean_CFG_H = np.mean(bin_df['CFG_H'])
                unc_CFG_tf[jws-1] = max(abs(mean_CFG_M - mean_CFG_L), abs(mean_CFG_H - mean_CFG_M))

                unc_CFL_sd[jws-1] = np.std(cf_lorentz_clean)
                unc_CFG_sd[jws-1] = np.std(cf_gauss_clean)

                lut_cf_temp['ws'].loc[jws-1,'unc_L_tf']=unc_CFL_tf[jws-1]
                lut_cf_temp['ws'].loc[jws-1,'unc_L_sd']=unc_CFL_sd[jws-1]
                lut_cf_temp['ws'].loc[jws-1,'unc_G_tf']=unc_CFG_tf[jws-1]
                lut_cf_temp['ws'].loc[jws-1,'unc_G_sd']=unc_CFG_sd[jws-1]

                all_cf_ws['ws_l'].append(ws_values.loc[cf_lorentz_clean.index])
                all_cf_ws['cf_l'].append(cf_lorentz_clean)
                all_cf_ws['ws_g'].append(ws_values.loc[cf_gauss_clean.index])
                all_cf_ws['cf_g'].append(cf_gauss_clean)
               
            all_cf_ws['ws_l'] = pd.concat(all_cf_ws['ws_l']).reset_index(drop=True) if len(all_cf_ws['ws_l']) else pd.Series(dtype=float)
            all_cf_ws['cf_l'] = pd.concat(all_cf_ws['cf_l']).reset_index(drop=True) if len(all_cf_ws['cf_l']) else pd.Series(dtype=float)
            all_cf_ws['ws_g'] = pd.concat(all_cf_ws['ws_g']).reset_index(drop=True) if len(all_cf_ws['ws_g']) else pd.Series(dtype=float)
            all_cf_ws['cf_g'] = pd.concat(all_cf_ws['cf_g']).reset_index(drop=True) if len(all_cf_ws['cf_g']) else pd.Series(dtype=float)

            all_cf_ws['mean_cf_l'].append(lut_cf_temp['ws'].loc[:,'CF_L'])
            all_cf_ws['mean_cf_g'].append(lut_cf_temp['ws'].loc[:,'CF_G'])
            all_cf_ws['mean_ws'].append(lut_cf_temp['ws'].loc[:,'ws_mean'])
            
            all_data_cf[jc]=all_cf_ws
            LUT_CF[jc]=lut_cf_temp
        
    run_suffix = run_tag

    # Plot CF as function of WS
    if plot[2]==1:
        file_tag = None
        if gss == 2:
            file_tag = f"5_CF_vs_ws__h2o__{stability[sts-1]}__{run_suffix}"
        else:
            file_tag = f"5_CF_vs_ws__{stability[sts-1]}__{run_suffix}"
        fqplt.plot_u_cf(all_data_cf, gss, plot, outputpath, stability[sts-1], file_tag=file_tag)
        
    # Average Kaimal and Massman cospectra computation and plotting
    if not(np.isnan(massman_coef).all()):
        A0 = massman_coef[0]
        kf0 = massman_coef[1]
        mu = massman_coef[2]
    else:
        A0 = np.nan
        kf0 = np.nan
        mu = np.nan

    if plot[2] == 1:
        cache_key = (run_suffix, gss)
        if sts == 1:
            av_cosp_k = av_Kaimal(Cos_i, freq, freqn, Zeta, A0, kf0, mu)
            _AV_IDEALS_CACHE[cache_key] = av_cosp_k
        elif sts == 2:
            av_unst = _AV_IDEALS_CACHE.get(cache_key, None)
            if av_unst is not None:
                av_cosp_k = av_Kaimal(Cos_i, freq, freqn, Zeta, A0, kf0, mu)
                fig, axs = plt.subplots(2, 1, figsize=(16, 12), sharex=True)

                title_fs = 16
                label_fs = 16
                tick_fs = 12
                legend_fs = 14

                for ax, df, title_suffix in (
                    (axs[0], av_unst, 'unstable'),
                    (axs[1], av_cosp_k, 'stable'),
                ):
                    ax.plot(df['Ideal'], '.', label='Ideal')
                    ax.plot(df['Kaimal'], label='Kaimal')
                    ax.plot(df['Massman'], label='Massman')
                    ax.set_yscale('log')
                    ax.set_xscale('log')
                    ax.set_ylim(1e-5, 1)
                    ax.set_xlim(1e-4, 1e3)
                    ax.grid(alpha=0.5)
                    ax.legend(fontsize=legend_fs)
                    ax.set_title(
                        f'Ideal and model cospectra. {title_suffix} conditions',
                        fontsize=title_fs,
                    )
                    ax.set_ylabel('Cosp/covar*fnat', fontsize=label_fs)
                    ax.tick_params(axis='both', which='both', labelsize=tick_fs)

                axs[1].set_xlabel('Freq norm [-]', fontsize=label_fs)
                plt.tight_layout()

                if plot[3] == 1:
                    file_tag = f"5_comp_reference_cospectra__{run_suffix}"
                    plt.savefig(outputpath + '/' + file_tag + '.png')
                    plt.close()
                else:
                    plt.show()

                _AV_IDEALS_CACHE.pop(cache_key, None)

        
   
        
    # Write stats to text file
    stats_path = outputpath + f"/7_stats__{run_tag}.txt"
    with open(stats_path, 'a') as file:
        file.write('\nNumber of CF wind classes: '+ str(classnumCF) + 
                       '\nNumber of (co)spectra used per wind class: ' +
                       str(classize_ws)+ '\n')

        # Write outlier removal details
        file.write('\nOutliers removed per wind class (Lorentz):\n')
        for class_id, info in outlier_info_L.items():
            file.write(f"Wind class {class_id} (range {info['wind_speed_range'][0]:.2f}-{info['wind_speed_range'][1]:.2f} m/s):\n")
            file.write(f"  Removed values: {info['removed_outliers']}\n")

        file.write('\nOutliers removed per wind class (Gauss):\n')
        for class_id, info in outlier_info_G.items():
            file.write(f"Wind class {class_id} (range {info['wind_speed_range'][0]:.2f}-{info['wind_speed_range'][1]:.2f} m/s):\n")
            file.write(f"  Removed values: {info['removed_outliers']}\n")
       
    
    # Delete intermediate vectors
    del matCF
    del CFLL_M
    del CFLL_H
    del CFLL_L
    del CFLG_M
    del CFLG_H
    del CFLG_L
    del Icoscf
    del RcoscfL_M
    del RcoscfL_H
    del RcoscfL_L
    del RcoscfG_M
    del RcoscfG_H
    del RcoscfG_L
    del RcoscfL_MS
    del RcoscfL_HS
    del RcoscfL_LS
    del RcoscfG_MS
    del RcoscfG_HS
    del RcoscfG_LS

    return LUT_CF, mat_meteo_cf, classize_ws



