# -*- coding: utf-8 -*-

"""Look-up table construction for cut-off frequencies (cof)."""

from .FREQCOR_functions import fun_Gauss, fun_Lorentz, fun_Lorentz_peltola
from . import  as fcof
from . import  as fqplt
from datetime import datetime
import numpy as np
import os
import pandas as pd

# Custom exception class for FREQCOR_LUT_cof errors
class FREQCORLUTCofError(Exception):
    """Exception raised for errors in `FREQCOR_LUT_cof`."""
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)

def FREQCOR_LUT_cof(sps, gss, classnum_first, classnum_ws, Icon_cof, Rcon_cof, WS, WD, Zeta, 
                    meteo_df, nfreq, jtmin_hz,jtmax_hz,
                    varclim,freq,freqn,tfmin_hz,tfmax_hz,plot,
                    outputpath, tf_peltola, run_tag):
    """
    Build the cut-off frequency (cof) look-up table.

    The routine groups (co)spectra by a first sorting variable (wind direction for CO2
    or relative humidity for H2O) and then by wind speed classes. For each class it
    estimates transfer functions and derives cut-off frequencies stored in a LUT.

    Parameters
    ----------
    sps : int
        Approach selector (cospectral or spectral).
    gss : int
        Gas species selector.
    classnum_first : int
        Number of classes for the first sorting variable.
    classnum_ws : int
        Number of wind speed classes.
    Icon_cof, Rcon_cof : pandas.DataFrame
        Ideal content (cof) and real content (cof) (co)spectral data (either cospectra or spectra; depends on `sps`).
    WS, WD, Zeta : pandas.Series
        Wind speed, wind direction and stability parameter time series.
    meteo_df : pandas.DataFrame
        Half-hourly meteorological and flux data (used for RH when `gss == 2`).
    nfreq : int
        Number of frequency bins.
    jtmin_hz, jtmax_hz : float
        Frequency range (Hz) used for transfer function normalization.
    varclim : float
        Threshold for coefficient of variation of the normalization parameter.
    freq : pandas.Series
        Natural frequencies (Hz).
    freqn : pandas.DataFrame
        Normalized frequencies per time step.
    tfmin_hz, tfmax_hz : float
        Frequency range (Hz) used for transfer function fitting.
    plot : list[int]
        Plotting options selector.
    outputpath : str
        Output directory.
    tf_peltola : int
        If 1, use the Peltola et al. (2021) fit form for the transfer function.
    run_tag : str
        Run identifier used to generate output filenames.

    Returns
    -------
    LUT_cof : dict
        Cut-off frequency look-up table.
    classize : int
        Number of (co)spectra used per wind speed class.

    Raises
    ------
    ValueError
        If `run_tag` is missing.
    FREQCORLUTCofError
        If there is not enough data to compute the transfer function for a class.
    """
    if not run_tag:
        raise ValueError("run_tag is required for output naming")

#%% Main LUT_cof routine    
    hh=0 # hard-coded option to not perform 1/2h computation of cof
    
    # Initialisation mean (co)spectra dataframe which will be saved to csv file
    mean_cosp=pd.DataFrame(index=freq)
    mean_cosp.index.name="Nat. freq."
    # column names basis for saved file (x being the general scalar)
    if sps == 1: cn = ['cosp(wt)','cosp(wx)'] 
    else: cn = ['sp(t)','sp(cx)'] 
    first_key = Rcon_cof.keys()[0]
    last_key = Rcon_cof.keys()[-1]
    start_date = first_key[:8]  
    end_date = last_key[:8]     
    
    # Keep date and time in row name for plotting
    if gss != 2:
        metcon = pd.concat([WS,WD],axis = 1).set_index(Rcon_cof.keys())    
    elif gss == 2:
        metcon = pd.concat([WS,WD,meteo_df.loc[:,'RH']],axis = 1).set_index(Rcon_cof.keys())
    
    valid_mask = np.isfinite(Rcon_cof).all(axis=0)
    metcon = metcon[valid_mask]
    Icon_cof = Icon_cof.loc[:,valid_mask]
    Rcon_cof = Rcon_cof.loc[:,valid_mask]
    zL = Zeta[valid_mask.reset_index(drop=True)]
    wd = metcon.wind_dir

    # Prepare column name for sorted data
    freq_labels = [f'f{i+1}' for i in range(nfreq)]
    ideal_cols = [f'Icos_{f}' for f in freq_labels]
    real_cols = [f'Rcos_{f}' for f in freq_labels]

    #%%% Individual half-hours
    if not(gss==2) and hh==1:
        freqn = freqn.loc[:,Rcon_cof.columns]
        ws = metcon.wind_speed
        [x_h,y_h], cofmat_h, fnmat_h, trafun_h, denoising = fcof.FREQCOR_cof(Icon_cof, Rcon_cof, jtmin_hz,
                                                             jtmax_hz,varclim,freq,
                                                             freqn,sps,zL,tfmin_hz,tfmax_hz,
                                                             plot,ws,outputpath, 
                                                             wd, tf_peltola, run_tag=run_tag)
        del ws

#%%% LUT approach    
    sizemat = np.isfinite(metcon.loc[:, 'wind_dir']).sum()
    all_cols = list(metcon.columns) + ideal_cols + real_cols
    matsort = pd.concat([metcon,Icon_cof.T,Rcon_cof.T],axis=1)
    matsort.columns = all_cols
    cof_cols = ['ws_mean','ws_max','cof_L','unc_L','cof_G','unc_G','fn_L','uncfn_L','fn_G','uncfn_G']
    
    if not(gss==2):
        li = list(np.linspace(1, classnum_first, classnum_first, dtype=int))
        LUT_cof = dict.fromkeys(li)
        all_individual_cosp = {
            'freq': freq,
            'Rcos': Rcon_cof.values if hasattr(Rcon_cof, 'values') else Rcon_cof,  
            'Icos': Icon_cof.values if hasattr(Icon_cof, 'values') else Icon_cof  
        }
        
        if plot[2] == 1:
            run_suffix = run_tag
            fqplt.plot_all_individual_cosp(all_individual_cosp['freq'], 
                                           all_individual_cosp['Rcos'], 
                                           all_individual_cosp['Icos'], 
                                           plot, outputpath, gas_type='co2',
                                           file_tag=f"2_all_individual_co2__{run_suffix}")
        
        # WD preprocessing (apply symmetry: map [180â€“360Â°] to [0â€“180Â°])
        matsort['wd_reduced'] = matsort['wind_dir'] % 360
        matsort['wd_reduced'] = matsort['wd_reduced'].apply(lambda wd: wd - 180 if wd > 180 else wd)
        bin_edges = np.linspace(0, 180, classnum_first + 1)
        wd_bin_labels = list(range(1, classnum_first + 1))
        matsort['wd_bin'] = pd.cut(matsort['wd_reduced'], bins=bin_edges, labels=wd_bin_labels, include_lowest=True)
        
        for jc in range(1, classnum_first + 1):
            all_tf_data = {
                'x': [],
                'y': [],
                'cofmat': [],
                'fnmat': [],
                'original_trafun': [],
                'original_freq': [],
                'ws': [],
                'num_classes': classnum_ws
            }
                
            all_cosp_data = {
                'freq': freq,
                'Rcosv': [],
                'Icosv': [],
                'ws': [],
                'zL': zL,
                'num_classes': classnum_ws
            }
            
            denoising = {}
            
            wd_class_df = matsort[matsort['wd_bin'] == jc].copy().reset_index(drop=True)
            sizemat_wd = len(wd_class_df)
            classizews = int(np.floor(sizemat_wd / classnum_ws))
            if classizews < 2:
                raise FREQCORLUTCofError("Error in LUT_cof: Not enough data to compute the transfer function")
            
        
            lut_cof_temp = {
                'wd': pd.DataFrame(),
                'ws': pd.DataFrame(np.zeros(shape=(classnum_ws, len(cof_cols))), columns=cof_cols)
            }
        
            lut_cof_temp['wd'].loc[0, 'wdmean'] = wd_class_df['wd_reduced'].mean()
            lut_cof_temp['wd'].loc[0, 'wdmax'] = wd_class_df['wd_reduced'].max()
            lut_cof_temp['wd'].loc[0, 'uncclass'] = wd_class_df['wd_reduced'].std()
        
            wd_class_df_sorted = wd_class_df.sort_values(by='wind_speed').reset_index(drop=True)
            
            matsortspI = wd_class_df_sorted[ideal_cols]
            matsortspR = wd_class_df_sorted[real_cols]
            matsortspI.columns = range(matsortspI.shape[1])
            matsortspR.columns = range(matsortspR.shape[1])
            
            ws_bins = np.array_split(wd_class_df_sorted.index, classnum_ws)
            for jws, bin_idx in enumerate(ws_bins, 1):
                bin_ws = wd_class_df_sorted.loc[bin_idx]
                Icosv = np.mean(matsortspI.iloc[bin_ws.index, :], axis=0).T
                Rcosv = np.mean(matsortspR.iloc[bin_ws.index, :], axis=0).T
                ws = bin_ws['wind_speed']

                [x, y], cofmat, fnmat, trafun, denoising_out = fcof.FREQCOR_cof(
                    Icosv, Rcosv, jtmin_hz, jtmax_hz, varclim, freq, freqn,
                    sps, zL, tfmin_hz, tfmax_hz, plot, ws, outputpath,
                    wd, tf_peltola, jc, run_tag=run_tag
                )

                lut_cof_temp['ws'].loc[jws - 1, 'ws_mean'] = ws.mean()
                lut_cof_temp['ws'].loc[jws - 1, 'ws_max'] = ws.max()
                lut_cof_temp['ws'].loc[jws - 1, 'cof_L'] = cofmat['cof_L']
                lut_cof_temp['ws'].loc[jws - 1, 'unc_L'] = cofmat['unc_L']
                lut_cof_temp['ws'].loc[jws - 1, 'cof_G'] = cofmat['cof_G']
                lut_cof_temp['ws'].loc[jws - 1, 'unc_G'] = cofmat['unc_G']
                lut_cof_temp['ws'].loc[jws - 1, 'fn_L'] = fnmat['fn_L']
                lut_cof_temp['ws'].loc[jws - 1, 'uncfn_L'] = fnmat['uncfn_L']
                lut_cof_temp['ws'].loc[jws - 1, 'fn_G'] = fnmat['fn_G']
                lut_cof_temp['ws'].loc[jws - 1, 'uncfn_G'] = fnmat['uncfn_G']

                denoising[jws-1]=denoising_out

                col_suffix = f'_wd{jc}_ws{jws}'
                mean_cosp.loc[:, cn[0] + col_suffix] = Icosv.values
                mean_cosp.loc[:, cn[1] + col_suffix] = Rcosv.values

                all_tf_data['x'].append(x)
                all_tf_data['y'].append(y)
                all_tf_data['cofmat'].append(cofmat)
                all_tf_data['fnmat'].append(fnmat)
                all_tf_data['ws'].append(ws.mean())
                all_tf_data['original_trafun'].append(trafun)
                all_tf_data['original_freq'].append(freq)
                all_cosp_data['Rcosv'].append(Rcosv)
                all_cosp_data['Icosv'].append(Icosv)
                all_cosp_data['ws'].append(ws.mean())
                
            LUT_cof[jc] = lut_cof_temp
        
            if plot[1] == 1:
                norm_range=[jtmin_hz, jtmax_hz]
                class_pad = len(str(classnum_ws))
                jc_pad = len(str(classnum_first))
                run_suffix = run_tag
                if sps == 1 and tf_peltola==1:
                    fqplt.plot_TF_unified(fun_Lorentz_peltola, fun_Gauss, None, None, None,
                                        denoising, plot, outputpath, norm_range, None, jc, 
                                        gas_type='co2', all_classes_data=all_tf_data, file_tag=f"4_mean_TF_all_classes__wd{jc:0{jc_pad}d}__{run_suffix}")
                else:
                    fqplt.plot_TF_unified(fun_Lorentz, fun_Gauss, None, None, None,
                                        denoising, plot, outputpath, norm_range, None, jc, 
                                        gas_type='co2', all_classes_data=all_tf_data, file_tag=f"4_mean_TF_all_classes__wd{jc:0{jc_pad}d}__{run_suffix}")
                
                fqplt.plot_cosp_unified(freq, None, None, None, denoising, zL, sps,
                                          None, jc, plot, outputpath, 'co2',
                                          all_classes_data=all_cosp_data, file_tag=f"4_mean_cosp_all_classes__wd{jc:0{jc_pad}d}__{run_suffix}")
               
        classize=classizews 

    elif gss==2:
        matsort = matsort.sort_values(by = ['RH']).reset_index(drop=True)
        li = list(np.linspace(1,classnum_first,classnum_first,dtype=int))
        LUT_cof = dict.fromkeys(li) 
        
        all_individual_cosp = {
            'freq': freq,
            'Rcos': Rcon_cof.values if hasattr(Rcon_cof, 'values') else Rcon_cof,  
            'Icos': Icon_cof.values if hasattr(Icon_cof, 'values') else Icon_cof  
        }
        
        if plot[2] == 1:
            run_suffix = run_tag
            fqplt.plot_all_individual_cosp(all_individual_cosp['freq'], 
                                           all_individual_cosp['Rcos'], 
                                           all_individual_cosp['Icos'], 
                                           plot, outputpath, gas_type='h2o',
                                           file_tag=f"2_all_individual_h2o__{run_suffix}")
        
        classize_rh=int(np.floor(sizemat/classnum_first)) 
        classizews=int(np.floor(classize_rh/classnum_ws)) 
        if classizews < 2:
            raise FREQCORLUTCofError("Error in LUT_cof: Not enough data to compute the transfer function")
        
        # Use np.array_split to divide sorted DataFrame into RH classes
        rh_bins = np.array_split(matsort.sort_values(by=['RH']).reset_index(drop=True), classnum_first)
        for jc, bin_rh in enumerate(rh_bins, 1):
            lut_cof_temp={
                'rh':pd.DataFrame(),
                'ws': pd.DataFrame(np.zeros(shape=(classnum_ws, len(cof_cols))), columns=cof_cols)
            }
            
            all_tf_data = {
                'x': [],
                'y': [],
                'cofmat': [],
                'fnmat': [],
                'original_trafun': [],
                'original_freq': [],
                'ws': [],
                'num_classes': classnum_ws
            }
            
            all_cosp_data = {
                'freq': freq,
                'Rcosv': [],
                'Icosv': [],
                'ws': [],
                'zL': zL,
                'num_classes': classnum_ws
            }
            
            denoising = {}
            
            lut_cof_temp['rh'].loc[0,'rhmean']=np.mean(bin_rh['RH'])             
            lut_cof_temp['rh'].loc[0,'rhmax']=np.max(bin_rh['RH'])             
            lut_cof_temp['rh'].loc[0,'uncclass']=np.std(bin_rh['RH']) 
            
            matsortws=bin_rh.sort_values(by = ['wind_speed']).reset_index(drop = True)
            
            matsortspI = matsortws[ideal_cols]
            matsortspR = matsortws[real_cols]
            matsortspI.columns = range(matsortspI.shape[1])
            matsortspR.columns = range(matsortspR.shape[1])
            # Use np.array_split to divide sorted DataFrame into wind speed classes
            ws_bins = np.array_split(matsortws.index, classnum_ws)
            for jws, bin_idx in enumerate(ws_bins, 1):
                bin_ws = matsortws.loc[bin_idx]
                Icosv = np.mean(matsortspI.iloc[bin_ws.index, :], axis=0).T
                Rcosv = np.mean(matsortspR.iloc[bin_ws.index, :], axis=0).T
                ws = bin_ws['wind_speed']

                [x,y], cofmat, fnmat, trafun, denoising_out = fcof.FREQCOR_cof(
                    Icosv, Rcosv, jtmin_hz,jtmax_hz,varclim,freq,freqn,
                    sps,zL,tfmin_hz,tfmax_hz,plot,ws,outputpath,wd,
                    tf_peltola, jc, run_tag=run_tag)

                lut_cof_temp['ws'].loc[jws - 1, 'ws_mean'] = ws.mean()
                lut_cof_temp['ws'].loc[jws - 1, 'ws_max'] = ws.max()
                lut_cof_temp['ws'].loc[jws - 1, 'cof_L'] = cofmat['cof_L']
                lut_cof_temp['ws'].loc[jws - 1, 'unc_L'] = cofmat['unc_L']
                lut_cof_temp['ws'].loc[jws - 1, 'cof_G'] = cofmat['cof_G']
                lut_cof_temp['ws'].loc[jws - 1, 'unc_G'] = cofmat['unc_G']
                lut_cof_temp['ws'].loc[jws - 1, 'fn_L'] = fnmat['fn_L']
                lut_cof_temp['ws'].loc[jws - 1, 'uncfn_L'] = fnmat['uncfn_L']
                lut_cof_temp['ws'].loc[jws - 1, 'fn_G'] = fnmat['fn_G']
                lut_cof_temp['ws'].loc[jws - 1, 'uncfn_G'] = fnmat['uncfn_G']

                denoising[jws-1]=denoising_out
                col_suffix = f'_rh{jc}_ws{jws}'

                mean_cosp.loc[:,cn[0]+col_suffix]=Icosv.values
                mean_cosp.loc[:,cn[1]+col_suffix]=Rcosv.values

                valid = np.isfinite(x) & np.isfinite(y)

                all_tf_data['x'].append(x[valid])
                all_tf_data['y'].append(y[valid])
                all_tf_data['cofmat'].append(cofmat)
                all_tf_data['fnmat'].append(fnmat)
                all_tf_data['ws'].append(ws.mean())
                all_tf_data['original_trafun'].append(trafun)
                all_tf_data['original_freq'].append(freq)
                all_cosp_data['Rcosv'].append(Rcosv)
                all_cosp_data['Icosv'].append(Icosv)
                all_cosp_data['ws'].append(ws.mean())
            
            LUT_cof[jc]=lut_cof_temp
            
            if plot[1] == 1:
                norm_range=[jtmin_hz, jtmax_hz]
                class_pad = len(str(classnum_ws))
                jc_pad = len(str(classnum_first))
                run_suffix = run_tag
                if sps == 1 and tf_peltola==1:
                    fqplt.plot_TF_unified(fun_Lorentz_peltola, fun_Gauss, None, None, None,
                                        denoising, plot, outputpath, norm_range, None, jc, 
                                        gas_type='h2o', all_classes_data=all_tf_data, file_tag=f"4_mean_TF_all_classes__rh{jc:0{jc_pad}d}__{run_suffix}")
                else:
                    fqplt.plot_TF_unified(fun_Lorentz, fun_Gauss, None, None, None,
                                        denoising, plot, outputpath, norm_range, None, jc, 
                                        gas_type='h2o', all_classes_data=all_tf_data, file_tag=f"4_mean_TF_all_classes__rh{jc:0{jc_pad}d}__{run_suffix}")
                
                fqplt.plot_cosp_unified(freq, None, None, None, denoising, zL, sps,
                                          None, jc, plot, outputpath, 'h2o',
                                          all_classes_data=all_cosp_data, file_tag=f"4_mean_cosp_all_classes__rh{jc:0{jc_pad}d}__{run_suffix}")
               
        classize=classizews

   
    run_suffix = run_tag
    csv_filename = f"4_mean_cosp__{run_suffix}.csv"
    mean_cosp.to_csv(os.path.join(outputpath, csv_filename))
    
#%% Write stats to txt file
    now = datetime.now()
    current_time = now.strftime("%d/%m/%Y %H:%M:%S")
    stats_path = outputpath + f"/7_stats__{run_tag}.txt"
    with open(stats_path, 'a') as file:
        file.write('\n' + current_time + 
                   '\nNumber of cof wind classes: '+ str(classnum_ws) + 
                   '\nNumber of (co)spectra used per wind class: ' + 
                   str(classize))
    del now, current_time
#%% Uncertainty estimation    
    # Estimate of the uncertainty due to cof class width    
    # Why only for Lorentz?? And using min wind speed?
    # slo_cof = pd.Series(np.zeros(classnum_ws))
    # for jsm in range(1,classnum_ws+1):
    #     if classnum_ws==1:
    #         LUT_cof.iloc[jsm-1,3]=matcof.iloc[jsm,1]
    #     else:
    #         if jsm==1:
    #             slo_cof.iloc[jsm-1]=((LUT_cof.iloc[1,2]-LUT_cof.iloc[0,2])/(LUT_cof.iloc[1,0]-LUT_cof.iloc[0,0]))
    #         elif jsm==classnum_ws:
    #             slo_cof[jsm-1]=((LUT_cof.iloc[classnum_ws,2]-LUT_cof.iloc[classnum_ws-1,2])/(LUT_cof.iloc[classnum_ws,0]-LUT_cof.iloc[classnum_ws-1,0]))
    #         else:
    #             slo_cof[jsm-1]=((LUT_cof.iloc[jsm+1,2]-LUT_cof.iloc[jsm-1,2])/(LUT_cof.iloc[jsm+1,0]-LUT_cof.iloc[jsm-1,0]))
                
    #     LUT_cof.iloc[jsm-1,3]=(matcof.iloc[jsm-1,2]**2+(uncclass[jsm-1]*slo_cof[jsm-1])**2)**0.5
    
    # Delete useless intermediate vectors
    # del cospcon
    del metcon
    del cofmat
    return LUT_cof, classize
