
# -*- coding: utf-8 -*-
"""
Created on Dec 28 2022

@author: Bernard Heinesch

This subroutine is the FREQCOR_Read.py routine, adapted for BVOC data coming from GEddySoft 

Inputs :
      - GEddySoft daily output files
Outputs : 
      -(co)spectra matrices: Icon_raw(nfreq, nspec), Rcon_raw(nfreq, nspec) (either cospectra or spectra; depends on `sps`);
      - meteo and flux: x vectors(nspec);
      - General Flags: x vectors(nspec);

Comments :
    - ! the mz choice is hard-coded, at the beginning of the routine !
    - some variables are missing in the GEddySoft output but are needed by FREQCOR for filtering steps. 
      They are therefore created, filled with values that will end in no-filtering, and returned, to mimic
      FREQCOR_Read
TO DO:
    - filename selection and target_mz should not be hard-coded but aim was to not modify the FREQCOR structure
"""

import glob
import numpy as np
from FREQCOR_VM_flag import FREQCOR_VM_flag
import pandas as pd
import h5py

def FREQCOR_Read_TOF(config, jsite):

    # %%  Initialisation

    # USER CHOICE !
    # mz to be selected (will choose the closest one)
    # 33.033 : proposed flux treshold of 0.1
    # 45.033 : proposed flux treshold of 0.01
    # 69.0699 : proposed flux treshold of 0.5
    # 137.132 : proposed flux treshold of 0.2
    target_mz = 69.0699

    # Map config names to variables used in the script
    dirname = config['IO']['input_path']
    massmanpath = config['IO']['massman_path']
    fcname = config['IO']['binned_cosp']
    foname = config['IO']['flx_meteo']
    sps = int(config['PROCEDURE_OPTIONS']['sps'])
    gss = int(config['PROCEDURE_OPTIONS']['gss'])
    zmeas = float(config['EC_SETUP']['zmeas']) 
    disph = float(config['EC_SETUP']['disph']) 
    
    total_to_process = glob.glob(dirname + fcname)

    # get total number of available (co)spectra 
    nspec = 0
    for filen in total_to_process:
        with h5py.File(filen, 'r') as hdf5_f:

            if 'TRACER' in hdf5_f:
                name_key = hdf5_f['TRACER'].keys()
                if len(name_key):  # file with TRACER values
                    ds = hdf5_f['time']
                    nspec = nspec + ds.shape[0]
    print(str(nspec) + ' nspec available')

    # initialize vectors and matrices
    timestamp = pd.Series(dtype='float64')
    freqm = pd.DataFrame(dtype='float64')
    Icon_raw = pd.DataFrame(dtype='float64')
    Rcon_raw = pd.DataFrame(dtype='float64')
    H = pd.DataFrame(dtype='float64')
    Fc = pd.DataFrame(dtype='float64')
    FcEP = pd.DataFrame(dtype='float64')
    LE = pd.DataFrame(dtype='float64')
    RH = pd.DataFrame(dtype='float64')
    meteo_df = pd.DataFrame(dtype='float64')
    WS = pd.Series(dtype='float64', name='wind_speed')
    WD = pd.Series(dtype='float64', name='wind_dir')
    Ustar = pd.Series(dtype='float64', name='u*')
    Zeta = pd.Series(dtype='float64', name='(z-d)/L')
    OOR_w = pd.Series(dtype='float64', name='(z-d)/L')
    num_spikes_w = pd.Series(dtype='float64', name='spikes_w')

    # %% loop over GEddySoft daily output files

    # loop over the daily files in the given directory matching the name pattern
    k = 0
    for filen in total_to_process:
    
        with h5py.File(filen, 'r') as hdf5_f:

            if 'TRACER' in hdf5_f:
                name_key = hdf5_f['TRACER'].keys()
                if len(name_key):  # file with TRACER values
    
                    # find the mz closest to target_mz 
                    for x in name_key:
                        if abs(float(hdf5_f['TRACER'][x]['mz'][()])-target_mz) < 0.01:
                            target_channel = x
                            break
                        
                    # find the number of 1/2h 
                    ds = hdf5_f['time']
                    n_hh = ds.shape[0]
        
                    # Create all the tables that will contain the 1/2 (co)spectra or frequency and the FLUX, METEO and QUALITY FLAGS
                    # Each colum will correspond to one specific 1/2 hour
                    # Dimensions: max nb of bins x nb of (co)spectra
            
                    # timestamp = timestamp.append(pd.Series(np.transpose(np.array(hdf5_f['time']).astype('U17'))),ignore_index=True)
                    timestamp = pd.concat([timestamp, pd.Series(np.array(hdf5_f['time']).astype('U17'))], ignore_index=True)                    

                    # Natural frequency of each bin
                    freqm_d = pd.DataFrame(np.transpose(np.array(hdf5_f['freq'])))
                    freqm_d = freqm_d[np.repeat(freqm_d.columns.values,n_hh)]
                    freqm = pd.concat([freqm,freqm_d],axis=1,ignore_index=True)
            
                    # Ideal (co)spectral data from file (sensible heat)
                    Icon_raw = pd.concat([Icon_raw,pd.DataFrame(np.transpose(np.array(hdf5_f['MET']['cospec_wT_scaled'])))],axis=1,ignore_index=True)
                    # Real (co)spectral data from file (CO2 or H2O)
                    Rcon_raw = pd.concat([Rcon_raw,pd.DataFrame(np.transpose(np.array(hdf5_f['TRACER'][target_channel]['cospec_scaled'])))],axis=1,ignore_index=True)
        
                    # # Creation of a flux and meteo dataframe with all the variables of interest (*1206 for wT because converted from m-1 K to Wm-2)
                    meteo_df_d = pd.concat([pd.DataFrame(np.nan_to_num(np.transpose(np.array(hdf5_f['MET']['wT'])*1206))),
                                            pd.DataFrame(np.nan_to_num(np.transpose(np.array(hdf5_f['TRACER'][target_channel]['flux'])))),
                                            pd.DataFrame(np.zeros(shape=(n_hh))),
                                            pd.DataFrame(np.zeros(shape=(n_hh))),
                                            pd.DataFrame(np.zeros(shape=(n_hh)))]
                                           ,axis=1,ignore_index=True)
                    meteo_df_d.iloc[:,2:6] = 999
                    meteo_df = pd.concat([meteo_df,meteo_df_d],axis=0,ignore_index=True)

                    WS = pd.concat([WS, pd.DataFrame(np.transpose(np.array(hdf5_f['MET']['wsh'])), columns=['wind_speed']).squeeze()], ignore_index=True)
                    WD = pd.concat([WD, pd.DataFrame(np.transpose(np.array(hdf5_f['MET']['wdir'])), columns=['wind_dir']).squeeze()], ignore_index=True)
                    Ustar = pd.concat([Ustar, pd.DataFrame(np.transpose(np.array(hdf5_f['MET']['ust'])), columns=['u*']).squeeze()], ignore_index=True)
                    Zeta = pd.concat([Zeta, pd.DataFrame(np.transpose(np.array(hdf5_f['MET']['zoL'])), columns=['(z-d)/L']).squeeze()], ignore_index=True)
                    OOR_w = pd.concat([OOR_w, pd.DataFrame(np.transpose(np.array(hdf5_f['MET']['qaqc']['IPT_w'][:,11])), columns=['OOR_w']).squeeze()], ignore_index=True)
                    num_spikes_w = pd.concat([num_spikes_w, pd.DataFrame(np.transpose(np.array(hdf5_f['MET']['qaqc']['num_spikes_w'])), columns=['num_spikes_w']).squeeze()], ignore_index=True)                    
                    k = k + n_hh
                
    # Rename the columns in the DataFrames with the info on the date and
    # time of the 1/2 hour to which it corresponds.
    timestamp = timestamp.replace({'-':''}, regex=True)
    timestamp = timestamp.replace({' ':'-'}, regex=True)
    freqm.columns = timestamp.iloc[:]
    Icon_raw.columns = timestamp.iloc[:]
    Rcon_raw.columns = timestamp.iloc[:]

    # %% Preparation of variables

    meteo_df.columns = ['H','Fc', 'FcEP', 'LE', 'RH']

    # Add debugging prints for dimensions
    print("\nDimension check:")
    print(f"freqm shape: {freqm.shape}")
    print(f"WS shape before transpose: {WS.shape}")
    print(f"WS shape after transpose: {np.transpose(WS).shape}")
    print(f"zmeas - disph: {zmeas - disph}")

    # Normalised frequency of each bin
    freqn = np.divide(freqm * (zmeas - disph), WS.values.reshape(1, -1))
        
    # Division of (co)spectral density by frequency (because GEddysoft output gives 
    # f_nat*cospec/cov or f_nat*spec/var) 
    Icon_raw=Icon_raw.div(freqm)
    Rcon_raw=Rcon_raw.div(freqm)

    # FlagVM_sp = np.zeros(shape=(nspec)); FlagVM_sp[:] = 999  # Vickers and Mahrt flag spikes
    # FlagVM_do = np.zeros(shape=(nspec)); FlagVM_do[:] = 999  # Vickers and Mahrt flag drop out
    # FlagVM_skh = np.zeros(shape=(nspec)); FlagVM_skh[:] = 999  # Vickers and Mahrt flag skewness hard
    # FlagVM_sks = np.zeros(shape=(nspec)); FlagVM_sks[:] = 999  # Vickers and Mahrt flag skewness soft
    # FlagVM_dh = np.zeros(shape=(nspec)); FlagVM_dh[:] = 999  # Vickers and Mahrt flag discontinuitites hard
    # FlagVM_ds = np.zeros(shape=(nspec)); FlagVM_ds[:] = 999  # Vickers and Mahrt flag discontinuitites soft
    # FlagF_H = np.zeros(shape=(nspec)); FlagF_H[:] = 999  # Foken sensible heat
    # FlagF_g = np.zeros(shape=(nspec)); FlagF_g[:] = 999  # Foken CO2
    FlagF_H = pd.Series(np.full(nspec,1.))#; FlagF_H[:] = 999  # Foken sensible heat
    FlagF_g = pd.Series(np.full(nspec,1.))#; FlagF_g[:] = 999  # Foken CO2

    # Remove freq bins fully filled with Nans
    idx = Icon_raw.notna().any(axis=1) & Rcon_raw.notna().any(axis=1)
    freqm = freqm[idx]
    freqn = freqn[idx]
    Icon_raw = Icon_raw[idx]
    Rcon_raw = Rcon_raw[idx]
    
    # Remove 1/2h fully filled with Nans + additional data filtering
    idx = (Icon_raw.notna().any(axis=0) & 
           Rcon_raw.notna().any(axis=0) & 
           ~WS.isnull().values &
           #(num_spikes_w > 0).values  &
           (OOR_w == 0).values)

    freqm = freqm.loc[:,idx]
    freqn = freqn.loc[:,idx]
    Icon_raw = Icon_raw.loc[:,idx]
    Rcon_raw = Rcon_raw.loc[:,idx]
    idx = idx.reset_index(drop=True)
    meteo_df = meteo_df.loc[idx,:]
    WS = WS[idx]
    WD = WD[idx]
    Ustar = Ustar[idx]
    Zeta = Zeta[idx]
    OOR_w = OOR_w[idx]
    num_spikes_w = num_spikes_w[idx]
    FlagF_H = FlagF_H[idx]
    FlagF_g = FlagF_g[idx]

    nspec = sum(idx)
    
    # # Reset the indexes so to have a continuous integer series of indexes
    freqm = freqm.reset_index(drop=True)
    freqn = freqn.reset_index(drop=True)
    Icon_raw = Icon_raw.reset_index(drop=True)    
    Rcon_raw = Rcon_raw.reset_index(drop=True)    
    meteo_df = meteo_df.reset_index(drop=True)
    WS = WS.reset_index(drop=True)
    WD = WD.reset_index(drop=True)
    Ustar = Ustar.reset_index(drop=True)
    Zeta = Zeta.reset_index(drop=True)
    
    # Series of natural frequencies, taken as the first column of the freqm
    # DataFrame where all the columns should contain the same values (see here above)
    freq = freqm.copy()
    freq = freq.iloc[:,0]
    # Same for normalised frequencies
    # freqn = freqn.iloc[:,0]

    # Cleaned (co)spectra vector length (number of actually present frequency bins)
    nfreq = freq.shape[0]
    
    # %% Creation of a global quality flag for VM test
    # Calls FREQCOR_VM_flag
    # FlagVM_w, FlagVM_T, FlagVM_g = FREQCOR_VM_flag(NumberDataFLAGF, FlagVM_sp, FlagVM_do, FlagVM_skh,FlagVM_dh, gss)
    FlagVM_w = np.full(nspec,0)
    FlagVM_T = np.full(nspec,0)
    FlagVM_g = np.full(nspec,0)
    
    #%% Read Massman coefficients from file
    massman_coef = pd.read_csv(massmanpath, header=0)
    massman_coef = massman_coef.set_index('var')
    massman_coef = massman_coef.loc[:,jsite]
    
    A0 = [float(x) for x in massman_coef['A0'].split(',')]
    kf0 = [float(x) for x in massman_coef['kf0'].split(',')]
    mu = [float(x) for x in massman_coef['mu'].split(',')]

    massman_coef = [A0, kf0, mu]

    #%%
    # Suppression of useless intermediate variables
    del freqm
    del freqm_d
    del meteo_df_d
    # del FlagVM_sp
    # del FlagVM_do
    # del FlagVM_skh
    # del FlagVM_sks
    # del FlagVM_ds
    # del FlagVM_dh

    # create a copy of Icon_raw for further use
    Icosp_raw = Icon_raw.copy()
    
    # import matplotlib.pyplot as plt
    # plt.plot(meteo_df.loc[:,'Fc'])

    return [nspec, Icon_raw, Rcon_raw, Icosp_raw, nfreq, freq, freqn, meteo_df, WS, WD,
            Ustar, Zeta, FlagF_H, FlagF_g, FlagVM_w, FlagVM_T, FlagVM_g, massman_coef]


