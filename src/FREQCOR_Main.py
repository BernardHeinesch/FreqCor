# -*- coding: utf-8 -*-

"""Main routine entry point for running the FREQCOR processing chain."""

# Custom exception class for FREQCOR_Main errors
from FREQCOR_Compute import FREQCOR_Compute
from FREQCOR_Flux import FREQCOR_Flux
from FREQCOR_Read_EP import FREQCOR_Read_EP
from FREQCOR_Read_TOF import FREQCOR_Read_TOF
from FREQCOR_Sel_general import FREQCOR_Sel_general
from FREQCOR_validate import print_validation_results, validate_config
from FREQCOR_write_outputs import FREQCOR_write_outputs
from datetime import datetime
import numpy as np
import os
import pandas as pd

class FREQCORMainError(Exception):
    """Exception raised for errors in the FREQCOR_Main function.
    
    Attributes:
        message -- explanation of the error
    """
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)

def FREQCOR_Main(config):
    """
    Run the full FREQCOR processing chain for one configuration.

    Parameters
    ----------
    config : configparser.ConfigParser
        Parsed configuration.

    Returns
    -------
    LUTcof : dict or None
        Look-up table of cut-off frequencies, with two-levels
        sorting.
    LUTCF_u : dict or None
        Look-up table of correction factors for unstable conditions, with two-levels
        sorting.
    LUTCF_s : dict or None
        Look-up table of correction factors for stable conditions, with two-levels
        sorting.

    Raises
    ------
    FREQCORMainError
        If an unrecoverable error occurs during processing.
    """

    
    # %% Validate configuration parameters
    valid, errors = validate_config(config)
    print_validation_results(valid, errors)
    
    if not valid:
        print("Configuration validation failed. Please fix the errors and try again.")
        return None, None, None
    
    # %% Variables extraction
    jsite = config['EC_SETUP']['jsite']
    gss = int(config['PROCEDURE_OPTIONS']['gss']) 
    sps = int(config['PROCEDURE_OPTIONS']['sps']) 
    zmeas = float(config['EC_SETUP']['zmeas']) 
    disph = float(config['EC_SETUP']['disph'])  
    rsz = float(config['EC_SETUP']['rsz'])  
    rseast = float(config['EC_SETUP']['rseast'])
    rsnorth = float(config['EC_SETUP']['rsnorth'])  
    hl_cof_min = float(config['USER_LIMITS']['hl_cof_min'])
    hl_cof_max = float(config['USER_LIMITS']['hl_cof_max'])
    fcl_cof_min = float(config['USER_LIMITS']['fcl_cof_min'])
    fcl_cof_max = float(config['USER_LIMITS']['fcl_cof_max'])
    lel_cof_min = float(config['USER_LIMITS']['lel_cof_min'])
    lel_cof_max = float(config['USER_LIMITS']['lel_cof_max'])
    cospec_abs_limit = float(config['USER_LIMITS']['cospec_abs_limit'])
    hl_cf_u_min = float(config['USER_LIMITS']['hl_cf_u_min'])
    hl_cf_s_max = float(config['USER_LIMITS']['hl_cf_s_max'])
    hl_cf_abs_max = float(config['USER_LIMITS']['hl_cf_abs_max'])
    FVM = int(config['USER_LIMITS']['fvm'])   
    FF = int(config['USER_LIMITS']['ff'])  
    jtmin_hz = float(config['USER_LIMITS']['jtmin'])  
    jtmax_hz = float(config['USER_LIMITS']['jtmax'])  
    varclim = float(config['USER_LIMITS']['varclim'])  
    tfmin_hz = float(config['USER_LIMITS']['tfmin'])  
    tfmax_hz = float(config['USER_LIMITS']['tfmax'])   
    plot = [
        int(config['PROCEDURE_OPTIONS']['plot_hh']),
        int(config['PROCEDURE_OPTIONS']['plot_main']),
        int(config['PROCEDURE_OPTIONS']['plot_aux']),
        int(config['PROCEDURE_OPTIONS']['plot_save']),
    ]
   
    outputpath = config['IO']['output_path']

    gas_map = {1: 'co2', 2: 'h2o', 3: 'o3', 4: 'ch4', 5: 'n2o'}
    gas_tag = gas_map.get(gss, f'gss{gss}')
    method_tag = 'cosp' if sps == 1 else 'sp'
    run_dt = datetime.now().strftime('%y%m%dT%H%M')
    window_tag = 'all'
    if 'TIME_WINDOW' in config and 'enable_time_window' in config['TIME_WINDOW']:
        try:
            if int(config['TIME_WINDOW']['enable_time_window']) == 1:
                start_d = datetime.strptime(config['TIME_WINDOW']['start_datetime'], '%Y-%m-%d %H:%M')
                end_d = datetime.strptime(config['TIME_WINDOW']['end_datetime'], '%Y-%m-%d %H:%M')
                window_tag = f"{start_d.strftime('%y%m%d')}-{end_d.strftime('%y%m%d')}"
        except Exception:
            window_tag = 'all'

    run_tag = f"{jsite}__{gas_tag}__{method_tag}__{window_tag}__{run_dt}"
    
    tf_sonic = int(config['PROCEDURE_OPTIONS']['tf_sonic'])
    tf_peltola = int(config['PROCEDURE_OPTIONS']['tf_peltola'])
    eq = int(config['PROCEDURE_OPTIONS']['eq'])

    read_routine = config['IO'].get('read_routine', None)
    if not read_routine:
        read_routine = config['PROCEDURE_OPTIONS'].get('read_routine', 'FREQCOR_Read_EP')
    
    
    if gss==2:
        classnum_first=int(config['PROCEDURE_OPTIONS']['classnumrh'])
    else:
        classnum_first=int(config['PROCEDURE_OPTIONS']['classnumwd'])
    
    if classnum_first==1:
        
        classnum=int(config['PROCEDURE_OPTIONS']['classnumcofws'])
        classnumCF_u=int(config['PROCEDURE_OPTIONS']['classnumcf_u'])
        classnumCF_s=int(config['PROCEDURE_OPTIONS']['classnumcf_s'])
    else:
        classnum = 4
        classnumCF_u = 4
        classnumCF_s = 4
    
    
    print(f"run_tag={run_tag}")
    print(f"read_routine={read_routine}")

    # %% Read datas
    if read_routine == 'FREQCOR_Read_EP':
        print("\n=== Reading ===")
        print("  Reading inputs")
        nspec, Icon_raw, Rcon_raw, Icosp_raw, nfreq, freq, freqn, freqcount, meteo_df, WS, \
            WD, Ustar, Zeta, FlagF_H, FlagF_g, FlagVM_w, FlagVM_T, FlagVM_g, flag_wd, \
            massman_coef = FREQCOR_Read_EP(config, jsite)
    elif read_routine == 'FREQCOR_Read_TOF':
        print("\n=== Reading ===")
        print("  Reading inputs")
        nspec, Icon_raw, Rcon_raw, Icosp_raw, nfreq, freq, freqn, meteo_df, WS, WD, \
            Ustar, Zeta, FlagF_H, FlagF_g, FlagVM_w, FlagVM_T, FlagVM_g, massman_coef = \
                FREQCOR_Read_TOF(config, jsite)
        freqcount = pd.Series(320 * np.asarray(freq), index=getattr(freq, 'index', None))
        flag_wd = pd.Series(np.nan, index=FlagF_H.index)
    else:
        raise FREQCORMainError(f"Unknown read_routine: {read_routine}")

    print(f"  Read complete: nspec={nspec}, nfreq={nfreq}")

    # %% Filter data
    print("\n=== Filtering ===")
    print("  Filtering data")
    Icon_sel, Rcon_sel, Icosp_sel, meteo_df, WS, WD, Ustar, Zeta, FlagF_H, FlagF_g, \
        FlagVM_T, FlagVM_w, FlagVM_g = FREQCOR_Sel_general(Icon_raw, Rcon_raw, Icosp_raw,
                                                           meteo_df, 
                                                           WS, WD, Ustar, Zeta, 
                                                           FlagF_H,FlagF_g,FlagVM_T,
                                                           FlagVM_w,FlagVM_g, flag_wd, 
                                                           config)
    
    # %%% Compute: cut-off frequency and correction factors
    print("\n=== Computing cof/CF ===")
    print("  Computing cut-off frequencies and correction factors")
    LUT_cof, classize, LUTCF_u, matsortCF_u, classizeCF_u, LUTCF_s, matsortCF_s, classizeCF_s=\
        FREQCOR_Compute(Icon_sel, Rcon_sel, Icosp_sel, nspec, nfreq, gss, sps, 
                        classnum_first, classnum, classnumCF_u, classnumCF_s,   
                        WS, WD, Zeta, meteo_df, hl_cof_min, fcl_cof_min, lel_cof_min,   
                        hl_cf_u_min, hl_cf_s_max, hl_cof_max, fcl_cof_max, lel_cof_max, hl_cf_abs_max,
                        FlagVM_T, FlagVM_g, FlagVM_w, FlagF_H, FlagF_g, FVM, FF,
                        jtmin_hz, jtmax_hz, varclim, freq, freqn, freqcount,
                        tfmin_hz, tfmax_hz, plot, outputpath,
                        massman_coef, tf_sonic, tf_peltola, cospec_abs_limit, run_tag)

    # %% Flux correction 
    print("\n=== Flux correction ===")
    print("  Computing corrected fluxes")
    FccorL, FccorG, CFHL, CF_L, CF_G  = \
        FREQCOR_Flux(nspec, sps, gss, meteo_df,WS,WD,Zeta,LUTCF_u, 
                     classnumCF_u, LUTCF_s, classnumCF_s, zmeas, disph,rsz,
                     rseast,rsnorth,eq)
    
    # %% Write outputs to files
    print("\n=== Writing outputs ===")
    print("  Writing outputs")
    FREQCOR_write_outputs(gss, LUT_cof, LUTCF_u, LUTCF_s, outputpath, meteo_df, FccorL, FccorG,
                          freqn, WS, WD, Zeta, sps, CF_L, 
                          CF_G, CFHL, config, run_tag)
    print("\n=== Done ===")
    print(f"Summary: run_tag={run_tag}; nspec={nspec}; nfreq={nfreq}")
    print(f"output_path: {outputpath}")
    print("Outputs: 4_LUT_cof__run_tag.csv")
    print("Outputs: 5_LUT_CF__run_tag.csv")
    print("Outputs: 6_flux_out__run_tag.csv")
    print("Outputs: 0_run_input__run_tag.ini")
    
    
    return LUT_cof, LUTCF_u, LUTCF_s

