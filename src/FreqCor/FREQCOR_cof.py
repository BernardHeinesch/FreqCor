# -*- coding: utf-8 -*-

"""Transfer function computation utilities for cut-off frequency estimation."""

import os
import numpy as np
import pandas as pd
import scipy.stats
from scipy.optimize import curve_fit
from . import  as fqfct
from .FREQCOR_functions import fun_Gauss, fun_Lorentz, fun_Lorentz_peltola, fit_Aslan21, test_and_apply_denoising
from . import  as fqplt

def FREQCOR_cof(Icon,Rcon,jtmin_hz,jtmax_hz,varclim,freq, freqn,
                sps, zL, tfmin_hz,tfmax_hz, plot, ws,outputpath, 
                wd, tf_peltola, jrh_class=None, run_tag=None):
    """
    Estimate transfer functions and cut-off frequencies from ideal and real (co)spectra.

    The routine computes the transfer function as the ratio between real and ideal
    (co)spectra and fits parametric forms (Lorentz and Gauss) to derive cut-off
    frequencies and normalization factors.

    Parameters
    ----------
    Icon, Rcon : pandas.Series, numpy.ndarray, or pandas.DataFrame
        Ideal (reference) and real (measured) (co)spectra.
        - If 1D, the routine treats the input as an ensemble-average spectrum.
        - If 2D, columns are treated as individual half-hours.
    jtmin_hz, jtmax_hz : float
        Frequency range (Hz) used for transfer function quality checks.
    varclim : float
        Maximum allowed coefficient of variation for the normalization factor.
    freq : pandas.Series or numpy.ndarray
        Natural frequencies (Hz).
    freqn : pandas.DataFrame or numpy.ndarray
        Normalized frequencies per time step.
    sps : int
        Approach selector (cospectral or spectral).
    zL : pandas.Series or numpy.ndarray
        Stability parameter time series.
    tfmin_hz, tfmax_hz : float
        Frequency range (Hz) used for parametric fitting.
    plot : list[int]
        Plotting options selector.
    ws, wd : pandas.Series
        Wind speed and wind direction time series.
    outputpath : str
        Output directory.
    tf_peltola : int
        Use Peltola et al. form of Lorentzian fit.
    jrh_class : int, optional
        First-variable class index (WD or RH class; used for plotting/output naming).
    run_tag : str, optional
        Run identifier used to generate output filenames.

    Returns
    -------
    xy : list[numpy.ndarray]
        ``[x, y]`` vectors used for the fit (frequency and transfer function samples).
    cofmat : pandas.Series
        Cut-off frequencies and uncertainties with keys:
        ``'cof_L', 'unc_L', 'cof_G', 'unc_G'``.
    fnmat : pandas.Series
        Normalization factors and uncertainties with keys:
        ``'fn_L', 'uncfn_L', 'fn_G', 'uncfn_G'``.
    trafun : pandas.Series, numpy.ndarray, or pandas.DataFrame
        Transfer function values (ratio real/ideal).
    denoising_out : dict
        Denoising results (when active; see Aslan et al., 2021).

    Raises
    ------
    ValueError
        If `run_tag` is missing.
    ValueError
        If plotting of PDF outputs is enabled and `jrh_class` is missing.
    """
    if not run_tag:
        raise ValueError("run_tag is required for output naming")

#%% Extraction of variables
    # find indexes jtmin, jtmax, tfmin, tfmax according to the input value in Hz
    differences = np.abs(np.array(freq) - jtmin_hz)
    jtmin = np.argmin(differences)
    differences = np.abs(np.array(freq) - jtmax_hz)
    jtmax = np.argmin(differences)
    differences = np.abs(np.array(freq) - tfmin_hz)
    tfmin = np.argmin(differences)
    differences = np.abs(np.array(freq) - tfmax_hz)
    tfmax = np.argmin(differences)
    
#%% Variables initialisation   
    cofmat = pd.Series(index=['cof_L', 'unc_L', 'cof_G', 'unc_G'], dtype=float)
    fnmat = pd.Series(index=['fn_L', 'uncfn_L', 'fn_G', 'uncfn_G'], dtype=float)
#%% Transfer function computation
    # Quality check of transfer function
    trafun = Rcon.div(Icon)
    trafun = fqfct.check_var_tf(trafun,jtmin,jtmax,varclim)
    
#%% Fit parameters [cof, Fn]
    p0 = [1.0,1.0]
    bounds = ([0,0.5],[np.inf,1.5])
# %% Denoising (Aslan et al. 2021)
    # For spectra and ensemble averaged data (not individual half-hours)
    # and in case the transfer function has passed the quality check: perform
    # denoising, if necessary, according to Aslan et al. 2021
    denoising_out = {'active': False}
    if sps == 2 and len(Rcon.shape) == 1 and trafun.all():
        denoising_out, x, y = test_and_apply_denoising(Rcon, Icon, freq, trafun, 
                                                    tfmin, tfmax, fit_Aslan21, 
                                                    bounds)

# %% Individual half-hours    
    if len(trafun.shape) > 1:
        ws.index = trafun.columns
        wd.index = trafun.columns
        zL.index = trafun.columns
        
        valid_mask = trafun.all(axis=0)
        invalid_cols = np.where(~valid_mask)[0]
        if len(invalid_cols)>0:
            print(f"{len(invalid_cols)} half-hours had the coefficient of variation of the normalisation factor above its maximum limit. Spectral data is discarded and the cut-off frequency is not computed for these half-hours.")
        
        trafun = trafun.loc[:,valid_mask]
        ws = ws.loc[valid_mask]
        wd = wd.loc[valid_mask]
        zL = zL.loc[valid_mask]
        
        hh_names_list = trafun.columns
        
        cofmat_full = pd.DataFrame(
            np.nan,
            index=hh_names_list,
            columns=['cof_L', 'unc_L', 'cof_G', 'unc_G']
        )
        fnmat_full = pd.DataFrame(
            np.nan,
            index=hh_names_list,
            columns=['fn_L', 'uncfn_L', 'fn_G', 'uncfn_G']
        )
        
        for i, hh_name in enumerate(trafun.columns):
            cofmat_i = pd.Series(index=['cof_L', 'unc_L', 'cof_G', 'unc_G'], dtype=float)
            fnmat_i = pd.Series(index=['fn_L', 'uncfn_L', 'fn_G', 'uncfn_G'], dtype=float)
            
            x=(freq.iloc[tfmin:tfmax+1]).reset_index(drop = True)
            y=(trafun.iloc[tfmin:tfmax+1,i]).reset_index(drop = True)
            
            # Remove nan values: they're not supported by curve_fit
            y[y[:]<0] = np.nan
            valid = ~(np.isnan(x)|np.isnan(y))
            
            if sps==1 and tf_peltola==1:
                try:
                    popt, pcov = curve_fit(fun_Lorentz_peltola,x[valid],y[valid],p0,bounds=bounds) 
                    cof, fn = popt
                    errs = np.sqrt(np.diag(pcov))
                    perr = errs[0]
                    fnerr = errs[1]
                    
                    cofmat_i[0]=cof
                    cofmat_i[1]=perr
                    fnmat_i[0]=fn
                    fnmat_i[1]=fnerr
                    
                except (RuntimeError, ValueError) as e:
                    print(hh_name + ' : ' + str(e))
           
                
            else:   
                try:
                    popt, pcov = curve_fit(fun_Lorentz,x[valid],y[valid],p0=p0,bounds=bounds)
                    cof, fn = popt
                    errs = np.sqrt(np.diag(pcov))
                    perr = errs[0]
                    fnerr = errs[1]
                    
                    cofmat_i[0]=cof
                    cofmat_i[1]=perr 
                    fnmat_i[0]=fn
                    fnmat_i[1]=fnerr
                    
                except (RuntimeError, ValueError) as e:
                    print(hh_name + ' : ' + str(e))
                    
            
            try:
                popt, pcov = curve_fit(fun_Gauss,x[valid],y[valid],p0=p0,bounds=bounds) 
                cof, fn = popt
                errs = np.sqrt(np.diag(pcov))
                perr = errs[0]
                fnerr = errs[1]
                cofmat_i[2]=cof
                cofmat_i[3]=perr
                fnmat_i[2]=fn
                fnmat_i[3]=fnerr
                
                
            except(RuntimeError, ValueError) as e:
                print(hh_name + ' : ' + str(e))
                x = np.zeros(tfmax-tfmin)
                x[:] = np.nan
                y = np.zeros(tfmax-tfmin)
                y[:] = np.nan    
                
            
            cofmat_full.loc[hh_name, :] = cofmat_i
            fnmat_full.loc[hh_name, :] = fnmat_i
            
            if plot[0] == 1:
                fqplt.plot_hh_cosp(freq,freqn,sps,zL,Rcon,Icon,ws,i,plot,outputpath) 
                fqplt.plot_TF(fun_Lorentz_peltola,fun_Gauss,x[valid],y[valid],
                              cofmat_i,plot,outputpath,hh_name=hh_name)    
                     
            
        # Remove bad quality fit results and missing values
        final_valid_mask = (
            (cofmat_full['unc_L'] < 1) &
            (cofmat_full['unc_G'] < 1) &
            cofmat_full.notna().all(axis=1)
            )
        
        cofmat_full = cofmat_full[final_valid_mask]
        fnmat_full = fnmat_full[final_valid_mask]
        ws = ws[final_valid_mask]
        wd = wd[final_valid_mask]
        zL = zL[final_valid_mask]
        
        cofmat_hh = pd.concat([ws,wd,zL,cofmat_full,fnmat_full],axis=1)
        cofmat_hh.index.name='Timestamp'
        cofmat_hh.to_csv(os.path.join(outputpath,hh_name[:6]+'_hh_cof.csv'))
        
        hist_L = np.histogram(cofmat_full.iloc[:,0], bins=20) 
        hist_G = np.histogram(cofmat_full.iloc[:,2], bins=20)
        hist_Ldist = scipy.stats.rv_histogram(hist_L)
        hist_Gdist = scipy.stats.rv_histogram(hist_G)
        # Get the latest frequency bin value to compute the pdf on the whole
        # set and not just up to a given threshold.
        X_bound = hist_L[1][-1]
        
        if plot[2]==1:
            if jrh_class is None:
                raise ValueError("jrh_class is required for output naming")
            file_tag = f"4_PDF_cof__class{int(jrh_class):02d}__{run_tag}"
            fqplt.hist_pdf_cof(cofmat_full, hist_Ldist, hist_Gdist, X_bound, plot, outputpath, file_tag=file_tag)

        # Extract cof: get median and std of the scipy distribution and put  
        cofmat[:] = [hist_Ldist.median(),hist_Ldist.std(),hist_Gdist.median(),
                  hist_Gdist.std()]
        fnmat[:] = np.nan 
        
# %% Ensemble average (co)spectra in LUT
    else: 
        if trafun.all(): 
            if not(denoising_out['active']):
                x=(freq.iloc[tfmin:tfmax+1]).reset_index(drop = True)
                y=(trafun.iloc[tfmin:tfmax+1]).reset_index(drop = True)
            
            # Remove negative TF values
            # Remove nan values (not supported by curve_fit)
            y[y[:]<0] = np.nan
            valid = ~(np.isnan(x)|np.isnan(y))
            
            if sps==1 and tf_peltola==1:
                popt, pcov = curve_fit(fun_Lorentz_peltola,x[valid],y[valid],p0=p0,bounds=bounds)
                cof, fn = popt
                errs = np.sqrt(np.diag(pcov))
                perr = errs[0]
                fnerr = errs[1]
                
            else:
                popt, pcov = curve_fit(fun_Lorentz,x[valid],y[valid],p0=p0,bounds=bounds) 
                cof, fn = popt
                errs = np.sqrt(np.diag(pcov))
                perr = errs[0]
                fnerr = errs[1]
                
            
            cofmat['cof_L'] = cof
            cofmat['unc_L'] = perr
            fnmat['fn_L'] = fn
            fnmat['uncfn_L'] = fnerr
            
            popt, pcov = curve_fit(fun_Gauss,x[valid],y[valid],p0=p0,bounds=bounds) 
            cof, fn = popt
            errs = np.sqrt(np.diag(pcov))
            perr = errs[0]
            fnerr = errs[1]
            
            cofmat['cof_G'] = cof  
            cofmat['unc_G'] = perr
            fnmat['fn_G'] = fn
            fnmat['uncfn_G'] = fnerr
        
        else:
            print('The coefficient of variation of the normalisation factor is above its maximum limit. Spectral data is discarded and the cut-off frequency is not computed.')
            x = np.zeros(tfmax-tfmin)
            x[:] = np.nan
            y = np.zeros(tfmax-tfmin)
            y[:] = np.nan
    
    
    return [x,y], cofmat, fnmat, trafun, denoising_out
