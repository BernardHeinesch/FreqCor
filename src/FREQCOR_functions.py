# -*- coding: utf-8 -*-

"""Utility functions used across the FREQCOR processing pipeline."""
import numpy as np
import pandas as pd
import os
from FREQCOR_Ref_cospectrum_for_plotting import Kaimal_cosp, Massman_cosp
import scipy.stats
from scipy.optimize import curve_fit

def check_var_tf(trafun,jtmin,jtmax,varclim):
    """
    Function that computes the coefficient of variation of the transfer
    function normalisation factor and compares it to the user-defined limit. 
    If the coefficient of variation exceeds its limit value, the transfer 
    function is discarded and the cut-off frequency won't be computed for that 
    specific half-hour or wind class.

    Parameters
    ----------
    trafun : TYPE
        Experimental transfer function computed as the ratio of real to ideal
        (co)spectra.
    jtmin : int
        Index for the lower threshold of the user-defined normalisation range.
    jtmax : int
        Index for the upper threshold of the user-defined normalisation range.
    varclim : float
        User-defined limit for the coefficient of variation of the 
        normalisation factor. If exceeded, the transfer function is discarded.

    Returns
    -------
    trafun : TYPE
        Experimental transfer function after potential filtering performed 
        according to the results of the check of the coefficient of variation
        against the set limit. Is "False" if it didn't pass the check.

    """
    mf=np.nanmean(trafun.iloc[jtmin:jtmax+1],axis=0) 
    varcof= np.nanstd(trafun.iloc[jtmin:jtmax+1],axis=0)/(np.nanmean(trafun.iloc[jtmin:jtmax+1],axis=0))
    
    if isinstance(mf,float): # if work with mean values (mf not dataframe but single float)
        if varcof > varclim:
            trafun.iloc[:] = False #Series as the method .all will be applied on it after
    else: # if half-hour data
        for i in range(0,len(mf)):
            if varcof[i] <= 0 or varcof[i] > varclim:
                trafun.iloc[:,i]=False 
        
    return trafun

def fun_Lorentz(x,a,c):
    """
    Non-linear regression according to Lorentz
    """

    # Lorentzian regression
    # The coefficient (a) is the cut-off frequency (fco)    
    # The coefficient (c) is the normalization factor (Fn)
    return c/(1+(x/a)**2)             


def fun_Gauss(x,a,c):
    """
    Non-linear regression according to Gauss
    """

    # Gaussian regression
    # The coefficient (a) is the cut-off frequency (fco)    
    # The coefficient (c) is the normalization factor (Fn)
    return c*np.exp(-np.log(2)*(x/a)**2)

def fun_Lorentz_peltola(x,a,c):
    """
    Non-linear regression according to Lorentz, applied with the square-root
    according to Peltola et al. 2021
    """

    # Lorentzian regression
    # The coefficient (a) is the cut-off frequency (fco)    
    # The coefficient (c) is the normalization factor (Fn)
    return c*np.sqrt(1/(1+(x/a)**2))             

def av_Kaimal(Icos,freq,freqn,zL,A0,kf0,mu):
    """
    Function to compute the average Kaimal, Massman and ideal cospectra 
    
    Parameters
    ----------
    Icos : pd.DataFrame
        [nfreq,nspec_f]: contains 1/2 hour (co)spectral data of sonic 
        temperature, already filtered for best data quality.
    freq : pd.Series
        [nfreq,]: natural frequency data, in log-spaced bins.
    freqn : pd.DataFrame
        [nfreq,nspec_f]: normalised frequency data, in log-spaced bins, for
        each available 1/2 hour also present in the "Icos" DataFrame.
    
    Returns
    -------
    grouped_data : pd.DataFrame
        [nfreq,3]: DataFrame with newly log-binned frequencies and averaged 
        values for the ideal, Kaimal and Massman cospectra

    """
    # Concatenate all columns of the dataframe to create a series for both the
    # ideal cospectral data and the normalised frequencies
    Icos_n = Icos.T.stack()
    freq_n = freqn.T.stack()
    
    # Create a series of the same size as the two above with the natural
    # frequencies by repeating the same list as many time as needed and 
    # concatenating them
    freq_list = pd.Series(freq.tolist()*freqn.shape[1])
    Icos_n = Icos_n.droplevel(level=1)
    Icos_n = Icos_n.reset_index(drop=True)
    Icos_nf = Icos_n*freq_list
    freq_n = freq_n.droplevel(level=[1])
    freq_n = freq_n.reset_index(drop=True)

    # # Create the full dataframe with normalised frequencies, ideal and Kaimal
    # # cospectra presented in the form cosp/cov*freq
    # # Sort the dataframe by the normalised frequencies increasing order
    # cosp_av = pd.concat([freq_n,Icos_nf,kaimal_tot_f],axis=1)
    cosp_av = pd.concat([freq_n,Icos_nf],axis=1)
    cosp_av_sort = cosp_av.sort_values(by=0)
    cosp_av_sort = cosp_av_sort.reset_index(drop=True)
    
    # Create the new vector of "averaged" normalised frequencies: log-space in
    # as many bins as the original ones (len(freq)) the vector with all the  
    # normalised frequencies.
    # Then add a "0" line so that the first bin goes from 0 to log_bins[1]
    log_bins = np.logspace(np.log10(cosp_av_sort.iloc[0,0]),np.log10(cosp_av_sort.iloc[-1,0]),len(freq))
    log_bins_groups = np.append(np.array(0),log_bins) 
    
    # Use pandas .cut and .groupby to first create the groups according to the
    # pre-defined bins and the normalised frequency column. Each group includes
    # the right-extremum and excludes the left one (x1, x2].
    # Groupby the created groups and take the mean of the values in the bin
    groups = pd.cut(cosp_av_sort.iloc[:,0], log_bins_groups)
    grouped_data = cosp_av_sort.groupby(groups, observed=False).mean()
    
    # Drop the averaged normalised frequency column as it makes no sense
    # Replace the group index by the log_bins values
    # Rename the columns so that they are easily recognised
    grouped_data = grouped_data.drop(0, axis=1) 
    
    # grouped_data.columns = ['Ideal','Kaimal']
    grouped_data.columns = ['Ideal']
    
    grouped_data = grouped_data.reset_index(drop=True)
    kaimal = pd.Series(np.zeros(shape=freq.shape))
    
    for k in range(len(log_bins)):
        kaimal[k] = Kaimal_cosp(freq[k],log_bins[k],zL[0])*freq[k]
    
    grouped_data = pd.concat([grouped_data,kaimal.rename('Kaimal')],axis=1)
    
    
    # Add Massman model cospectrum
    massman = pd.Series(np.zeros(shape=freq.shape))
    if not(np.isnan(A0).all()) and not(np.isnan(kf0).all()) and not(np.isnan(mu).all()):
        for m in range(len(log_bins)):
            massman[m] = Massman_cosp(freq[m],log_bins[m],zL[0],A0,kf0,mu)*freq[m]
    else:
        massman.loc[:]=np.nan
    grouped_data = pd.concat([grouped_data,massman.rename('Massman')],axis=1)
    
    grouped_data.index = log_bins
    
    return(grouped_data)


def ReadLUT(Var, LUT, mxcl, gss, var0=None):
    """
    Read a LUT class index from a look-up table.

    This helper supports two cases:
    - 1D LUT access when `var0 is None`.
    - 2-level LUT access (first by WD/RH, then by wind speed) when `var0` is given.

    Parameters
    ----------
    Var : float
        Query value for the second-level variable (typically wind speed).
    LUT : pandas.DataFrame or dict
        Look-up table.
    mxcl : int
        Number of classes for 1D LUT access (ignored for 2-level access).
    gss : int
        Gas species selector.
    var0 : float, optional
        First-level query value (wind direction for CO2 or relative humidity for H2O).

    Returns
    -------
    VarI : int or float
        Index of the selected wind-speed class (0-based). Returns NaN if inputs
        are not finite.
    clrh : int or float
        Selected first-level class key (WD/RH class). Returns NaN when `var0` is
        None or not finite.
    """
    if var0 is None:
        clrh = np.nan
        if np.isfinite(Var):
            VarI = np.searchsorted(LUT.iloc[:, 1], Var, side='left')
            if VarI == mxcl:
                VarI = mxcl - 1
        else:
            VarI = np.nan
    else:
        if np.isfinite(var0):
            if gss != 2 and var0 > 180:
                var0 = var0 - 180
            key = 'rh' if gss == 2 else 'wd'
            first_level_keys = sorted(LUT.keys())
            first_level_bounds = []
            for k in first_level_keys:
                first_level_bounds.append(LUT[k][key][f'{key}max'].iloc[0])
            clrh_idx = np.searchsorted(first_level_bounds, var0, side='left')
            if clrh_idx == len(first_level_keys):
                clrh_idx = len(first_level_keys) - 1
            selected_key = first_level_keys[clrh_idx]
            clrh = selected_key
            ws_df = LUT[selected_key]['ws']
            mxclws = len(ws_df)
            ws_bounds = ws_df.iloc[:mxclws,1].values
            VarI = np.searchsorted(ws_bounds, Var, side='left')
            if VarI == mxclws:
                VarI = mxclws - 1
        else:
            VarI = np.nan
            clrh = np.nan

    return VarI, clrh


def Simpson(nonweighted):
    """
    Apply Simpson weighting along a 1D array.

    Parameters
    ----------
    nonweighted : numpy.ndarray or pandas.Series
        Input array.

    Returns
    -------
    weighted : numpy.ndarray
        Weighted array computed as ``(f[i-1] + f[i+1] + 4*f[i]) / 6``.
    """

    weightedfw=np.roll(nonweighted,1)
    weightedbk=np.roll(nonweighted,-1)
    weighted=(weightedfw+weightedbk + 4*nonweighted)/6
    return weighted

def TFsonic(freq,ws):
    """Compute the sonic transfer function correction.

    Parameters
    ----------
    freq : numpy.ndarray or pandas.Series
        Natural frequencies (Hz).
    ws : float
        Mean wind speed (m s-1).

    Returns
    -------
    TFs : numpy.ndarray
        Transfer function values.
    """
    sonic_len = np.sqrt(0.11**2 + 0.125**2)# path length (m). Where do the values come from?
    sonic_freq = 20 # acquisition frequency (Hz)
    fp_sonic = freq*abs(sonic_len/ws) # Normalised frequency defined for the path averaging TF definition

    # sonic dynamic time response
    TFFdsonic = 1 / np.sqrt(1 + (2 * np.pi * freq[:] / sonic_freq)**2 )    
    # sonic path-averaging
    TFFwsonic = (2 / (np.pi * fp_sonic)) * (1 + np.exp(-2 * np.pi * fp_sonic) / 2 - 3 * (1 - np.exp(-2 * np.pi * fp_sonic)) / (4 * np.pi * fp_sonic))
        
    TFs= TFFdsonic * np.sqrt(TFFwsonic) # sqrt for the term of path averaging since estimated through spectra
    return TFs


def spectral_outlier_mask_iqr(
    df,
    freq,
    fmin=0.005,
    fmax=0.1,
    k=2.5,
    min_fraction=0.15,
):
    """
    Outlier detection using pointwise IQR in a frequency band.
    """

    freq_iqr = np.asarray(freq)
    band_mask = (freq_iqr >= fmin) & (freq_iqr <= fmax)
    eps = 1e-12
    df_or = df.copy()
    df_pos = df.where((df + eps) > 0)
    df_log = np.log10(df_pos + eps)
    df_band = df_log.loc[band_mask]
    # df_band = df.loc[band_mask]

    arr = df_band.values

    q1 = np.nanpercentile(arr, 25, axis=1, keepdims=True)
    q3 = np.nanpercentile(arr, 75, axis=1, keepdims=True)
    iqr = q3 - q1

    iqr[iqr == 0] = 1e-12

    lower = q1 - k * iqr
    upper = q3 + k * iqr

    outlier_points = (arr < lower) | (arr > upper)
    valid_points = ~np.isnan(arr)

    numer = np.sum(outlier_points & valid_points, axis=0)
    denom = np.sum(valid_points, axis=0)
    fraction_out = np.ones_like(numer, dtype=float)
    np.divide(numer, denom, out=fraction_out, where=denom > 0)

    median = np.nanmedian(arr, axis=1)
    keep_mask = fraction_out < min_fraction
    stats = {
        'median': median,
        'lower': lower,
        'upper': upper,
        'band_mask': band_mask,
    }

    return pd.Series(keep_mask, index=df.columns), stats

def remove_outliers(data, class_id, ws_values, outlier_info=None):
    """
    Detect and remove outliers from a data series using the IQR method.
    Optionally stores information about removed outliers in outlier_info dict.
    Returns filtered data (same type as input), and optionally updates outlier_info.
    """
    if len(data) <= 3:  # Not enough data to reliably detect outliers
        return data
    q1 = np.percentile(data, 25)
    q3 = np.percentile(data, 75)
    iqr = q3 - q1
    lower_bound = q1 - 1.5 * iqr
    upper_bound = q3 + 1.5 * iqr
    mask = (data >= lower_bound) & (data <= upper_bound)
    # Return original data in case all values are filtered
    if np.sum(mask) == 0:
        return data
    ws_range = (np.min(ws_values), np.max(ws_values))
    if outlier_info is not None and np.sum(~mask) > 0:
        removed_values = data[~mask].values
        outlier_info[class_id] = {
            "wind_speed_range": ws_range,
            "removed_outliers": removed_values.tolist()
        }
    return data[mask]

def write_coflut_to_csv(data_dict, gss, outputpath, run_tag=None):
    """Write cut-off frequency LUT to CSV.

    Parameters
    ----------
    data_dict : dict
        LUT structure as returned by `FREQCOR_LUT_cof`.
    gss : int
        Gas species selector.
    outputpath : str
        Output directory.
    run_tag : str, optional
        Run identifier used to generate output filenames.

    Raises
    ------
    ValueError
        If `run_tag` is missing.

    Returns
    -------
    None
        Writes ``4_LUT_cof__{run_tag}.csv``.
    """
    if not run_tag:
        raise ValueError("run_tag is required for output naming")
    cof_file = os.path.join(outputpath, f"4_LUT_cof__{run_tag}.csv")
    if gss == 2:
        concat_keys = ['rh', 'ws']
        var0 = 'rh'
        var0_unit = '%'
    else:
        concat_keys = ['wd', 'ws']
        var0 = 'wd'
        var0_unit = 'deg'

    base_columns = [
        f'{var0}_mean ({var0_unit})',
        f'{var0}_max ({var0_unit})',
        f'unc_{var0}class ({var0_unit})',
        'ws_mean (m s-1)',
        'ws_max (m s-1)',
        'cof_L (Hz)',
        'unc_L_tf (Hz)',
        'cof_G (Hz)',
        'unc_G_tf (Hz)',
        'fn_L (Hz)',
        'uncfn_L (Hz)',
        'fn_G (Hz)',
        'uncfn_G (Hz)',
    ]
    with open(cof_file, "w", newline="") as f:
        first_entry = True
        for key, dfs in data_dict.items():
            merged_df = pd.concat([dfs[concat_keys[0]], dfs[concat_keys[1]]], axis=1)

            columns = list(base_columns)
            if 'ws_n' in dfs:
                n_series = dfs['ws_n']
                n_values = n_series.to_numpy() if isinstance(n_series, pd.Series) else n_series
                merged_df.insert(5, 'N', n_values)
                columns.insert(5, 'N (-)')

            merged_df.columns = columns
            try:
                if first_entry:
                    merged_df.to_csv(f, index=False, lineterminator="\n")
                    first_entry = False
                else:
                    f.write("\n")
                    merged_df.to_csv(f, index=False, header=False, lineterminator="\n")
            except TypeError as e:
                print(e)
                if first_entry:
                    merged_df.to_csv(f, index=False, lineterminator="\n")
                    first_entry = False
                else:
                    f.write("\n")
                    merged_df.to_csv(f, index=False, header=False, lineterminator="\n")

def write_cflut_to_csv(outfile,data_dict,gss):
    """Write correction factor LUT to an open CSV file handle.

    Parameters
    ----------
    outfile : io.TextIOBase
        Open file handle in text write mode.
    data_dict : dict
        LUT structure as returned by `FREQCOR_LUT_CF`.
    gss : int
        Gas species selector.

    Returns
    -------
    None
        Writes CSV blocks to `outfile`.
    """
    first_entry = True
    # Write the different dictionary entries
    var0 ='rh' if gss == 2 else 'wd'
    for key, dfs in data_dict.items():
        # Concatenate DataFrames a and b horizontally
        merged_df = pd.concat([dfs[var0], dfs['ws']], axis=1)
        if 'ws_n' in dfs:
            n_series = dfs['ws_n']
            merged_df['N (-)'] = n_series.to_numpy() if isinstance(n_series, pd.Series) else n_series
        if var0 == 'wd':
            var0_unit = 'deg'
        else:
            var0_unit = '%'
        base_cols=[
            f'{var0}_mean ({var0_unit})',
            f'{var0}_max ({var0_unit})',
            f'unc_{var0}class ({var0_unit})',
            'ws_mean (m s-1)',
            'ws_max (m s-1)',
            'CF_L (-)',
            'unc_L_tf (-)',
            'unc_L_sd (-)',
            'CF_G (-)',
            'unc_G_tf (-)',
            'unc_G_sd (-)',
            'unc_wsCFclass (m s-1)',
        ]
        if 'ws_n' in dfs:
            base_cols.append('N (-)')
        merged_df.columns = base_cols
        # Write to CSV
        try:
            if first_entry:
                merged_df.to_csv(outfile, index=False, lineterminator="\n")  # Write with header
                first_entry = False
            else:
                outfile.write("\n")  # Add blank line before next block
                merged_df.to_csv(outfile, index=False, header=False, lineterminator="\n")  # Append without header
        except TypeError as e:
            print(e)
            if first_entry:
                merged_df.to_csv(outfile, index=False, lineterminator="\n")  # Write with header
                first_entry = False
            else:
                outfile.write("\n")  # Add blank line before next block
                merged_df.to_csv(outfile, index=False, header=False, lineterminator="\n")  # Append without header


#%% Aslan et al. 2021 - DENOISING PROCEDURE
# Define non-linear function to fit
# perform regression
# plot results
# substract noise term from signal to compute noise-free attenuated spectrum


# Define the model to fit: PS_original * b / (1 + (2 * np.pi * freq_binned * a)^2) + freq_binned * d
# If PS_original is passed explicitely as a input, then the curve_fit method
# cannot be directly used because PS_original is not the independent variable 
# nor one of the parameters that have to be optimised. 
# So to still perform the non-linear regression, use of partial function application 
def fit_func(freq_binned, a, b, d, PS_original):
    """
    @author: Aslan et al. 2021
    Translated by A. FaurÃ¨s from Matlab to Python on 21/10/2024
    
    Function to fit power spectra. The function is defined as the sum of an 
    attenuated turbulent signal component and an independent white noise 
    compontent.
    

    Parameters
    ----------
    freq_binned : numpy array
        Natural frequency, binned.
    a : float64
        Low-pass filtering time constant (s).
    b : float64
        Transfer function normalisation factor.
    d : float64
        Intercept of white noise line with the y-axis.
    PS_original : numpy array
        Binned power spectra of sonic temperature (unattenuated, noise-free).
        Is in the form: f*sp(T)/var(T)

    Returns
    -------
    Numpy array
        Attenuated and noise-contaminated power spectra.

    """
    return PS_original * b / (1 + (2 * np.pi * freq_binned * a) ** 2) + freq_binned * d


# Perform fit

def fit_Aslan21(f, ideal_sp, real_sp, par_bounds, ini_par):
    """
    @author: Aslan et al. 2021
    Translated by A. Faurès from Matlab to Python on 21/10/2024 
    
    Function to fit the equation developed by Aslan et al. 2021 to any spectral
    case and dataset.
    
    Parameters
    ----------
    f : numpy array
        Natural frequency, binned.
    ideal_sp : numpy array
        Binned power spectra of sonic temperature (unattenuated, noise-free).
        Is in the form: f*sp(T)/var(T)
    real_sp : numpy array
        Binned power spectra of tracer of interest (attenuated, noisy).
        Is in the form: f*sp(c)/var(c)
    lim : list of float64
        Thresholds for fitting [Hz]. In the form [Low_lim, Up_lim]
    ini_par : list of float64
        Starting points for a, b, d parameters. In the form [a0, b0, d0]
        Reminder: a: filter time constant (s); b: normalisation factor; d: intercept of white noise line with y-axis.
    Returns
    -------
    popt : numpy array
        a, b, d parameters optimised (least-squares).
    PSA_A21 : numpy array
        Binned power spectra modelled with optimised parameters.
    r_squared : float64
        Goodness of fit (calculated from residuals).
    noise : numpy array
        Calculation of the noise term as f*d.

    """
    
    # Input parameters:
    # Version adapted from Aslan et al. 2021 code to correctly set the
    # parameter bounds and not use the fitting thresholds to do so (as it
    # was the case)
    
    # Low_lim = lim[0]  # lower threshold for fitting (Hz)
    # Up_lim = lim[1]       # upper threshold for fitting (Hz)

    # Bounds for a, b, d (from Low_lim to Up_lim)
    # param_bounds = ([Low_lim, Low_lim, -5], [Up_lim, Up_lim, 5])
    
    
    # Perform the non-linear curve fitting, passing ideal_sp via a lambda function
    # See above (befor def function) why the lamba function has to be used
    popt, pcov = curve_fit(lambda f, a, b, d: fit_func(f, a, b, d, ideal_sp), 
                           f, real_sp, p0=ini_par, bounds=par_bounds)
    
    # Extract optimized parameters
    Tc, Fn, Int_point = popt
    # Tc = a : time constant for low-pass filtering
    # Fn = b : normalisation factor for the transfer function
    # Int_point = d : intercept of noise line with y axis; noise contribution
    
    # Print optimized parameters
    # print(f"Optimized parameters: Tc = {Tc}, Fn = {Fn}, Int_point = {Int_point}")
    
    # Calculate PSA_A21 based on the fitted parameters
    PSA_A21 = ideal_sp * Fn / (1 + (2 * np.pi * f * Tc) ** 2) + f * Int_point
    
    
    # Calculate R^2 (goodness of fit)
    residuals = real_sp - fit_func(f, *popt, ideal_sp)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((real_sp - np.mean(real_sp))**2)
    r_squared = 1 - (ss_res / ss_tot)
    
    # Output the results
    print(f"R^2 (Goodness of fit): {r_squared}")
    
    noise = f*Int_point
    
    return popt, pcov, PSA_A21, r_squared, noise

def test_and_apply_denoising(Rcos, Icos, freq, trafun, tfmin, tfmax, fit_Aslan21, bounds):
    """
    Conditionally applies the denoising procedure described by Aslan et al. 2021 to spectral data.

    The function:
      1. Determines a frequency threshold (3 Hz if max(freq) <= 5 Hz, else 5 Hz).
      2. Performs a linear regression on log-log(f * spectrum/variance) above the threshold.
      3. If the slope of the regression is > 0.2, applies the denoising fit using fit_Aslan21.
      4. Returns the denoising results and the (possibly denoised) transfer function.

    Parameters
    ----------
    Rcos : array-like or pd.Series
        Real (measured) (co)spectrum.
    Icos : array-like or pd.Series
        Ideal (reference) (co)spectrum.
    freq : array-like or pd.Series
        Frequency vector (Hz).
    trafun : array-like or pd.Series
        Transfer function (real/ideal cospectra ratio).
    tfmin : int
        Index for the lower bound of the fit frequency range.
    tfmax : int
        Index for the upper bound of the fit frequency range.
    fit_Aslan21 : callable
        Function to perform the Aslan et al. 2021 denoising fit.
    bounds : tuple of lists
        Bounds for the fit parameters.

    Returns
    -------
    denoising_out : dict
        Dictionary with denoising status and results. Contains keys:
            - 'active': bool, True if denoising was performed.
            - 'sp_re_denoised': array, denoised spectrum (if denoising performed).
            - 'x': array, frequency vector used in denoising (if performed).
    y_result : array or None
        The denoised transfer function if denoising was performed, otherwise None.

    Notes
    -----
    - Denoising is only performed if the slope of the log-log regression is > 0.2.
    - If denoising fails, denoising_out['active'] is False and y_result is None.
    """

    denoising_out = {'active': False}
    # Step 1: Threshold selection
    threshold = 3 if freq.iloc[-1] <= 5 else 5
    indices = np.where(freq >= threshold)[0]
    x = freq[indices]
    y = (Rcos * freq)[indices]
    valid = np.logical_not(np.isnan(x) | np.isnan(y))
    x_c = x[valid]
    y_c = y[valid]
    x_l = np.log(x_c)
    y_l = np.log(y_c)

    try:
        m1, b1, r_value1, p_value1, std_err1 = scipy.stats.linregress(x_l, y_l)
    except ValueError as e:
        print(e)
        m1 = np.nan

    y_result = None
    x_fit = None
    if m1 > 0.2:
        print('Performing denoising')
        x_fit = (freq.iloc[tfmin:tfmax+1]).reset_index(drop=True)
        sp_id = (Icos.iloc[tfmin:tfmax+1]).reset_index(drop=True)
        sp_re = (Rcos.iloc[tfmin:tfmax+1]).reset_index(drop=True)
        try:
            bounds_den = (bounds[0] + [-5], bounds[1] + [5])
            p0_den = [0.10, 1, -3]
            popt, pcov, PSA_A21, r_squared, noise = fit_Aslan21(
                np.array(x_fit), np.array(x_fit * sp_id), np.array(x_fit * sp_re), bounds_den, p0_den)
            sp_re_denoised = (x_fit * sp_re - noise) / x_fit
            y_result = sp_re_denoised / sp_id
            denoising_out['active'] = True
            denoising_out['sp_re_denoised'] = sp_re_denoised
            denoising_out['x'] = x_fit
        except Exception as e:
            print('Could not perform denoising:', e)
            y_result = None
            x_fit = None
    return denoising_out, x_fit, y_result
