# -*- coding: utf-8 -*-
"""Plotting utilities used across the FREQCOR processing pipeline."""
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from theor_cosp_Kaimal import Kaimal_cosp


_FREQCOR_PLOT_COLORS = {
    'lorentz': '#009E73',
    'gauss': '#CC79A7',
    'cf_points': '#6E6E6E',
    'tf_original': '#9E9E9E',
    'tf_data': '#0072B2',
}


def _single_plot_fontsizes():
    return {
        'title_fs': 14,
        'label_fs': 12,
        'tick_fs': 11,
        'legend_fs': 11,
    }

def plot_TF_unified(fun1, fun2, x, y, cofmat, denoising, plot, outputpath, norm_range, j_class=None, jrh_class=None, hh_name=None, pause=3, gas_type='co2', all_classes_data=None, file_tag=None):
    """
    Plot experimental and fitted transfer functions.

    The function supports:
    - single-class plots (using `x`, `y`, `cofmat`),
    - multi-class plots (using `all_classes_data`).

    Parameters
    ----------
    fun1, fun2 : callable
        Transfer function models used for fitting (Lorentz and Gauss).
    x, y : array-like
        Frequency (Hz) and transfer function samples for a single class.
    cofmat : array-like
        Fit parameters for a single class.
    denoising : dict or list[dict]
        Denoising metadata (e.g., Aslan et al. 2021), for single- or multi-class mode.
    plot : list[int]
        Plot configuration.
    outputpath : str
        Output directory.
    norm_range : list[float]
        Normalization range bounds in Hz.
    j_class : int, optional
        Wind speed class index.
    jrh_class : int, optional
        First-level class index (WD or RH class).
    hh_name : str, optional
        Half-hour identifier.
    pause : int, optional
        Pause duration (seconds) for interactive use.
    gas_type : str, optional
        Gas label used for titles.
    all_classes_data : dict, optional
        Multi-class container holding lists of `x`, `y`, `cofmat`, `fnmat`, etc.
    file_tag : str, optional
        Output filename tag (required when saving is enabled).

    Returns
    -------
    None
        A figure is shown or saved depending on `plot` options.

    Raises
    ------
    ValueError
        If saving is enabled and `file_tag` is missing.
    """

    def _pos(v, idx):
        return v.iloc[idx] if hasattr(v, 'iloc') else v[idx]

    if all_classes_data is not None:
        # Multi-class mode
        x_list = all_classes_data['x']
        y_list = all_classes_data['y']
        cofmat_list = all_classes_data['cofmat']
        fn_list = all_classes_data['fnmat']
        num_classes = len(x_list)
        
        # Get original transfer function data if available
        original_trafun_data = all_classes_data.get('original_trafun', None)
        original_freq_data = all_classes_data.get('original_freq', None)
        ws_list = all_classes_data.get('ws', [None] * num_classes)
        
        # Determine subplot layout based on number of classes
        if num_classes <= 3:
            nrows, ncols = 1, num_classes
        elif num_classes <= 6:
            nrows, ncols = 2, 3
        else:
            nrows, ncols = 3, 3
        
        # Create figure and subplots
        fig, axes = plt.subplots(nrows, ncols, figsize=(16, 8), sharex=True, sharey=True)
        if nrows == 1 and ncols == 1:
            axes = np.array([axes])  # Make it indexable for single subplot
        axes = axes.flatten()  # Flatten for easy indexing
        
        # Hide unused subplots
        for i in range(num_classes, len(axes)):
            axes[i].set_visible(False)
        
        for i in range(num_classes):
            if i < len(x_list) and i < len(y_list) and i < len(cofmat_list):
                if len(x_list[i]) > 0 and not np.isnan(_pos(cofmat_list[i], 0)) and not np.isnan(_pos(cofmat_list[i], 2)):
                    # Plot original transfer function data if available
                    if original_trafun_data is not None and original_freq_data is not None:
                        if i < len(original_trafun_data) and i < len(original_freq_data):
                            axes[i].plot(
                                original_freq_data[i],
                                original_trafun_data[i],
                                'o',
                                alpha=0.3,
                                color=_FREQCOR_PLOT_COLORS['tf_original'],
                                label='Original Data',
                            )
                    
                    # Plot data used for fitting
                    axes[i].plot(
                        x_list[i],
                        y_list[i],
                        'o',
                        color=_FREQCOR_PLOT_COLORS['tf_data'],
                        label='Fitting Data',
                    )
                    
                    # Generate fitted curves using the original x values
                    yfit1 = fun1(x_list[i], _pos(cofmat_list[i], 0), _pos(fn_list[i], 0))
                    yfit2 = fun2(x_list[i], _pos(cofmat_list[i], 2), _pos(fn_list[i], 2))
                    if denoising[i]['active']:
                        axes[i].plot(
                            x_list[i],
                            yfit1,
                            '-',
                            color=_FREQCOR_PLOT_COLORS['lorentz'],
                            label='Lorentz '
                            f'(fco={_pos(cofmat_list[i], 0):.2f}, Fn={_pos(fn_list[i], 0):.2f})\n A21 denoising',
                        )
                    else:
                        axes[i].plot(
                            x_list[i],
                            yfit1,
                            '-',
                            color=_FREQCOR_PLOT_COLORS['lorentz'],
                            label='Lorentz '
                            f'(fco={_pos(cofmat_list[i], 0):.2f}, Fn={_pos(fn_list[i], 0):.2f})',
                        )
                        
                    axes[i].plot(
                        x_list[i],
                        yfit2,
                        '-',
                        color=_FREQCOR_PLOT_COLORS['gauss'],
                        label='Gauss '
                        f'(fco={_pos(cofmat_list[i], 2):.2f}, Fn={_pos(fn_list[i], 2):.2f})',
                    )
                    axes[i].set_xscale('log')
                    
                    # Add vertical lines for normalization range:
                        # not applicable if normalization is done in the fit.
                    # axes[i].axvline(norm_range[0], color='red', linestyle='--', linewidth=1, label='Norm range')
                    # axes[i].axvline(norm_range[1], color='red', linestyle='--', linewidth=1)

                    # Set y-axis limit to minimum of 10^-8 if needed
                    ymin, ymax = axes[i].get_ylim()
                    if ymin < 1e-8:
                        axes[i].set_ylim(bottom=1e-8)
                    
                    # Add coefficient info
                    # ymin, ymax = axes[i].get_ylim()
                    # xmin, xmax = axes[i].get_xlim()
                    # axes[i].text(xmin*2, ymin*2, f'cof_L={cofmat_list[i][0]:.2f}\ncof_G={cofmat_list[i][2]:.2f}')
                    
                    # Add labels and grid
                    if i % ncols == 0:  # First column
                        axes[i].set_ylabel('TF')
                    if i >= (nrows-1) * ncols:  # Last row
                        axes[i].set_xlabel('Frequency [Hz]')
                    
                    axes[i].grid(alpha=0.5)
                    axes[i].legend()
                    
                    # Add class title with wind speed if available
                    if ws_list[i] is not None:
                        axes[i].set_title(f'Class {i+1}: WS = {ws_list[i]:.2f} m/s')
                    else:
                        axes[i].set_title(f'Class {i+1}')
                else:
                    axes[i].text(0.5, 0.5, 'Insufficient data for fitting', 
                               horizontalalignment='center', verticalalignment='center',
                               transform=axes[i].transAxes)
            else:
                axes[i].text(0.5, 0.5, 'No data available', 
                           horizontalalignment='center', verticalalignment='center',
                           transform=axes[i].transAxes)
        
        # Add overall title
        if gas_type.lower() == 'co2':
            plt.suptitle(f'Experimental and fitted TF for different wind speed classes - CO2 (WD Class {jrh_class})')
        else:  # h2o
            plt.suptitle(f'Experimental and fitted TF for different wind speed classes - H2O (RH Class {jrh_class})')
        
        plt.tight_layout()
        
        if plot and plot[3] == 1:
            if file_tag is None:
                raise ValueError('file_tag is required for saving plots')
            plt.savefig(f'{outputpath}/{file_tag}.png')
            plt.close()
        else:
            plt.show()
            
    else:
        plt.figure()
        fs = _single_plot_fontsizes()
        yfit1 = fun1(x, _pos(cofmat, 0))
        yfit2 = fun2(x, _pos(cofmat, 2))
        plt.plot(x, y, 'o', label='Data')
        plt.plot(x, yfit1, '-', color=_FREQCOR_PLOT_COLORS['lorentz'], label='Fit Lorentz')
        plt.plot(x, yfit2, '-', color=_FREQCOR_PLOT_COLORS['gauss'], label='Fit Gauss')
        # Add vertical lines for normalization range
        plt.axvline(norm_range[0], color='red', linestyle='--', linewidth=1, label='Norm range')
        plt.axvline(norm_range[1], color='red', linestyle='--', linewidth=1)
        plt.xscale('log')
        plt.xlabel('Frequency [Hz]', fontsize=fs['label_fs'])
        plt.ylabel('TF', fontsize=fs['label_fs'])
        plt.grid(alpha=0.5)
        plt.tick_params(axis='both', which='both', labelsize=fs['tick_fs'])
        plt.legend(fontsize=fs['legend_fs'])
        
        if hh_name is None:
            plt.title('Experimental and fitted TF. Wind class: ' + str(j_class), fontsize=fs['title_fs'])
        else:
            plt.title(hh_name + ' TF fit', fontsize=fs['title_fs'])
            
        if plot and plot[3] == 1:
            if file_tag is None:
                raise ValueError('file_tag is required for saving plots')
            plt.savefig(outputpath + '/' + file_tag + '.png')
            plt.close()
        else:
            plt.show()



def plot_cosp_unified(freq, Rcosv, Icosv, ws, denoising, zL=None, sps=None, j_class=None, jrh_class=None, plot=None, outputpath=None, gas_type='co2', all_classes_data=None, file_tag=None):
    """
    Plot mean (co)spectra for one class or for all classes.

    Parameters
    ----------
    freq : array-like
        Natural frequencies (Hz).
    Rcosv, Icosv : array-like
        Real and ideal (co)spectra for a single class.
    ws : float
        Mean wind speed for current class.
    denoising : dict or list[dict]
        Denoising metadata (e.g., Aslan et al. 2021), for single- or multi-class mode.
    zL : float, optional
        Stability parameter. The default is None.
    sps : int, optional
        Sampling frequency. The default is None.
    j_class : int, optional
        Wind speed class number. The default is None.
    jrh_class : int, optional
        First-level class index (WD or RH class). The default is None.
    plot : list, optional
        Plot configuration. The default is None.
    outputpath : str, optional
        Output directory. The default is None.
    gas_type : str, optional
        Gas label used for titles. The default is 'co2'.
    all_classes_data : dict, optional
        Multi-class container holding lists of `x`, `y`, `cofmat`, `fnmat`, etc. The default is None.
    file_tag : str, optional
        Output filename tag (required when saving is enabled). The default is None.

    Returns
    -------
    None

    Raises
    ------
    ValueError
        If saving is enabled and `file_tag` is missing.
    """

    
    if all_classes_data is not None:
        # Plot all classes in a single figure with multiple subplots
        num_classes = all_classes_data['num_classes']
        freq = all_classes_data['freq']
        Rcosv_list = all_classes_data['Rcosv']
        Icosv_list = all_classes_data['Icosv']
        ws_list = all_classes_data['ws']
        zL = all_classes_data['zL']
        
        # Determine subplot layout based on number of classes
        if num_classes <= 3:
            nrows, ncols = 1, num_classes
        elif num_classes <= 6:
            nrows, ncols = 2, 3
        else:
            nrows, ncols = 3, 3
            
        fig, axes = plt.subplots(nrows, ncols, figsize=(16, 8), sharex=True, sharey=True)
        if nrows == 1 and ncols == 1:
            axes = np.array([axes])  # Make it indexable for single subplot
        axes = axes.flatten()  # Flatten for easy indexing
        
        # Hide unused subplots
        for i in range(num_classes, len(axes)):
            axes[i].set_visible(False)
        
        for i in range(num_classes):
            # Check if we have valid data for this class
            if i < len(Rcosv_list) and i < len(Icosv_list) and i < len(ws_list):
                Rcosv = Rcosv_list[i]
                Icosv = Icosv_list[i]
                ws = ws_list[i]
                
                # Check if data is valid
                if len(Rcosv) > 0 and len(Icosv) > 0:
                    # Plot data for each class
                    axes[i].loglog(freq, Icosv, '-', color='blue', label='Ideal')
                    axes[i].loglog(freq, Rcosv, 'o', markersize=2, color='red', label='Real')
                    if denoising[i]['active']:
                        axes[i].loglog(denoising[i]['x'], denoising[i]['sp_re_denoised'], '-.', color='red', label='Real denoised')
                    # Set y-axis limit to minimum of 10^-8 if needed
                    ymin, ymax = axes[i].get_ylim()
                    if ymin < 1e-8:
                        axes[i].set_ylim(bottom=1e-8)
                    if ymax > 1e4:
                        axes[i].set_ylim(top=1e4)
                    
                    # Add labels and grid
                    if i % ncols == 0:  # First column
                        axes[i].set_ylabel('(co)sp/(co)var')
                    if i >= (nrows-1) * ncols:  # Last row
                        axes[i].set_xlabel('Frequency [Hz]')
                    
                    axes[i].grid(alpha=0.5)
                    axes[i].legend()
                    
                    # Add class title
                    # Add class title with wind speed if available
                    if ws_list[i] is not None:
                        axes[i].set_title(f'Class {i+1}: WS = {ws_list[i]:.2f} m/s')
                    else:
                        axes[i].set_title(f'Class {i+1}')
                else:
                    axes[i].text(0.5, 0.5, 'Insufficient data', 
                               horizontalalignment='center', verticalalignment='center',
                               transform=axes[i].transAxes)
            else:
                axes[i].text(0.5, 0.5, 'No data available', 
                           horizontalalignment='center', verticalalignment='center',
                           transform=axes[i].transAxes)
        
        # Add overall title
        if gas_type.lower() == 'co2':
            plt.suptitle(f'Mean (co)spectra for different wind speed classes - CO2 (WD Class {jrh_class})')
        else:  # h2o
            plt.suptitle(f'Mean (co)spectra for different wind speed classes - H2O (RH Class {jrh_class})')
        
        plt.tight_layout()
        
        if plot is not None and plot[3] == 1 and outputpath is not None:
            if file_tag is None:
                raise ValueError('file_tag is required for saving plots')
            plt.savefig(f'{outputpath}/{file_tag}.png')
            plt.close()
        else:
            plt.show()
    else:
        # Single class plotting (legacy mode)
        plt.figure(figsize=(10, 6))
        fs = _single_plot_fontsizes()
        plt.loglog(freq, Icosv, '-', color='blue', label='Ideal')
        plt.loglog(freq, Rcosv, 'o', markersize=2, color='red', label='Real')
        plt.xlabel('Frequency [Hz]', fontsize=fs['label_fs'])
        plt.ylabel('(co)sp/(co)var', fontsize=fs['label_fs'])
        plt.grid(alpha=0.5)
        plt.tick_params(axis='both', which='both', labelsize=fs['tick_fs'])
        plt.legend(fontsize=fs['legend_fs'])
        
        if gas_type.lower() == 'co2':
            plt.title(f'Mean (co)spectra - Class {j_class}: WS = {ws:.2f} m/s', fontsize=fs['title_fs'])
        else:  # h2o
            plt.title(f'Mean (co)spectra - RH Class {jrh_class}, WS Class {j_class}: WS = {ws:.2f} m/s', fontsize=fs['title_fs'])
        
        if plot is not None and plot[3] == 1 and outputpath is not None:
            if file_tag is None:
                raise ValueError('file_tag is required for saving plots')
            plt.savefig(f'{outputpath}/{file_tag}.png')
            plt.close()
        else:
            plt.show()
    

        
def plot_hh_cosp(freq,freqn,sps,zL,cosR,cosI,ws,i,plot,outputpath,pause=3, file_tag=None):
    """
    Plot an individual half-hour (co)spectrum.

    Parameters
    ----------
    freq : array-like
        Natural frequencies (Hz).
    freqn : pandas.DataFrame
        Normalized frequencies per half-hour.
    sps : int
        Approach selector (cospectral or spectral).
    zL : pandas.Series
        Stability parameter time series.
    cosR, cosI : pandas.DataFrame
        Real and ideal (co)spectral matrices.
    ws : pandas.Series
        Wind speed time series.
    i : int
        Column index of the half-hour to plot.
    plot : list[int]
        Plot configuration.
    outputpath : str
        Output directory.
    pause : int, optional
        Pause duration for interactive use.
    file_tag : str, optional
        Output filename tag (required when saving is enabled).

    Returns
    -------
    None

    Raises
    ------
    ValueError
        If saving is enabled and `file_tag` is missing.
    """
    fs = _single_plot_fontsizes()
    f, axs = plt.subplots(1,2,figsize=(16,8))
    
    # (Co)sp/(co)var as a function of normalised frequency
    axs[0].plot(freqn.iloc[:,i],cosR.iloc[:,i],'.',label='Real')
    axs[0].plot(freqn.iloc[:,i],cosI.iloc[:,i],'.',label='Ideal')
    if sps == 1: # If cospectra plot Kaimal as well, if spectra don't
        kaimal = np.zeros(len(freq))
        for x in range(0,len(freq)):
            kaimal[x] = Kaimal_cosp(freq[x],freqn.iloc[x,i],zL[i])    
        axs[0].plot(freqn.iloc[:,i],kaimal, label = 'Kaimal')
        axs[0].set_ylabel('Cospectrum/covar')
    else: axs[0].set_ylabel('Spectrum/var')

    axs[0].set_yscale('log')
    axs[0].set_xscale('log')
    axs[0].set_xlabel('Frequency norm. [-]', fontsize=fs['label_fs'])
    axs[0].set_title(cosI.iloc[:,i].name, fontsize=fs['title_fs'])
    axs[0].tick_params(axis='both', which='both', labelsize=fs['tick_fs'])
    axs[0].legend(fontsize=fs['legend_fs'])
    
    ymin, ymax = axs[0].get_ylim()
    xmin, xmax = axs[0].get_xlim()
    axs[0].text(xmin*2,ymin*2,'u='+str("%.2f" % ws[i]) + ' ms-1')
    
    # (Co)sp/(co)var*f_nat as a function of normalised frequency
    axs[1].plot(freqn.iloc[:,i],cosR.iloc[:,i]*freq.iloc[:],'.',label='Real')
    axs[1].plot(freqn.iloc[:,i],cosI.iloc[:,i]*freq.iloc[:],'.',label='Ideal')
    if sps == 1: # If cospectra plot Kaimal as well, if spectra don't
        kaimal = np.zeros(len(freq))
        for x in range(0,len(freq)):
            kaimal[x] = Kaimal_cosp(freq[x],freqn.iloc[x,i],zL[i])    
        axs[1].plot(freqn.iloc[:,i],kaimal*freq.iloc[:], label = 'Kaimal')
        axs[1].set_ylabel('Cospectrum/covar*f_nat')
    else: axs[1].set_ylabel('Spectrum/var*f_nat')
    
    axs[1].set_yscale('log')
    axs[1].set_xscale('log')
    axs[1].set_xlabel('Frequency norm. [-]', fontsize=fs['label_fs'])
    axs[1].set_title(cosI.iloc[:,i].name, fontsize=fs['title_fs'])
    axs[1].tick_params(axis='both', which='both', labelsize=fs['tick_fs'])
    axs[1].legend(fontsize=fs['legend_fs'])
    
    ymin, ymax = axs[1].get_ylim()
    xmin, xmax = axs[1].get_xlim()
    axs[1].text(xmin*2,ymin*2,'u='+str("%.2f" % ws[i]) + ' ms-1')
    
    plt.tight_layout()
    plt.show()
    
    # If requested: save plot to figure and close it in python
    hh_name = cosI.iloc[:,i].name
    if plot[3]==1:
        if file_tag is None:
            raise ValueError('file_tag is required for saving plots')
        plt.savefig(outputpath + '/' + file_tag + '.png')
        plt.close()


def plot_av_ideals(av_cosp_k,plot,outputpath,stability, file_tag=None):
    """
    Function to plot averaged, ideal and model cospectra.    

    Parameters
    ----------
    av_cosp_k : pd.DataFrame
        [nspec,3]: DataFrame with log-binned frequencies and corresponding 
        averaged values for the ideal, Kaimal and Massman cospectra.
    plot : list
        User inputs on plotting preferences (choice of plots, saving options)
    outputpath : str
        Name (date and hour) of the plotted half-hour
    stability : list
        List of two strings identifying the stability regime of the data 
        (unstable or stable).

    Returns
    -------
    None.

    """
    plt.figure(figsize=(16,8))
    fs = _single_plot_fontsizes()
    plt.plot(av_cosp_k['Ideal'], '.',label = 'Ideal')       
    plt.plot(av_cosp_k['Kaimal'],label = 'Kaimal')
    plt.plot(av_cosp_k['Massman'],label = 'Massman')
    plt.yscale('log')
    plt.xscale('log') 
    plt.ylim(1e-5,1)
    plt.xlim(1e-4,1e3)
    plt.grid(alpha=0.5)
    plt.tick_params(axis='both', which='both', labelsize=fs['tick_fs'])
    plt.legend(fontsize=fs['legend_fs'])
    plt.title('Ideal and model cospectra. '+stability+'able conditions', fontsize=fs['title_fs'])
    plt.ylabel('Cosp/covar*fnat', fontsize=fs['label_fs'])
    plt.xlabel('Freq norm [-]', fontsize=fs['label_fs'])
    if plot[3]==1:
        if file_tag is None:
            raise ValueError('file_tag is required for saving plots')
        plt.savefig(outputpath + '/' + file_tag + '.png')
        # plt.savefig(outputpath + '/av_kaimal'+'_' + stability + '.png')
        plt.close()
    

    
def hist_pdf_cof(cofmat_full, hist_Ldist, hist_Gdist, X_bound, plot, outputpath, file_tag=None):
    """
    Plots the histograms of the cut-off frequencies found with two different 
    fitting methods (Lorentz and Gauss functions) and their corresponding 
    probability density function and cumulative distribution function.

    Parameters
    ----------
    cofmat_full : pd.DataFrame
        Cut-off frequencies for all the half-hours, for both fitting methods.
    hist_Ldist : scipy.stats.rv_histogram
        Probability distribution of cut-off frequencies resulting from fitting
        Lorentz transfer function.
    hist_Gdist : scipy.stats.rv_histogram
        Probability distribution of cut-off frequencies resulting from fitting
        Gauss transfer function.
    X_bound : float
        Upper bound for the x-axis used to compute PDF/CDF.
    plot : list[int]
        Plot configuration.
    outputpath : str
        Output directory.
    file_tag : str, optional
        Output filename tag (required when saving is enabled).
        
    Returns
    -------
    None

    Raises
    ------
    ValueError
        If saving is enabled and `file_tag` is missing.

    """

    fs = _single_plot_fontsizes()
    f, axs = plt.subplots(2, 2,figsize=(16,8))
    axs[0,0].hist(cofmat_full.iloc[:,0], color = 'blue', edgecolor = 'black',
             bins = 100)
    axs[0,0].set_title('Histogram of 1/2 h cof (Lorentz fit)', fontsize=fs['title_fs'])
    axs[1,0].hist(cofmat_full.iloc[:,2], color = 'green', edgecolor = 'black',
             bins = 100)
    axs[1,0].set_title('Histogram of 1/2 h cof (Gauss fit)', fontsize=fs['title_fs'])

    axs[0,0].set_ylabel('nÂ° of 1/2 hours', fontsize=fs['label_fs'])
    axs[1,0].set_ylabel('nÂ° of 1/2 hours', fontsize=fs['label_fs'])
    axs[1,0].set_xlabel('Cut-off Frequency [Hz]', fontsize=fs['label_fs'])
    
    axs[1, 0].sharex(axs[0, 0])
    axs[1, 0].sharey(axs[0, 0])
    # Distribution fct from scipy
    X = np.linspace(0, np.round(X_bound), 100)
    
    axs[0,1].hist(cofmat_full.iloc[:,0], density=True, bins=100)
    axs[0,1].plot(X, hist_Ldist.pdf(X), label='PDF')
    axs[0,1].plot(X, hist_Ldist.cdf(X), label='CDF')
    axs[0,1].set_ylabel('Probability', fontsize=fs['label_fs'])
    axs[0,1].tick_params(axis='both', which='both', labelsize=fs['tick_fs'])
    axs[0,1].legend(fontsize=fs['legend_fs'])
    axs[0,1].set_title("PDF and CDF of cof (Lorentz)", fontsize=fs['title_fs'])
    
    axs[1,1].hist(cofmat_full.iloc[:,2], density=True, bins=100)
    axs[1,1].plot(X, hist_Gdist.pdf(X), label='PDF')
    axs[1,1].plot(X, hist_Gdist.cdf(X), label='CDF')
    axs[1,1].set_ylabel('Probability', fontsize=fs['label_fs'])
    axs[1,1].set_xlabel('Cut-off Frequency [Hz]', fontsize=fs['label_fs'])
    axs[1,1].tick_params(axis='both', which='both', labelsize=fs['tick_fs'])
    axs[1,1].legend(fontsize=fs['legend_fs'])
    axs[1,1].set_title("PDF and CDF of cof (Gauss)", fontsize=fs['title_fs'])
    
    axs[1, 1].sharex(axs[0, 1])
    axs[1, 1].sharey(axs[0, 1])
    
    plt.tight_layout()
    if plot[3] == 1:
        if file_tag is None:
            raise ValueError('file_tag is required for saving plots')
        plt.savefig(outputpath + '/' + file_tag + '.png')
        plt.close()
    
   # For the errors on the cof
    # f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, sharey=True)
    # ax1.hist(cofmat_full[:,1], color = 'blue', edgecolor = 'black',
    #          bins = 100)
    # ax1.set_title('Histogram for 1/2 h cof (Lorentz fit)')
    # ax2.hist(cofmat_full[:,3], color = 'green', edgecolor = 'black',
    #          bins = 100)
    # ax2.set_title('Histogram for 1/2 h cof (Gauss fit)')


def plot_u_cf(all_data_cf, gss, plot, outputpath, stability, file_tag=None):
    """
    Plot correction factors as a function of wind speed class.

    Parameters
    ----------
    all_data_cf : dict
        Per-class CF and wind speed data used for plotting.
    gss : int
        Gas species selector.
    stability : str
        Stability regime label.
    plot : list[int]
        Plot configuration.
    outputpath : str
        Output directory.
    file_tag : str, optional
        Output filename tag (required when saving is enabled).

    Returns
    -------
    None

    Raises
    ------
    ValueError
        If saving is enabled and `file_tag` is missing.

    """
    fs = _single_plot_fontsizes()
    title_fs = fs['title_fs'] + 3
    label_fs = fs['label_fs'] + 4
    tick_fs = fs['tick_fs'] + 3
    legend_fs = fs['legend_fs'] + 3
    suptitle_fs = fs['title_fs'] + 4
    
    num_classes = len(all_data_cf)
    mainvar = 'rh' if gss == 2 else 'wd'
    # Determine subplot layout based on number of classes
    if num_classes <= 3:
        nrows, ncols = 1, num_classes
    elif num_classes <= 6:
        nrows, ncols = 2, 3
    elif num_classes <=9:
        nrows, ncols = 3, 3
    else:
        nrows, ncols= 4, 3
        
    fig, axes = plt.subplots(nrows, ncols, figsize=(16, 8), sharex=True, sharey=True)
    if nrows == 1 and ncols == 1:
        axes = np.array([axes])  # Make it indexable for single subplot
    axes = axes.flatten()  # Flatten for easy indexing
    
    # Hide unused subplots
    for i in range(num_classes, len(axes)):
        axes[i].set_visible(False)
    
    for i in range(num_classes):
        ws = all_data_cf[i+1]['ws']
        cf_l = all_data_cf[i+1]['cf_l']
        cf_g = all_data_cf[i+1]['cf_g']
        av_ws = all_data_cf[i+1]['mean_ws']
        av_cf_l = all_data_cf[i+1]['mean_cf_l']
        av_cf_g = all_data_cf[i+1]['mean_cf_g']
        av_rh = all_data_cf[i+1][mainvar]
        # Check if data is valid
        if len(ws) > 0 and len(av_ws) > 0:
            try:
                # Plot data for each class
                axes[i].plot(
                    ws,
                    cf_l,
                    '.',
                    color=_FREQCOR_PLOT_COLORS['cf_points'],
                    label='_nolegend_',
                )
                axes[i].plot(
                    ws,
                    cf_g,
                    '.',
                    color=_FREQCOR_PLOT_COLORS['cf_points'],
                    label='_nolegend_',
                )
                axes[i].plot(
                    av_ws,
                    av_cf_l,
                    'o',
                    color=_FREQCOR_PLOT_COLORS['lorentz'],
                    markersize=10,
                    label='_nolegend_',
                )
                axes[i].plot(
                    av_ws,
                    av_cf_g,
                    'o',
                    color=_FREQCOR_PLOT_COLORS['gauss'],
                    markersize=10,
                    label='_nolegend_',
                )
                
                # Add labels and grid
                if i % ncols == 0:  # First column
                    axes[i].set_ylabel('Correction Factor [-]', fontsize=label_fs)
                if i >= (nrows-1) * ncols:  # Last row
                    axes[i].set_xlabel('Wind speed [ms-1]', fontsize=label_fs)
                
                axes[i].grid(alpha=0.5)
                axes[i].tick_params(axis='both', which='both', labelsize=tick_fs)
                if gss == 2:
                    axes[i].set_title(f'RH Class {i+1}: rh_mean = {av_rh:.2f} %', fontsize=title_fs)
                else:
                    axes[i].set_title(f'WD Class {i+1}: wd_mean = {av_rh:.2f} Â°', fontsize=title_fs)
            except: 
                print('Error in plotting CF vs ws data')    
                continue
            
    
    fig.suptitle('CF vs WS. ' + stability + 'able conditions.', fontsize=suptitle_fs, y=0.99)

    if num_classes > 0:
        legend_elements = [
            Line2D(
                [0],
                [0],
                marker='.',
                linestyle='None',
                color=_FREQCOR_PLOT_COLORS['cf_points'],
                label='Individual CF (Lorentz & Gauss)',
                markersize=10,
            ),
            Line2D(
                [0],
                [0],
                marker='o',
                linestyle='None',
                color=_FREQCOR_PLOT_COLORS['lorentz'],
                label='LUT CF Lorentz (used for flux computation)',
                markersize=8,
            ),
            Line2D(
                [0],
                [0],
                marker='o',
                linestyle='None',
                color=_FREQCOR_PLOT_COLORS['gauss'],
                label='LUT CF Gauss (used for flux computation)',
                markersize=8,
            ),
        ]
        axes[0].legend(
            handles=legend_elements,
            fontsize=legend_fs,
            loc='upper right',
            frameon=True,
        )

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    
    if plot is not None and plot[3] == 1 and outputpath is not None:
        if file_tag is None:
            raise ValueError('file_tag is required for saving plots')
        plt.savefig(f'{outputpath}/{file_tag}.png')
        plt.close()
    else:
        plt.show()
   

def plot_all_individual_cosp(freq, Rcos, Icos, plot, outputpath, gas_type='co2', file_tag=None):
    """
    Plot all individual half-hour (co)spectra on a single figure.
    
    Parameters
    ----------
    freq : numpy array
        Frequency array with shape (n_freq,)
    Rcos : numpy array
        Real part of the cospectrum with shape (n_spectra, n_freq)
    Icos : numpy array
        Imaginary part of the cospectrum with shape (n_spectra, n_freq)
    plot : list
        Plot configuration
    outputpath : str
        Path to save the plot
    gas_type : str, optional
        Gas type ('co2' or 'h2o'). The default is 'co2'.
        
    Returns
    -------
    None.
    """
    plt.figure(figsize=(10, 6))
    fs = _single_plot_fontsizes()
    
    # Define better colors for visibility
    real_color = '#404040'  # Darker grey
    ideal_color = '#1E88E5'  # Brighter blue
    
    # Plot each individual spectrum - transpose the arrays if needed
    if Rcos.shape[1] == len(freq):
        # Rcos shape is (n_spectra, n_freq)
        plt.plot(freq, Rcos.T, '-', color=real_color, alpha=0.3)
        plt.plot(freq, Icos.T, '-', color=ideal_color, alpha=0.3)
    else:
        # Rcos shape might be (n_freq, n_spectra)
        plt.plot(freq, Rcos, '-', color=real_color, alpha=0.3)
        plt.plot(freq, Icos, '-', color=ideal_color, alpha=0.3)

    # Set scales and labels
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Frequency [Hz]', fontsize=fs['label_fs'])
    plt.ylabel('(co)sp/(co)var', fontsize=fs['label_fs'])
    plt.title(f'All individual half-hour {gas_type.upper()} (co)spectra', fontsize=fs['title_fs'])
    plt.grid(alpha=0.3)
    plt.tick_params(axis='both', which='both', labelsize=fs['tick_fs'])
    
    # Add legend with the updated labels and colors
    plt.plot([], [], '-', color=ideal_color, label='Ideal')
    plt.plot([], [], '-', color=real_color, label='Real')
    plt.legend(fontsize=fs['legend_fs'])
    
    # Save or show the plot
    if plot[3] == 1:
        if file_tag is None:
            raise ValueError('file_tag is required for saving plots')
        plt.savefig(f'{outputpath}/{file_tag}.png')
        plt.close()
    else:
        plt.show()
        
def plot_stepwise_filtering(
    *,
    mask_history,
    step_names,
    fluxes_df,
    flux_col,
    cospectra_df,
    freq,
    outputpath,
    plot,
    init_valid_mask=None,
    step_plot_info=None,
    filename="filtering_steps.png",
):
    """
    Plot cumulative step-wise filtering of fluxes and (co)spectra.

    Each row shows one filtering step:
      - blue: data retained up to this step
      - red: data removed at this step only

    The last row shows the final retained data after all filters and reports
    availability relative to the initially valid dataset.

    Parameters
    ----------
    mask_history : list of pandas.Series
        Boolean masks (indexed by time) defining each filtering step.
        True = retained. All masks must align with fluxes_df.index and
        cospectra_df.columns.

    step_names : list of str
        Names of the filtering steps (same length as mask_history).

    fluxes_df : pandas.DataFrame
        Flux time series (indexed by time).

    flux_col : str
        Column of fluxes_df to plot.

    cospectra_df : pandas.DataFrame
        (Co)spectra with frequency as rows and time as columns.

    freq : array-like
        Frequency vector for the (co)spectra.

    outputpath : str
        Directory where the figure is saved.
    
    plot : list of int
        Plot configuration .
        
    init_valid_mask : pandas.Series, optional
        Boolean mask of initially available data.
        Defaults to columns of cospectra_df that are not fully NaN.

    step_plot_info : dict, optional
        Optional step-specific plotting info (e.g. flux thresholds,
        absolute (co)spectral limits, IQR statistics).

    filename : str
        Output figure name.

    Notes
    -----
    Filtering is cumulative: each step operates only on data that passed all
    previous steps. Percentages shown per step refer to data entering that step.
    """

    fs = _single_plot_fontsizes()
    axis_label_fs = fs['label_fs'] + 2

    flux_col_labels = {
        'H': 'H flux [W m-2]',
        'Fc': 'CO2 flux [umol m-2 s-1]',
        'LE': 'Latent heat flux (LE) [W m-2]',
    }
    flux_ylabel = flux_col_labels.get(flux_col, f"{flux_col} flux")

    if init_valid_mask is None:
        init_valid_mask = ~cospectra_df.isna().all()

    if step_plot_info is None:
        step_plot_info = {}
    
    time_index = pd.Series(cospectra_df.columns)
    time_index = pd.to_datetime(time_index, format="%Y%m%d-%H%M")
    time_index.index=cospectra_df.columns
    n_steps = len(mask_history)

    fig, axes = plt.subplots(
        nrows=n_steps + 1,
        ncols=2,
        figsize=(16, 8),
        sharex='col',
        sharey=False
    )

    if n_steps == 1:
        axes = axes.reshape(1, 2)

    prev_keep = init_valid_mask.copy()

    # =========================
    # Stepwise filtering plots
    # =========================
    for i in range(n_steps):
        ax_flux, ax_spec = axes[i]
        step_name = step_names[i]

        removed_now = prev_keep & ~mask_history[i]

        n_prev = prev_keep.sum()
        pct_removed = 100 * removed_now.sum() / n_prev if n_prev > 0 else 0

        # ---- Fluxes
        ax_flux.scatter(
            time_index[prev_keep],
            fluxes_df.loc[prev_keep, flux_col],
            color='steelblue', alpha=0.25, s=20
        )

        if removed_now.any():
            ax_flux.scatter(
                time_index[removed_now],
                fluxes_df.loc[removed_now, flux_col],
                color='crimson', alpha=0.6, s=5
            )

        ax_flux.set_title(f"Step {i+1}: {step_name}")
        ax_flux.grid(True, ls=':')

        # Step-specific flux info
        info = step_plot_info.get(step_name, {})
        if "flux_limits" in info:
            for y in info["flux_limits"]:
                ax_flux.axhline(y, color='k', ls='--', lw=1)
                ax_flux.axhline(-y, color='k', ls='--', lw=1)

        # ---- (Co)spectra
        ax_spec.loglog(
            freq,
            cospectra_df.loc[:, prev_keep].values,
            color='steelblue', alpha=0.25
        )

        if removed_now.any():
            ax_spec.loglog(
                freq,
                cospectra_df.loc[:, removed_now].values,
                color='crimson', alpha=0.6, lw=1
            )

        # Step-specific cospectra info
        if "cospec_abs_limit" in info:
            ax_spec.axhline(info["cospec_abs_limit"], color='k', ls='--', lw=1)

        if "iqr_stats" in info:
            s = info["iqr_stats"]
            ax_spec.plot(freq[s["band_mask"]], 10**s["lower"], 'k--', lw=1)
            ax_spec.plot(freq[s["band_mask"]], 10**s["upper"], 'k--', lw=1)
            ax_spec.plot(freq[s["band_mask"]], 10**s["median"], 'k', lw=1)
            ax_spec.axvspan(
                freq[s["band_mask"]].iloc[0],
                freq[s["band_mask"]].iloc[-1],
                color='grey', alpha=0.15
            )

        ax_spec.set_title(f"Step {i+1}: {step_name}")
        ax_spec.grid(True, which='both', ls=':')

        ax_spec.text(
            0.98, 0.95,
            f"Removed: {pct_removed:.1f}%",
            transform=ax_spec.transAxes,
            ha='right', va='top',
            fontsize=9,
            bbox=dict(facecolor='white', alpha=0.7, edgecolor='none')
        )

        prev_keep &= mask_history[i]

    # =========================
    # Final retained data
    # =========================
    final_keep = prev_keep
    n_init = init_valid_mask.sum()
    n_final = final_keep.sum()
    pct_available = 100 * n_final / n_init if n_init > 0 else np.nan

    ax_flux_f, ax_spec_f = axes[-1]

    ax_flux_f.scatter(
        time_index[final_keep],
        fluxes_df.loc[final_keep, flux_col],
        color='steelblue', alpha=0.5, s=10
    )
    ax_flux_f.set_title("Final retained fluxes")
    ax_flux_f.grid(True, ls=':')

    ax_spec_f.loglog(
        freq,
        cospectra_df.loc[:, final_keep].values,
        color='steelblue', alpha=0.5
    )

    ax_spec_f.text(
        0.98, 0.95,
        f"Available: {n_final}/{n_init}\n({pct_available:.1f}%)",
        transform=ax_spec_f.transAxes,
        ha='right', va='top',
        fontsize=9,
        bbox=dict(facecolor='white', alpha=0.7, edgecolor='none')
    )

    ax_spec_f.set_title("Final retained (co)spectra")
    ax_spec_f.grid(True, which='both', ls=':')

    axes[-1, 0].set_xlabel("Time", fontsize=axis_label_fs)
    axes[-1, 1].set_xlabel("Frequency [Hz]", fontsize=axis_label_fs)
    axes[n_steps // 2, 0].set_ylabel(flux_ylabel, fontsize=axis_label_fs)
    axes[n_steps // 2, 1].set_ylabel("Spectral density", fontsize=axis_label_fs)

    # ---- Global legend
    legend_elements = [
        Line2D([0], [0], marker='o', color='w',
               label='Retained data', markerfacecolor='steelblue', markersize=8),
        Line2D([0], [0], marker='o', color='w',
               label='Removed at this step', markerfacecolor='crimson', markersize=8),
    ]

    fig.legend(handles=legend_elements, loc='lower center', ncol=2, frameon=False)

    plt.tight_layout()

    if plot is not None and plot[3] == 1 and outputpath is not None:
        plt.savefig(f"{outputpath}/{filename}", dpi=300)
        plt.close(fig)
    else: 
        plt.show()
