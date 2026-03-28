# -*- coding: utf-8 -*-
"""
Script to write an .ini file in a python-dictionary format.

Created on Thu Mar 2 16:11:33 2023

@author: Ariane Faurès
"""
import configparser
import os

os.chdir(r'C:\Ariane_FaurÃ¨s_Archives\Outils_traitement\Python\GitLab\freqcor')

config = configparser.ConfigParser()

config['IO'] = {
    '#PATHS_1': 'Path to input files: the folder should contain a full_output file',
    'input_path': r'C:\ICOS_Treatment\Database_Lonzee\2_Eddy_covariance\3_LI7200\1_FLUX\EP\202206',
    '#PATHS_2': 'Path to input files: the folder should contain a eddypro_binned_cospectra folder',
    'input_path_sp': r'C:\ICOS_Treatment\Database_Lonzee\2_Eddy_covariance\3_LI7200\1_FLUX\EP\202206',
    '#PATHS_2b': 'Read routine selection. FREQCOR_Read_EP (standard EddyPro binned cospectra) or FREQCOR_Read_TOF (GEddySoft/TOF)',
    'read_routine': 'FREQCOR_Read_EP',
    '#PATHS_3': 'Path to output folder',
    'output_path': r'C:\ICOS\BE-Lon\1_Analyses\Spectral_analysis_Aubinet\outputs_code\202206\co2\cospec\peltola\windsurf',
    '#PATHS_4': 'Path (full, including file name) to file containing Massman coefficients for the site of interest. Optional. Used for plotting of the reference Massman cospectra only.',
    'massman_path': r'input\data\Massman_coef.csv',
    '#PATHS_5': 'Path to input file (csv) containing the ICOS-ETC global flux flag computed according to Vitale et al. 2020. Is only used when the "vitale_qc_flags" option is set to 1.',
    'vitale_path': r'C:\Users\Ariane Faurès\Documents\Multi-site\spectral_corrections_paper\data\icos_sites\c_outputs\0_from_archive_L2_cp\vitale_qc_files',
    '#NAMES_1': 'Name of binned cospectra files (and folder)',
    'binned_cosp': '_binned_cospectra_',
    '#NAMES_2': 'Name of file containing flux, meteo and QC data.',
    'flx_meteo': 'full_output',
    '#NAMES_3': 'Name of output flux file with corrected fluxes',
    'output_file': 'Lonzee_June_2022_out_cosp',
}

config['TIME_WINDOW'] = {
    '#1': 'Enable time window filtering (0: disabled, 1: enabled)',
    'enable_time_window': 0,
    '#2': 'Start date and time of the processing window. Format: YYYY-MM-DD HH:MM',
    'start_datetime': '2022-06-01 00:00',
    '#3': 'End date and time of the processing window. Format: YYYY-MM-DD HH:MM',
    'end_datetime': '2022-06-30 23:30',
    '#4': 'Enable exclusion of sub-periods within the processing window (0: disabled, 1: enabled)',
    'enable_exclusion_windows': 0,
    '#5': 'Start date(s) of exclusion window(s). Format: YYYY-MM-DD HH:MM. Multiple windows: comma-separated (e.g. 2022-06-05 00:00, 2022-06-20 12:00)',
    'date_exclusion_start': '2022-06-05 00:00',
    '#6': 'End date(s) of exclusion window(s). Format: YYYY-MM-DD HH:MM. Must match the number of entries in date_exclusion_start',
    'date_exclusion_end': '2022-06-10 23:30',
}

config['PROCEDURE_OPTIONS'] = {
    '#1': 'Procedure. Cospectra = 1; Spectra = 2',
    'sps': 1,
    '#2': 'Gas. CO2 = 1; H2O = 2; O3 = 3; CH4 = 4; N2O = 5',
    'gss': 1,
    '#5': 'Number of classes for cut-off frequency look-up tables first sorting level. The first sorting variable is wind direction (wd) for non-h2o gases and relative humidity (rh) for h2o. Used both for cof and CF computation. If set to 1, the number of classes for the second sorting level is considered.',
    'classnumwd': 1,
    'classnumrh': 8,
    '#6': 'Number of classes for cut-off frequency look-up tables second sorting level (wind speed). The number of classes of "ws" defined by the user is considered only if number of classes for the first level of sorting for the selected gas is set to 1 (#5). Otherwise, by default 4 ws classes are used (hard-coded).',
    'classnumcofws': 5,
    '#7': 'Number of classes for correction factor look-up tables second sorting level, unstable (u) and stable (s). This user input is considered only if the number of classes at  #5 is set to 1, otherwise their value is hard-coded.',
    'classnumcf_u': 8,
    'classnumcf_s': 8,
    '#10': 'Horst and Lenschow (2009) equation to use in sensor separation correction factor computation. Can be either 13 or 16. It is only used when #1 = 2 (procedure : Spectra)',
    'eq': 16,
    '#11': 'Anemometer losses for path averaging and time response. Applies the theoretical correction on reference cospectra before computing the correction factors.',
    'tf_sonic': 1,
    '#12': 'Transfer function fit according to Peltola et al. 2021 (square-root of H). Activated = 1; Inactivated = 0. If activated, will be performed for cospectral approach only (#1 = 1)',
    'tf_peltola': 1,
    '#13': 'Plot options (0 = disabled, 1 = enabled). plot_hh: individual half-hours; plot_main: main results (averaged/unified (co)spectra, TF, CF vs WS); plot_aux: auxiliary/diagnostic (all individual, filtering, additional diagnostics); plot_save: save to output_path',
    'plot_hh': 0,
    'plot_main': 1,
    'plot_aux': 1,
    'plot_save': 1,
    '#14': 'Option to use the ICOS-ETC global data flux flag computed according to Vitale et al. 2020. When activated (= 1), the flags stored in the dedicated csv file accessible at "vitale_path" is used instead of the VM and F flags defined above.', 
    'vitale_qc_flags': 1,
}

config['EC_SETUP'] = {
    '#1': 'Site name',
    'jsite': 'BE-Lon',
    '#2': 'Measurement height (m). It is only used when [PROCEDURE_OPTIONS] #1 = 2 (procedure : Spectra), for the sensor separation term.',
    'zmeas': 2,
    '#3': 'Displacement height (m). It is only used when [PROCEDURE_OPTIONS] #1 = 2 (procedure : Spectra), for the sensor separation term.',
    'disph': 0.536,
    '#4': 'Vertical sensor separation (m). It is only used when [PROCEDURE_OPTIONS] #1 = 2 (procedure : Spectra), for the sensor separation term.',
    'rsz': -0.04,
    '#5': 'Eastward sensor separation (m). It is only used when [PROCEDURE_OPTIONS] #1 = 2 (procedure : Spectra), for the sensor separation term.',
    'rseast': -0.1,
    '#6': 'Northward sensor separation (m). It is only used when [PROCEDURE_OPTIONS] #1 = 2 (procedure : Spectra), for the sensor separation term.',
    'rsnorth': 0.17,
}

config['USER_LIMITS'] = {
    '#1': 'Wind sector exclusion range (Â° from north). Data included in this range will be discarded from the correction procedure.',
    'wdmin': 0,
    'wdmax': 0,
    '#2': 'Limits for (co)spectra selection for cof computation. h : sensible heat [Wm-2]; fc : main tracer [umolm-2s-1]; le : latent heat [Wm-2]. The flux and limit are used in absolute values',
    'hlim': 50,
    'fclim': 10,
    'lelim': 30,
    '#4': 'Limits for cospectra selection for CF computation, unstable conditions. On sensible heat only [Wm-2]. Only fluxes greater than this threshold are kept.',
    'hlcfu': 30,
    '#5': 'Limits for cospectra selection for CF computation, stable conditions. On sensible heat only [Wm-2]. Only fluxes smaller than this threshold are kept.',
    'hlcfs': 0,
    '#6': 'Quality flags for Vickers and Mahrt (VM) tests. The VM tests used, when available, are: spikes, drop-out, skewness hard, discontinuities hard. 0: limit value. If > 0, at least a test failed, and the data is discarded.',
    'fvm': 0,
    '#7': 'Quality flags for Foken (F) tests (FS and ITCÏƒ combined): 0 if both < 30%, 1 if any is 30-100%, 2 if any > 100%. If value larger than the threshold, the data is discarded.',
    'ff': 0,
    '#8': 'Frequency range for transfer function (TF) quality check [Hz]. The coefficient of variation of the experimental TF is checked in this range and compared to #9 (varclim).',
    'jtmin': 0.021,
    'jtmax': 0.34,
    '#9': 'Limit for the coefficient of variation (std dev/mean) of the TF quality check range. If it is > varclim, the TF is discarded.',
    'varclim': 1,
    '#10': 'Frequency range for TF fitting [Hz]',
    'tfmin': 0.021,
    'tfmax': 10,
    '#11': 'Frequency range for local similarity procedure: factor of cof between low and high frequencies',
    'nulim': 0.05,
}

config['PROCESSED_DATA'] = {
    '#1': 'Save Read routine output data to files for later use (0: disabled, 1: enabled)',
    'enable_saving': 1,
    '#2': 'Inclusion of timestamp in saved intermediate file names (0: use site name only, 1: include timestamp). Applies only if enable_saving = 1',
    'use_timestamp': 0,
    '#3': 'Load previously saved data to speed up reading time (0: disabled, 1: enabled)',
    'enable_loading': 0,
    '#4': 'Prefix for the saved files ("site-name" or "site-name_timestamp")',
    'file_prefix': 'BE-Lon',
}

with open('FREQCOR_input.ini', 'w') as configfile:
    config.write(configfile)

