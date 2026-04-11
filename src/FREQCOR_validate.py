# -*- coding: utf-8 -*-
"""
Created on Tue Apr 8 2025

Validation functions for FREQCOR configuration parameters
"""
import os
import glob
import datetime

def validate_config(config):
    """
    Validate FREQCOR configuration parameters to avoid later errors.

    Parameters
    ----------
    config : configparser.ConfigParser
        Configuration object with all parameters, including PROCEDURE_OPTIONS and USER_LIMITS sections.

    Returns
    -------
    bool
        True if all validations pass, False otherwise.
    list
        List of error messages if any validation fails.
    """
    valid = True
    errors = []

    def _warn_windows_path_length(*, base_dir, sample_filenames, limit=250, margin=10):
        if os.name != 'nt':
            return
        if not isinstance(base_dir, str) or not base_dir:
            return
        for name in sample_filenames:
            full_path = os.path.join(base_dir, name)
            n = len(full_path)
            if n >= (limit - margin):
                print(
                    "WARNING: output_path may be too long for Windows when combined with output filenames. "
                    "This can cause save failures (e.g. plots/CSVs not written) and lead to incomplete outputs. "
                    "Reduce/shorten output_path to avoid this. "
                    f"Example path length={n} (limit~{limit}): {full_path}"
                )
                return
    
    # Extract parameters for validation
    try:
        # PROCEDURE_OPTIONS
        gss = int(config['PROCEDURE_OPTIONS']['gss'])
        sps = int(config['PROCEDURE_OPTIONS']['sps'])
        tf_sonic = int(config['PROCEDURE_OPTIONS']['tf_sonic'])
        tf_peltola = int(config['PROCEDURE_OPTIONS']['tf_peltola'])
        eq = int(config['PROCEDURE_OPTIONS']['eq'])
        classnum = int(config['PROCEDURE_OPTIONS']['classnumcofws'])
        classnumCF_u = int(config['PROCEDURE_OPTIONS']['classnumcf_u'])
        classnumCF_s = int(config['PROCEDURE_OPTIONS']['classnumcf_s'])
        if gss == 2 : classnum_first = int(config['PROCEDURE_OPTIONS']['classnumrh']) 
        else : classnum_first = int(config['PROCEDURE_OPTIONS']['classnumwd'])

        read_routine = config['IO'].get('read_routine', None)

        # USER_LIMITS
        jtmin_hz = float(config['USER_LIMITS']['jtmin'])
        jtmax_hz = float(config['USER_LIMITS']['jtmax'])
        tfmin_hz = float(config['USER_LIMITS']['tfmin'])
        tfmax_hz = float(config['USER_LIMITS']['tfmax'])
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
        
        # Plot options
        plot_hh = int(config['PROCEDURE_OPTIONS']['plot_hh'])
        plot_save = int(config['PROCEDURE_OPTIONS']['plot_save'])

        plot_main = int(config['PROCEDURE_OPTIONS']['plot_main'])
        plot_aux = int(config['PROCEDURE_OPTIONS']['plot_aux'])
        vitale_qc_flags = int(config['PROCEDURE_OPTIONS']['vitale_qc_flags'])
        
        # Paths
        outputpath = config['IO']['output_path']
        input_path = config['IO']['input_path']
        input_path_sp = config['IO']['input_path_sp']
        flx_meteo = config['IO']['flx_meteo']
        binned_cosp = config['IO']['binned_cosp']
        vitale_path = config['IO']['vitale_path']
        
    except KeyError as e:
        valid = False
        errors.append(f"Missing configuration parameter: {e}")
        return valid, errors
    except ValueError as e:
        valid = False
        errors.append(f"Invalid value type in configuration: {e}")
        return valid, errors
    
    # Validate gss (must be between 1 and 5)
    if not 1 <= gss <= 5:
        valid = False
        errors.append(f"gss must be between 1 and 5, got {gss}")
    
    # Validate sps (must be 1 or 2)
    if sps not in [1, 2]:
        valid = False
        errors.append(f"sps must be 1 or 2, got {sps}")

    # Validate jtmax > jtmin
    if jtmax_hz <= jtmin_hz:
        valid = False
        errors.append(f"jtmax ({jtmax_hz}) must be greater than jtmin ({jtmin_hz})")
    
    # Validate tfmax > tfmin
    if tfmax_hz <= tfmin_hz:
        valid = False
        errors.append(f"tfmax ({tfmax_hz}) must be greater than tfmin ({tfmin_hz})")
    
    # Validate plot options (each must be 0 or 1)
    for name, value in zip(['plot_hh', 'plot_main', 'plot_aux', 'plot_save'], 
                          [plot_hh, plot_main, plot_aux, plot_save]):
        if value not in [0, 1]:
            valid = False
            errors.append(f"{name} must be 0 or 1, got {value}")
    
    # Validate tf_sonic and tf_peltola (must be 0 or 1)
    if tf_sonic not in [0, 1]:
        valid = False
        errors.append(f"tf_sonic must be 0 or 1, got {tf_sonic}")
    if tf_peltola not in [0, 1]:
        valid = False
        errors.append(f"tf_peltola must be 0 or 1, got {tf_peltola}")
    
    # Validate eq (must be 13 or 16)
    if eq not in [13, 16]:
        valid = False
        errors.append(f"eq must be 13 or 16, got {eq}")
    
    # Validate all class numbers (must not be 0)
    if classnum_first == 0:
        valid = False
        errors.append("classnum cannot be 0")
    if classnum == 0:
        valid = False
        errors.append("classnumws cannot be 0")
    if classnumCF_u == 0:
        valid = False
        errors.append("classnumCF_u cannot be 0")
    if classnumCF_s == 0:
        valid = False
        errors.append("classnumCF_s cannot be 0")
    # Validate outputpath (check if folder exists, create if not)
    if not os.path.exists(outputpath):
        try:
            os.makedirs(outputpath)
            print(f"Created output directory: {outputpath}")
        except Exception as e:
            valid = False
            errors.append(f"Failed to create output directory {outputpath}: {e}")

    gas_map = {1: 'co2', 2: 'h2o', 3: 'o3', 4: 'ch4', 5: 'n2o'}
    gas_tag = gas_map.get(gss, f'gss{gss}')
    method_tag = 'cosp' if sps == 1 else 'sp'
    window_tag = 'all'
    try:
        if 'TIME_WINDOW' in config and int(config['TIME_WINDOW'].get('enable_time_window', 0)) == 1:
            start_d = datetime.datetime.strptime(config['TIME_WINDOW']['start_datetime'], '%Y-%m-%d %H:%M')
            end_d = datetime.datetime.strptime(config['TIME_WINDOW']['end_datetime'], '%Y-%m-%d %H:%M')
            window_tag = f"{start_d.strftime('%y%m%d')}-{end_d.strftime('%y%m%d')}"
    except Exception:
        window_tag = 'all'

    _warn_windows_path_length(
        base_dir=outputpath,
        sample_filenames=[
            f"4_mean_TF_all_classes__wd1__SITE__{gas_tag}__{method_tag}__{window_tag}__YYMMDDTHHMM.png",
            f"5_CF_vs_ws__unst__SITE__{gas_tag}__{method_tag}__{window_tag}__YYMMDDTHHMM.png",
            f"5_LUT_CF__SITE__{gas_tag}__{method_tag}__{window_tag}__YYMMDDTHHMM.csv",
            f"7_stats__SITE__{gas_tag}__{method_tag}__{window_tag}__YYMMDDTHHMM.txt",
        ],
    )

    # Validate input paths and expected files
    if not os.path.isdir(input_path):
        valid = False
        errors.append(f"input_path does not exist or is not a directory: {input_path}")
    else:
        expected_flux_files = glob.glob(os.path.join(input_path, f"*{flx_meteo}*.csv"))
        if len(expected_flux_files) == 0:
            valid = False
            errors.append(
                f"No input flux/meteo CSV found in input_path='{input_path}' matching pattern '*{flx_meteo}*.csv'"
            )

    if not os.path.isdir(input_path_sp):
        valid = False
        errors.append(f"input_path_sp does not exist or is not a directory: {input_path_sp}")
    else:
        expected_binned_files = glob.glob(os.path.join(input_path_sp, f"*{binned_cosp}*"))
        if len(expected_binned_files) == 0:
            expected_binned_files = glob.glob(
                os.path.join(input_path_sp, 'eddypro_binned_cospectra', f"*{binned_cosp}*")
            )
        if len(expected_binned_files) == 0:
            valid = False
            errors.append(
                f"No binned cospectra file found in input_path_sp='{input_path_sp}' matching pattern '*{binned_cosp}*' (also checked 'eddypro_binned_cospectra/' subfolder)"
            )

    if vitale_qc_flags == 1:
        if not os.path.isdir(vitale_path):
            valid = False
            errors.append(f"vitale_qc_flags=1 but vitale_path does not exist or is not a directory: {vitale_path}")
        else:
            vitale_csv_files = glob.glob(os.path.join(vitale_path, "*.csv"))
            if len(vitale_csv_files) == 0:
                valid = False
                errors.append(f"vitale_qc_flags=1 but no CSV file found in vitale_path='{vitale_path}'")

    # Validate TIME_WINDOW
    if 'TIME_WINDOW' in config:
        try:
            enable_time_window = int(config['TIME_WINDOW'].get('enable_time_window', 0))
        except ValueError:
            valid = False
            errors.append("TIME_WINDOW.enable_time_window must be 0 or 1")
            enable_time_window = 0

        if enable_time_window not in [0, 1]:
            valid = False
            errors.append("TIME_WINDOW.enable_time_window must be 0 or 1")

        if enable_time_window == 1:
            if 'start_datetime' not in config['TIME_WINDOW'] or 'end_datetime' not in config['TIME_WINDOW']:
                valid = False
                errors.append("TIME_WINDOW.enable_time_window = 1 but start_datetime and/or end_datetime are missing")
            else:
                start_datetime_str = config['TIME_WINDOW']['start_datetime']
                end_datetime_str = config['TIME_WINDOW']['end_datetime']
                try:
                    start_datetime = datetime.datetime.strptime(start_datetime_str, '%Y-%m-%d %H:%M')
                except ValueError:
                    valid = False
                    errors.append(
                        f"TIME_WINDOW.start_datetime '{start_datetime_str}' is not parseable. Expected format: YYYY-MM-DD HH:MM"
                    )
                    start_datetime = None
                try:
                    end_datetime = datetime.datetime.strptime(end_datetime_str, '%Y-%m-%d %H:%M')
                except ValueError:
                    valid = False
                    errors.append(
                        f"TIME_WINDOW.end_datetime '{end_datetime_str}' is not parseable. Expected format: YYYY-MM-DD HH:MM"
                    )
                    end_datetime = None

                if start_datetime is not None and end_datetime is not None:
                    if start_datetime >= end_datetime:
                        valid = False
                        errors.append(
                            f"Invalid TIME_WINDOW: start_datetime ({start_datetime}) must be strictly before end_datetime ({end_datetime})"
                        )

        # Optional exclusion windows
        try:
            enable_exclusion = int(config['TIME_WINDOW'].get('enable_exclusion_windows', 0))
        except ValueError:
            valid = False
            errors.append("TIME_WINDOW.enable_exclusion_windows must be 0 or 1")
            enable_exclusion = 0

        if enable_exclusion not in [0, 1]:
            valid = False
            errors.append("TIME_WINDOW.enable_exclusion_windows must be 0 or 1")

        if enable_exclusion == 1:
            if 'date_exclusion_start' not in config['TIME_WINDOW'] or 'date_exclusion_end' not in config['TIME_WINDOW']:
                valid = False
                errors.append(
                    "TIME_WINDOW.enable_exclusion_windows = 1 but date_exclusion_start and/or date_exclusion_end are missing"
                )
            else:
                excl_starts_raw = [s.strip() for s in config['TIME_WINDOW']['date_exclusion_start'].split(',') if s.strip()]
                excl_ends_raw = [s.strip() for s in config['TIME_WINDOW']['date_exclusion_end'].split(',') if s.strip()]
                if len(excl_starts_raw) != len(excl_ends_raw):
                    valid = False
                    errors.append(
                        f"TIME_WINDOW exclusion windows mismatch: date_exclusion_start has {len(excl_starts_raw)} entries but date_exclusion_end has {len(excl_ends_raw)}"
                    )
                else:
                    for i, (es_raw, ee_raw) in enumerate(zip(excl_starts_raw, excl_ends_raw), start=1):
                        try:
                            es = datetime.datetime.strptime(es_raw, '%Y-%m-%d %H:%M')
                        except ValueError:
                            valid = False
                            errors.append(
                                f"TIME_WINDOW exclusion {i}: date_exclusion_start '{es_raw}' is not parseable. Expected format: YYYY-MM-DD HH:MM"
                            )
                            continue
                        try:
                            ee = datetime.datetime.strptime(ee_raw, '%Y-%m-%d %H:%M')
                        except ValueError:
                            valid = False
                            errors.append(
                                f"TIME_WINDOW exclusion {i}: date_exclusion_end '{ee_raw}' is not parseable. Expected format: YYYY-MM-DD HH:MM"
                            )
                            continue

                        if es >= ee:
                            valid = False
                            errors.append(
                                f"TIME_WINDOW exclusion {i}: start ({es}) must be strictly before end ({ee})"
                            )

    # Validate read_routine
    if read_routine not in ['FREQCOR_Read_EP', 'FREQCOR_Read_TOF']:
        valid = False
        errors.append("read_routine must be one of: FREQCOR_Read_EP, FREQCOR_Read_TOF")
    
    return valid, errors


def print_validation_results(valid, errors):
    """
    Print validation results in a formatted way.
    
    Parameters
    ----------
    valid : bool
        Whether validation passed.
    errors : list
        List of error messages.
    """
    if valid:
        print("Configuration validation: PASSED")
    else:
        print("Configuration validation: FAILED")
        print("Errors found:")
        for i, error in enumerate(errors, 1):
            print(f"  {i}. {error}")

