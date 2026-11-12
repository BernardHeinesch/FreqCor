# -*- coding: utf-8 -*-
"""
Created on Tue Apr 8 2025

Validation functions for FREQCOR configuration parameters
"""
import os
import glob
import datetime
import pandas as pd

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
        plot_hh = int(config['PLOTS']['plot_hh'])
        plot_save = int(config['PLOTS']['plot_save'])
        plot_main = int(config['PLOTS']['plot_main'])
        plot_aux = int(config['PLOTS']['plot_aux'])
        vitale_qc_flags = int(config['PROCEDURE_OPTIONS']['vitale_qc_flags'])
        rh_from_file = int(config['PROCEDURE_OPTIONS'].get('rh_from_file', 0))
        
        # Paths
        outputpath = config['IO']['output_path']
        flx_meteo = config['IO'].get('flx_meteo', '')
        binned_cosp = config['IO']['binned_cosp']
        # vitale_path and rh_path are only used by FREQCOR_Read_EP
        if read_routine == 'FREQCOR_Read_EP':
            if vitale_qc_flags == 1:
                vitale_path = config['IO']['vitale_path']
            if rh_from_file == 1:
                rh_path = config['IO']['rh_path']
        
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
    if read_routine == 'FREQCOR_Read_EP':
        expected_flux_files = glob.glob(flx_meteo)
        if len(expected_flux_files) == 0:
            valid = False
            errors.append(f"No input flux/meteo CSV found matching flx_meteo pattern: {flx_meteo}")
        elif len(expected_flux_files) > 1:
            valid = False
            errors.append(
                f"Multiple input flux/meteo CSV files found matching flx_meteo pattern: {flx_meteo} ({len(expected_flux_files)} matches)"
            )

    expected_binned_files = glob.glob(binned_cosp)
    if len(expected_binned_files) == 0:
        valid = False
        errors.append(f"No binned cospectra file found matching binned_cosp pattern: {binned_cosp}")

    if read_routine == 'FREQCOR_Read_EP' and vitale_qc_flags == 1:
        vitale_csv_files = glob.glob(vitale_path)
        if len(vitale_csv_files) == 0:
            valid = False
            errors.append(f"vitale_qc_flags=1 but no Vitale QC CSV found matching vitale_path pattern: {vitale_path}")
        elif len(vitale_csv_files) > 1:
            valid = False
            errors.append(
                f"vitale_qc_flags=1 but multiple Vitale QC CSV files found matching vitale_path pattern: {vitale_path} ({len(vitale_csv_files)} matches)"
            )
        else:
            try:
                vitale_cols = pd.read_csv(vitale_csv_files[0], nrows=0).columns
                missing_cols = []
                required_cols = [
                    'TIMESTAMP_END',
                    'H_SSITC_TEST',
                    'H_DATA_FLAG',
                    'WSECT_FLAG',
                ]
                if gss == 1:
                    required_cols.extend(['FC_SSITC_TEST', 'NEE_DATA_FLAG'])
                elif gss == 2:
                    required_cols.extend(['LE_SSITC_TEST', 'LE_DATA_FLAG'])
                for col in required_cols:
                    if col not in vitale_cols:
                        missing_cols.append(col)
                if missing_cols:
                    valid = False
                    errors.append(
                        "vitale_qc_flags=1 but Vitale/ICOS FLUX_L2 CSV is missing required column(s): "
                        + ", ".join(missing_cols)
                    )
            except Exception as e:
                valid = False
                errors.append(f"Could not read Vitale/ICOS FLUX_L2 CSV header to validate columns: {e}")

    if rh_from_file not in [0, 1]:
        valid = False
        errors.append("PROCEDURE_OPTIONS.rh_from_file must be 0 or 1")
    elif rh_from_file == 1 and read_routine == 'FREQCOR_Read_EP':
        rh_csv_files = glob.glob(rh_path)
        if len(rh_csv_files) == 0:
            valid = False
            errors.append(f"rh_from_file=1 but no RH CSV found matching rh_path pattern: {rh_path}")
        elif len(rh_csv_files) > 1:
            valid = False
            errors.append(
                f"rh_from_file=1 but multiple RH CSV files found matching rh_path pattern: {rh_path} ({len(rh_csv_files)} matches)"
            )
        else:
            try:
                rh_cols = pd.read_csv(rh_csv_files[0], nrows=0).columns
                missing_cols = []
                if 'TIMESTAMP_END' not in rh_cols:
                    missing_cols.append('TIMESTAMP_END')
                if 'RH' not in rh_cols:
                    missing_cols.append('RH')
                if missing_cols:
                    valid = False
                    errors.append(
                        "rh_from_file=1 but RH CSV is missing required column(s): " + ", ".join(missing_cols)
                    )
            except Exception as e:
                valid = False
                errors.append(f"Could not read RH CSV header to validate columns: {e}")

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
    if read_routine not in ['FREQCOR_Read_EP', 'FREQCOR_Read_GEddySoft']:
        valid = False
        errors.append("read_routine must be one of: FREQCOR_Read_EP, FREQCOR_Read_GEddySoft")

    # GEddySoft-specific checks
    if read_routine == 'FREQCOR_Read_GEddySoft':
        src_type = config['IO'].get('source_type', '').strip().lower()
        if src_type not in ('irga', 'tof', 'tracer'):
            valid = False
            errors.append("IO.source_type must be 'irga' or 'tof' "
                          "when read_routine = FREQCOR_Read_GEddySoft")
        target_mz_str = config['IO'].get('target_mz', '').strip()
        if src_type in ('tof', 'tracer') and not target_mz_str:
            valid = False
            errors.append("IO.target_mz is required when source_type = 'tof'")
        elif target_mz_str:
            try:
                float(target_mz_str)
            except ValueError:
                valid = False
                errors.append(f"IO.target_mz must be a number (got '{target_mz_str}')")
    
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

