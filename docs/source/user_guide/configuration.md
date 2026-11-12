# Configuration

FREQCOR is configured using `.ini` files (read with Python `configparser`). A configuration
file defines:

- where to find input files (binned (co)spectra, EddyPro `full_output`, optional QC files),
- some fundamental processing choices (cospectra vs spectra mode, gas selection, fitting and QC options),
- thresholds used for data selection and quality screening,
- where to write outputs.

Example configurations are provided in `examples/metadata/`:

- `FREQCOR_config_example_EddyPro.ini` (EddyPro inputs)
- `FREQCOR_config_example_GEddySoft.ini` (GEddySoft inputs)

## Why a graphical configurator?

FREQCOR exposes **about 60 configuration keys**, of which **~40 are true processing and filtering options**
(algorithm switches, LUT classing, quality thresholds, time windows, transfer-function checks). This richness is
intentional: FREQCOR is designed both as an **operational tool** — routine flux-correction runs by non-experts —
and as a **research tool** where experienced users explore alternative processing choices.

To help non-expert users navigate this space without editing the INI by hand, FREQCOR ships with a graphical
configurator:

```bash
python tools/FREQCOR_ini_gui.py
```

`tools/FREQCOR_ini_gui.py` is a standalone wizard that:

- groups all options in tabs mirroring the INI sections, with cascading enable/disable of dependent options
  (e.g. spectral-only settings are greyed out in cospectral mode),
- shows a syntax-colored INI preview, with hover-to-highlight of the key being edited,
- provides per-parameter help (`?` buttons) including synthetic explanatory schematics for the main choices,
- **validates the resulting INI with the same `FREQCOR_validate.validate_config` routine used at run time**, so
  configuration errors are caught before the run rather than during it.

The tool starts blank: load a starting INI (e.g. one of the example files), edit, and save. It writes a commented
INI preserving the template structure.

## File structure

The INI file is organized into sections:

- `[IO]`: input/output locations and read routine.
- `[TIME_WINDOW]`: optional time window filtering.
- `[PROCEDURE_OPTIONS]`: main algorithm switches (mode, gas, fitting, LUT classes, QC options).
- `[EC_SETUP]`: site metadata and sensor separation (spectral mode only).
- `[USER_LIMITS]`: thresholds for filtering and transfer-function checks.
- `[PROCESSED_DATA]`: optional caching of intermediate read outputs.
- `[PLOTS]`: plotting switches.

## [IO]

**Purpose:** select the reader and point to all required input files.

File locations are given as **full glob patterns** (path + wildcard + extension). Paths may be absolute or
relative to the repository root.

- `read_routine`
  - `FREQCOR_Read_EP`: standard workflow using EddyPro binned cospectra exports.
  - `FREQCOR_Read_GEddySoft`: workflow for GEddySoft daily HDF5 outputs.

- `source_type` (`FREQCOR_Read_GEddySoft` only)
  - `irga`: reads the `IRGA` group (channel selected by name/`gss`).
  - `tof`: reads the `TRACER` group (PTR-TOF-MS); channel selected by `target_mz`.

- `target_mz` (`source_type = tof` only)
  - Target m/z used to pick the TRACER channel (tolerance 0.01).

- `binned_cosp`
  - Glob pattern for the binned cospectra files (EddyPro `*binned_cospectra*_adv.csv`, or GEddySoft `*.hdf5`).

- `flx_meteo`
  - Glob pattern for the EddyPro flux/meteo/QC file (typically `full_output`).
    Not used by `FREQCOR_Read_GEddySoft`.

- `rh_path`
  - Optional glob pattern for an alternative RH CSV (columns `TIMESTAMP_END`, `RH`),
    used only when `rh_from_file = 1`.

- `vitale_path`
  - Glob pattern for the ICOS-ETC global flux flag file (Vitale et al., 2020), used only
    when `vitale_qc_flags = 1`.

- `massman_path`
  - Optional Massman coefficient file (used for plotting reference cospectra).

- `output_path`, `output_file`
  - Output directory and base name for intermediate/cache files.

## [TIME_WINDOW]

**Purpose:** restrict computations to a time period.

- `enable_time_window`
  - `0`: process all available data.
  - `1`: process only the window defined by `start_datetime` and `end_datetime`.

- `start_datetime`, `end_datetime`
  - Format: `YYYY-MM-DD HH:MM`.

- `enable_exclusion_windows`
  - `0`: no exclusion.
  - `1`: exclude sub-periods defined by `date_exclusion_start` and `date_exclusion_end`.

- `date_exclusion_start`, `date_exclusion_end`
  - Comma-separated list of datetimes (same format as above). The number of start and end
    entries must match.

## [PROCEDURE_OPTIONS]

**Purpose:** select the processing mode, gas, sorting strategy, fitting options, and plotting.

- `sps`
  - `1`: cospectral approach.
  - `2`: spectral approach (enables sensor-separation corrections).

- `gss`
  - Gas selector:
    - `1`: CO2
    - `2`: H2O
    - `3`: O3
    - `4`: CH4
    - `5`: N2O

- LUT classing parameters
  - `classnumwd` / `classnumrh`: number of bins for the first sorting variable.
    - for H2O (`gss = 2`) the first sorting variable is RH;
    - otherwise it is wind direction (WD).
  - `classnumcofws`: number of wind-speed classes for cut-off frequency LUTs.
  - `classnumcf_u`, `classnumcf_s`: number of wind-speed classes for correction-factor LUTs
    (unstable and stable conditions).

- Spectral-mode (sensor separation) only
  - `eq`: Horst & Lenschow equation selector (e.g. 13 or 16).

- Transfer-function options
  - `tf_sonic`: apply theoretical correction to sonic cospectra (`0`/`1`).
  - `tf_peltola`: enable Peltola et al. (2021) fit form where applicable (`0`/`1`) (affects the Lorentz fit form only; Gaussian TF is unchanged).

- Vitale QC integration
  - `vitale_qc_flags`
    - `0`: use internal VM/F quality thresholds.
    - `1`: use ICOS-ETC global flux flags (Vitale et al., 2020) from `vitale_path`.

- Alternative RH source
  - `rh_from_file`
    - `0`: RH read from `flx_meteo` (EddyPro `full_output`).
    - `1`: RH read from the alternative CSV pointed to by `rh_path` (aligned on `TIMESTAMP_END`).
    Useful when EddyPro runs lack external RH, in which case RH derived from sonic temperature may be biased.

## [EC_SETUP]

**Purpose:** site metadata and sensor separation geometry.

- `jsite`
  - Site identifier (e.g. ICOS-style `BE-Lon`).

The following are **used only when `sps = 2` (spectral approach)**:

- `zmeas`: measurement height (m)
- `disph`: displacement height (m)
- `rsz`: vertical sensor separation (m)
- `rseast`: eastward sensor separation (m)
- `rsnorth`: northward sensor separation (m)

## [USER_LIMITS]

**Purpose:** numerical thresholds controlling selection and TF checks.

- Wind-sector exclusion
  - `wdmin`, `wdmax`: range (deg from North) to discard.

- Flux thresholds for LUT construction
  - `hl_cof_min`/`hl_cof_max`, `fcl_cof_min`/`fcl_cof_max`, `lel_cof_min`/`lel_cof_max`:
    absolute-value ranges for sensible heat, tracer flux, and latent heat used to select
    (co)spectra for cut-off frequency estimation.
  - `cospec_abs_limit`: absolute limit on individual normalized (co)spectrum points; a single
    exceedance discards the whole half-hour.
  - `hl_cf_u_min`: sensible-heat threshold for selecting cospectra used for CF (unstable).
  - `hl_cf_s_max`: sensible-heat threshold for selecting cospectra used for CF (stable).
  - `hl_cf_abs_max`: additional absolute cap on sensible heat for CF selection (both regimes).

- Quality-flag thresholds
  - `fvm`: Vickers & Mahrt (1997) threshold (discard if VM flag > `fvm`).
  - `ff`: Mauder & Foken (2004) threshold (discard if F flag > `ff`).

- Transfer-function quality check and fitting ranges (Hz)
  - `jtmin`, `jtmax`: frequency range used for TF quality checks.
  - `varclim`: coefficient of variation limit for the TF quality check window.
  - `tfmin`, `tfmax`: frequency range used for TF fitting.

## [PROCESSED_DATA]

**Purpose:** optional saving/loading of intermediate read outputs.

- `enable_saving`: save intermediate read outputs (`0`/`1`).
- `use_timestamp`: include timestamp in filenames (`0`/`1`).
- `enable_loading`: load previously saved outputs to speed up reading (`0`/`1`).
- `file_prefix`: prefix for saved files.

## [PLOTS]

**Purpose:** plotting switches (`0`/`1`).

- `plot_hh`: individual half-hour (co)spectra plots (very verbose).
- `plot_main`: main summary plots (averaged/unified (co)spectra, transfer functions, CF vs wind speed).
- `plot_aux`: auxiliary/diagnostic plots (all individual cospectra, filtering steps, additional diagnostics).
- `plot_save`: save figures to `output_path` (1) or display interactively (0).

```{note}
Plot options must live in `[PLOTS]`. Older INIs that keep them in `[PROCEDURE_OPTIONS]` will fail
validation and must be updated.
```

## Validation

At runtime, FREQCOR validates the configuration with `FREQCOR_validate.validate_config` to catch
common issues early (invalid enumerations, inconsistent ranges, missing input files, output path
creation, time window parsing). The INI configurator (`tools/FREQCOR_ini_gui.py`) runs the same
validation when saving, so a GUI-generated INI is guaranteed to pass the run-time checks.
