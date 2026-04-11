# Configuration

FREQCOR is configured using `.ini` files (read with Python `configparser`). A configuration
file defines:

- where to find input files (binned (co)spectra, EddyPro `full_output`, optional QC files),
- which processing mode to run (cospectra vs spectra, gas selection, plotting options),
- thresholds used for data selection and quality screening,
- where to write outputs.

An example configuration is provided in:

`examples/metadata/FREQCOR_config_example.ini`

## File structure

The INI file is organized into sections:

- `[IO]`: input/output locations and read routine.
- `[TIME_WINDOW]`: optional time window filtering.
- `[PROCEDURE_OPTIONS]`: main algorithm switches (mode, gas, fitting, plotting, LUT classes).
- `[EC_SETUP]`: site metadata and sensor separation (spectral mode only).
- `[USER_LIMITS]`: thresholds for filtering and transfer-function checks.
- `[PROCESSED_DATA]`: optional caching of intermediate read outputs.

## [IO]

**Purpose:** select the reader and point to all required input files.

- `read_routine`
  - `FREQCOR_Read_EP`: standard workflow using EddyPro binned cospectra exports.
  - `FREQCOR_Read_TOF`: alternative workflow for GEddySoft/TOF-type datasets.

- `input_path`
  - Folder containing the EddyPro `full_output` file.

- `input_path_sp`
  - Folder containing the binned (co)spectra (used by some workflows).

- `binned_cosp`
  - Substring used to identify binned cospectra files/folders.

- `flx_meteo`
  - Base name of the EddyPro flux/meteo/QC file (typically `full_output`).

- `vitale_path`
  - Folder containing the ICOS-ETC global flux flag file (Vitale et al., 2020), used only
    when `vitale_qc_flags = 1`.

- `massman_path`
  - Optional Massman coefficient file (used for plotting reference cospectra).

- `output_path`, `output_file`
  - Output directory and base name for the corrected flux output file.

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

- Plotting options (`0`/`1`)
  - `plot_hh`: individual half-hour plots (very verbose).
  - `plot_main`: main summary plots.
  - `plot_aux`: additional diagnostic plots.
  - `plot_save`: save figures to `output_path` instead of interactive display.

- Vitale QC integration
  - `vitale_qc_flags`
    - `0`: use internal VM/F quality thresholds.
    - `1`: use ICOS-ETC global flux flags (Vitale et al., 2020) from `vitale_path`.

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
  - `hlim`, `fclim`, `lelim`: absolute-value thresholds for selecting (co)spectra used for
    cut-off frequency estimation.
  - `hlcfu`: sensible-heat threshold for selecting cospectra used for CF (unstable).
  - `hlcfs`: sensible-heat threshold for selecting cospectra used for CF (stable).

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

## Validation

At runtime, FREQCOR validates the configuration to catch common issues early (e.g. invalid
enumerations, inconsistent ranges, output path creation, and time window parsing).
