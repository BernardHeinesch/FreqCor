# Welcome to FREQCOR Documentation

FREQCOR (**FREQ**uency **COR**rection) is a Python toolbox to compute and apply spectral/cospectral correction factors to eddy covariance fluxes.

It produces:
- cut-off frequency look-up tables (LUTs), and
- correction factor LUTs (stable/unstable), which can then be applied to half-hourly flux time series.

Two types of inputs are currently supported: (i) outputs from an EddyPro run (binned (co)spectra and the full_output file), or (ii) outputs from the GEddySoft software, a flux computation toolbox developed by the same research team (https://github.com/BernardHeinesch/GEddySoft). Handling of other input formats can be easily implemented by replacing the FREQCOR_Read routine by your own.

The repository ships with **two example datasets covering the same time window** at the ICOS Lonzée station (BE-Lon): one processed with **EddyPro** and one processed with **GEddySoft**, each with its own example configuration file (`FREQCOR_config_example_EddyPro.ini` and `FREQCOR_config_example_GEddySoft.ini` in `examples/metadata/`). This makes it easy to compare the two input formats on identical data.

The workflow implemented in FREQCOR follows a classical frequency-domain approach:
1) read and quality-screen binned (co)spectra and associated flags,
2) estimate transfer functions and cut-off frequencies from ideal vs real (co)spectra,
3) compute correction factors by integrating degraded vs ideal reference (co)spectra,
4) apply the resulting correction factors.

This tool has been written by Ariane Faurès, from the BioDynE research axis, Gembloux Agro-Bio Tech, University of Liège in Belgium and its version 1.3 has been used in the study "Spectral versus co-spectral correction of eddy-covariance fluxes: a multi-site assessment across the ICOS network", published in the Atmospheric Measurements and Technologies journal ([doi.org/10.5194/amt-19-6037-2026](https://doi.org/10.5194/amt-19-6037-2026)).

## Key capabilities

- **Cospectral and spectral modes**
  - `sps = 1`: cospectra-based approach
  - `sps = 2`: spectra-based approach (includes sensor separation correction following Horst & Lenschow, 2009 — equation 13 or 16 selectable via `eq`)

- **Look-up tables capturing meteorological influences**
  - Wind speed, wind direction, and (for H2O) relative humidity dependencies are represented through LUTs for cut-off frequencies and correction factors

- **Transfer function fitting**
  - Lorentzian and Gaussian forms
  - Optional Peltola et al. (2021) fit form for the cospectral approach
  - (Co-)Spectra denoising adapted from Aslan et al. 2021

- **Data-quality screening and selection**
  - Vickers & Mahrt (1997) global flags aggregation produced by eddyPro can be replaced by Vitale (2020) flagging system
  - Selection of (co)spectra for LUT construction via flux thresholds, wind-sector exclusion, quality flags and spectral-value caps

- **Graphical INI configurator**
  - `tools/FREQCOR_ini_gui.py`: a standalone wizard to build configuration files — options grouped in tabs, dependent options enabled/disabled automatically, contextual help with explanatory schematics, and INI validation on save using the same routine as run-time checks

- **Batch processing**
  - CLI-driven multi-site runs from a directory of configuration files
  - Optional flat mode to process a single directory of `.ini` files

- **Plotting and diagnostics**
  - Optional plots for (co)spectra, transfer functions, and LUT diagnostics for cut-off frequencies and correction factors
  - Plot saving controlled by configuration switches

- **Outputs**
  - LUTs for cut-off frequencies and correction factors
  - corrected fluxes (Lorentz/Gauss) and associated diagnostics

## Contents

```{toctree}
:maxdepth: 2
:caption: User Guide

user_guide/installation
user_guide/quickstart
user_guide/batch_mode
user_guide/configuration
user_guide/running_example
user_guide/outputs
contact
```

```{toctree}
:maxdepth: 2
:caption: Project Info

version
```

```{toctree}
:maxdepth: 2
:caption: Implementation Details

theory/software_architecture
```

```{toctree}
:caption: API Reference
:maxdepth: 2

api/index
```
