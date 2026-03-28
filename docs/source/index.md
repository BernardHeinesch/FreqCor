# Welcome to FREQCOR Documentation

FREQCOR (**FREQ**uency **COR**rection) is a Python toolbox to compute and apply spectral/cospectral correction factors to eddy covariance fluxes.

It produces:
- cut-off frequency look-up tables (LUTs), and
- correction factor LUTs (stable/unstable), which can then be applied to half-hourly flux time series.

Two types of inputs are currently supported: (i) outputs from an EddyPro run (binned (co)spectra and the full_output file), or (ii) outputs from the GEddySoft software, a flux computation toolbox developed by the same research team (https://github.com/BernardHeinesch/GEddySoft). Handling of other input formats can be easily implemented by replacing the FREQCOR_Read routine by your own.

The workflow implemented in FREQCOR follows a classical frequency-domain approach:
1) read and quality-screen binned (co)spectra and associated flags,
2) estimate transfer functions and cut-off frequencies from ideal vs real (co)spectra,
3) compute correction factors by integrating degraded vs ideal reference (co)spectra,
4) apply the resulting correction factors.

This tool has been written by Ariane Faurès, from the BioDynE research axis, Gembloux Agro-Bio Tech, University of Liège in Belgium and its version 1.3 has been used in the study "A multi-site comparison of spectral and co-spectral approaches for correction of turbulent gas fluxes with ICOS set-up", to appear soon in the Atmospheric Measurements and Technologies journal peer-review process.

## Key capabilities

- **Cospectral and spectral modes**
  - `sps = 1`: cospectra-based approach
  - `sps = 2`: spectra-based approach (optionally with sensor separation correction)

- **Look-up tables capturing meteorological influences**
  - Wind speed, wind direction, and (for H2O) relative humidity dependencies are represented through LUTs for cut-off frequencies and correction factors

- **Transfer function fitting**
  - Lorentzian and Gaussian forms
  - Optional Peltola et al. (2021) fit form for the cospectral approach
  - Optional denoising support (Aslan et al. 2021)

- **Batch processing**
  - CLI-driven multi-site runs from a directory of configuration files
  - Optional flat mode to process a single directory of `.ini` files

- **Plotting and diagnostics**
  - Optional plots for (co)spectra, transfer functions, and LUT diagnostics for cut-off frequencies and correction factors
  - Plot saving controlled by configuration switches

- **Data-quality screening and selection**
  - Vickers & Mahrt (1997) global flags aggregation produced by eddyPro can be replaced by Vitale (2020) flagging system
  - Flux-based thresholds for selecting (co)spectra used in LUT construction

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
