# Software Architecture

FREQCOR follows a modular architecture designed to compute spectral/cospectral correction factors from
binned (co)spectra (e.g., EddyPro exports). The codebase is organized around a small set of orchestration
modules that drive the processing chain, and a collection of scientific/processing modules that implement
selection, LUT construction, correction-factor computation, and plotting.

## Processing flow

```
FREQCOR_Start
    |
    v
FREQCOR_Main
    |
    |-> Configuration validation
    |     |- FREQCOR_validate
    |
    |-> Input reading (EddyPro / GEddySoft)
    |     |- FREQCOR_Read_EP
    |     |- FREQCOR_Read_GEddySoft
    |     |- FREQCOR_Read_TOF
    |
    |-> Computation pipeline
          |
          |-> Select (co)spectra used to build LUTs
          |     |- FREQCOR_Sel_general
          |     |- FREQCOR_Sel_stunst
          |
          |-> Cut-off frequency estimation and LUT construction
          |     |- FREQCOR_Sel_cof
          |     |- FREQCOR_cof
          |     |- FREQCOR_LUT_cof
          |
          |-> Correction factor estimation and LUT construction
          |     |- FREQCOR_Sel_CF
          |     |- FREQCOR_LUT_CF
          |
          |-> Apply correction factors to flux time series
          |     |- FREQCOR_Flux
          |
          |-> Output generation
                |- FREQCOR_write_outputs
                |- write_ini
                |- FREQCOR_plot
```

## Module categories

### Core orchestration

- `FREQCOR_Start`
  - Script-like entry point (manual mode and batch patterns).
- `FREQCOR_Main`
  - High-level dispatcher that selects the read routine and orchestrates the run.
- `FREQCOR_Compute`
  - Main computation chain once inputs are loaded.

### Input and validation

- `FREQCOR_validate`
  - Configuration checks and defaults.
- `FREQCOR_Read_EP`
  - Reader for EddyPro exports (binned (co)spectra and associated files).
- `FREQCOR_Read_GEddySoft`, `FREQCOR_Read_TOF`
  - Alternative readers for GEddySoft/TOF-like datasets.

### Processing and scientific core

- **Selection**
  - `FREQCOR_Sel_general`: baseline filtering and preparation.
  - `FREQCOR_Sel_stunst`: stability split (stable/unstable classes).
  - `FREQCOR_Sel_cof`: selection for cut-off frequency estimation.
  - `FREQCOR_Sel_CF`: selection for correction-factor estimation.

- **Cut-off frequencies and transfer functions**
  - `FREQCOR_cof`: cut-off frequency estimation.
  - `FREQCOR_LUT_cof`: cut-off frequency LUT construction.

- **Correction factors**
  - `FREQCOR_LUT_CF`: correction factor LUT construction.
  - `theor_cosp_Kaimal`: theoretical reference (co)spectra.
  - `FREQCOR_Sensor_Separation`: sensor-separation correction (spectral workflow).

- **Quality screening**
  - `FREQCOR_VM_flag`: Vickers & Mahrt-type flagging utilities.

### Output and diagnostics

- `FREQCOR_write_outputs`, `write_ini`
  - Writing LUTs, flux outputs, and configuration snapshots.
- `FREQCOR_plot`
  - Diagnostic plots for selection, transfer functions, LUTs, and correction factors.

## Data flow

1. **Inputs**
   - A configuration (`.ini`) file selects the reader and points to the binned (co)spectra and associated
     EddyPro (or GEddySoft) outputs.

2. **Selection and quality screening**
   - (Co)spectra are filtered using quality flags and threshold criteria.

3. **LUT construction**
   - Cut-off frequencies and correction factors are estimated and aggregated into LUTs (typically split
     into stable/unstable regimes).

4. **Application to fluxes**
   - The correction factor LUTs are applied to half-hourly flux time series.

5. **Outputs**
   - LUTs, corrected fluxes, and diagnostic plots are written to the output folder.
