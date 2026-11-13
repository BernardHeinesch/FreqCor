# FREQCOR

FREQCOR (**FREQ**uency **COR**rection) is a Python toolbox to compute and apply spectral/cospectral correction factors to eddy covariance fluxes.

It produces:

- cut-off frequency look-up tables (LUTs), and
- correction factor LUTs (stable/unstable), which can then be applied to half-hourly flux time series.

Two types of inputs are currently supported: (i) outputs from an EddyPro run (binned (co)spectra and the full_output file), or (ii) outputs from the GEddySoft software (https://github.com/BernardHeinesch/GEddySoft).

Handling of other input formats can be easily implemented by replacing the read routine by your own.

## Documentation

The HTML documentation is available online via GitHub Pages:
[https://bernardheinesch.github.io/FreqCor/](https://bernardheinesch.github.io/FreqCor/)

## Installation

### Option 1: pip + venv

From the repository root:

```bash
python -m venv .venv
# Windows PowerShell:
.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

### Option 2: conda

```bash
conda create -n freqcor python=3.10
conda activate freqcor
pip install -r requirements.txt
```

## Quick start

1. Use one of the example configurations (same dataset, two input formats):

   - `examples/metadata/FREQCOR_config_example_EddyPro.ini` (EddyPro inputs)
   - `examples/metadata/FREQCOR_config_example_GEddySoft.ini` (GEddySoft inputs)

2. Run the entry point:

   - Open `src/FREQCOR_Start.py`.
   - In the *MANUAL MODE* section, set `ini_file` to the example INI path.
   - Run `FREQCOR_Start.py`.

Outputs are written under `examples/output/` (or under the `output_path` configured in the INI file).

## Example datasets

The repository ships with two example datasets from the ICOS Lonzée station (BE-Lon), covering the **same time window** — two spring months of 2024 over a winter wheat crop — processed with two different flux software: **EddyPro** and **GEddySoft**.

Both example INI files are configured for a **CO2 cospectral** run, making it easy to compare the two input formats on identical data.

## Description

FreqCor (Frequency Correction) is a Python toolbox to compute and apply spectral/cospectral correction factors to eddy covariance fluxes.

It produces:

- cut-off frequency look-up tables (LUTs)
- correction factor LUTs (stable/unstable), which can then be applied to half-hourly flux time series

The main entry point is `src/FREQCOR_Start.py`, which reads an `.ini` configuration file and runs the workflow through `FREQCOR_Main`.

## Authors

This tool has been written by Ariane Faurès (BioDynE research axis, Gembloux Agro-Bio Tech, University of Liege, Belgium).

## Citation

If you use FREQCOR in your work, please cite the associated manuscript when available.

## License

Apache-2.0