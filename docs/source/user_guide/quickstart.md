# Quick Start

Get started with FREQCOR in a few simple steps.

## Basic usage

1. Use one of the example configurations shipped with the repository (same dataset, two input formats):
   - `examples/metadata/FREQCOR_config_example_EddyPro.ini` (EddyPro inputs)
   - `examples/metadata/FREQCOR_config_example_GEddySoft.ini` (GEddySoft inputs)

2. Edit the INI file paths — either directly, or with the graphical configurator `python tools/FREQCOR_ini_gui.py` (recommended for a first configuration: it groups options in tabs, disables irrelevant fields automatically, and validates the INI on save):
   - In `[IO]`, update the input file patterns (e.g. `binned_cosp`, `flx_meteo`) to point to your own outputs — or keep the shipped example data paths.
   - `output_path` to an output folder you can write to.
   - Review your processing switches in `[PROCEDURE_OPTIONS]` (cospectra vs spectra, gas, plotting, etc.). The proposed ini is well adapted for an ICOS set-up.

3. Run the processing:
   - Open `src/FREQCOR_Start.py` in your Python editor.
   - Locate the *MANUAL MODE* section near the top of the file and set the `ini_file` path to your INI.
   - Run `FREQCOR_Start.py` from the `src` directory.

4. Inspect outputs:
   - LUT files and corrected flux outputs are written to `output_path`.
   - If plotting is enabled, figures are saved/displayed according to `plot_save`.

For details on all configuration keys, see the **Configuration** page.
