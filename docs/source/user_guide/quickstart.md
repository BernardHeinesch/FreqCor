# Quick Start

Get started with FREQCOR in a few simple steps.

## Basic usage

1. Use the example configuration: `examples/metadata/FREQCOR_config_example.ini`.

2. Edit the INI file paths:
   - In `[IO]`, update:
     - `input_path` and `input_path_sp` to point to your EddyPro outputs (or example data if provided).
     - `output_path` to an output folder you can write to.
   - Review your processing switches in `[PROCEDURE_OPTIONS]` (cospectra vs spectra, gas, plotting).

3. Run the processing:
   - Open `src/FREQCOR_Start.py` in your Python editor.
   - Locate the *MANUAL MODE* section near the top of the file and set the `ini_file` path to your INI.
   - Run `FREQCOR_Start.py` from the `src` directory.

4. Inspect outputs:
   - LUT files and corrected flux outputs are written to `output_path`.
   - If plotting is enabled, figures are saved/displayed according to `plot_save`.

For details on all configuration keys, see the **Configuration** page.
