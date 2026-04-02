# Running the Example

FREQCOR is shipped with an example dataset from the ICOS Lonzée station (BE-Lon). The repository includes
EddyPro-generated inputs covering two spring months of 2024 over a winter wheat crop.

The provided example configuration file:

`examples/metadata/FREQCOR_config_example.ini`

is set up to run a **cospectral** workflow for **CO2** on this example dataset.

1. Run the main entry point:

   - Open `src/FREQCOR_Start.py` and run it.

2. In case of path problems, check the *MANUAL MODE* section in FREQCOR_Start.py for .ini file path, and the file paths in the `examples/metadata/FREQCOR_config_example.ini`.

Outputs are written under `examples/output/` (or under the `output_path` configured in the INI file). For a description of the generated results, see the **Outputs** section.
