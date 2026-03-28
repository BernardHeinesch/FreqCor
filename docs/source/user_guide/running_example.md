# Running the Example

FREQCOR is shipped with an example dataset from the ICOS Lonzée station (BE-Lon). The repository includes
EddyPro-generated inputs covering two spring months of 2024 over a winter wheat crop.

The provided example configuration file:

`examples/metadata/FREQCOR_config_example.ini`

is set up to run a **cospectral** workflow for **CO2** on this example dataset.

1. Open the example configuration file and verify that paths still point to your local repository.

2. Run the main entry point:

   - Open `src/FREQCOR_Start.py`.
   - In the *MANUAL MODE* section, set `ini_file` to `examples/metadata/FREQCOR_config_example.ini`.
   - Run `FREQCOR_Start.py` from the `src` directory.

Outputs are written under `examples/output/` (or under the `output_path` configured in the INI file). For a description of the generated results, see the **Outputs** section.
