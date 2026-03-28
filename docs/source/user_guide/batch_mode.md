# Batch mode (CLI)

A single FREQCOR run (i.e., a single `.ini` configuration processed once) is designed to produce one
consistent set of outputs.

This means that within one run you cannot:

- process several gases at once,
- mix spectral and cospectral workflows,
- process several sites,
- run two independent analyses on two distinct time periods.

For time windows, you may exclude several sub-periods within one run (exclusion windows), but you cannot
define multiple independent processing windows within a single run.

Batch mode addresses this by allowing you to launch many runs conveniently (multiple sites and/or multiple
`.ini` files), while keeping each run internally consistent. It also makes it easier to use the available CPU
resources by running processing repeatedly without manual intervention.

## Batch mode: multi-site runs

Batch mode processes sites discovered under an input directory, where each site has its own folder.
Each site is expected to have a corresponding folder of `.ini` files.

### Folder structure

A typical layout is:

```text
inputs/
  BE-Lon/
    ... EddyPro outputs ...
  BE-Dor/
    ... EddyPro outputs ...

ini/
  BE-Lon/
    run_co2_cosp.ini
    run_h2o_cosp.ini
  BE-Dor/
    run_co2_cosp.ini

results/
```

### List available sites

```bash
python FREQCOR_Start.py --list-sites
```

### Run a subset of sites

```bash
python FREQCOR_Start.py --sites BE-Lon BE-Dor
```

### Override base folders

```bash
python FREQCOR_Start.py --input-dir inputs --ini-dir ini --results-dir results
```

## Flat mode: run many INI files from one directory

Flat mode processes every `.ini` file in a single directory, without any site-name pattern or subfolder
structure.

```bash
python FREQCOR_Start.py --flat-ini-dir path/to/ini_files --results-dir path/to/results
```

If `--results-dir` is omitted in flat mode, the summary CSV is written to the `--flat-ini-dir` directory.

## Outputs

Each `.ini` file produces its own outputs in the folder configured by `output_path`.

Additionally, `FREQCOR_Start.py` writes a summary CSV in the selected results directory to help track
success/failures across runs.
