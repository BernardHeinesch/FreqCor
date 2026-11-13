# Installation

## Requirements

FREQCOR runs on Python 3 and requires a standard scientific Python stack.

- Python 3.9+ (recommended)
- Dependencies listed in `requirements.txt`

## Install

You can install FREQCOR dependencies either with **conda** (recommended on Windows) or with a
standard Python virtual environment.

### Option A: conda (recommended)

```bash
conda create -n freqcor python=3.11
conda activate freqcor
pip install -r requirements.txt
```

### Option B: venv

```bash
python -m venv .venv
# Windows PowerShell:
.venv\Scripts\Activate.ps1
# macOS/Linux:
# source .venv/bin/activate
pip install -r requirements.txt
```

## Quick sanity check

After installation, you can perform the example runs shipped with the repository. Two example runs are available on the same dataset (ICOS Lonzée, BE-Lon), processed from two different input formats:

- `examples/metadata/FREQCOR_config_example_EddyPro.ini` — EddyPro inputs
- `examples/metadata/FREQCOR_config_example_GEddySoft.ini` — GEddySoft inputs

Each can be run independently; together they let you compare the two input formats on identical data.

The main entry point for batch/CLI runs is `src/FREQCOR_Start.py`.
