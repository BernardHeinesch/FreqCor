# Version history

## FREQCOR v1.5 (2026-11-11)

### Configuration (flexible data sources)
- File locations are now configured using **full glob patterns** (path + wildcard + extension) directly in the INI, rather than splitting directory + filename.

### Optional ICOS-ETC / Vitale global flags
- More flexibility provided for the Vitale-flags input file (filename was hard-coded).

### Optional RH data source
- Instead of relying on relative humidity from the EddyPro full_output file, the user can now provide their own relative humidity data file using `PROCEDURE_OPTIONS.rh_from_file` which controls whether RH is read from `IO.rh_path` (alternative CSV) instead of `IO.flx_meteo`.

### INI configurator tool
- New standalone graphical tool `tools/FREQCOR_ini_gui.py` to fill FREQCOR configuration files: fields grouped in tabs, cascading enable/disable of dependent options, syntax-colored INI preview with hover-to-highlight of the key being edited, grouped parameter help (`?` buttons) with synthetic explanatory schematics, and validation through `FREQCOR_validate.validate_config` on save.
- The tool writes a commented INI that preserves the template structure; it starts blank and requires the user to load a starting INI.

### New `[PLOTS]` INI section
- Plot options (`plot_hh`, `plot_main`, `plot_aux`, `plot_save`) were moved out of `[PROCEDURE_OPTIONS]` into a dedicated `[PLOTS]` section at the end of the configuration file.
- `FREQCOR_Main` and `FREQCOR_validate` read the plot options exclusively from `[PLOTS]`; older INIs keeping them in `[PROCEDURE_OPTIONS]` will fail and must be updated.

### Read routines
- The former `FREQCOR_Read_IRGA` and `FREQCOR_Read_TOF` routines were merged into a unified `FREQCOR_Read_GEddySoft` reader for GEddySoft daily HDF5 outputs.
- The source is selected via `IO.source_type` (`irga` reads the `IRGA` group, `tof` reads the `TRACER` group); in `tof` mode the channel is selected by its `mz` value closest to the new `IO.target_mz` setting (tolerance 0.01).
- The active reader is chosen in the INI via `IO.read_routine` (`FREQCOR_Read_EP` for EddyPro output, `FREQCOR_Read_GEddySoft` for GEddySoft HDF5).
- In `tof` mode `PROCEDURE_OPTIONS.gss` is ignored: the non-H2O branch is forced and run outputs are tagged with the `target_mz` value.

### Removed dead parameter
- `USER_LIMITS.nulim` was removed from the example INI (and the configurator): it has not been read by the code since v1.3 (it was only passed through `FREQCOR_LUT_CF` in v1.2).

### Example INI comment clarifications
- Several INI comments were reworded for clarity and to document applicability conditions and recommended values.

## FREQCOR v1.4 (2026-04-11)

### Uncertainties (transfer function vs variability)
- Some uncertainty variables were renamed for clarity:
  - `unc_L` -> `unc_L_tf`
  - `unc_G` -> `unc_G_tf`
- The component uncertainties are now also written explicitly to outputs:
  - `unc_L_tf`, `unc_G_tf` (transfer-function-related component)
  - `unc_L_sd`, `unc_G_sd` (within-class CF spread, standard deviation after outlier removal)
- Total uncertainties are no longer written to outputs; users can recompute them as needed from the component uncertainties.
- The transfer-function uncertainty component is now computed uniformly as the maximum deviation from the central CF ("M") value.
- Output CSV headers were updated accordingly (notably `5_LUT_CF__*.csv`, as well as LUT cof uncertainty columns).

### Plot titles (clarity + consistency)
- Plot titles now include stability + spectra/cospectra where relevant, and are more consistent across steps (3/4/5).
- Peltola runs are annotated (`*` on Lorentz fco + brief note on TF interpretation).

### CSV headers (units)
- Added units to `4_LUT_cof__*.csv` and `4_mean_cosp__*.csv` headers (incl. `Nat. freq. (Hz)` and `(-)`), and fixed mojibake characters in outputs (e.g. `Âµmol` -> `µmol`, `Â°` -> `°`).

### Peltola interpretation note
- Added a concise note in the example INI next to `tf_peltola` clarifying the half-power vs half-amplitude interpretation.

### LUT CF (wind-speed statistics)
- In `FREQCOR_LUT_CF`, `ws` statistics written to the CF LUT (mean/max/std per WS class) are now computed from the points retained after CF outlier filtering (`remove_outliers`), rather than from all data in the class.

### LUT CF (sample counts)
- The number of cospectra retained per wind-speed class is now tracked as `ws_n` (N per WS class) alongside the CF LUT structure.

### Reference cospectra module rename
- Reference cospectrum functions used for plotting (Kaimal / Massman) have been consolidated under `FREQCOR_Ref_cospectrum_for_plotting`.
- Documentation has been updated accordingly.

### Plot filename (reference cospectra comparison)
- The saved comparison plot is now named `5_comp_reference_cospectra__*.png` (previously `5_av_kaimal_massman__*.png`).

---

## FREQCOR v1.3 (2026-03-28, first public GitHub release)

Version 1.3 is the first version shared publicly on GitHub.

No version history is available for earlier internal development versions.

v1.3 was used in the study “Spectral versus co-spectral correction of eddy-covariance fluxes: a multi-site assessment across the ICOS network”, published in Atmospheric Measurement Techniques ([doi.org/10.5194/amt-19-6037-2026](https://doi.org/10.5194/amt-19-6037-2026)).
