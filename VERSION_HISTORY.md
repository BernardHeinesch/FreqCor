# Version history

## FREQCOR v1.4

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

## FREQCOR v1.3 (first public GitHub release)

Version 1.3 is the first version shared publicly on GitHub.

No version history is available for earlier internal development versions.

v1.3 was used in the study “A multi-site comparison of spectral and co-spectral approaches for correction of turbulent gas fluxes with ICOS set-up” (to appear soon in the review process in Atmospheric Measurement Techniques).
