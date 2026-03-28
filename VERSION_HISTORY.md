# Version history

## FREQCOR v1.3 (local modifications)

### Output naming
- Output and plot filenames now use double underscores as separators and include a run tag in the form:
  - `site__gas__method__window__run_dt`
- Dates in filename tags now use 2-digit years:
  - window: `yymmdd-yymmdd`
  - run datetime: `yymmddTHHMM`
- Step-number output prefixes no longer use leading zero padding (e.g. `07_...` -> `7_...`).
- Step-numbered outputs and per-run stats:
  - `0_run_input__... .ini`
  - `4_LUT_cof__... .csv`
  - `5_LUT_CF__... .csv`
  - `6_flux_out__... .csv`
  - `7_stats__... .txt`
- CF-related plot filenames now use capital `CF` consistently:
  - `5_CF_vs_ws__... .png`

### Plot outputs
- Plot saving now requires an explicit `file_tag` for saved figures.
- Unified plots updated to use run-aware filenames (including class index zero-padding where used).
- Combined unstable/stable Kaimal+Massman averaged cospectra into a single vertically stacked figure:
  - `05_av_kaimal_massman__{run_tag}.png`
- Renamed mean all-classes cospectra plot output:
  - `04_mean_cospectra_all_classes__...` -> `04_mean_cosp_all_classes__...`

### Read routine selection
- Added INI option `read_routine` under `[FILES_PATHS]` to select the input reader:
  - `FREQCOR_Read`
  - `FREQCOR_Read_TOF`
- `FREQCOR_Main` dispatches to the selected read routine.
- Configuration validation updated to validate `read_routine`.

### Warning fixes / forward-compatibility
- Fixed divide-by-zero warning in IQR outlier filtering (`FREQCOR_Sel_cof`) using `np.divide(..., where=denom>0)`.
- Removed pandas `DataFrame.swapaxes` deprecation warnings by avoiding passing DataFrames to `np.array_split` (split indices, then subset).
- Fixed pandas positional indexing FutureWarning in TF plotting by using `.iloc` for Series positional access.
- Fixed pandas `groupby(observed=...)` FutureWarning by setting `observed=False` explicitly.

### CSV writing compatibility
- Fixed pandas `to_csv` argument incompatibility by using `lineterminator=` instead of `line_terminator=`.

### Signature cleanup (remove passed-but-unused arguments)
- Removed `zmeas` and `disph` from `FREQCOR_cof` signature and updated its call sites.
- Removed unused `zmeas`/`disph` forwarding through plotting/LUT_cof plumbing:
  - `plot_cosp_unified` signature simplified.
  - `FREQCOR_LUT_cof` signature simplified.
  - `FREQCOR_Compute` signature updated accordingly.
- Removed unused arguments from `FREQCOR_Sel_general` (`nspec`, `nfreq`, `jsite`) and updated its call site.
