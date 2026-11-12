# -*- coding: utf-8 -*-
"""
FREQCOR INI configurator.

Standalone graphical wizard to fill a FREQCOR configuration (.ini) file.
The user loads a starting .ini file (e.g.
``examples/metadata/FREQCOR_config_example_EddyPro.ini``), edits the fields
through forms grouped in tabs, shows each INI comment as help (with
illustrative figures: synthetic schematics drawn with matplotlib, or PNGs
stored in ``tools/figures`` when available),
validates the result with ``FREQCOR_validate.validate_config`` and
writes a commented .ini file.

Usage:
    python tools\\FREQCOR_ini_gui.py

This script is self-contained: it does not modify any other source file.
Optional (used if present, silently skipped otherwise):
    - src/FREQCOR_validate.py  -> validates the configuration on save
    - tools/figures/*.png      -> illustrative figures in the help popups
"""

import configparser
import glob
import os
import re
import sys
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

# --------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------
GUI_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(GUI_DIR)          # FreqCor_v1.5
SRC_DIR = os.path.join(PROJECT_ROOT, 'src')
FIG_DIR = os.path.join(GUI_DIR, 'figures')

if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)
try:
    from FREQCOR_validate import validate_config
except Exception:
    validate_config = None


# --------------------------------------------------------------------------
# Synthetic schematic figures for the help popups
# --------------------------------------------------------------------------
def _schem_coslim(ax):
    """cospec_abs_limit: reject a half-hour if any point of the
    normalised (co)spectral density leaves the +/- limit band."""
    import numpy as np
    rng = np.random.default_rng(4)
    f = np.logspace(-3, 1.3, 300)
    fm = 0.1
    # f-weighted cospectrum: low-f plateau, inertial-range f^(-4/3) decay
    y = (135.0 / (1.0 + (f / fm) ** (4.0 / 3.0))
         * np.exp(rng.normal(0, 0.35, f.size)))
    lim = 250.0
    ax.loglog(f, y, 'k-', lw=1.4, label='f Co(f) / flux')
    over = np.where(y > lim)[0]          # actual points above the limit
    ax.loglog(f[over], y[over], 'o', ms=6, c='tab:red')
    ax.axhline(lim, ls='--', c='tab:red', lw=1.3)
    ax.text(8.0, lim * 0.78, 'limit', fontsize=8, color='tab:red',
            ha='center')
    if over.size:
        ax.annotate('any point outside +/-limit -> half-hour rejected',
                    xy=(f[over[0]], y[over[0]]), xytext=(0.006, 430),
                    fontsize=8, color='tab:red', ha='left', va='top',
                    arrowprops=dict(arrowstyle='->', lw=0.8))
    ax.set_xlabel('frequency [Hz]')
    ax.set_ylabel('f Co(f) / flux')
    ax.set_ylim(1, 500)
    ax.set_title('Absolute (co)spectrum limit (cospec_abs_limit)',
                 fontsize=9)
    ax.grid(True, which='both', alpha=0.3)


def _schem_tf(ax):
    import numpy as np
    rng = np.random.default_rng(5)
    f = np.logspace(-2, 1.3, 300)
    tau = 0.4
    tf2 = 1.0 / (1.0 + (2 * np.pi * f * tau) ** 2)
    fb = np.logspace(-2, 1.2, 24)
    pts = np.clip(tf2[np.searchsorted(f, fb)]
                  + rng.normal(0, 0.04, fb.size), 0, 1.15)
    ax.semilogx(fb, pts, 'o', ms=5, c='0.55', label='experimental TF2')
    ax.semilogx(f, tf2, 'k-', lw=1.8, label='fitted TF2')
    ax.axvspan(0.04, 2.0, color='tab:green', alpha=0.15)
    ax.axvspan(0.10, 1.5, color='tab:orange', alpha=0.20)
    ax.text(0.05, 0.95, 'fit band\n[tfmin, tfmax]', fontsize=8,
            color='tab:green', va='top')
    ax.text(0.8, 0.08, 'QC band [jtmin, jtmax]', fontsize=8,
            color='tab:orange', ha='center')
    cof = 1.0 / (2 * np.pi * tau)
    ax.plot([cof], [0.5], 'o', c='tab:red')
    ax.annotate('cof (half-power)', xy=(cof, 0.5), xytext=(0.8, 0.65),
                fontsize=8, color='tab:red',
                arrowprops=dict(arrowstyle='->', lw=0.8))
    ax.axhline(0.5, ls=':', c='tab:red', lw=0.8)
    ax.set_xlim(1e-2, 20)
    ax.set_xlabel('frequency [Hz]')
    ax.set_ylabel('TF2')
    ax.set_ylim(0, 1.05)
    ax.legend(fontsize=7.5, loc='lower left')
    ax.tick_params(labelsize=8)
    ax.set_title('Transfer function fit', fontsize=9)
    ax.grid(True, which='major', alpha=0.4)
    ax.grid(True, which='minor', alpha=0.15)


_schem_tf.figsize = (5.2, 5.6)   # larger: avoids a compressed look


def _schem_lut(ax):
    """Two-level class sorting used to build the cof / CF look-up tables."""
    ax.axis('off')

    def box(x, y, s, fs=8, fc='#eaf1fb'):
        ax.text(x, y, s, ha='center', va='center', fontsize=fs,
                bbox=dict(boxstyle='round,pad=0.35', fc=fc, ec='#4c72b0'))

    def link(x0, y0, x1, y1):
        ax.annotate('', xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle='-', lw=0.9, color='0.5'))

    box(0.5, 0.93, 'all selected half-hours', 9)
    for i, x in enumerate((0.30, 0.50, 0.70)):
        link(0.5, 0.88, x, 0.75)
        box(x, 0.68, f'wd {i + 1}' + (' ...' if i == 2 else ''), 7)
    ax.text(0.5, 0.80, '1st level : wd (non-H2O)\n'
            'or rh (H2O)', ha='center', fontsize=7, color='0.4')
    for i, x in enumerate((0.18, 0.28, 0.38, 0.48)):
        link(0.30, 0.62, x, 0.52)
        box(x, 0.45, f'ws{i + 1}', 7)
    ax.text(0.72, 0.45, '2nd level :\n'
            'wind-speed\nclasses\n'
            '(cof / CF)',
            fontsize=7, color='0.4', va='center')
    ax.text(0.5, 0.25, 'CF LUTs are built twice :\n'
            'unstable (z/L < 0) -> classnumcf_u\n'
            'stable (z/L > 0) -> classnumcf_s',
            fontsize=7, ha='center')
    ax.text(0.5, 0.06,
            'if the 1st-level class number is 1,\n'
            'only the wind-speed split applies',
            fontsize=7, ha='center', style='italic', color='0.4')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_title('Two-level LUT sorting', fontsize=9)


def _schem_cof_sel(ax):
    import numpy as np
    rng = np.random.default_rng(3)
    h = np.concatenate([rng.normal(-90, 35, 200), rng.normal(90, 35, 200),
                        rng.normal(0, 12, 60)])
    lo, hi = 15, 300
    ax.hist(h, bins=40, color='tab:blue', alpha=0.6)
    ax.axvspan(-hi, -lo, color='tab:green', alpha=0.15)
    ax.axvspan(lo, hi, color='tab:green', alpha=0.15)
    for x in (-hi, -lo, lo, hi):
        ax.axvline(x, c='tab:red', ls='--', lw=1)
    ax.set_xlim(-380, 380)
    ax.text(0.5, 0.96, 'kept where  min <= |flux| <= max',
            transform=ax.transAxes, color='tab:green', fontsize=8,
            ha='center', va='top')
    for x, lab, yy in ((-hi, '-max', 0.80), (-lo, '-min', 0.80),
                       (lo, 'min', 0.62), (hi, 'max', 0.80)):
        ax.annotate(lab, xy=(x, yy), xycoords=('data', 'axes fraction'),
                    fontsize=8, color='tab:red', ha='center',
                    xytext=(0, 8), textcoords='offset points')
    ax.set_xlabel('flux (H, Fc or LE), selection on |flux|')
    ax.set_ylabel('count of half-hours')
    ax.set_title('Flux-based (co)spectrum selection (absolute values)',
                 fontsize=9)
    ax.grid(True, alpha=0.3)


_schem_cof_sel.figsize = (5.2, 5.6)   # taller: avoids a squashed look


def _schem_massman(ax):
    """Model cospectrum fit with the 3 Massman coefficients a, b, c."""
    import numpy as np
    rng = np.random.default_rng(4)
    a, b, c = 1.0, 0.08, 0.6        # norm. factor, peak freq., broadness
    model = lambda f: (a * (f / b) /
                       (1.0 + (f / b) ** (2 * c)) ** (1.1667 / c))
    fb = np.logspace(-3, 0.8, 40)
    pts = np.clip(model(fb) + rng.normal(0, 0.06, fb.size), 0, None)
    f = np.logspace(-3, 1.2, 300)
    ax.semilogx(fb, pts, 'o', ms=4, c='0.55', label='binned cospectra')
    ax.semilogx(f, model(f), 'tab:blue', lw=1.8, label='Massman model')
    # maximum of x/(1+x^(2c))^(1.1667/c) is at x = (1/1.3334)^(1/(2c))
    f_peak = b * (1.0 / (2.0 * 1.1667 - 1.0)) ** (1.0 / (2 * c))
    ax.plot([f_peak], [model(f_peak)], 'v', c='tab:red', ms=8)
    ax.text(1.5e-3, 0.68,
            r"$\frac{f\,Co(f)}{\overline{w'c'}} ="
            r"\,\frac{a\,(f/b)}{[1+(f/b)^{2c}]^{1.1667/c}}$",
            fontsize=11)
    ax.text(1.5e-3, 0.44,
            'a : normalisation\nb : peak frequency\nc : broadness',
            fontsize=8)
    ax.annotate('b', xy=(f_peak, model(f_peak)), xytext=(0.4, 1.02),
                fontsize=10, color='tab:red', ha='center',
                arrowprops=dict(arrowstyle='->', lw=0.8))
    ax.text(0.02, -0.28, 'Massman & Clement, 2004: "Uncertainty in eddy\n'
            'covariance flux estimates resulting from spectral\n'
            'attenuation" (Eq. 4.2). In: Handbook of Micrometeorology\n'
            '(Lee, Massman & Law eds), Kluwer, ch. 4, pp. 67-99.',
            transform=ax.transAxes, fontsize=7, style='italic',
            color='0.35', va='top')
    ax.set_xlabel('frequency [Hz]')
    ax.set_ylabel('normalised cospectrum')
    ax.set_ylim(0, 1.15)
    ax.legend(fontsize=8, loc='lower right')
    ax.set_title('Massman coefficient fit (massman_path)', fontsize=9)
    ax.grid(True, which='both', alpha=0.3)


def _schem_peltola(ax):
    """Peltola TF: Lorentz vs sqrt(Lorentz) fits, cof at half-power."""
    import numpy as np
    f = np.logspace(-3, 1.3, 300)
    cof = 0.35
    lorentz = 1.0 / (1.0 + (f / cof) ** 2)
    peltola = np.sqrt(lorentz)
    rng = np.random.default_rng(7)
    noise = rng.normal(0, 0.03, 60)
    fd = np.logspace(-3, 0.6, 60)
    data = np.clip(1.0 / (1.0 + (fd / cof) ** 2) + noise, 0, 1.05)
    ax.semilogx(fd, data, 'o', ms=3, c='0.6', label='experimental TF2')
    ax.semilogx(f, lorentz, 'tab:blue', lw=1.6, label='Lorentz fit')
    ax.semilogx(f, peltola, 'tab:red', lw=1.6, ls='--',
                label='sqrt(Lorentz) - Peltola')
    ax.axhline(0.5, ls=':', c='0.4', lw=0.8)
    ax.annotate('cof = TF2 = 0.5', xy=(cof, 0.5), xytext=(1.0, 0.72),
                fontsize=8, color='tab:blue',
                arrowprops=dict(arrowstyle='->', lw=0.8))
    ax.plot([cof], [0.5], 'o', c='tab:blue')
    ax.set_xlabel('frequency [Hz]')
    ax.set_ylabel('TF')
    ax.set_ylim(0, 1.1)
    ax.legend(fontsize=7, loc='lower left')
    ax.set_title('Peltola et al. 2021 TF fit (tf_peltola)', fontsize=9)
    ax.grid(True, which='both', alpha=0.3)


def _schem_flags(ax):
    import numpy as np
    rng = np.random.default_rng(1)
    t = np.arange(60)
    flags = np.zeros(60, dtype=int)
    flags[rng.choice(60, 12, replace=False)] = 1
    flags[rng.choice(60, 4, replace=False)] = 2
    ok = flags == 0
    ax.scatter(t[ok], flags[ok], c='tab:green', s=18, label='kept')
    ax.scatter(t[~ok], flags[~ok], c='tab:red', s=18, label='rejected')
    ax.axhline(0.5, ls='--', c='tab:red')
    ax.text(1, 0.62, 'threshold (fvm / ff)', fontsize=8, color='tab:red')
    ax.set_xlabel('half-hour index')
    ax.set_ylabel('QC flag value')
    ax.set_yticks([0, 1, 2])
    ax.set_ylim(-0.4, 2.6)
    ax.legend(fontsize=7, loc='upper right')
    ax.set_title('QC flag filtering over the run', fontsize=9)
    ax.grid(True, alpha=0.3)


def _schem_sector(ax):
    import numpy as np
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    wdmin, wdmax = np.deg2rad([200, 250])
    ax.fill_between(np.linspace(wdmin, wdmax, 60), 0, 1,
                    color='tab:red', alpha=0.5)
    ax.annotate('excluded sector\n[wdmin, wdmax]', xy=(np.deg2rad(225), 0.7),
                fontsize=8, ha='center', color='black')
    ax.set_rticks([])
    ax.set_title('Wind-direction exclusion', fontsize=9)


def _schem_timewin(ax):
    import numpy as np
    import matplotlib.dates as mdates
    t0 = np.datetime64('2024-01-01')
    t1 = np.datetime64('2024-12-31')
    # top track: full data + processing window
    ax.broken_barh([(t0, (t1 - t0).astype(int))], (0.4, 0.8),
                   facecolors='lightgrey')
    ax.broken_barh([(np.datetime64('2024-04-01'), 150)], (0.4, 0.8),
                   facecolors='tab:green', alpha=0.6)
    # bottom track: excluded sub-periods (same date axis)
    ax.broken_barh([(np.datetime64('2024-04-20'), 15),
                    (np.datetime64('2024-07-15'), 20)], (-0.5, 0.7),
                   facecolors='tab:red', alpha=0.75)
    ax.text(np.datetime64('2024-05-15'), 1.5,
            'processing window\n[start, end]',
            fontsize=8, color='tab:green', ha='center')
    ax.text(np.datetime64('2024-07-25'), -0.95,
            'excluded\nsub-periods',
            fontsize=8, color='tab:red', ha='center')
    ax.set_xlim(t0, t1)
    ax.set_ylim(-1.5, 2.1)
    ax.set_yticks([])
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %y'))
    for lab in ax.get_xticklabels():
        lab.set_rotation(30)
        lab.set_fontsize(7)
        lab.set_ha('right')
    ax.set_xlabel('date')
    ax.set_title('Time window and exclusions', fontsize=9)
    ax.grid(True, axis='x', alpha=0.3)


SCHEMATICS = {
    'coslim': _schem_coslim,
    'tf': _schem_tf,
    'lut': _schem_lut,
    'massman': _schem_massman,
    'peltola': _schem_peltola,
    'cof_sel': _schem_cof_sel,
    'flags': _schem_flags,
    'sector': _schem_sector,
    'timewin': _schem_timewin,
}

# --------------------------------------------------------------------------
# Field specification
#   (section, key, label, kind, choices/browse, help, fig_keyword, tab)
#   kind : 'entry' | 'choice' | 'check' | 'glob' | 'dir' | 'file'
# --------------------------------------------------------------------------
E = 'Essentials'
P = 'Paths & options'
A = 'Advanced procedure'
L = 'Limits & QC'
S = 'Site & cache'
G = 'Plots'

FIELDS = [
    # ---- IO / essentials --------------------------------------------------
    ('IO', 'read_routine', 'Read routine', 'choice',
     ['FREQCOR_Read_EP', 'FREQCOR_Read_GEddySoft'],
     "FREQCOR_Read_EP : standard EddyPro run output.\n"
     "FREQCOR_Read_GEddySoft : GEddySoft HDF5 files (IRGA or TRACER/TOF).",
     None, E),
    ('IO', 'source_type', "Source GEddySoft", 'choice', ['irga', 'tof'],
     "irga : 'IRGA' group of the HDF5 files, channel selected by name/gss.\n"
     "tof  : 'TRACER' group (PTR-TOF-MS), channel selected by target_mz.\n"
     "Ignored when read_routine = FREQCOR_Read_EP.",
     None, E),
    ('IO', 'target_mz', "Target m/z (TOF)", 'entry', None,
     "Tracer target mass (tolerance 0.01), REQUIRED when source_type = tof,\n"
     "ignored otherwise. Examples: 33.033 (methanol), 69.0699 (isoprene),\n"
     "137.132 (monoterpenes).",
     None, E),
    ('IO', 'binned_cosp', 'Binned cospectra (glob)', 'glob', None,
     "Full pattern of the input files.\n"
     "  EP        : *binned_cospectra*_adv.csv\n"
     "  GEddySoft : *.hdf5",
     None, E),
    ('IO', 'output_path', 'Output folder', 'dir', None,
     "Folder where results and figures are written.",
     None, E),
    ('PROCEDURE_OPTIONS', 'sps', 'Procedure (sps)', 'choice', ['1', '2'],
     "1 = co-spectral approach (w-gas): recommended.\n"
     "2 = spectral approach (gas alone): requires the sensor separation\n"
     "    parameters of [EC_SETUP].",
     None, E),
    ('PROCEDURE_OPTIONS', 'gss', 'Corrected gas (gss)', 'choice',
     ['1', '2', '3', '4', '5'],
     "Variable corrected by FREQCOR:\n"
     "  1=CO2  2=H2O  3=O3  4=CH4  5=N2O.\n"
     "Unused when source_type = tof: the species is selected by\n"
     "target_mz and the non-H2O branch is forced.",
     None, E),

    # ---- IO / advanced paths ----------------------------------------------
    ('IO', 'flx_meteo', 'Flux & meteo CSV (glob)', 'glob', None,
     "EddyPro 'full_output' file containing fluxes, meteo and QC.\n"
     "Required for FREQCOR_Read_EP only.",
     None, P),
    ('PROCEDURE_OPTIONS', 'rh_from_file', 'RH from file', 'check', None,
     "1 : read air relative humidity from rh_path (cols TIMESTAMP_END,\n"
     "RH) instead of flx_meteo, aligned on TIMESTAMP_END.\n"
     "Useful because with EddyPro runs without external RH, RH is\n"
     "computed from Tsonic and might be biased.",
     None, P),
    ('IO', 'rh_path', 'Alternative RH CSV (glob)', 'glob', None,
     "Pattern of the alternative RH CSV (TIMESTAMP_END, RH columns).\n"
     "Used only when rh_from_file = 1.\n"
     "See ICOSETC_BE-Lon_METEO_L2.csv in examples/data for the format.",
     None, P),
    ('PROCEDURE_OPTIONS', 'vitale_qc_flags', 'Vitale flags', 'check', None,
     "1 : use the global ICOS-ETC flags (Vitale et al. 2020) read from\n"
     "vitale_path instead of the locally computed VM and F flags.",
     'schem:flags', P),
    ('IO', 'vitale_path', 'Vitale flags CSV (glob)', 'glob', None,
     "Pattern of the ICOS-ETC flags CSV. Used when vitale_qc_flags = 1.\n"
     "See ICOSETC_BE-Lon_FLUXES_L2.csv in examples/data for the format.",
     None, P),
    ('IO', 'massman_path', 'Massman coefficients', 'file', None,
     "Full path (including file name) to the file containing the\n"
     "Massman coefficients of modeled local cospectra.\n"
     "Optional. When available, used only for plotting the\n"
     "reference Massman cospectra.\n"
     "See Massman_coef.csv in examples/data for the format.",
     'schem:massman', P),
    ('IO', 'output_file', 'Base name (cache)', 'entry', None,
     "Base name of the intermediate outputs (PROCESSED_DATA).\n"
     "Does not affect the main outputs (run_tag).",
     None, P),

    # ---- TIME_WINDOW -------------------------------------------------------
    ('TIME_WINDOW', 'enable_time_window', 'Time window', 'check', None,
     "1 : only process files/half-hours between start and end.",
     'schem:timewin', E),
    ('TIME_WINDOW', 'start_datetime', 'Start (YYYY-MM-DD HH:MM)', 'entry', None,
     "Start of the processing window.",
     'schem:timewin', E),
    ('TIME_WINDOW', 'end_datetime', 'End (YYYY-MM-DD HH:MM)', 'entry', None,
     "End of the processing window.",
     'schem:timewin', E),
    ('TIME_WINDOW', 'enable_exclusion_windows', 'Exclusions', 'check', None,
     "1 : exclude the sub-periods listed below.",
     'schem:timewin', P),
    ('TIME_WINDOW', 'date_exclusion_start', 'Excl. start(s)', 'entry', None,
     "Exclusion start date(s), comma-separated.",
     'schem:timewin', P),
    ('TIME_WINDOW', 'date_exclusion_end', 'Excl. end(s)', 'entry', None,
     "Exclusion end date(s) (same count as the starts).",
     'schem:timewin', P),

    # ---- PROCEDURE_OPTIONS / advanced --------------------------------------
    ('PROCEDURE_OPTIONS', 'classnumwd', 'wd classes (cof+CF)', 'entry', None,
     "Number of classes of the 1st split (wind direction) for gases != H2O.\n"
     "1 = no wd split -> the 2nd level (ws) applies.",
     'schem:lut', A),
    ('PROCEDURE_OPTIONS', 'classnumrh', 'rh classes (cof+CF)', 'entry', None,
     "Number of classes of the 1st split (relative humidity) for H2O (gss=2).\n"
     "1 = no rh split -> the 2nd level (ws) applies.",
     'schem:lut', A),
    ('PROCEDURE_OPTIONS', 'classnumcofws', 'ws classes (cof)', 'entry', None,
     "Wind-speed classes for the cutoff-frequency LUTs, used only when\n"
     "the 1st-level field is 1; otherwise hard-coded to 4.",
     'schem:lut', A),
    ('PROCEDURE_OPTIONS', 'classnumcf_u', 'ws classes (CF unstable)', 'entry', None,
     "Wind-speed classes of the correction-factor LUT built on the\n"
     "unstable subset (z/L < 0). Used only when the 1st-level class\n"
     "number is 1; otherwise hard-coded to 4.",
     'schem:lut', A),
    ('PROCEDURE_OPTIONS', 'classnumcf_s', 'ws classes (CF stable)', 'entry', None,
     "Wind-speed classes of the correction-factor LUT built on the\n"
     "stable subset (z/L > 0). Used only when the 1st-level class\n"
     "number is 1; otherwise hard-coded to 4.",
     'schem:lut', A),
    ('PROCEDURE_OPTIONS', 'eq', 'HL equation (eq)', 'choice', ['13', '16'],
     "Horst & Lenschow (2009) equation for the sensor-separation\n"
     "correction: 13 or 16. Equation 16 is recommended.\n"
     "Used only when sps = 2.",
     None, A),
    ('PROCEDURE_OPTIONS', 'tf_sonic', 'Anemo losses (tf_sonic)', 'check', None,
     "1 : apply the theoretical correction (path averaging, time\n"
     "response) to reference cospectra before computing CFs\n"
     "(recommended).",
     None, A),
    ('PROCEDURE_OPTIONS', 'tf_peltola', 'TF Peltola 2021', 'check', None,
     "TF fit according to Peltola et al. 2021 (square-root of H).\n"
     "Considered only for the cospectral approach (sps = 1); the\n"
     "spectral approach always applies the square-root.\n"
     "Set to 1 if the timelag search through covariance maximisation\n"
     "was used before (co-)spectra computation.\n"
     "Applies to the Lorentz fit form only (Gaussian TF is unchanged).\n"
     "Warning: cof then corresponds to the half-power point\n"
     "(TF2 = 0.5), not the half-amplitude point of TF.",
     'schem:peltola', A),
    # ---- USER_LIMITS / general & QC filtering -------------------------------
    ('USER_LIMITS', 'wdmin', 'Excluded wind sector min', 'entry', None,
     "Lower bound of the excluded wind sector (degrees from North).\n"
     "0/0 = no exclusion.",
     'schem:sector', L),
    ('USER_LIMITS', 'wdmax', 'Excluded wind sector max', 'entry', None,
     "Upper bound of the excluded wind sector (degrees).",
     'schem:sector', L),
    ('USER_LIMITS', 'fvm', 'VM flag threshold (fvm)', 'entry', None,
     "Vickers & Mahrt flags (spikes, drop-outs, skewness, discontinuities).\n"
     "0: limit; any value > threshold -> half-hour rejected.",
     'schem:flags', L),
    ('USER_LIMITS', 'ff', 'Foken flag threshold (ff)', 'entry', None,
     "Foken flag: 0 (<30%), 1 (30-100%), 2 (>100%).\n"
     "Value > threshold -> half-hour rejected. 0 is strictest.",
     'schem:flags', L),
    ('USER_LIMITS', 'cospec_abs_limit', 'Abs (co)spectrum limit', 'entry', None,
     "Applied to the normalised (co)spectral density (Co(f)/cov):\n"
     "if any point leaves the band [-limit, +limit], the whole\n"
     "half-hour is rejected. 250 is a reasonable choice.",
     'schem:coslim', L),

    # ---- USER_LIMITS / cof selection (cutoff frequency) ---------------------
    ('USER_LIMITS', 'hl_cof_min', 'cof: |H| min [W/m2]', 'entry', None,
     "Cospectrum selection for cof computation: keep |H| >= min.\n"
     "Limits apply to the absolute value of the flux.",
     'schem:cof_sel', L),
    ('USER_LIMITS', 'hl_cof_max', 'cof: |H| max [W/m2]', 'entry', None,
     "Keep |H| <= max (absolute value).",
     'schem:cof_sel', L),
    ('USER_LIMITS', 'fcl_cof_min', 'cof: |Fc| min', 'entry', None,
     "Main tracer flux [umol/m2/s], absolute value (non-H2O).",
     'schem:cof_sel', L),
    ('USER_LIMITS', 'fcl_cof_max', 'cof: |Fc| max', 'entry', None,
     "Upper bound of the tracer flux for cof selection.",
     'schem:cof_sel', L),
    ('USER_LIMITS', 'lel_cof_min', 'cof: |LE| min [W/m2]', 'entry', None,
     "Latent heat flux, absolute value; used when gss = 2 (H2O).",
     'schem:cof_sel', L),
    ('USER_LIMITS', 'lel_cof_max', 'cof: |LE| max [W/m2]', 'entry', None,
     "Upper bound of the latent heat flux for cof selection.",
     'schem:cof_sel', L),

    # ---- USER_LIMITS / CF selection (correction factors) --------------------
    ('USER_LIMITS', 'hl_cf_u_min', 'CF unstable: H min', 'entry', None,
     "CF selection under unstable conditions: keep H > threshold [W/m2].",
     'schem:cof_sel', L),
    ('USER_LIMITS', 'hl_cf_s_max', 'CF stable: H max', 'entry', None,
     "CF selection under stable conditions: keep H < threshold [W/m2]\n"
     "(0 = neutral).",
     'schem:cof_sel', L),
    ('USER_LIMITS', 'hl_cf_abs_max', 'CF: |H| cap', 'entry', None,
     "Absolute cap on H for CF selection, unstable and stable.",
     'schem:cof_sel', L),

    # ---- USER_LIMITS / TF quality & fitting bands ---------------------------
    ('USER_LIMITS', 'jtmin', 'TF check f min [Hz]', 'entry', None,
     "Lower bound of the QC band of the experimental TF (coefficient of\n"
     "variation compared to varclim).",
     'schem:tf', L),
    ('USER_LIMITS', 'jtmax', 'TF check f max [Hz]', 'entry', None,
     "Upper bound of the TF quality-control band.",
     'schem:tf', L),
    ('USER_LIMITS', 'varclim', 'varclim', 'entry', None,
     "If the TF coefficient of variation over [jtmin, jtmax] exceeds\n"
     "varclim, the TF is rejected.",
     'schem:tf', L),
    ('USER_LIMITS', 'tfmin', 'TF fit f min [Hz]', 'entry', None,
     "Lower bound of the transfer-function fitting band.",
     'schem:tf', L),
    ('USER_LIMITS', 'tfmax', 'TF fit f max [Hz]', 'entry', None,
     "Upper bound of the transfer-function fitting band.",
     'schem:tf', L),

    # ---- EC_SETUP -----------------------------------------------------------
    ('EC_SETUP', 'jsite', 'Site name', 'entry', None,
     "Site identifier (used in output names).",
     None, S),
    ('EC_SETUP', 'zmeas', 'Measurement height z [m]', 'entry', None,
     "Required only when sps = 2 (HL separation formulas).",
     None, S),
    ('EC_SETUP', 'disph', 'Displacement height d [m]', 'entry', None,
     "Required only when sps = 2.",
     None, S),
    ('EC_SETUP', 'rsz', 'Vertical separation [m]', 'entry', None,
     "Vertical separation between the gas-sampling inlet and the\n"
     "centre of the anemometer measurement volume; sps = 2 only.",
     None, S),
    ('EC_SETUP', 'rseast', 'East separation [m]', 'entry', None,
     "Eastward separation between the gas-sampling inlet and the\n"
     "centre of the anemometer measurement volume; sps = 2 only.",
     None, S),
    ('EC_SETUP', 'rsnorth', 'North separation [m]', 'entry', None,
     "Northward separation between the gas-sampling inlet and the\n"
     "centre of the anemometer measurement volume; sps = 2 only.",
     None, S),

    # ---- PROCESSED_DATA -----------------------------------------------------
    ('PROCESSED_DATA', 'enable_saving', 'Save intermediates', 'check', None,
     "1 : save the reader outputs for later reuse.",
     None, S),
    ('PROCESSED_DATA', 'use_timestamp', 'Timestamp files', 'check', None,
     "1 : include a timestamp in intermediate file names.",
     None, S),
    ('PROCESSED_DATA', 'enable_loading', 'Load intermediates', 'check', None,
     "1 : reload previously saved data (fast read).",
     None, S),
    ('PROCESSED_DATA', 'file_prefix', 'File prefix', 'entry', None,
     "Prefix of the intermediate files (site or site_timestamp).",
     None, S),

    # ---- PLOTS ----------------------------------------------------------------
    ('PLOTS', 'plot_hh', 'Individual 1/2h plots', 'check', None,
     "1 : plot every half-hourly (co)spectrum (very verbose).",
     None, G),
    ('PLOTS', 'plot_main', 'Main plots', 'check', None,
     "1 : main figures (means, transfer functions, CF vs WS).",
     None, G),
    ('PLOTS', 'plot_aux', 'Diagnostic plots', 'check', None,
     "1 : auxiliary figures (filtering, detailed diagnostics).",
     None, G),
    ('PLOTS', 'plot_save', 'Save figures', 'check', None,
     "1 : PNGs in output_path ; 0 : interactive display.",
     None, G),
]


# Logical help groups: one '?' button per group, shown on its first key.
HELP_GROUPS = [
    {'title': 'GEddySoft source selection',
     'keys': ['source_type', 'target_mz'], 'fig': None},
    {'title': 'Processing & exclusion time windows',
     'keys': ['enable_time_window', 'start_datetime', 'end_datetime',
              'enable_exclusion_windows', 'date_exclusion_start',
              'date_exclusion_end'],
     'fig': 'schem:timewin'},
    {'title': 'RH input',
     'keys': ['rh_from_file', 'rh_path'], 'fig': None},
    {'title': 'Vitale QC flags',
     'keys': ['vitale_qc_flags', 'vitale_path'], 'fig': 'schem:flags'},
    {'title': 'Class numbers (LUT levels)',
     'keys': ['classnumwd', 'classnumrh', 'classnumcofws',
              'classnumcf_u', 'classnumcf_s'],
     'fig': 'schem:lut'},
    {'title': 'Plot options',
     'keys': ['plot_hh', 'plot_main', 'plot_aux', 'plot_save'],
     'fig': None},
    {'title': 'Wind-sector exclusion',
     'keys': ['wdmin', 'wdmax'], 'fig': 'schem:sector'},
    {'title': 'cof flux limits',
     'keys': ['hl_cof_min', 'hl_cof_max', 'fcl_cof_min', 'fcl_cof_max',
              'lel_cof_min', 'lel_cof_max'],
     'fig': 'schem:cof_sel'},
    {'title': 'CF flux limits',
     'keys': ['hl_cf_u_min', 'hl_cf_s_max', 'hl_cf_abs_max'],
     'fig': 'schem:cof_sel'},
    {'title': 'QC flag thresholds',
     'keys': ['fvm', 'ff'], 'fig': 'schem:flags'},
    {'title': 'TF frequency bands',
     'keys': ['jtmin', 'jtmax', 'varclim', 'tfmin', 'tfmax'],
     'fig': 'schem:tf'},
    {'title': 'EC sensor geometry',
     'keys': ['zmeas', 'disph', 'rsz', 'rseast', 'rsnorth'], 'fig': None},
    {'title': 'Intermediate (processed) data',
     'keys': ['enable_saving', 'use_timestamp', 'enable_loading',
              'file_prefix'], 'fig': None},
]

# Bold section headers inserted before these fields inside their tab.
ROW_HEADERS = {
    'wdmin': 'Half-hour QC & general filtering',
    'hl_cof_min': 'Filtering for cof establishment (cutoff frequency)',
    'hl_cf_u_min': 'Filtering for CF computation (correction factors)',
    'jtmin': 'Transfer-function QC and fitting bands',
}

# key -> group dict
KEY2GROUP = {k: g for g in HELP_GROUPS for k in g['keys']}
# key -> (help text, fig keyword)
FIELD_META = {key: (help_txt, fig)
              for (_, key, _, _, _, help_txt, fig, _) in FIELDS}


class IniGui(tk.Tk):
    """Main window: tabbed form + INI preview + load/save/validate."""

    def __init__(self):
        super().__init__()
        self.title('FREQCOR - INI configurator')
        self.geometry('980x720')

        self.vars = {}            # (section, key) -> tk.Variable
        self.widgets = {}         # (section, key) -> widget
        self._fig_cache = {}
        self._help_windows = []

        self.template_text = ''      # text of the loaded .ini (or empty)
        self._loaded_values = None   # values of the loaded .ini

        self._build_menu()
        self._build_body()
        self._update_preview()
        self._apply_dependencies()

    # ------------------------------------------------------------------ UI
    def _build_menu(self):
        top = ttk.Frame(self, padding=4)
        top.pack(side='top', fill='x')
        ttk.Button(top, text='Open .ini...',
                   command=self.load_ini).pack(side='left', padx=2)
        ttk.Button(top, text='Reset fields',
                   command=self._reset_values).pack(side='left', padx=2)
        ttk.Button(top, text='Validate',
                   command=self._validate_only).pack(side='left', padx=10)
        ttk.Button(top, text='Save .ini...',
                   command=self.save_ini).pack(side='right', padx=2)

    def _build_body(self):
        paned = ttk.PanedWindow(self, orient='horizontal')
        paned.pack(fill='both', expand=True, padx=4, pady=4)

        # left: notebook of field tabs
        style = ttk.Style(self)
        style.configure('TNotebook.Tab', padding=(12, 5),
                        font=('Segoe UI', 9, 'bold'))
        style.map('TNotebook.Tab',
                  background=[('selected', '#b7d1f2')],
                  foreground=[('selected', '#0b3d91')])

        tab_colors = {'Essentials': '#2e7d32',
                      'Paths & options': '#1565c0',
                      'Advanced procedure': '#e65100',
                      'Limits & QC': '#b71c1c',
                      'Site & cache': '#6a1b9a',
                      'Plots': '#00838f'}
        self._tab_dots = []

        nb = ttk.Notebook(paned)
        paned.add(nb, weight=3)

        tabs = {}
        for _, _, _, _, _, _, _, tabname in FIELDS:
            if tabname not in tabs:
                frame = ttk.Frame(nb)
                canvas = tk.Canvas(frame, highlightthickness=0)
                sb = ttk.Scrollbar(frame, orient='vertical',
                                   command=canvas.yview)
                inner = ttk.Frame(canvas)
                inner.bind('<Configure>',
                           lambda e, c=canvas: c.configure(
                               scrollregion=c.bbox('all')))
                canvas.create_window((0, 0), window=inner, anchor='nw')
                canvas.configure(yscrollcommand=sb.set)
                canvas.bind('<Enter>',
                            lambda e, c=canvas: c.bind_all(
                                '<MouseWheel>',
                                lambda ev: c.yview_scroll(
                                    int(-ev.delta / 120), 'units')))
                canvas.bind('<Leave>',
                            lambda e, c=canvas: c.unbind_all('<MouseWheel>'))
                canvas.pack(side='left', fill='both', expand=True)
                sb.pack(side='right', fill='y')
                tabs[tabname] = inner
                dot = tk.PhotoImage(width=10, height=10)
                dot.put(tab_colors.get(tabname, '#555555'),
                        (0, 0, 10, 10))
                self._tab_dots.append(dot)
                nb.add(frame, text=f'  {tabname}', image=dot,
                       compound='left')

        current_grp = None
        current_rows = None
        for (sec, key, label, kind, aux, help_txt, fig, tabname) in FIELDS:
            header = ROW_HEADERS.get(key)
            if header:
                ttk.Separator(tabs[tabname]).pack(fill='x', pady=(10, 2),
                                                  padx=4)
                ttk.Label(tabs[tabname], text=header,
                          font=('Segoe UI', 9, 'bold')).pack(anchor='w',
                                                             padx=6)
            parent = tabs[tabname]
            grp = KEY2GROUP.get(key)
            if grp is not None:
                if current_grp is not grp:
                    # container: rows on the left, bracket + '?' on the right
                    container = ttk.Frame(parent)
                    container.pack(fill='x', pady=(0, 8))
                    current_rows = ttk.Frame(container)
                    current_rows.pack(side='left', fill='x', expand=True)
                    strip = tk.Frame(container)
                    strip.pack(side='right', fill='y', padx=(0, 4))
                    bracket = tk.Canvas(strip, width=9, height=10,
                                      highlightthickness=0)
                    bracket.pack(side='left', fill='y')
                    bracket.bind('<Configure>',
                                 lambda e, c=bracket: self._draw_bracket(c))
                    ttk.Button(strip, text='?', width=3,
                               command=lambda g=grp:
                                   self._show_group_help(g)
                               ).pack(side='left', expand=True)
                    current_grp = grp
                parent = current_rows
            else:
                current_grp = None
            self._add_row(parent, sec, key, label, kind, aux,
                          help_txt, fig)

        # right: INI preview
        right = ttk.Frame(paned)
        paned.add(right, weight=2)
        ttk.Label(right, text='INI preview').pack(anchor='w')
        self.preview = tk.Text(right, wrap='none', width=50,
                               font=('Consolas', 9))
        # editor-like syntax highlighting
        self.preview.tag_configure('section', foreground='#0000cc',
                                   font=('Consolas', 9, 'bold'))
        self.preview.tag_configure('comment', foreground='#008000')
        self.preview.tag_configure('key', foreground='#7a3e9d')
        self.preview.tag_configure('number', foreground='#b35c00')
        self.preview.tag_configure('string', foreground='#a31515')
        self.preview.tag_configure('hover', background='#fff2a8')
        psb = ttk.Scrollbar(right, orient='vertical',
                            command=self.preview.yview)
        self.preview.configure(yscrollcommand=psb.set)
        self.preview.pack(side='left', fill='both', expand=True)
        psb.pack(side='right', fill='y')

    def _add_row(self, parent, sec, key, label, kind, aux, help_txt, fig):
        row = ttk.Frame(parent, padding=(4, 2))
        row.pack(fill='x')

        ttk.Label(row, text=label, width=26).pack(side='left')

        var = tk.StringVar()
        self.vars[(sec, key)] = var

        if kind == 'choice':
            w = ttk.Combobox(row, textvariable=var, values=aux,
                             state='readonly', width=30)
        elif kind == 'check':
            w = ttk.Checkbutton(row, variable=var, onvalue='1', offvalue='0')
        else:
            w = ttk.Entry(row, textvariable=var, width=34)

        w.pack(side='left', padx=2)
        self.widgets[(sec, key)] = w

        for widget in (row, w):
            widget.bind('<Enter>',
                        lambda e, k=key: self._hover_ini_line(k))
            widget.bind('<Leave>', lambda e: self._hover_ini_line(None))

        if kind == 'glob':
            ttk.Button(row, text='...', width=3,
                       command=lambda v=var: self._browse_glob(v)
                       ).pack(side='left')
        elif kind == 'dir':
            ttk.Button(row, text='...', width=3,
                       command=lambda v=var: self._browse_dir(v)
                       ).pack(side='left')
        elif kind == 'file':
            ttk.Button(row, text='...', width=3,
                       command=lambda v=var: self._browse_file(v)
                       ).pack(side='left')

        if key not in KEY2GROUP:
            ttk.Button(row, text='?', width=3,
                       command=lambda k=key, h=help_txt, f=fig:
                           self._show_help(k, h, f)).pack(side='left', padx=2)

        var.trace_add('write', lambda *_: self._on_change())

    def _draw_bracket(self, canvas):
        """Redraw the group bracket as ']' : bar near the '?' button and
        end caps pointing left toward the covered fields."""
        canvas.delete('all')
        h = canvas.winfo_height()
        w = canvas.winfo_width()
        x = w - 3
        canvas.create_line(x, 0, x, h, fill='#7f9fcf', width=3)
        canvas.create_line(0, 1, x, 1, fill='#7f9fcf', width=3)
        canvas.create_line(0, h - 1, x, h - 1, fill='#7f9fcf', width=3)

    def _hover_ini_line(self, key):
        """Highlight the ``key = ...`` line in the preview while hovering."""
        if not hasattr(self, 'preview'):
            return
        self.preview.tag_remove('hover', '1.0', 'end')
        line = getattr(self, '_key_line', {}).get(key) if key else None
        if line:
            self.preview.tag_add('hover', f'{line}.0', f'{line}.end')
            self.preview.see(f'{line}.0')

    # -------------------------------------------------------------- helpers
    def _browse_glob(self, var):
        path = filedialog.askopenfilename(title='Pick a template file')
        if path:
            path = os.path.normpath(path)
            base, ext = os.path.splitext(path)
            var.set(f"{os.path.dirname(path)}{os.sep}*{ext}")

    def _browse_dir(self, var):
        path = filedialog.askdirectory(title='Pick a folder')
        if path:
            var.set(os.path.normpath(path))

    def _browse_file(self, var):
        path = filedialog.askopenfilename(title='Pick a file')
        if path:
            var.set(os.path.normpath(path))

    def _on_change(self):
        self._apply_dependencies()
        self._update_preview()

    def _apply_dependencies(self):
        """Enable/disable child options according to parent choices."""
        v = self.vars

        def val(sec, key):
            return v[(sec, key)].get().strip()

        rr = val('IO', 'read_routine')
        is_ep = rr == 'FREQCOR_Read_EP'
        is_ges = rr == 'FREQCOR_Read_GEddySoft'
        is_tof = val('IO', 'source_type') == 'tof'
        tof_ges = is_ges and is_tof          # gss forced to non-H2O branch
        is_h2o = (val('PROCEDURE_OPTIONS', 'gss') == '2'
                  and not tof_ges)
        sps2 = val('PROCEDURE_OPTIONS', 'sps') == '2'
        tw = val('TIME_WINDOW', 'enable_time_window') == '1'
        excl = val('TIME_WINDOW', 'enable_exclusion_windows') == '1'
        rh_ff = val('PROCEDURE_OPTIONS', 'rh_from_file') == '1'
        vit = val('PROCEDURE_OPTIONS', 'vitale_qc_flags') == '1'
        lvl1 = (val('PROCEDURE_OPTIONS', 'classnumrh') if is_h2o
                else val('PROCEDURE_OPTIONS', 'classnumwd'))

        rules = [
            # reader-driven
            ('IO', 'source_type', is_ges),
            ('IO', 'target_mz', is_ges and is_tof),
            ('PROCEDURE_OPTIONS', 'gss', not tof_ges),
            ('IO', 'flx_meteo', is_ep),
            ('PROCEDURE_OPTIONS', 'rh_from_file', is_ep),
            ('IO', 'rh_path', is_ep and rh_ff),
            ('PROCEDURE_OPTIONS', 'vitale_qc_flags', is_ep),
            ('IO', 'vitale_path', is_ep and vit),
            # time window
            ('TIME_WINDOW', 'start_datetime', tw),
            ('TIME_WINDOW', 'end_datetime', tw),
            ('TIME_WINDOW', 'date_exclusion_start', excl),
            ('TIME_WINDOW', 'date_exclusion_end', excl),
            # gas-driven class fields and cof flux bounds
            ('PROCEDURE_OPTIONS', 'classnumwd', not is_h2o),
            ('PROCEDURE_OPTIONS', 'classnumrh', is_h2o),
            ('PROCEDURE_OPTIONS', 'classnumcofws', lvl1 in ('', '1')),
            ('USER_LIMITS', 'fcl_cof_min', not is_h2o),
            ('USER_LIMITS', 'fcl_cof_max', not is_h2o),
            ('USER_LIMITS', 'lel_cof_min', is_h2o),
            ('USER_LIMITS', 'lel_cof_max', is_h2o),
            # spectral-procedure-only fields
            ('PROCEDURE_OPTIONS', 'eq', sps2),
            ('EC_SETUP', 'zmeas', sps2),
            ('EC_SETUP', 'disph', sps2),
            ('EC_SETUP', 'rsz', sps2),
            ('EC_SETUP', 'rseast', sps2),
            ('EC_SETUP', 'rsnorth', sps2),
        ]
        for sec, key, enabled in rules:
            w = self.widgets.get((sec, key))
            if w is not None:
                w.configure(state='normal' if enabled else 'disabled')

    # ---------------------------------------------------------- load/reset
    def load_ini(self):
        path = filedialog.askopenfilename(
            title='Open a starting INI file',
            filetypes=[('INI files', '*.ini'), ('All files', '*.*')])
        if not path:
            return
        try:
            with open(path, 'r', encoding='utf-8') as fh:
                self.template_text = fh.read()
        except OSError as exc:
            messagebox.showerror('Error', f'Cannot read file:\n{exc}')
            return
        cfg = configparser.ConfigParser(inline_comment_prefixes=('#', ';'))
        cfg.read_string(self.template_text)
        self._loaded_values = {}
        for (sec, key, _, _, _, _, _, _) in FIELDS:
            value = ''
            if cfg.has_option(sec, key):
                value = cfg[sec][key]
            else:
                # fallback: key may live in another (legacy) section
                for alt in cfg.sections():
                    if cfg.has_option(alt, key):
                        value = cfg[alt][key]
                        break
            self._loaded_values[(sec, key)] = value
            self.vars[(sec, key)].set(value)
        self.title(f'FREQCOR - INI configurator - {os.path.basename(path)}')
        self._update_preview()

    def _reset_values(self):
        """Restore the values of the loaded .ini (or clear all fields)."""
        for (sec, key, *_rest) in FIELDS:
            value = (self._loaded_values or {}).get((sec, key), '')
            self.vars[(sec, key)].set(value)
        self._update_preview()

    # ------------------------------------------------------------- generate
    def _generate_ini_text(self):
        """Write managed values back into the template text (comments kept)."""
        text = self.template_text
        if not text:
            # fallback: plain configparser output
            cfg = configparser.ConfigParser()
            for (sec, key, *_rest) in FIELDS:
                if not cfg.has_section(sec):
                    cfg.add_section(sec)
                cfg[sec][key] = self.vars[(sec, key)].get()
            import io
            buf = io.StringIO()
            cfg.write(buf)
            return buf.getvalue()

        for (sec, key, *_rest) in FIELDS:
            val = self.vars[(sec, key)].get()
            pattern = re.compile(rf'(?m)^({re.escape(key)}\s*=).*$')
            # replace every occurrence for keys that moved to their own
            # section (e.g. plot_* under both PLOTS and PROCEDURE_OPTIONS)
            cnt = 0 if sec == 'PLOTS' else 1
            text = pattern.sub(lambda m, v=val: f"{m.group(1)} {v}",
                               text, count=cnt)

        # make sure the output contains a [PLOTS] section
        if '[PLOTS]' not in text:
            block = '\n[PLOTS]\n' + '\n'.join(
                f'{key} = {self.vars[(sec, key)].get()}'
                for (sec, key, *_r) in FIELDS if sec == 'PLOTS') + '\n'
            text = text.rstrip() + '\n' + block
        return text

    def _update_preview(self):
        if not hasattr(self, 'preview'):
            return
        text = self._generate_ini_text()
        self.preview.configure(state='normal')
        self.preview.delete('1.0', 'end')
        self.preview.insert('1.0', text)
        self._highlight_ini()
        self.preview.configure(state='disabled')

    def _highlight_ini(self):
        """Apply syntax colouring to the preview widget."""
        self._key_line = {}          # INI key -> preview line number
        n_lines = int(self.preview.index('end-1c').split('.')[0])
        for i in range(1, n_lines + 1):
            line = self.preview.get(f'{i}.0', f'{i}.end')
            stripped = line.lstrip()
            indent = len(line) - len(stripped)
            if not stripped:
                continue
            if stripped.startswith(('#', ';')):
                self.preview.tag_add('comment', f'{i}.{indent}', f'{i}.end')
            elif stripped.startswith('['):
                self.preview.tag_add('section', f'{i}.{indent}', f'{i}.end')
            elif '=' in line:
                eq = line.index('=')
                self._key_line[line[:eq].strip()] = i
                self.preview.tag_add('key', f'{i}.{indent}',
                                     f'{i}.{eq}')
                value = line[eq + 1:].strip()
                if not value:
                    continue
                vstart = eq + 1 + len(line[eq + 1:]) - len(line[eq + 1:].lstrip())
                tag = ('number' if re.fullmatch(r'[-+]?[\d.,eE\s]+', value)
                       else 'string')
                self.preview.tag_add(tag, f'{i}.{vstart}', f'{i}.end')

    def _build_config(self):
        cfg = configparser.ConfigParser(inline_comment_prefixes=('#', ';'))
        cfg.read_string(self._generate_ini_text())
        return cfg

    # ------------------------------------------------------------ validation
    def _validate_only(self):
        errors = self._validate()
        if errors is None:
            messagebox.showinfo('Validation',
                                'Validator unavailable\n'
                                '(src/FREQCOR_validate.py not found).')
        elif not errors:
            messagebox.showinfo('Validation', 'Configuration is valid.')
        else:
            messagebox.showwarning('Validation',
                                   'Errors:\n- ' + '\n- '.join(errors))

    def _validate(self):
        if validate_config is None:
            return None
        try:
            valid, errors = validate_config(self._build_config())
        except Exception as exc:
            return [f'validation failed: {exc}']
        return [] if valid else errors

    def save_ini(self):
        errors = self._validate()
        if errors:
            go = messagebox.askyesno(
                'Validation',
                'The configuration contains errors:\n- '
                + '\n- '.join(errors)
                + '\n\nSave anyway?')
            if not go:
                return
        path = filedialog.asksaveasfilename(
            title='Save the INI file',
            defaultextension='.ini',
            initialfile='FREQCOR_config.ini',
            filetypes=[('INI files', '*.ini'), ('All files', '*.*')])
        if not path:
            return
        try:
            with open(path, 'w', encoding='utf-8') as fh:
                fh.write(self._generate_ini_text())
        except OSError as exc:
            messagebox.showerror('Error', f'Cannot write file:\n{exc}')
            return
        messagebox.showinfo('Saved', f'File written:\n{path}')

    # ------------------------------------------------------------------ help
    def _find_figure(self, keyword):
        """Return the path of the newest PNG in tools/figures matching keyword."""
        if not keyword or not os.path.isdir(FIG_DIR):
            return None
        if keyword in self._fig_cache:
            return self._fig_cache[keyword]
        matches = sorted(glob.glob(os.path.join(FIG_DIR, f'*{keyword}*.png')),
                         key=os.path.getmtime)
        fig = matches[-1] if matches else None
        self._fig_cache[keyword] = fig
        return fig

    def _show_group_help(self, grp):
        """Open one help popup aggregating all fields of a logical group."""
        parts = []
        for k in grp['keys']:
            meta = FIELD_META.get(k)
            if meta and meta[0]:
                parts.append(f'{k}:\n{meta[0]}')
        self._show_help(grp['title'], '\n\n'.join(parts), grp['fig'])

    def _show_help(self, key, help_txt, fig_keyword):
        win = tk.Toplevel(self)
        win.title(f'Help: {key}')
        self._help_windows.append(win)

        n_lines = (help_txt or '').count('\n') + 1
        txt = tk.Text(win, wrap='word', height=min(24, n_lines + 1),
                      font=('Segoe UI', 9))
        txt.pack(fill='x', padx=6, pady=6)
        txt.insert('1.0', help_txt or 'No documentation for this field.')
        txt.configure(state='disabled')

        if fig_keyword and fig_keyword.startswith('schem:'):
            self._embed_schematic(win, fig_keyword[6:])
            return

        fig = self._find_figure(fig_keyword)
        if fig:
            try:
                img = tk.PhotoImage(file=fig)
                factor = max(1, -(-img.width() // 540))   # ceil to fit 540 px
                if factor > 1:
                    img = img.subsample(factor, factor)
                lbl = ttk.Label(win, image=img)
                lbl.image = img          # keep reference
                lbl.pack(padx=6, pady=4)
                ttk.Label(win, text=os.path.basename(fig),
                          font=('Segoe UI', 8)).pack()
            except tk.TclError:
                ttk.Label(win,
                          text=f'(unreadable figure: {os.path.basename(fig)})'
                          ).pack(padx=6, pady=4)
        elif fig_keyword:
            ttk.Label(win, text='(no figure available)'
                      ).pack(padx=6, pady=4)

    def _embed_schematic(self, win, name):
        """Draw a synthetic explanatory figure into the help popup."""
        try:
            from matplotlib.figure import Figure
            from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
        except ImportError:
            ttk.Label(win, text='(matplotlib unavailable: cannot draw '
                                'schematic)').pack(padx=6, pady=4)
            return
        drawer = SCHEMATICS.get(name)
        if drawer is None:
            ttk.Label(win, text='(no schematic for this field)'
                      ).pack(padx=6, pady=4)
            return
        fig = Figure(figsize=getattr(drawer, 'figsize', (4.8, 3.2)),
                     dpi=100)
        ax = fig.add_subplot(111, polar=(name == 'sector'))
        drawer(ax)
        fig.tight_layout()
        canvas = FigureCanvasTkAgg(fig, master=win)
        canvas.draw()
        canvas.get_tk_widget().pack(padx=6, pady=4)


def main():
    app = IniGui()
    app.mainloop()


if __name__ == '__main__':
    main()
