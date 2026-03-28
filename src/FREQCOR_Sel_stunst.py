# -*- coding: utf-8 -*-
"""Selection utilities for splitting (co)spectra by stability class."""
import numpy as np
import pandas as pd

def FREQCOR_Sel_stunst(Icosp_sel, sts, Zeta):
    """
    Select stable or unstable conditions based on the sign of `Zeta`.

    Parameters
    ----------
    Icosp_sel : pd.DataFrame
        Ideal cospectra for CF (sel) [nfreq, nspec].
    sts : int
        Stability class selector (1: unstable, 2: stable).
    Zeta : pandas.Series
        Stability parameter time series.

    Returns
    -------
    Icos_CF : pandas.DataFrame
        Stability-filtered cospectra used for CF computation.

    """
    # Unstable: select Zeta <= 0
    # Stable: select Zeta > 0
    mask = (Zeta <= 0).to_numpy() if sts == 1 else (Zeta > 0).to_numpy()
    mask_s = pd.Series(mask, index=Icosp_sel.columns)

    Icos_CF = Icosp_sel.copy()
    Icos_CF.loc[:, ~mask_s] = np.nan

    return Icos_CF
