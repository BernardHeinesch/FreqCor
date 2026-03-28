# -*- coding: utf-8 -*-
"""Reference cospectra models used in FREQCOR.

This module provides:
- Kaimal (1972)-type theoretical cospectra parameterizations.
- A Massman et al. (2004) model cospectrum.
"""

def Kaimal_cosp (nf, kf, zL):
    """
    Compute the Kaimal theoretical cospectrum.

    Parameters
    ----------
    nf : float
        Natural frequency (Hz).
    kf : float
        Normalized frequency (dimensionless).
    zL : float
        Stability parameter (z/L).

    Returns
    -------
    kaimal_cosp : float
        Theoretical cospectrum value.

    """
    if zL > 0: # Stable conditions
        kf0 = 0.23*(1 + 6.4*zL)**0.75
        kaimal_cosp=0.81*kf/kf0/(1 + 1.5*(kf/kf0)**(2.1))/nf
        # reynolds stresses
        Au = 0.124 * ((1 + 7.9 * zL)**0.75)
        Bu = 2.34 * (Au**(-1.1))
        Cospwu = kf / (nf * (Au + Bu * kf**2.1))
    else : # Unstable conditions
        if kf <= 1:
            kaimal_cosp=11 *kf/(1 + 13.3*kf)**(7 /4 )/nf
        else:
            kaimal_cosp=4 *kf/(1 +3.8*kf)**(7 /3 )/nf
            
    return kaimal_cosp

def Kaimal_cosp_EP(fnorm,zL): # According to the EddyPro source code (G. Fratini)
    """Compute the Kaimal cospectrum shape used in EddyPro.

    Parameters
    ----------
    fnorm : float
        Normalized frequency (dimensionless).
    zL : float
        Stability parameter (z/L).

    Returns
    -------
    kaimal : float
        Theoretical cospectrum shape value.
    """
    if zL > 0: # Stable conditions
       Ak = 0.284 * ((1 + 6.4 * zL)**0.75)
       Bk = 2.34  * (Ak**(-1.1))
       kaimal = fnorm / (Ak + Bk * fnorm**2.1)
    else : # Unstable conditions
        if fnorm <= 0.54:
            kaimal = 12.92 * fnorm / ((1 + 26.7 * fnorm)**1.375)
        else:
            kaimal = 4.378 * fnorm / ((1 +  3.8 * fnorm)**2.4)
            
    
    return kaimal


def theor_Massman(nf, kf, zL, A0, kf0, mu):
    """
    Compute a model cospectrum according to Massman et al. (2004, Eq. 4.2).

    Parameters
    ----------
    nf : float
        Natural frequency (Hz).
    kf : float
        Normalized frequency (dimensionless).
    zL : float
        Stability parameter (z/L).
    A0 : array-like
        Normalisation factor, [stable, unstable].
    kf0 : array-like
        Frequency at which cospectrum attains highest value (fpeak), 
        [stable, unstable].
    mu : array-like
        Broadness factor, [stable, unstable].

    Returns
    -------
    Cospws : float
        Modelled cospectrum value.
    """

    if zL > 0: # Stable
      Cospws=A0[0]*(kf/kf0[0])/((1+(kf/kf0[0])**(2*mu[0]))**(1.1667/mu[0]))/nf
    else: # Unstable
      Cospws=A0[1]*(kf/kf0[1])/((1+(kf/kf0[1])**(2*mu[1]))**(1.1667/mu[1]))/nf

    return(Cospws)

