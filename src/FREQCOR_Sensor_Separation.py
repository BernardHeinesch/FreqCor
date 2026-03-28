# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 12:05:39 2023

This module implements the combined crosswind and vertical separation
correction factor used when applying spectral corrections.
"""

import numpy as np

# Computation of Horst and Lenshow correction factors
# for crosswind and vertical separation
# Horst and Lenshow BLM 2009
# Parameters needed
# rseast and rsnorth: Eastward and Northward crosswind distances between instruments 
# rsz: vertical separation of instruments 
# zmeas : measurement height 
# disph : zero plane displacement height
# WD : wind direction
###########################################################################
# Crosswind separation:
    # Application Eq. (16) from Horst and Lenshow, parameters from Eq. (18)
# Vertical separation:
    # Application of Eq. (28) from HL, parameters from Eq. (29) and (30)


def FREQCOR_Sensor_Separation(zmeas,disph,rsz,WD,rseast,rsnorth,Zeta,eq):
    """Compute the Horst and Lenschow (2009) sensor separation correction factor.

    The correction accounts for:
    - crosswind separation between instruments (Eq. 16 or Eq. 13), and
    - vertical separation (Eq. 28).

    Parameters
    ----------
    zmeas : float
        Measurement height (m).
    disph : float
        Zero-plane displacement height (m).
    rsz : float
        Vertical separation between instruments (m).
    WD : pandas.Series or array-like
        Wind direction (degrees).
    rseast : float
        Eastward separation distance between instruments (m).
    rsnorth : float
        Northward separation distance between instruments (m).
    Zeta : array-like
        Stability parameter series (z/L).
    eq : int
        Horst and Lenschow equation selector for the crosswind term (13 or 16).

    Returns
    -------
    CFHL : numpy.ndarray
        Multiplicative correction factor (same length as `Zeta`).

    Notes
    -----
    If `eq` is not 13 or 16, the function prints an error message and returns
    an array that will typically contain NaNs.
    """
    HLcw = np.zeros(len(Zeta))
    HLcw[:]=np.nan
    HLvert = np.zeros(len(Zeta))
    HLvert[:] = np.nan
    CFHL = np.zeros(len(Zeta))
    CFHL[:] = np.nan
    
    zinst=zmeas-disph
    zson=zinst
    zanal=zinst+rsz
    
    # Calculate per stability regime
    nmy = np.zeros(len(Zeta))
    nmy[:] = np.nan
    nmy[Zeta <= -0.05] = 0.15
    nmy[Zeta > -0.05] = 2.43-2.28/(1.01+0.2*Zeta[Zeta>-0.05])**2
    
    bee = np.zeros(len(Zeta))
    bee[:] = np.nan
    if eq == 16:
        bee[Zeta<-0.01] = 1.09
        bee[Zeta>0.01] = 1.22
        bee[(Zeta>=-0.01) & (Zeta <= 0.01)] = 1.14
    elif eq == 13:
        bee[Zeta<-0.01] = 1
        bee[Zeta>0.01] = 1
        bee[(Zeta>=-0.01) & (Zeta <= 0.01)] = 1
    else:
        print('Impossible value of "eq": must be either 13 or 16')
    
    # cross wind distances
    rscw=abs(-rseast*np.cos(np.radians(WD.values))+rsnorth*np.sin(np.radians(WD.values))) # transform the angles in radians 
    
    # Rem: here we calculate directly the correction factor, so F0/F, thus in
    # practice doing 1/Eq.HL, which explains why there is no negative sign
    # (1/exp(1) = exp(-1))
    
    # Cross wind HL: eq. (16) from Horst and Lenschow
    HLcw=np.exp((2*np.pi*nmy*rscw/zinst)**bee)
    
    # Vertical separation
    # Application formules (29) et (30) from Horst and Lenschow 
    nmz = np.zeros(len(Zeta))    
    if zson<zanal:
         nmz[zson*Zeta/zinst <=0.03]= 0.1
         nmz[zson*Zeta/zinst >0.03]=0.43-0.33/(0.964+1.2*zson*Zeta[zson*Zeta/zinst >0.03]/zinst)**2
         # Vertical HL : eq. (28)
         HLvert=np.exp(2*np.pi*nmz*abs(rsz)/zson) # the formula uses zmin
    else:
        nmz[zanal*Zeta/zinst <=-0.03]= 0.013
        nmz[zanal*Zeta/zinst >-0.03]=0.3-0.287/(1.051+1.7*zanal*Zeta[zanal*Zeta/zinst >-0.03]/zinst)**2
        # Vertical HL : eq. (28)
        HLvert=np.exp(2*np.pi*nmz*abs(rsz)/zanal)
        # take abs(rz) as per publication - the asymmetry is already accounted 
        # for in the definition of nmz
    
    # Total HL
    CFHL=HLcw*HLvert
    
    return CFHL
