# -*- coding: utf-8 -*-
"""Vickers & Mahrt (1997) quality flags aggregation.

This module aggregates the per-variable Vickers & Mahrt (VM) test flags
exported by EddyPro into a compact set of global flags (ONE FOR EACH
VARIABLE: w, T, CO2, H2O), used downstream for screening (co)spectra and fluxes.

This global flag should be zero (all flags zero)
The considered flag types are : spikes,drop out,skewness/kurtosis(hard),
discontinuity(hard) 

"""
import numpy as np

# Define function to spilt a chain of characters into single characters
def split(word): 
    """Split a string into a list of single-character strings.

    Parameters
    ----------
    word : str
        Input string.

    Returns
    -------
    list[str]
        List of characters.
    """
    return [char for char in word]  


def FREQCOR_VM_flag(NumberDataFLAGF, FlagVM_sp, FlagVM_do, FlagVM_skh, FlagVM_dh, gss):
    """Aggregate VM test flags into one global flag per variable.

    EddyPro stores VM flags as a multi-digit code per half-hour, where each
    digit corresponds to a variable (u, v, w, T, CO2, H2O, CH4, ...). This
    function decodes those codes for several VM tests (spikes, dropout,
    skewness/kurtosis, discontinuities) and sums the digits for each target
    variable to produce a global flag.

    Parameters
    ----------
    NumberDataFLAGF : int
        Number of half-hours.
    FlagVM_sp, FlagVM_do, FlagVM_skh, FlagVM_dh : array-like
        EddyPro VM flag codes (one per half-hour) for spikes, dropout,
        skewness/kurtosis, and discontinuities.
    gss : int
        Gas species selector (1: CO2, 2: H2O, 3: O3, 4: CH4, 5: N2O).

    Returns
    -------
    FlagVM_w : numpy.ndarray
        Global flag for vertical wind speed (w).
    FlagVM_T : numpy.ndarray
        Global flag for temperature.
    FlagVM_g : numpy.ndarray
        Global flag for the selected gas.

    Notes
    -----
    The function discards a VM test column if it contains only ``8`` and ``9``
    digits (a typical pattern indicating the test was not selected in EddyPro).
    """
# Initialise flag matrix

# Decoding Eddypro flag: matflag(js, jtest, jvar)
# jtest: spikes, dropout, skeness, discontinuities 
# jvar: 8, u, v, w, T, CO2, H2O, CH4, 9
# Grouping VM flags by variables
    matflag = [[ [0] * 4 ] * NumberDataFLAGF]*9
    matflag = np.asarray(matflag)
        
    # Warning !!! Indexing of a 3D Matrix is [x, y, z] where
    # x : level of the matrix (3rd dimension, depth)
    # y : rows
    # z : columns

    # Check if any of the flag columns contain only 8 and 9 values
    # This would indicate they haven't been selected and should be discarded
    flag_columns_valid = [True, True, True, True]  # [sp, do, skh, dh]
    
    # Function to check if a flag value contains only 8 and 9
    def contains_only_8_and_9(flag_value):
        # Convert to string and check each character
        flag_str = str(int(flag_value))
        return all(char in ['8', '9'] for char in flag_str)
    
    # Check FlagVM_sp
    valid_flags_sp = [flag for flag in FlagVM_sp if flag > 800000000]
    if len(valid_flags_sp) > 0 and all(contains_only_8_and_9(flag) for flag in valid_flags_sp):
        flag_columns_valid[0] = False
        print("Warning: Spike flag column contains only 8 and 9 values and will be discarded.")
    
    # Check FlagVM_do
    valid_flags_do = [flag for flag in FlagVM_do if flag > 800000000]
    if len(valid_flags_do) > 0 and all(contains_only_8_and_9(flag) for flag in valid_flags_do):
        flag_columns_valid[1] = False
        print("Warning: Dropout flag column contains only 8 and 9 values and will be discarded.")
    
    # Check FlagVM_skh
    valid_flags_skh = [flag for flag in FlagVM_skh if flag > 800000000]
    if len(valid_flags_skh) > 0 and all(contains_only_8_and_9(flag) for flag in valid_flags_skh):
        flag_columns_valid[2] = False
        print("Warning: Skewness/kurtosis flag column contains only 8 and 9 values and will be discarded.")
    
    # Check FlagVM_dh
    valid_flags_dh = [flag for flag in FlagVM_dh if flag > 800000000]
    if len(valid_flags_dh) > 0 and all(contains_only_8_and_9(flag) for flag in valid_flags_dh):
        flag_columns_valid[3] = False
        print("Warning: Discontinuities flag column contains only 8 and 9 values and will be discarded.")

    for js in range(0,NumberDataFLAGF):
        if FlagVM_sp[js]>800000000:           # Elimination of missing data : why sp flag as reference?
            # Split the flag code into the 9 variables (3rd dimension)
            # This is done for 4 flags of interest (4 columns) over the whole
            # dataset (nspec or NumberDataFLAGF 1/2 hours)
            if flag_columns_valid[0]:
                matflag[:,js,0]=np.asarray(split(str(int(FlagVM_sp[js])))).astype('int32')
            else:
                matflag[:,js,0]=np.asarray([0]*9)  # Use zeros if column is invalid
                
            if flag_columns_valid[1]:
                matflag[:,js,1]=split(str(int(FlagVM_do[js])))
            else:
                matflag[:,js,1]=[0]*9  # Use zeros if column is invalid
                
            if flag_columns_valid[2]:
                matflag[:,js,2]=split(str(int(FlagVM_skh[js])))
            else:
                matflag[:,js,2]=[0]*9  # Use zeros if column is invalid
                
            if flag_columns_valid[3]:
                matflag[:,js,3]=split(str(int(FlagVM_dh[js])))
            else:
                matflag[:,js,3]=[0]*9  # Use zeros if column is invalid
        # Why this else statement??    
        else:
            matflag[:,js,0]=np.asarray([9]*9)
            matflag[:,js,1]=[9]*9
            matflag[:,js,2]=[9]*9
            matflag[:,js,3]=[9]*9
    
    
    # Creating global flags for each (co)spectra: w, T, gas (CO2 or H2O)
    FlagVM_w=(matflag[3,:,:]).sum(axis=1)
    FlagVM_T=(matflag[4,:,:]).sum(axis=1)
    
    if gss==1:
        FlagVM_g=(matflag[5,:,:]).sum(axis=1)        #### Flag for CO2 
    elif gss==2:
        FlagVM_g=(matflag[6,:,:]).sum(axis=1)       #### Flag for H2O
    elif gss==3:
        FlagVM_g=(matflag[8,:,:]).sum(axis=1)       #### Flag for O3
    elif gss==4:
        FlagVM_g=(matflag[7,:,:]).sum(axis=1)    #### Flag for CH4
    elif gss==5:
        FlagVM_g=(matflag[8,:,:]).sum(axis=1)   #### Flag for n2o
    # Deleting useless intermediate files
    del matflag
    
    return FlagVM_w, FlagVM_T, FlagVM_g
