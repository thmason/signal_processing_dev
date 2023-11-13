# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 15:55:44 2023

@author: timhm
"""

import numpy as np
from scipy import signal as sig
import sys


def fourierFFT(beamformed_interp, gridNum):
    """
    this appears to only procses the first row

    Parameters
    ----------
    beamformed_interp : TYPE
        DESCRIPTION.
    gridNum : TYPE
        DESCRIPTION.

    Returns
    -------
    P1 : TYPE
        DESCRIPTION.

    """
    Y = np.fft.fft(beamformed_interp,axis=0)
    P2 = np.abs(Y/gridNum)    
    P1 = P2[:int(np.floor(gridNum/2.0)+1)]
    P1[1:-2] = 2*P1[1:-2]
    return P1


def filterBPF(N,Fc1,Fc2,Fs):
    """
    FILTER1 Returns a discrete-time filter object.
    
    FIR Window Bandpass filter designed using the FIR1 function.
    
    All frequency values are in MHz.
    
    Parameters
    ----------
    N : TYPE
        DESCRIPTION.
    Fc1 : TYPE
        DESCRIPTION.
    Fc2 : TYPE
        DESCRIPTION.
    Fs : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """    

    # Calculate the coefficients using the firwin function.
    b = sig.firwin(N+1,  np.array([Fc1, Fc2])/(Fs/2), window='hann', pass_zero=False, scale = True)
    return(b)
    
def SGF(FFT):
    SFG = sig.savgol_filter(FFT, window_length=101, polyorder=5, axis=0)
    return SFG
    

