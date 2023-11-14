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
    
def compress(dB_Range,beamformed_env):
    
    log_env=20*np.log10(beamformed_env/100)
    Reject_level = 0.1
    log_env=1/dB_Range*(log_env+dB_Range)
    #set rejected values to 0 (0.1??)
    log_env[log_env<Reject_level] = 0
    return log_env

def image(gain,log_env,Lateral_Line_Count,depth,gridNum,x_ele):
    #depth = 20, "double";
    ImgSize_z=depth
    gridSize = depth/gridNum

    pre_z=np.arange(0,gridNum)*gridSize;

    # Image size x is in mm
 
    ImgSize_x=(x_ele(end)-x_ele(1));
    # --- scan conversion (place holder)

    X_RESOLUTION = 1024;     # [pixels]
    # Temporary scale fix mm,cm
    Z_RESOLUTION = 10*(ImgSize_z/ImgSize_x)*(X_RESOLUTION)     # [pixels]

    pre_x=linspace(x_ele(1),x_ele(end),Lateral_Line_Count)

    # use factor of 10 for L7-4
    # aspect_ratio_correction = (10*ImgSize_z)/ImgSize_x;

    aspect_ratio_correction = (ImgSize_z)/ImgSize_x
    Z_RESOLUTION=Z_RESOLUTION*aspect_ratio_correction

    # x_pxlsize=ImgSize_x/X_RESOLUTION
    # z_pxlsize=ImgSize_z/Z_RESOLUTION

    pos_vec_x_new = (0:1 / (X_RESOLUTION-1):1) .* ImgSize_x - ImgSize_x / 2;  #lateral
    pos_vec_z_new = (0:1 / (Z_RESOLUTION-1):1) .* ImgSize_z; % axial
    [pos_x,pos_z]=meshgrid(pos_vec_x_new,pos_vec_z_new)
    scale_log_env=interp2(pre_x,pre_z,log_env,pos_x,pos_z)

    # --- output image normalization (Do we want this???)
    # env_disp=uint8(255*scale_log_env/max(max(scale_log_env)));

    env_disp=uint8(gain*scale_log_env)
    
    return env_disp

    