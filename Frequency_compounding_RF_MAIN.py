# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 08:46:23 2023

@author: Tim Mason
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
from scipy import signal as sig
from scipy import interpolate
import math_utilities as mu

import sys
import os
import h5py
import mat73



fundamental_filename = './data/C4-2_Fundamental.mat'
harmonic_filename = './data/C4-2_Harmonic.mat'
filterParams_filename = './data/FilterParam_C4-2_Harmonic.mat'
comboParams_filename = './data/ComboParams.mat'

#get beamfored data
# get filter parameters
fundamental_data = sio.loadmat(fundamental_filename)
harmonic_data = sio.loadmat(harmonic_filename)
filterParams_C4_2_Harmonic = mat73.loadmat(filterParams_filename)
comboParams = mat73.loadmat(comboParams_filename)


plt.close('all')    
# gamma is used to creat a "poor man's" grayscale map" using the matlab function, imadjust
gamma = 1.6


# get filter parameters
filterParams = filterParams_C4_2_Harmonic['FilterParam']

# C4_2 Harmonic Data
filter_order = int(float(str(filterParams['Filter_Order'])))
# parameters for non-compounded image
# Filter0
FcLo_0=float(filterParams['FcLo_0'])
FcHi_0=float(filterParams['FcHi_0'])
Weighting_0=float(filterParams['Weighting_0'])
dB_Range=float(filterParams['dB_Range'])

# Compounding Channels, using filters 1, 2 and 3 (where number of channels
# to be compounded is variable from 1 to 3)

Num_Compounding_Channels=filterParams['Num_Compounding_Channels']

# Filter1
FcLo_1=float(filterParams['FcLo_1'])
FcHi_1=float(filterParams['FcHi_1'])
Weighting_1=float(filterParams['Weighting_1'])
dB_Range1=float(filterParams['dB_Range1'])

# Filter2
FcLo_2=float(filterParams['FcLo_2'])
FcHi_2=float(filterParams['FcHi_2'])
Weighting_2=float(filterParams['Weighting_2'])
dB_Range2=float(filterParams['dB_Range2'])

# Filter3
FcLo_3=float(filterParams['FcLo_3'])
FcHi_3=float(filterParams['FcHi_3'])
Weighting_3=float(filterParams['Weighting_3'])
dB_Range3=float(filterParams['dB_Range3'])

#print(comboParams)
param = comboParams['ComboParam']
beamformed = harmonic_data['Beamformed']




# Edit Values here as desired


# Num_Compounding_Channels=3
# Filter_Order=Filter_Order


# FcLo_0=FcLo_0
# FcHi_0=FcHi_0
# Weighting_0=Weighting_0
# dB_Range=dB_Range

# FcLo_1=FcLo_1
# FcHi_1=FcHi_1
# Weighting_1=Weighting_1
# dB_Range1=dB_Range1

# FcLo_2=FcLo_2
# Weighting_2=Weighting_2
# dB_Range2=dB_Range2

# FcLo_3=FcLo_3
# FcHi_3=FcHi_3
# Weighting_3=Weighting_3
# dB_Range3=dB_Range3




# %  Find Maximum value in Beamformed
Q=np.max(np.max(beamformed))

# FUDGE factor: apply digital gain to overcome weak signal from FrontEnd
DGain = (32767/(4*Q))
# DGain = 1
beamformed = DGain*beamformed

#param.gridNum = size(Beamformed,1);
param['gridNum'] = beamformed.shape[0]

Fs = float(param['fs'])


# Post beamformer processing
# --- lateral interpolation, 4x128 512


# Compute interpolation filter coefficients
interpFactor=4;
x=param['x_ele'] # To be reconstructed pixel coordinates
step_size = 1.0/interpFactor
interp_points = np.arange(1,len(x) + step_size, step_size)
interpolator = interpolate.interp1d(np.arange(1,len(x)+1), x)
x_interp = interpolator(interp_points)


# % Low pass filtering for post-interpolation
lat_fs=1/(x_interp[1]-x_interp[0])
lat_cutoff=0.5*1/(x[1]-x[0])

b, a  = sig.butter(12, lat_cutoff/(lat_fs/2))
sos = sig.butter(12, lat_cutoff/(lat_fs/2), output='sos')


beamformed_interp = np.zeros((beamformed.shape[0],len(interp_points)))

for nrow in range(beamformed.shape[0]):
    
    # set interpolater
    interpolator = interpolate.interp1d(x, beamformed[nrow,:].flatten())
    interpRow = interpolator(x_interp)    
    #interpRow_filt = sig.sosfilt(sos, interpRow)
    interpRow_filt = sig.filtfilt(b, a, interpRow) # for some reason this isn't being stored in an array??    
    beamformed_interp[nrow,:]=(interpRow)



Lateral_Line_Count = beamformed_interp.shape[1]


# **********   below function only seems to process the first row ********
FFT = mu.fourierFFT(beamformed_interp,param['gridNum'])
beamformedSpectralContent=mu.SGF(FFT);

# --- Filter0 (Used for non compounded image)
filter_taps = mu.filterBPF(filter_order,FcLo_0,FcHi_0,Fs)
beamformed_Filtered0 = Weighting_0*sig.lfilter(filter_taps, 1, beamformed_interp, axis=0)

# %% --- Filter1
filter_taps1 = mu.filterBPF(filter_order,FcLo_1,FcHi_1,Fs)
beamformed_Filtered1 = Weighting_1*sig.lfilter(filter_taps1, 1, beamformed_interp, axis=0)

# %% --- Filter2
filter_taps2 = mu.filterBPF(filter_order,FcLo_2,FcHi_2,Fs)
beamformed_Filtered2 = Weighting_2*sig.lfilter(filter_taps2, 1, beamformed_interp, axis=0)

# %% --- Filter3
filter_taps3 = mu.filterBPF(filter_order,FcLo_3,FcHi_3,Fs)
beamformed_Filtered3 = Weighting_3*sig.lfilter(filter_taps3, 1, beamformed_interp, axis=0)


FFT1 = mu.fourierFFT(beamformed_Filtered1,param['gridNum'])
beamformedSpectralContent1 = mu.SGF(FFT1)

FFT2 = mu.fourierFFT(beamformed_Filtered2,param['gridNum'])
beamformedSpectralContent2 = mu.SGF(FFT2)

FFT3 = mu.fourierFFT(beamformed_Filtered3,param['gridNum'])
beamformedSpectralContent3 = mu.SGF(FFT3)


fig, ax = plt.subplots(2,1)

f = Fs/param['gridNum']*np.arange(0,param['gridNum']/2 + 1)
ax[0].plot(f, FFT[:,0])
ax[0].set_title('Single-Sided Amplitude Spectrum')
ax[0].set_xlabel('f (MHz)')
ax[0].set_ylabel('|P1(f)|')

ax2 = ax[0].twinx()
ax2.set_ylabel('db')

if Num_Compounding_Channels == 1 :
    pass
elif Num_Compounding_Channels == 2 :
    pass
elif Num_Compounding_Channels == 3 :
    w, h = sig.freqz(filter_taps1,  fs=Fs)
    ax2.plot(w, 20 * np.log10(np.abs(h)), label='filter 1')

    w, h = sig.freqz(filter_taps2,  fs=Fs)
    ax2.plot(w, 20 * np.log10(np.abs(h)), label = 'filter 2')

    w, h = sig.freqz(filter_taps3,  fs=Fs)
    ax2.plot(w, 20 * np.log10(np.abs(h)), label = 'filter 3')    
    
    ax2.legend()



ax[1].plot(f, FFT[:,0], label = 'FFT')
ax[1].plot(f, beamformedSpectralContent[:,0], label='FFTSG')
ax[1].set_title('Single-Sided Amplitude Spectrum')
ax[1].set_xlabel('f (MHz)')
ax[1].set_ylabel('|P1(f)|')

#ax_1_left = ax[1].twinx()

if Num_Compounding_Channels == 1 :
    pass
elif Num_Compounding_Channels == 2 :
    pass
elif Num_Compounding_Channels == 3 :
    ax[1].plot(f, beamformedSpectralContent1[:,0], label = 'Filter 1')
    ax[1].plot(f, beamformedSpectralContent2[:,0], label = 'Filter 3')
    ax[1].plot(f, beamformedSpectralContent3[:,0], label = 'Filter 4')
    ax[1].legend()


# Envelop detection (abs and Hilbert returns amplitude detected)
beamformed_env_0=np.abs(sig.hilbert(beamformed_Filtered0, axis=0))
beamformed_env_1=np.abs(sig.hilbert(beamformed_Filtered1, axis=0))
beamformed_env_2=np.abs(sig.hilbert(beamformed_Filtered2, axis=0))
beamformed_env_3=np.abs(sig.hilbert(beamformed_Filtered3, axis=0))

# log compression
log_env_0 = Weighting_0*mu.compress(dB_Range1,beamformed_env_0)
log_env_1 = Weighting_0*mu.compress(dB_Range1,beamformed_env_1)
log_env_2 = Weighting_0*mu.compress(dB_Range2,beamformed_env_2)
log_env_3 = Weighting_0*mu.compress(dB_Range3,beamformed_env_3)

if Num_Compounding_Channels == 1:
    channel_1 = 1
    channel_2 = 0
    channel_3 = 0;
elif Num_Compounding_Channels == 2:
    channel_1 = 1
    channel_2 = 1
    channel_3 = 0
elif Num_Compounding_Channels == 3:
    channel_1 = 1
    channel_2 = 1
    channel_3 = 0
else:
    channel_1 = 1
    channel_2 = 1
    channel_3 = 1

tempParams={}
tempParams['channel_1'] = channel_1
tempParams['channel_2'] = channel_2
tempParams['channel_3'] = channel_3

# Initialize Gain levels
gain=50;
# envelope gain for NonCompounded
tempParams['gain1'] = gain;
tempParams['gain2'] = gain;
tempParams['Num_Compounding_Channels'] = Num_Compounding_Channels
tempParams['gamma'] = gamma;
tempParams['lateral_line_count'] = Lateral_Line_Count

# path='C:\Users\DavidRoundhill\DavidFiles\SimulationTool\AE_FC_Line_Processor\';


# save([path,'ImageDataTemp.mat'],'TempParams','log_env_0','log_env_1','log_env_2','log_env_3','param','-v7.3');

# figure(4)

# set(gcf,'units','pixels','position', [1000 1000 800 80]);
# set(gcf, 'name', 'non-Compounded Image');

# % Create a slider
# slider1 = uicontrol('Style', 'slider', 'Position', [50, 30, 200, 20],...
#     'Min', 0, 'Max', 100, 'Value', Gain, 'Callback', @sliderCallback1);

# % Create another slider
# slider2 = uicontrol('Style', 'slider', 'Position', [300, 30, 200, 20],...
#     'Min', 0, 'Max', 100, 'Value', Gain, 'Callback', @sliderCallback2);


# % Create a text box to display the slider value 

# txt1 = uicontrol('Style', 'text','Position',[300, 50, 200, 20]);
# txt2 = uicontrol('Style', 'text','Position',[50, 50, 200, 20]);

   
# % Initialize the variable

# set(txt1,'String',num2str(Gain));
# set(txt2,'String',num2str(Gain));

env_disp = mu.image(3*tempParams['gain1'],
                    log_env_0,
                    tempParams['lateral_line_count'],
                    float(param['depth']),
                    param['gridNum'],
                    param['x_ele'])

 

# gamma=1.6;

# A = 0:255;
# A=A/255;
# figure(5)

# grayscale = imadjust(A,[0 1],[0 1],gamma);
# plot(grayscale,"LineWidth",1, "color", [0.75,0.75,0.75]) 
# hold off

# figure(1)
# set(gcf,'units','pixels','position', [1024 1200 900 600]);
# set(gcf, 'name', 'non-Compounded Image');

# env_disp_NonCompounded = imadjust(env_disp,[0 1],[0 1],gamma);

# imshow(env_disp_NonCompounded, gray)

# Gain=TempParams.Gain2;

# log_env_Compounded = (TempParams.Channel_1*log_env_1+TempParams.Channel_2*log_env_2+TempParams.Channel_3*log_env_3);
# env_disp_Compounded = Image(TempParams.Gain2,log_env_Compounded,TempParams.Lateral_Line_Count,param.depth,param.gridNum,param.x_ele);


# figure(2)
# set(gcf,'units','pixels','position', [1400 1200 900 600])
# set(gcf, 'name', 'Compounded Image');


# env_disp_Compounded = imadjust(env_disp_Compounded,[0 1],[0 1],gamma);


# imshow(env_disp_Compounded, gray)

