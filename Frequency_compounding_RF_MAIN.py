# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 08:46:23 2023

@author: Tim Mason
"""

import numpy as np
import scipy.io as sio
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


    
# gamma is used to creat a "poor man's" grayscale map" using the matlab function, imadjust
gamme = 1.6


# get filter parameters
filterParams = filterParams_C4_2_Harmonic['FilterParam']

# C4_2 Harmonic Data
filter_order = filterParams['Filter_Order']
# parameters for non-compounded image
# Filter0
FcLo_0=filterParams['FcLo_0']
FcHi_0=filterParams['FcHi_0']
Weighting_0=filterParams['Weighting_0']
dB_Range=filterParams['dB_Range']

# Compounding Channels, using filters 1, 2 and 3 (where number of channels
# to be compounded is variable from 1 to 3)

Num_Compounding_Channels=filterParams['Num_Compounding_Channels']

# Filter1
FcLo_1=filterParams['FcLo_1']
FcHi_1=filterParams['FcHi_1']
Weighting_1=filterParams['Weighting_1'];
dB_Range1=filterParams['dB_Range1'];

# Filter2
FcLo_2=filterParams['FcLo_2']
FcHi_2=filterParams['FcHi_2']
Weighting_2=filterParams['Weighting_2'];
dB_Range2=filterParams['dB_Range2'];

# Filter3
FcLo_3=filterParams['FcLo_3']
FcHi_3=filterParams['FcHi_3']
Weighting_3=filterParams['Weighting_3'];
dB_Range3=filterParams['dB_Range3'];

#print(comboParams)
param = comboParams['ComboParam']
beamformed = harmonic_data['Beamformed']
print(param)

# %  Find Maximum value in Beamformed
Q=np.max(np.max(beamformed))

print(Q)



# % FUDGE factor: apply digital gain to overcome weak signal from FrontEnd
# DGain = (32767/(4*Q));
# %DGain=1;
# Beamformed = DGain*Beamformed;
# param.gridNum = size(Beamformed,1);

# Fs=param.fs


print(filter_order)
print(dB_Range)
sys.exit(1)



# path='C:\Users\DavidRoundhill\DavidFiles\SimulationTool\L7-4_Data\';
# %load([path,'L7-4_Beamformed_data_A.mat']);
# %load([path,'L7-4_Beamformed_data.mat']);

# path='C:\Users\DavidRoundhill\DavidFiles\SimulationTool\C4-2_Data\';
# % load([path,'C4-2_Fundamental.mat']);

# load([path,'C4-2_Harmonic.mat']);
# load([path,'ComboParams.mat']);


# %% GET FILTER PARAMETERS %%
# load([path,'FilterParam_C4-2_Harmonic.mat'])
# %load([path,'FilterParam_C4-2_Fundamental.mat'])
# %load([path,'FilterParam_L7-4.mat'])

# % C4_2 Harmonic Data
# Filter_Order=FilterParam.Filter_Order;
# %parameters for non-compounded image
# %Filter0
# FcLo_0=FilterParam.FcLo_0;
# FcHi_0=FilterParam.FcHi_0;
# Weighting_0=FilterParam.Weighting_0;
# dB_Range=FilterParam.dB_Range;



# % Compounding Channels, using filters 1, 2 and 3 (where number of channels
# % to be compounded is variable from 1 to 3)
# % 
# Num_Compounding_Channels=FilterParam.Num_Compounding_Channels;

# %Filter1
# FcLo_1=FilterParam.FcLo_1;
# FcHi_1=FilterParam.FcHi_1;
# Weighting_1=FilterParam.Weighting_1;
# dB_Range1=FilterParam.dB_Range1;

# %Filter2
# FcLo_2=FilterParam.FcLo_2;
# FcHi_2=FilterParam.FcHi_2;
# Weighting_2=FilterParam.Weighting_2;
# dB_Range2=FilterParam.dB_Range2;

# %Filter3
# FcLo_3=FilterParam.FcLo_3;
# FcHi_3=FilterParam.FcHi_3;
# Weighting_3=FilterParam.Weighting_3;
# dB_Range3=FilterParam.dB_Range3;

# %Edit Values here as desired


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

# param = ComboParam;

# %  Find Maximum value in Beamformed
# Q=max(max(Beamformed));

# % FUDGE factor: apply digital gain to overcome weak signal from FrontEnd
# DGain = (32767/(4*Q));
# %DGain=1;
# Beamformed = DGain*Beamformed;
# param.gridNum = size(Beamformed,1);

# Fs=param.fs

# %define location for figures

# p = get(0, "MonitorPositions");
# Position = p(2, :); % second display
# Position = [Position(1)+100,Position(2)+300, Position(3)/4, Position(4)/2 ];


# %% Post beamformer processing
# % --- lateral interpolation, 4x128 512


# % Compute interpolation filter coefficients
# InterpFactor=4;
# x=param.x_ele(1,:); % To be reconstructed pixel coordinates
# x_interp=interp1(1:length(x),x,[1:1/InterpFactor:length(x)],'linear');
# % Low pass filtering for post-interpolation
# lat_fs=1/(x_interp(2)-x_interp(1));
# lat_cutoff=0.5*1/(x(2)-x(1));
# [b,a] = butter(12,lat_cutoff/(lat_fs/2));



# for nrow=1:size(Beamformed,1)
#     currentRow=Beamformed(nrow,:);
#     interpRow=interp1(x,currentRow,x_interp,'linear');
#     interpRow_filt=filtfilt(b,a,interpRow);  
#     Beamformed_Interp(nrow,:)=interpRow;
# end

# Lateral_Line_Count = size(Beamformed_Interp,2);

# FFT=FourierFFT(Beamformed_Interp,param.gridNum);

# BeamformedSpectralContent=SGF(FFT);

# %% --- Filter0 (Used for non compounded image)
# Hd0=FilterBPF(Filter_Order,FcLo_0,FcHi_0,Fs);
# Beamformed_Filtered0 = Weighting_0*filter(Hd0,Beamformed_Interp);

# %% --- Filter1
# Hd1=FilterBPF(Filter_Order,FcLo_1,FcHi_1,Fs);
# Beamformed_Filtered1 = Weighting_1*filter(Hd1,Beamformed_Interp);

# %% --- Filter2
# Hd2=FilterBPF(Filter_Order,FcLo_2,FcHi_2,Fs);
# Beamformed_Filtered2 = Weighting_2*filter(Hd2,Beamformed_Interp);

# %% --- Filter3
# Hd3=FilterBPF(Filter_Order,FcLo_3,FcHi_3,Fs);
# Beamformed_Filtered3 = Weighting_3*filter(Hd3,Beamformed_Interp);


# FFT1=FourierFFT(Beamformed_Filtered1,param.gridNum);
# BeamformedSpectralContent1=SGF(FFT1);

# FFT2=FourierFFT(Beamformed_Filtered2,param.gridNum);
# BeamformedSpectralContent2=SGF(FFT2);

# FFT3=FourierFFT(Beamformed_Filtered3,param.gridNum);
# BeamformedSpectralContent3=SGF(FFT3);

# figure(3);

# set(gcf,'units','pixels','position', Position);

# clf;

# subplot(2,1,1)
# hold on

# yyaxis right
# f = Fs/param.gridNum*(0:(param.gridNum/2));
# plot(f,FFT,"LineWidth",1, "color", [0.75,0.75,0.75]) 
# title("Single-Sided Amplitude Spectrum")
# xlabel("f (MHz)")
# ylabel("|P1(f)|") 

# yyaxis left
# ylabel("dB")
# f = Fs/param.gridNum*(0:(param.gridNum/2));

# if Num_Compounding_Channels == 1
#     [h,w] = freqz(Hd1,nrow,Fs);
     
#     plot(w,mag2db(abs(h)),"LineWidth",1,"color", [0,1,0])
    
#     legend("Filter1","FFT")


# elseif Num_Compounding_Channels==2
#     [h,w] = freqz(Hd1,nrow,Fs);
        
#     plot(w,mag2db(abs(h)),"LineWidth",1,"color", [0,1,0])
    
#     [h,w] = freqz(Hd2,nrow,Fs);
    
#     plot(w,mag2db(abs(h)),"LineWidth",1,"color", [0,0,1])
#     legend("Filter1","Filter2","FFT")


# else
#     [h,w] = freqz(Hd1,nrow,Fs);   
#     plot(w,mag2db(abs(h)),"LineWidth",1,"color", [0,1,0])
    
#     [h,w] = freqz(Hd2,nrow,Fs);
#     plot(w,mag2db(abs(h)),"LineWidth",1,"color", [0,0,1])
  
#     [h,w] = freqz(Hd3,nrow,Fs);
#     plot(w,mag2db(abs(h)),"LineWidth",1,"color", [0.25,0.5,0.5])
    
#     legend("Filter1","Filter2","Filter3","FFT")

# end

# yyaxis left
# ylabel("dB")


# hold off
# subplot(2,1,2)
# xlabel("f (MHz)")

# hold on

# plot(f,FFT,"LineWidth",1, "color", [0.75,0.75,0.75]) 

# plot(f,BeamformedSpectralContent,"LineWidth",1.5, "color",[1,0,0])

# if Num_Compounding_Channels == 1
# plot(f,BeamformedSpectralContent1,"LineWidth",1.5, "color",[0,1,0])
# legend('FFT','FFTSG',"Filtered1")


# elseif Num_Compounding_Channels==2
# plot(f,BeamformedSpectralContent1,"LineWidth",1.5, "color",[0,1,0])
# plot(f,BeamformedSpectralContent2,"LineWidth",1.5, "color",[0,0,1])
# legend('FFT','FFTSG',"Filtered1","Filtered2")

# else
# plot(f,BeamformedSpectralContent1,"LineWidth",1.5, "color",[0,1,0])
# plot(f,BeamformedSpectralContent2,"LineWidth",1.5, "color",[0,0,1])
# plot(f,BeamformedSpectralContent3,"LineWidth",1.5, "color",[.25,.5,.5])
# legend('FFT','FFTSG',"Filtered1","Filtered2","Filtered3")


# end

# hold off

# % --- Envelop detection (abs and Hilbert returns amplitude detected)
# Beamformed_env_0=abs(hilbert(Beamformed_Filtered0));
# Beamformed_env_1=abs(hilbert(Beamformed_Filtered1));
# Beamformed_env_2=abs(hilbert(Beamformed_Filtered2));
# Beamformed_env_3=abs(hilbert(Beamformed_Filtered3));

# % --- log compression0
# log_env_0 = Weighting_0*Compress(dB_Range1,Beamformed_env_0);

# % --- log compression1
# log_env_1 = Weighting_1*Compress(dB_Range1,Beamformed_env_1);

# % --- log compression2
# log_env_2 = Weighting_2*Compress(dB_Range2,Beamformed_env_2);

# % --- log compression3
# log_env_3 = Weighting_3*Compress(dB_Range3,Beamformed_env_3);

# if Num_Compounding_Channels == 1
#     Channel_1 = 1; Channel_2=0; Channel_3=0;

# elseif Num_Compounding_Channels==2
#     Channel_1 = 1; Channel_2=1; Channel_3=0;

# else
#     Channel_1 = 1; Channel_2=1; Channel_3=1;

# end

# TempParams.Channel_1 = Channel_1;
# TempParams.Channel_2 = Channel_2;
# TempParams.Channel_3 = Channel_3;

# %Initialize Gain levels
# Gain=50;
# %envelope gain for NonCompounded
# TempParams.Gain1 = Gain;
# %envelope gain for Compounded
# TempParams.Gain2 = Gain;
# TempParams.Num_Compounding_Channels = Num_Compounding_Channels;
# TempParams.gamma = gamma;

# TempParams.Lateral_Line_Count = Lateral_Line_Count;

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

# env_disp = Image(3*TempParams.Gain1,log_env_0,TempParams.Lateral_Line_Count,param.depth,param.gridNum,param.x_ele);

 

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

