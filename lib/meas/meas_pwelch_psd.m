function [pxx,freq] = meas_pwelch_psd(ieeg_data,window,overlap,freq,fs)
% MEAS_PWELCH_PSD Compute the power spectral densities (PSDs) of each
% channel of one iEEG segment using Welch's method.
%
%   [pxx,freq] = MEAS_PWELCH_PSD(ieeg_data,window,overlap,freq,fs) returns
%   a 2D array of the channel PSDs, pxx, and a vector of the corresponding
%   frequencies, freq, of the iEEG data ieeg_data with sampling rate (in 
%   Hz) fs. Each row in pxx and ieeg_data corresponds to one channel. To 
%   compute the PSD, the window size (window) and amount of overlap between 
%   windows (overlap) for Welch's method must also be specified in seconds. 
%   The desired frequencies at which to compute the PSD is provided as the 
%   vector freq (which will also be returned as an output). If freq is
%   empty ([]), pwelch's default frequency bins will be used and returned
%   as an output.
%
%   EXAMPLES:
%
%   EXAMPLE 1:
%   
%   % load preprocessed iEEG data (eeg_data) and its sampling frequency
%   % (eeg_fs)
%   load('preproc_interictal70s_segm1.mat','eeg_data','eeg_fs')
%
%   % PSD settings
%   window = 2;         % 2s window
%   overlap = 1;        % 1s overlap between windows
%   freq = 1:0.5:80;    % 1-80 Hz in 0.5 Hz wide bins
%
%   % compute PSD
%   [pxx,freq] = meas_pwelch_psd(eeg_data,window,overlap,freq,eeg_fs);
%
% See also PWELCH, MEAS_ALL_PWELCH_PSD, MEAS_ALL_FROM_PSD, MEAS_REL_BP
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% December 2022

% compute PSD
[pxx,freq] = pwelch(ieeg_data',ceil(fs)*window,ceil(fs)*overlap,freq,fs);

% transpose PSD so rows correspond to channels
pxx = pxx';

