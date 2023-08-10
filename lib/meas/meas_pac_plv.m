function [pac_plv, pac_plv_settings] = meas_pac_plv(ieeg_data,fs,freq_phase,freq_amplitude,options)
% MEAS_PAC_PLV Computes the phase-amplitude coupling, PAC, between
% phase-modulating frequency bands and amplitude-modulated frequency bands.
% PAC is defined as the phase locking value, PLV, between the time series
% of the phase-modulating frequency and the amplitude envelope of the
% amplitude-modulated frequency band. 
%
% This functio uses notation/terminology from Tort et al. 2010, Measuring
% phase-amplitude coupling between neuronal oscillations of different
% frequencies.
%
%   INPUTS:
%
%       ieeg_data: iEEG array, channels x time
% 
%       fs: sampling fequency in Hz
%
%       freq_phase and freq_amplitude: n x 2 arrays of the phase-modulating
%       and amplitude modulated frequency bands, respectively (first 
%       column provides lower bounds, second column provides upper bounds).
%       The PLV between the time series of frequency freq_phase(i,:) and  
%       the amplitude envelope of frequency freq_amplitude(i,:) is computed
%       for each i from 1 to n.
%
%       Optional arguments:
%           -'FilterOrder' scalar of the filter order to use for bandpass
%           filtering the iEEG time series in each frequency range; default
%           is 10. 
%
%   OUPUTS
%
%       pac_plv: array of PAC PLVs, size n channels x n frequency bands
%
%       pac_plv_settings: structure containing settings used to compute
%       PAC PLVs
%
%   See also MEAS_ALL_FROM_PREPROC_IEEG, HILBERT
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% February 2023

arguments
    ieeg_data double {mustBeNumeric}
    fs (1,1) {mustBeNumeric}
    freq_phase (:,2) {mustBeNumeric}
    freq_amplitude (:,2) {mustBeNumeric}
    options.FilterOrder (1,1) {mustBeNumeric} = 10;
end


% check number of bands match in freq_phase and freq_amplitude
if size(freq_phase,1) ~= size(freq_amplitude,1)
    error('Same number of bands must be provided for phase frequencies and amplitude frequencies')
else
    n_bands = size(freq_phase,1); % number of bands
end

% number of channels
n_chan = size(ieeg_data,1);

% array for holding PAC (in terms of PLV) of each channel in each frequency band pair
pac_plv = zeros(n_chan,n_bands);

% compute PLV
for i=1:n_bands
    
    % filter signal to frequencies of interest; transpose so columns = channel signals
    x_fp = pre_butter_filt(ieeg_data,fs,'bandpass',freq_phase(i,:),...      % phase frequency
        options.FilterOrder)';     
    x_fa = pre_butter_filt(ieeg_data,fs,'bandpass',freq_amplitude(i,:),...  % amplitude frequency
        options.FilterOrder)'; 
    
    % compute phase time series of x_fp, phi_fp
    phi_fp = angle(hilbert(x_fp));
    
    % compute amplitude time series of x_fa, A_fa
    A_fa = abs(hilbert(x_fa));
    
    % compute phase time series of the amplitude of x_fa, phi_A_fa
    phi_A_fa = angle(hilbert(A_fa));
    
    % compute PLV
    pac_plv_i = abs(mean(exp(1i * (phi_fp - phi_A_fa))));
    
    % store
    pac_plv(:,i) = pac_plv_i';
    
end

% save settings
pac_plv_settings = options;
pac_plv_settings.freq_phase = freq_phase;
pac_plv_settings.freq_amplitude = freq_amplitude;




