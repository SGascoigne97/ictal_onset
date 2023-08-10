function ieeg_data_filt = pre_butter_filt(ieeg_data, fs, filter_type, cutoff, final_order)
% PRE_BUTTER_FILT  Filter iEEG data with a zero-phase Butterworth filter.
% Uses z,p,k filter specification output from BUTTER for numerical stability.
%
%   ieeg_data_filt = PRE_BUTTER_FILT(ieeg_data, fs, filter_type, cutoff, final_order)
%   applies a Butterworth filter to a channels x time array of iEEG data, 
%   ieeg_data, that has sampling frequency fs. 
%
%   INPUTS:
%       ieeg_data: matrix of iEEG signals, each row = one channel (size: #
%       channels x # time points)
%
%       fs: sampling frequency of ieeg_data (scalar)
%
%       filter_type: type of filter (same options as MATLAB function butter
%       ('low','high','stop','bandpass')
%
%       cutoff: filter boundary/boundaries (scalar for lowpass or highpass,
%       vector of length 2 for bandstop or bandpass)
%
%       final_order: final order of filter after forwards and backwards
%       filtering is performed by filtfilt (scalar, must be divisible by 2)
%
%   OUTPUTS:
%       ieeg_data_filt: filtered iEEG matrix, same size as ieeg_data
%   
%   See also BUTTER, ZP2SOS, FILTFILT
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% November 2022
% Modified from preproc_butter_filt (GMS, 19 July 2021)

arguments
    ieeg_data double {mustBeNumeric}
    fs (1,1) double 
    filter_type char {mustBeMember(filter_type,{'low','high','stop','bandpass'})}
    cutoff (1,:) {mustBeNumeric,cutoffSizeCheck(cutoff,filter_type)}
    final_order (1,1) double
end

% Nyquist frequency
ny=fs/2; 

% design filter
divide_by = 2; % filter order will be doubled by filtfilt

if mod(final_order,divide_by)==0                        % check desired filter order is divisible by 2
    order = final_order/divide_by;
    [z,p,k] = butter(order, cutoff/ny, filter_type);    % create filter
    [sos,g] = zp2sos(z,p,k);
else
    msg = ['For ' filter_type ' filter, final filter order must be divisible by ' num2str(divide_by)];
    error(msg)
end

% filter iEEG
warning('off','signal:filtfilt:ParseSOS');   % suppress warning about how filter inputs are interpretted
ieeg_data_filt = filtfilt(sos,g,ieeg_data'); 
ieeg_data_filt = ieeg_data_filt';           % transpose back to channels x time
warning('on','signal:filtfilt:ParseSOS');   % turn warning back on so does not affect user's downstream code

end

% validation function for cutoff
function cutoffSizeCheck(cutoff,filter_order)

x = length(cutoff);

switch filter_order
    case {'low','high'}
        if x ~= 1
            error('Filter cutoff must be a scalar for lowpass or highpass filter.')
        end
    case {'stop','bandpass'}
        if x ~= 2
            error('Filter cutoff must be a vector of length 2 for bandstop or bandpass filter.')
        end
end
end
    