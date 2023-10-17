function bad_chan = pre_detect_bad_chan(ieeg_data,chan_name,fs,opts)
% PRE_DETECT_BAD_CHAN Algorithmically detect "bad" (noisy) iEEG channels
% by finding channels with outlier range and/or variance relative to the
% other channels. Two rounds of detection are perfomed; the first round is
% before preprocessing with (by default) less stringent detection
% thresholds, and the second is after basic preprocessing.
%
%   bad_chan = PRE_DETECT_BAD_CHAN(ieeg_data,chan_name,fs) returns a cell
%   array of "bad" channel names, bad_chan, which are a subset of the full 
%   cell array of  channel names, chan_name. Each channel name is a label
%   for a row in the iEEG data matrix, ieeg_data (size n channels x n time
%   points). The sampling frequency of ieeg_data is given by the scalar fs.
%   
%   bad_chan = PRE_DETECT_BAD_CHAN(...Name,Value) specifies settings for
%   outlier detection and preprocessing:
%       - 'VarThresh1' and 'RangeThresh1' (scalars) specify the first  
%       thresholds for outlier detection based on sigal variance and range,
%       respectively. Default for both is 16.
%       - 'VarThresh2' and 'RangeThresh2' (scalars) specify the second  
%       thresholds for outlier detection based on sigal variance and range,
%       respectively. Default for both is 12.
%       - 'BandpassFreq' (vector of length two) and 'BandpassOrder' 
%       (scalar, even) specify the settings for the band pass filter 
%       applied before the second round of outlier detection (see
%       PRE_BUTTER_FILT). Defaults are [1 100] and 4, respectively.
%       - 'NotchType' (character) and 'NotchOrder' (scalar, even) specify 
%       the settings for a notch (band stop) filter applied before the 
%       second round of outlier detection. NotchType can be 'UK' to remove 
%       UK line noise ([49 51] band stop filter), 'US' to remove US line 
%       noise ([59 61] band stop filter), or 'none' to skip the notch 
%       filtering step. See PRE_BUTTER_FILT fo details. Defaults are 'none' 
%       and 4, respectively. 
%
%
%   See also PRE_BUTTER_FILT, ISOUTLIER
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% November 2022
% Based on preproc_noisy_chan_detection.m and Chris P's RAM data pipeline; 
% modified to exclude NaNs from variance and range calculations
%
% TODO: filtering will not work with NaNs, so only first round will detect
% bad channels if signals contain NaNs - could use PRE_INTERP_SHORT_NANS to
% handle.

arguments
    ieeg_data {mustBeNumeric}
    chan_name cell 
    fs (1,1) {mustBeNumeric}
    opts.VarThresh1 (1,1) {mustBeNumeric} = 16
    opts.VarThresh2 (1,1) {mustBeNumeric} = 12
    opts.RangeThresh1 (1,1) {mustBeNumeric} = 16
    opts.RangeThresh2 (1,1) {mustBeNumeric} = 12;
    opts.BandpassFreq (1,2) {mustBeNumeric} = [1 100]
    opts.BandpassOrder (1,1) {mustBeNumeric} = 4
    opts.NotchType char {mustBeMember(opts.NotchType,{'UK','USA','none'})} = 'none'
    opts.NotchOrder (1,1) {mustBeNumeric} = 4
end

% warning if second set of thresholds are not more stringent than the first
if opts.VarThresh1 <= opts.VarThresh2
    warning('Second variance threshold is not stricter than first variance threshold.')
end
if opts.RangeThresh1 <= opts.RangeThresh2
    warning('Second range threshold is not stricter than first range threshold.')
end


%%% first round of exclusions - lenient thresholds

% exclude outliers in terms of variance (lenient exclusion first)
var_eeg = var(ieeg_data,0,2,'omitnan');
exclV = isoutlier(var_eeg, 'median','ThresholdFactor',opts.VarThresh1);

% exclude outliers in terms of min-max range (lenient exclusion first)
range_eeg = max(ieeg_data,[],2,'omitnan') -  min(ieeg_data,[],2,'omitnan');
exclR = isoutlier(range_eeg, 'median','ThresholdFactor',opts.RangeThresh1);

% first set of channels to exclude
excl1= exclV | exclR ;


%%% preprocessing

% warning if signal contains any nans
if sum(isnan(ieeg_data(:))) > 0
    warning('iEEG data contains NaNs; 2nd round of outlier detection will not detect bad channels that contain NaNs')
end

% common average excluding first set
carmedian=median(ieeg_data(~excl1,:));

% apply common average to all
ieeg_data=ieeg_data-carmedian;

% bandpass filter
ieeg_data = pre_butter_filt(ieeg_data, fs, 'bandpass', opts.BandpassFreq, opts.BandpassOrder);

% notch filter
switch opts.NotchType
    case 'UK'   % 50 Hz notch
       ieeg_data = pre_butter_filt(ieeg_data, fs, 'stop', [49 51], opts.NotchOrder); 
        
    case 'US'   % 60 Hz notch
        ieeg_data = pre_butter_filt(ieeg_data, fs, 'stop', [59 61], opts.NotchOrder);
        
    case 'none' % no notch filter
        
end

% replace the exclusions with nans
ieeg_data(excl1,:) = NaN;


%%% second round of exlusions - can use more stringent thresholds

% variance outliers
var_eeg = var(ieeg_data,0,2,'omitnan');
exclV = isoutlier(var_eeg, 'median', 'ThresholdFactor',opts.VarThresh2);

% range outliers
range_eeg = max(ieeg_data,[],2,'omitnan') -  min(ieeg_data,[],2,'omitnan');
exclR = isoutlier(range_eeg, 'median', 'ThresholdFactor',opts.RangeThresh2);    

% finalise the set of channels to exlude
excl2 = excl1 | exclV | exclR;

% output the names of the channels to exclude
bad_chan = chan_name(excl2);
