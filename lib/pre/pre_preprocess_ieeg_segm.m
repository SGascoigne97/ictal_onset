function [ieeg_data_proc,chan_names_proc,fs_proc] = ...
    pre_preprocess_ieeg_segm(ieeg_data,chan_names,fs,opts)
% PRE_PREPROCESS_IEEG_SEGM Preprocess one iEEG segment according to the
% specified options. Returns the preprocessed iEEG data (channels x time
% array), a cell array of the corresponding channel names, and the data's
% sampling frequency in Hz.
%
%   [ieeg_data_proc,chan_names_proc,fs_proc] = ...
%   pre_preprocess_ieeg_segm(ieeg_data,chan_names,fs,opts)
%
%   INPUTS:
%
%       ieeg_data: iEEG time series as a channels x time numeric array.
%
%       chan_names: cell array of channel names that correspond to the rows
%       in ieeg_data.
%
%       fs: sampling frequency of ieeg_data in Hz.
%
%       Name, Value pair arguments:
%
%           'BadChan': cell array of bad channels to remove (default: {},
%           in which case no channels are removed).
%
%           'InterpNaNsThresh': longest allowed length of NaNs in the data,
%           specified in seconds; <= InterpNaNsThresh will be linearly
%           interpolated. An error will be thrown if longer NaN segments 
%           are present. (default: 0, no NaNs allowed).
%
%           'Reference': how to re-reference the data; (default: 'none', 
%           no referencing; see pre_ref for other options).
%
%           'PassFilterType': string, type of pass filter to perform
%           (default: 'none', no filter; other options are 'low','high', or
%           'bandpass' - see pre_butter_filt for details).
%
%           'PassFilterCutoff': filter cutoff for pass filter; scalar for
%           'low' or 'high', vector of length two for 'bandpass' (see 
%           pre_butter_filt for details). Must be specified if
%           'PassFilterType' is not 'none'.
%
%           'PassFilterOrder': order of pass filter. Must be specified if
%           'PassFilterType' is not 'none'.
%
%           'NotchFilterCutoff': array of n x 2 specifying the filter
%           cutoffs for n bandstop filter (each row corresponding to one
%           filter). If empty (default), no notch filtering is performed.
%
%           'NotchFilterOrder': order of notch filter(s); if a scalar and
%           multiple notch filters are requested, all filters will have the
%           same order. Must be specified if 'Notch FilterCutoff' is not
%           empty. 
%
%           'DownsampleFreq': frequency to downsample to, in Hz
% 
%           'DownsampleMethod': method to use for downsampling (see
%           pre_downsample for options; default 'Interp')
%
%   See also PRE_REMOVE_BAD_CHAN, PRE_INTERP_SHORT_NANS, PRE_REREF, 
%   PRE_BUTTER_FILT, and PRE_DOWNSAMPLE
% 
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% December 2022

arguments
    ieeg_data double {mustBeNumeric}
    chan_names cell
    fs (1,1) {mustBeNumeric}
    opts.BadChan cell = {}
    opts.InterpNaNsThresh (1,1) {mustBeNumeric} = 0 % in seconds
    opts.Reference char = 'none'
    opts.PassFilterType char = 'none'
    opts.PassFilterCutoff = [];
    opts.PassFilterOrder = [];
    opts.NotchFilterCutoff = [];
    opts.NotchFilterOrder = [];
    opts.DownsampleFreq = [];
    opts.DownsampleMethod = 'Interp';
    
end

% Check that necessary options are specified
if ~isequal(opts.PassFilterType,'none')
    if isempty(opts.PassFilterCutoff)
        error('PassFilterCutoff must be specified if PassFilterType is not none')
    end
    if isempty(opts.PassFilterOrder)
        error('PassFilterOrder must be specified if PassFilterType is not none')
    end
end
if ~isempty(opts.NotchFilterCutoff)
    if isempty(opts.NotchFilterOrder)
        error('NotchFilterOrder must be specified if NotchFilterCutoff is not empty')
    end
end    

% 1) Remove bad channels
[ieeg_data_proc,chan_names_proc] = ...
    pre_remove_bad_chan(ieeg_data,chan_names,opts.BadChan);

% 2) Interpolate nans <= threshold length
ieeg_data_proc = ...
    pre_interp_short_nans(ieeg_data_proc,fs,opts.InterpNaNsThresh);

% 3) Re-reference
switch opts.Reference
    case 'none'
        disp('No re-referencing')
    otherwise
        disp(['Re-referencing to ' opts.Reference])
        ieeg_data_proc = pre_reref(ieeg_data_proc,opts.Reference);
end

% 4) Pass Filter (i.e., lowpass, highpass, or bandpass)
switch opts.PassFilterType
    case 'none'
        disp('No pass filter')
    otherwise
        disp([strrep(opts.PassFilterType,'pass','') 'pass filter, '...
            num2str(opts.PassFilterCutoff) ' Hz, ' ...
            num2str(opts.PassFilterOrder) 'th order'])
        ieeg_data_proc = pre_butter_filt(ieeg_data_proc,fs,...
            opts.PassFilterType,opts.PassFilterCutoff,opts.PassFilterOrder);
end

% 5) Notch Filter
n_notch = size(opts.NotchFilterCutoff,1);
if n_notch == 0
    disp('No notch filtering')
else
    % if same order used for multiple notch, extend scalar to vector
    if n_notch > 1 && length(opts.NotchFilterOrder)==1
        opts.NotchFilterOrder = opts.NotchFilterOrder*ones(1,n_notch);
    end
    % notch filter at each frequency
    for i=1:n_notch
        disp(['Notch (bandstop) filter, '...
            num2str(opts.NotchFilterCutoff(i,:)) ' Hz, ' ...
            num2str(opts.NotchFilterOrder(i)) 'th order'])
        ieeg_data_proc = pre_butter_filt(ieeg_data_proc,fs,...
            'stop',opts.NotchFilterCutoff(i,:),opts.NotchFilterOrder(i));
    end
end

% 6) Downsample
if isempty(opts.DownsampleFreq)
    disp('No downsampling performed')
    fs_proc = fs;
elseif opts.DownsampleFreq == fs
    disp(['Already at downsampling frequency (' num2str(opts.DownsampleFreq) ...
        ' Hz); no downsampling required'])
    fs_proc = fs;
else
    disp(['Downsampling from ' num2str(fs) ' to ' num2str(opts.DownsampleFreq) ...
        ' Hz using ' opts.DownsampleMethod ' method'])
    ieeg_data_proc = pre_downsample(ieeg_data_proc,fs,opts.DownsampleFreq,...
        'Method',opts.DownsampleMethod);
    fs_proc = opts.DownsampleFreq;
end
