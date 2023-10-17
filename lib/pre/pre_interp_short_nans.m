function ieeg_data_proc = pre_interp_short_nans(ieeg_data,fs,nan_thresh,opts)
% PRE_INTERP_SHORT_NANS Interpolate short sections of missing data in an
% iEEG segment.
%
%   ieeg_data_proc = PRE_INTERP_SHORT_NANS(ieeg_data,fs,nan_thresh)
%   interpolates NaN segments in the iEEG segment ieeg_data (channels x
%   time numeric array) that has sampling frequency fs Hz if the longest 
%   segment is less than or equal to a threshold length (nan_thresh, 
%   specified in seconds). By default, an error will be thrown if any
%   NaN segments exceed the threshold length.
%
%   ieeg_data_proc = PRE_INTERP_SHORT_NANS(...,Name,Value) can be used to  
%   specify the following options:
%       - 'ThrowError' (boolean) specifies whether to thrown an error if
%       any NaN segments are longer than the threshold length. The default
%       (true) is to throw an error. If false, a warning is produced
%       instead and ALL NaNs are interpolated, even if they exceed the
%       threshold length.
%   
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% December 2022
%
% Modified from preproc_handle_nans (GMS, 20 July 2021), which has more
% options for when to interpolate NaNs. 

arguments
    ieeg_data {mustBeNumeric}
    fs (1,1) {mustBeNumeric}
    nan_thresh {mustBeNumeric}
    opts.ThrowError = true
end

% find largest length of NaN in each channel
n_chan = size(ieeg_data,1);
max_nan_by_chan = zeros(n_chan,1);
for i=1:n_chan
    nan_counts = count_nans(ieeg_data(i,:));
    if ~isempty(nan_counts)
        max_nan_by_chan(i) = max(nan_counts);
    end
end

% max across all channels
max_nan = max(max_nan_by_chan);

% convert threshold from seconds to samples, rounding down
thresh_n_samp = floor(nan_thresh*fs);   

% check that largest length of nans is under threshold
if max_nan > thresh_n_samp
    msg = ['Largest length of NaNs is ' num2str(max_nan)...
            ' samples, which exceeds threshold (' num2str(thresh_n_samp) ' samples)'];
    if opts.ThrowError
        % by default, throw error if exceeds threshold
        error(msg)
    else
        warning([msg ' - will still interpolate nans'])
    end
end    

% interpolate
ieeg_data_proc = ieeg_data;
n = size(ieeg_data,2);
t = 1:n;

for i=1:n_chan
    x = ieeg_data(i,:);
    nanx = isnan(x);

    % interpolate
    ieeg_data_proc(i,nanx) = ...
        interp1(t(~nanx),x(~nanx),t(nanx),'linear','extrap');
end

end

function [counts,count_idx] = count_nans(x)
% [counts,count_idx] = count_nans(x,i)
% Counts number of consecutive occurences of NaNs in the vector x.
%
% Counts gives the length of each occurrence of NaNs, in the order that strings 
% of NaNs occur in x
% The length of counts gives the number of strings of NaNs that occur.
% count_idx gives the index of the first occurance of NaN in each string.

% boolean of NaN elements
nanx = isnan(x);

% indices where vec does not contain NaN
idx=find(~nanx);

% pad with zero and the length of the vector + 1 so strings at the
% beginning and end of the vector are captured
idx=[0 idx length(x)+1];

% difference between non-NaN indices
d=diff(idx);
counts=d(d>1) - 1;

bool = [0 nanx];
count_idx=find(diff(bool)==1);

end