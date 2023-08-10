function [ieeg_data_proc,chan_names_proc] = ...
    pre_remove_bad_chan(ieeg_data,chan_names,bad_chan)
% PRE_REMOVE_BAD_CHAN Removes list of "bad" channels from iEEG segment and
% the associated channel names.
%
%   [ieeg_data_proc, chan_names_proc] = PRE_REMOVE_BAD_CHAN(ieeg_data,
%   chan_names, bad_chan) removes "bad" channels (bad_chan, cell array)
%   from the iEEG segment ieeg_data (channels x time numeric array) and the
%   associated cell array of channel names, chan_names. The processed iEEG
%   data and channel names are returned as ieeg_data_proc and
%   chan_names_proc, respectively. Not all bad channels need to be present
%   in the array of channel names.
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% December 2022

arguments
    ieeg_data {mustBeNumeric}
    chan_names cell
    bad_chan cell
end

if size(ieeg_data,1) ~= length(chan_names)
    error('The length of the original list of channels (chan) must be the same as the number of rows in ieeg_data')
end

% find channels that should be removed
[rm_chan_bool] = ismember(strrep(chan_names, ' ', ''), strrep(bad_chan, ' ', ''));

% keep channels that are not in list of bad channels
disp(['Removing ' num2str(sum(rm_chan_bool)) ' channels'])
ieeg_data_proc = ieeg_data(~rm_chan_bool,:);
chan_names_proc = chan_names(~rm_chan_bool);
