function json_data = pre_find_miss_chan(json_data)
% PRE_FIND_MISS_CHAN Add list of channels not found in JSON channel details
% to JSON structure.
%
%   json_data = PRE_FIND_MISS_CHAN(json_data,mat_path) returns a
%   modified JSON structure, json_data, with a new field 'bad_chan_miss'
%   that contains a cell array of channels that are in the iEEG segment 
%   MATLAB workspace but not the segment's JSON channel details for each
%   iEEG segment. The variable mat_path specifies the folder that
%   contains the iEEG MATLAB workspaces.
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% January 2023

arguments
    json_data struct
end

disp('FINDING CHANNELS MISSING FROM CHANNEL DETAILS')

% number of iEEG segments
n_segm = length(json_data);

% find missing channels and add names to json structure
for i=1:n_segm

    eeg_channels = json_data(i).eeg_channels;
    eeg_channels = strrep(eeg_channels, ' ','');
    
    % check for missing channels
    chan_det_names = json_data(i).channel_details.chan_name;    % channel details names
    chan_det_names = strrep(chan_det_names, ' ','');
    is_member = ismember(eeg_channels,chan_det_names);          % which channels are in channel details
    miss_chan = eeg_channels(is_member==0);                     % missing channel names
    
    % save
    json_data(i).bad_chan_miss = miss_chan;
    
    % number of missing channels
    disp(['segment ' num2str(i) ': ' num2str(length(json_data(i).bad_chan_miss)) '/'...
        num2str(length(eeg_channels)) ' channels missing from channel details'])
end
