function [has_dupl_chan,n_missing_chan,n_missing_bad_chan, json_data] = db_check_mat_channels(json_data)
% DB_CHECK_MAT_CHANNELS Quality check for channel names and bad channel 
% names in iEEG segment MATLAB workspaces. All files listed in JSON data 
% structure are checked.
%
% Checks
%   1) Whether any channel names are duplicated.
%   2) Whether any channel names in the MATLAB workspaces are not present
%   in the segment's channel details table.
%   3) Whether any bad channel names in the MATLAB workspaces are not
%   present in the segment's channel details table.
%
%   [has_dupl_chan,n_missing_chan,n_missing_bad_chan] = 
%   DB_CHECK_MAT_CHANNELS(json_data, mat_path) returns a logical vector of 
%   whether each MATLAB workspace has duplicate channel names 
%   (has_dupl_chan) and numeric vectors of the number of channel names and 
%   bad channel names that are not in the segment's channel details 
%   (n_missing_chan and n_missing_bad_chan, respectively) for each iEEG 
%   segment in the JSON structure json_data. The path to the MATLAB 
%   workspaces is provided in mat_path.
%
% See also DB_ALL_JSON_EEG_FN
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% November 2022

arguments
    json_data struct
end

% number of iEEG segments
n_segm = length(json_data);

% initialise arrays
has_dupl_chan = zeros(n_segm,1);
n_missing_chan = zeros(n_segm,1);
n_missing_bad_chan = zeros(n_segm,1);

% load each workspace's channel names/bad channel names and perform checks
for i=1:n_segm
    
    eeg_channels = json_data(i).eeg_channels;
    eeg_bad_chan = json_data(i).eeg_bad_chan;
    
    % check for duplicates
    eeg_channels = strrep(eeg_channels,' ','');                 % remove spaces
    n_chan = length(eeg_channels);
    has_dupl_chan(i) = length(unique(eeg_channels)) < n_chan;
    
    % check for missing channels
    chan_det_names = json_data(i).channel_details.chan_name;    % channel details names
    chan_det_names = strrep(chan_det_names,' ','');             % remove spaces
    n_in_tbl = sum(ismember(eeg_channels,chan_det_names));      % number of mat channels in channel details table
    n_missing_chan(i) = n_chan - n_in_tbl;                      % number missing from table
    
    % Replace channel details to only include recorded channels
    json_data(i).channel_details = json_data(i).channel_details(ismember(chan_det_names,eeg_channels),:);
    
    % check for missing bad channels
    n_bad_chan = length(eeg_bad_chan);
    n_bad_in_tbl = sum(ismember(eeg_bad_chan,chan_det_names));  % still works for empty cell arrays {} (i.e., no bad chan) - returns 0
    n_missing_bad_chan(i) = n_bad_chan - n_bad_in_tbl; 
    
end

% as logical
has_dupl_chan = logical(has_dupl_chan);

% warn if duplicate channels found
n_dupl = sum(has_dupl_chan);
if n_dupl>0
    warning([num2str(n_dupl) ' iEEG segment MATLAB workspaces have duplicated channel names'])
else
    disp('No iEEG segments with duplicated channel names in their MATLAB workspaces')
end

% warn if missing channels found
n_segm_w_missing = sum(n_missing_chan>0);
if n_segm_w_missing>0
    warning([num2str(n_segm_w_missing) ' iEEG segments have MATLAB workspace channel names not found in JSON channel details'])
else
    disp('All MATLAB workspace channels were found in JSON channel details')
end

% warn if missing bad channels found
n_segm_w_bad_missing = sum(n_missing_bad_chan>0);
if n_segm_w_bad_missing>0
    warning([num2str(n_segm_w_bad_missing) ' iEEG segments have MATLAB workspace bad channel names not found in JSON channel details'])
else
    disp('All MATLAB workspace bad channels were found in JSON channel details')
end
