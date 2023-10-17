function [json_data,only_bad_chan,only_bad_chan_fn] = ...
    pre_rm_segm_with_only_bad_chan(json_data,mat_path,bad_chan_field)
% PRE_RM_SEGM_WITH_ONLY_BAD_CHAN Checks whether each segment in the JSON
% structure has any "good" channels (not listed in bad channel field); if
% not, removes segment from the structure and any downstream analysis.
%
%   [json_data,only_bad_chan,only_bad_chan_fn] = 
%   PRE_RM_SEGM_WITH_ONLY_BAD_CHAN(json_data, mat_path,bad_chan_field) 
%   modifies a JSON structure, json_data, so that it only contains iEEG 
%   segments with 1+ "good" channels.  
%   Each segment's channel names are pulled from the MATLAB workspaces
%   found in the directory mat_path. The variable bad_chan_field specifies 
%   the field in json_data that contains the full list of bad channels in 
%   each segment. In addition to a modified json_data, the function also
%   returns a boolean indicating which segments were removed and a cell
%   array of the correspoding MATLAB filenames of those segments.
%
% See also PRE_CONSOLIDATE_BAD_CHAN
% 
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% January 2023 


arguments
    json_data struct
    mat_path char
    bad_chan_field char
end

disp('REMOVING SEGMENTS WITH NO GOOD CHANNELS')
disp(['List of bad channels obtained from field ' bad_chan_field])

% number of iEEG segments
n_segm = length(json_data);

% get file names of all iEEG segments
json_eeg_fn = db_all_json_eeg_fn(json_data);

% initialise array for storing which segments to remove
only_bad_chan = zeros(n_segm,1);

% add bad channels to json structure
for i=1:n_segm
    
    % load
    load([mat_path '/' json_eeg_fn{i}],'eeg_channels');
    
    % check if all eeg_channels are listed in bad channels
    n_bad = sum(ismember(eeg_channels,json_data(i).(bad_chan_field)));
    only_bad_chan(i) = (n_bad == length(eeg_channels)); % 1 if all bad channels
    
end

% ensure logical
only_bad_chan = logical(only_bad_chan);

% filenames
only_bad_chan_fn = json_eeg_fn(only_bad_chan);

% modify json_data 
json_data = json_data(~only_bad_chan);
    
disp(['Removed ' num2str(sum(only_bad_chan)) ' segments'])