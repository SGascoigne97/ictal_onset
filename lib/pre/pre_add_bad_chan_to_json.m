function json_data = pre_add_bad_chan_to_json(json_data,mat_path)
% PRE_ADD_BAD_CHAN_TO_JSON Add list of bad channels in MATLAB workspaces to
% JSON structure.
%
%   json_data = PRE_ADD_BAD_CHAN_TO_JSON(json_data,mat_path) returns a
%   modified JSON structure, json_data, with a new field 'bad_chan_mat'
%   that contains a cell array of previously marked "bad" channels for each
%   iEEG segment. Marked "bad" channels come from the database iEEG MATLAB
%   workspaces, where they are stored as the variable "eeg_bad_chan" in
%   each workspace. The variable mat_path specifies the folder that
%   contains the iEEG MATLAB workspaces.
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% December 2022

arguments
    json_data struct
    mat_path char
end

disp('GETTING PREVIOUSLY MARKED BAD CHANNELS')

% number of iEEG segments
n_segm = length(json_data);

% get file names of all iEEG segments
json_eeg_fn = db_all_json_eeg_fn(json_data);

% add bad channels to json structure
for i=1:n_segm
    
    % load
    load([mat_path '/' json_eeg_fn{i}],'eeg_bad_chan','eeg_channels');
    
    % ensure output is column vector
    if size(eeg_bad_chan,2) > 1
        eeg_bad_chan = eeg_bad_chan';
    end
    
    % save
    json_data(i).bad_chan_mat = eeg_bad_chan;
    
    % number of previously labelled bad channels
    disp(['segment ' num2str(i) ': ' num2str(length(json_data(i).bad_chan_mat)) '/'...
        num2str(length(eeg_channels)) ' channels previously labelled bad'])
end

