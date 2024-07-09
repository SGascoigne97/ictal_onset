function json_data = pre_all_detect_bad_chan(json_data,mat_path,varargin)
% PRE_ALL_DETECT_BAD_CHAN Algorithmically detect "bad" (noisy) iEEG
% channels in all iEEG segments in JSON structure.
%
%   json_data = PRE_ALL_DETECT_BAD_CHAN(json_data, mat_path, varargin)
%   modifies the JSON structure json_data by adding a field,
%   'bad_chan_detected', that contains algorithmically detected "bad" iEEG
%   channels, identified using the function pre_detect_bad_chan. The
%   location of the iEEG MATLAB workspaces is provided by mat_path.
%   Any additional arguments must be Name/Value pair arguments that will be
%   pass to pre_detect_bad_chan to specify algorithm settings (see function
%   for details).
%
%   See also PRE_DETECT_BAD_CHAN
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% December 2022

% no argument validation - not recommended for functions with varargin

disp ('ALGORITHMICALLY DETECTING BAD CHANNELS')
% number of iEEG segments
n_segm = length(json_data);

% get file names of all iEEG segments
json_eeg_fn = db_all_json_eeg_fn(json_data);

% load each workspace's iEEG segments and detect bad channels
for i=1:n_segm
    
    % load
    load([mat_path '/' json_eeg_fn{i}],'eeg_data','eeg_channels');
    
    % detect bad channels and add to json structure
    json_data(i).bad_chan_detected = ...
        pre_detect_bad_chan(eeg_data,eeg_channels,json_data(i).eeg_fs,varargin{:});
    
    % number of bad channels found
    disp(['segment ' num2str(i) ': ' num2str(length(json_data(i).bad_chan_detected)) '/'...
        num2str(length(eeg_channels)) ' channels identified as noisy'])
    
end     
    