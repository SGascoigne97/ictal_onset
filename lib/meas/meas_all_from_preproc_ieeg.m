function json_data = meas_all_from_preproc_ieeg(json_data,meas_func,mat_path,...
    preproc_label,meas_name,feat_names,varargin)
% MEAS_ALL_FROM_PREPROC_IEEG Compute a specified measure from all 
% preprocessed iEEG segments whose filenames are stored in the JSON 
% structure. Filenames of  the new measure workspaces are also added to 
% the JSON structure.
%
%   INPUTS:
%
%       json_data: JSON structure; must have preprocessesd EEG filenames
%       stored in json_data.fn_analysis.(preproc_label).preproc (stored by
%       pre_all_preprocess).
%
%       meas_func: function handle to a function that computes the desired 
%       measure for one iEEG segment. First inputs of the function must
%       be the segment's iEEG array (channels x time) and the segment's 
%       sampling frequency in Hz. Function should then output the measure
%       array (channels x features) and a structure of the measure
%       settings. See meas_pac_plv as an example. 
%
%       mat_path: character specifying path to top level folder of the iEEG
%       workspaces (i.e., mat files). The output workspaces will also be
%       saved in this folder.
%
%       preproc_label: character specifying the options used to preprocess
%       the iEEG segments.
%
%       meas_name: character specifying a name for the computed measure;
%       will be used to set folder and structure field names for the
%       exported workspaces.
%
%       feat_names: cell array of names for each dimension of the computed
%       measure; will be stored in the exported workspaces.
%
%       Remaining arguments are passed to meas_func after the iEEG segment
%       time series and sampling frequency; they should specify any 
%       remaining parameters for the computed measure.
%
%   See also PRE_ALL_PREPROCESS, MEAS_ALL_FROM_PSD, MEAS_PAC_PLV
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% February 2023

% number of iEEG segments
n_segm = length(json_data);

% get filenames of all preprocessed iEEG segments
json_preproc_fn = db_all_json_preproc_fn(json_data,preproc_label);

% get folder structure and workspace names from original file names
[folder_structure,workspace_names] = db_decompose_filenames(json_data);

disp('COMPUTING MEASURE FROM PREPROCESSED IEEG')

% compute measure for each segment and save
for i=1:n_segm
    disp(json_preproc_fn{i})

    % load
    load([mat_path '/' json_preproc_fn{i}],...
        'eeg_data','eeg_channels','eeg_fs','x_oid');
    
    % compute measure
    [meas,meas_settings] = meas_func(eeg_data,eeg_fs,varargin{:});
    
    % folder and filename for saving workspace
    save_dir = [folder_structure{i} preproc_label '/' meas_name '/'];
    save_fn = [save_dir meas_name '_' workspace_names{i}];
    if ~exist([mat_path '/' save_dir],'dir')
        mkdir([mat_path '/' save_dir])
    end
    
    % save filename in json structure
    json_data(i).fn_analysis.(preproc_label).(meas_name) = save_fn;
    
    % save measure
    save([mat_path '/' save_fn],'meas','meas_settings','meas_name','eeg_channels','x_oid','feat_names');
end