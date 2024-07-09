function json_data = meas_all_from_psd(json_data,meas_func,mat_path,...
    preproc_label,psd_label,meas_name,feat_names,varargin)
% MEAS_ALL_FROM_PSD Compute a specified measure from all iEEG segment PSDs
% whose filenames are stored in the JSON structure. Filenames of the new
% measure workspaces are also added to the JSON structure
%
%   INPUTS:
%
%       json_data: JSON structure; must have PSD filenames stored in
%       json_data.fn_analysis.(preproc_label).(psd_label).psd (stored by
%       meas_all_pwelch_psd).
%
%       meas_func: function handle to a function that computes the desired 
%       measure for one iEEG segment. First two inputs of the function must
%       be the segment's PSD array (size n channels x n frequencies) and
%       a vector of the corresponding frequencies. See meas_rel_bp as an
%       example. 
%
%       mat_path: character specifying path to top level folder of the PSD
%       workspaces (i.e., mat files). The output workspaces will also be
%       saved in this folder.
%
%       preproc_label: character specifying the preprocessing options used
%       for the iEEG segments that the PSDs were computed from.
%       
%       psd_label: character specifying the label for the PSDs that will be
%       used to compute the measure. 
%
%       meas_name: character specifying a name for the computed measure;
%       will be used to set folder and structure field names for the
%       exported workspaces.
%
%       feat_names: cell array of names for each dimension of the computed
%       measure; will be stored in the exported workspaces.
%
%       Remaining arguments are passed to meas_func after each PSD array and
%       frequecy vector and should specify any parameters for the computed
%       measure.
%
%   See also MEAS_PWELCH_PSD, MEAS_ALL_PWELCH_PSD, MEAS_REL_BP
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% December 2022

% number of iEEG segments
n_segm = length(json_data);

% get filenames of all PSDs iEEG segments
json_psd_fn = db_all_json_psd_fn(json_data,preproc_label,psd_label);

% get folder structure and workspace names from original file names
[folder_structure,workspace_names] = db_decompose_filenames(json_data);

disp('COMPUTING MEASURE FROM PSDs')

% compute measure for each segment and save
for i=1:n_segm
    disp(json_psd_fn{i})

    % load
    load([mat_path '/' json_psd_fn{i}],...
        'pxx','freq','eeg_channels','x_oid');
    
    % compute measure
    [meas,meas_settings] = meas_func(pxx,freq,varargin{:});
    
    % folder and filename for saving workspace
    save_dir = [folder_structure{i} preproc_label '/' psd_label '/' meas_name '/'];
    save_fn = [save_dir meas_name '_' workspace_names{i}];
    if ~exist([mat_path '/' save_dir],'dir')
        mkdir([mat_path '/' save_dir])
    end
    
    % save filename in json structure
    json_data(i).fn_analysis.(preproc_label).(psd_label).(meas_name) = save_fn;
    
    % save measure
    save([mat_path '/' save_fn],'meas','meas_settings','meas_name','eeg_channels','x_oid','feat_names');
end

