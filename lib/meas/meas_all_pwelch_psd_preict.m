function json_data = meas_all_pwelch_psd_preict(json_data,mat_path,preproc_label,psd_label,...
    window,overlap,freq, opts)
% MEAS_ALL_PWELCH_PSD Uses Welch's method to computes and saves power 
% spectral densities for the preprocessed iEEG segments specified by the 
% JSON structure. Filenames of the resulting workspaces will also be stored 
% in the JSON structure.
%
%   INPUTS:
%
%       json_data: JSON structure; must have preprocessed iEEG filenames
%       stored in json_data.fn_analysis.(preproc_label).preproc (stored by
%       pre_all_preprocess).
%
%       mat_path: character specifying path to top level folder of the
%       preprocessed iEEG workspaces (i.e., the mat files). The output 
%       workspaces will also be saved in this folder.
%
%       preproc_label: character specifying the preprocessing options used
%       for the iEEG segments that the PSDs will be computed from.
%       
%       psd_label: character specifying a label for the PSDs. Will be used 
%       to set folder and structure field names for the exported
%       workspaces.
%
%       window: size of window used to compute PSDs, in seconds.
%
%       overlap: size of overlap between consecutive windows, in seconds
%
%       freq: vector of frequencies at which to compute the PSD.
%
%   See also PWELCH, MEAS_PWELCH_PSD, MEAS_REL_BP, MEAS_ALL_FROM_PSD
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% December 2022

arguments
    json_data struct
    mat_path char
    preproc_label char
    psd_label char
    window (1,1) {mustBeNumeric}
    overlap (1,1) {mustBeNumeric}
    freq {mustBeNumeric}
     % optional
    opts.pre_buffer {mustBeNumeric} = 10
end
pre_buffer = opts.pre_buffer;
% number of iEEG segments
n_segm = length(json_data);

% get filenames of all preprocessed iEEG segments
json_preproc_fn = db_all_json_preproc_fn(json_data,preproc_label);

% get folder structure and workspace names from original file names
[folder_structure,workspace_names] = db_decompose_filenames(json_data);

disp('COMPUTING PSDs')

% compute PSD of each segment and save
for i=1:n_segm
    disp(json_preproc_fn{i})

    % load
    load([mat_path '/' json_preproc_fn{i}],...
        'eeg_data','eeg_channels','eeg_fs','x_oid');

    % extract pre-ictal segment (excl 10 seconds closest to seizure onset)
    eeg_pre_dur = json_data(i).eeg_pre;
    eeg_pre_buffer = eeg_pre_dur - pre_buffer;
    pre_data = eeg_data(:,1:(eeg_pre_buffer*eeg_fs));
    
    % compute psd
    [pxx,freq] = meas_pwelch_psd(pre_data,window,overlap,freq,eeg_fs);
    
    % folder and filename for saving workspace
    save_dir = [folder_structure{i} preproc_label '/' psd_label '/psd/'];
    save_fn = [save_dir 'psd_' workspace_names{i}];
    if ~exist([mat_path '/' save_dir],'dir')
        mkdir([mat_path '/' save_dir])
    end
    
    % save filename in json structure
    json_data(i).fn_analysis.(preproc_label).(psd_label).psd = save_fn;
    
    % save settings in function
    meas_settings = struct();
    meas_settings.window = window;
    meas_settings.overlap = overlap;
    meas_settings.freq = freq;
    
    % save psd
    % Save freq in meas_settings so it can easily be saved in more generic
    % functions and also save as separate variable so easily accessible if
    % imported into python
    save([mat_path '/' save_fn],'pxx','freq','eeg_channels','x_oid','meas_settings');
end