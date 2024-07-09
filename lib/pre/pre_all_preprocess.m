function json_data = pre_all_preprocess(json_data,eeg_mat_path,pre_mat_path,preproc_label,opts)
% PRE_ALL_PREPROCESS Preprocess all iEEG segment specified in JSON
% structure.
%
% Note: if "bad" channels are removed, the modified JSON structure may have
% fewer segments than the original one. 
%
%   json_data = PRE_ALL_PREPROCESS(json_data,eeg_mat_path,pre_mat_path,
%   preproc_label,opts)
%
%   INPUTS:
%
%       json_data: JSON structure.
%
%       eeg_mat_path: folder containing MATLAB workspaces with iEEG data to
%       be preprocessed.
%
%       pre_mat_path: folder in which to save MATLAB workspaces with 
%       preprocessed iEEG data.
%
%       preproc_label: label for preprocessing settings; will be used for
%       folder names when workspaces are saved.
%
%       Name, Value pair arguments:
%
%           'RemoveBad': boolean, whether to remove "bad" channels
%           (default: true). If true, will also check for segments that
%           will not have any remaining channels after bad channel removal
%           and remove them from the analysis.
%
%           'BadChanField': string, field in json_data that specifies which
%           channels are "bad" in each segment (default, 'bad_chan_all',
%           created by pre_consolidate_bad_chan)
%
%       Name, Value pair arguments that are the same as in
%       pre_preprocess_ieeg_segm:
%
%           'InterpNaNsThresh': longest allowed length of NaNs in the data,
%           specified in seconds; <= InterpNaNsThresh will be linearly
%           interpolated. An error will be thrown if longer NaN segments 
%           are present. (default: 0, no NaNs allowed).
%
%           'Reference': how to re-reference the data; (default: 'none', 
%           no referencing; see pre_ref for other options).
%
%           'PassFilterType': string, type of pass filter to perform
%           (default: 'none', no filter; other options are 'low','high', or
%           'bandpass' - see pre_butter_filt for details).
%
%           'PassFilterCutoff': filter cutoff for pass filter; scalar for
%           'low' or 'high', vector of length two for 'bandpass' (see 
%           pre_butter_filt for details). Must be specified if
%           'PassFilterType' is not 'none'.
%
%           'PassFilterOrder': order of pass filter. Must be specified if
%           'PassFilterType' is not 'none'.
%
%           'NotchFilterCutoff': array of n x 2 specifying the filter
%           cutoffs for n bandstop filter (each row corresponding to one
%           filter). If empty (default), no notch filtering is performed.
%
%           'NotchFilterOrder': order of notch filter(s); if a scalar and
%           multiple notch filters are requested, all filters will have the
%           same order. Must be specified if 'Notch FilterCutoff' is not
%           empty. 
%
%           'DownsampleFreq': frequency to downsample to, in Hz
% 
%           'DownsampleMethod': method to use for downsampling (see
%           function pre_downsample for options; default 'Interp')
%
%   See also PRE_PREPROCESS_IEEG_SEGM, PRE_REMOVE_BAD_CHAN, 
%   PRE_INTERP_SHORT_NANS, PRE_REREF, PRE_BUTTER_FILT, PRE_DOWNSAMPLE, and
%   PRE_CONSOLIDATE_BAD_CHAN, PRE_RM_SEGM_WITH_ONLY_BAD_CHAN 
% 
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% December 2022

arguments
    json_data struct
    eeg_mat_path char
    pre_mat_path char
    preproc_label char
    opts.RemoveBad = true
    opts.BadChanField = 'bad_chan_all'
    % options below are same as pre_preprocess_ieeg_segm
    opts.InterpNaNsThresh = 0 % in seconds
    opts.Reference char = 'none'
    opts.PassFilterType char = 'none'
    opts.PassFilterCutoff = [];
    opts.PassFilterOrder = [];
    opts.NotchFilterCutoff = [];
    opts.NotchFilterOrder = [];
    opts.DownsampleFreq = [];
    opts.DownsampleMethod = 'Interp';
end

% if bad channels will be removed, first remove segments that will not have 
% any good channels
if opts.RemoveBad
    [json_data,~,~] = pre_rm_segm_with_only_bad_chan(...
        json_data,eeg_mat_path,opts.BadChanField);
end

% number of iEEG segments
n_segm = length(json_data);

% get filenames and specific parts of filenames of all iEEG segments
json_eeg_fn = db_all_json_eeg_fn(json_data);
[folder_structure,workspace_names] = db_decompose_filenames(json_data);

disp('PREPROCESSING')

% preprocess each segment and save
for i=1:n_segm
    disp(json_eeg_fn{i})
    
    eeg_channels = json_data(i).eeg_channels;
    eeg_fs = json_data(i).eeg_fs;
    % load
    load([eeg_mat_path '/' json_eeg_fn{i}], 'eeg_data');
        
    % get bad channels for removal if requested
    if opts.RemoveBad
        bad_chan = json_data(i).(opts.BadChanField);
    else
        bad_chan = {};
    end
    
    % preprocess
    [eeg_data,eeg_channels,eeg_fs] = ...
        pre_preprocess_ieeg_segm(eeg_data,eeg_channels,eeg_fs,...
        'BadChan',bad_chan,...
        'InterpNaNsThresh',opts.InterpNaNsThresh,...
        'Reference',opts.Reference,...
        'PassFilterType',opts.PassFilterType,...
        'PassFilterCutoff',opts.PassFilterCutoff,...
        'PassFilterOrder',opts.PassFilterOrder,...
        'NotchFilterCutoff',opts.NotchFilterCutoff,...
        'NotchFilterOrder',opts.NotchFilterOrder,...
        'DownsampleFreq',opts.DownsampleFreq,...
        'DownsampleMethod',opts.DownsampleMethod);
    
    % save preprocessing options as struct
    pre_opts = opts;
    pre_opts.BadChan = bad_chan;
    
    % segment id to save (anything else? patient?) 
    x_oid = json_data(i).x_id.x_oid;
    
    % folder and filename for saving workspace
    save_dir = [folder_structure{i} preproc_label '/preproc/']; 
    save_fn = [save_dir 'preproc_' workspace_names{i}];
    if ~exist([pre_mat_path '/' save_dir],'dir')
        mkdir([pre_mat_path '/' save_dir])
    end
    
    % save filename in json structure
    json_data(i).fn_analysis.(preproc_label).preproc = save_fn;

    json_data(i).pre_eeg_channels = eeg_channels;
    json_data(i).pre_eeg_fs = eeg_fs;
    
    % save
    save([pre_mat_path '/' save_fn],'eeg_data','eeg_channels','eeg_fs','x_oid','pre_opts');
end
    