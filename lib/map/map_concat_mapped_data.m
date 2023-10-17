function [roi_all, meta_all, json_all] = map_concat_mapped_data(data_dirs,mats_roi)
% MAP_CONCAT_MAPPED_DATA Concatenates together mapped iEEG data from
% different workspaces. Each set of mapped data must have the same atlas
% info (atlas_tbl) and measure settings (meas_name, meas_settings, and
% feat_names). Can be used, for example, to easily analyse data with
% slightly different preprocessing options (e.g., different notch filters).
%
% ROI data arrays as well as the corresponding metadata and JSON data are
% combined. 
%
% Recommend NOT combining data with different channel removal criteria
% (e.g., electrode types).
%
%   [roi_all, meta_all, json_all] = MAP_CONCAT_MAPPED_DATA(data_dirs,mats_roi)
%
%   INPUTS
%
%       data_dirs: cell array; each entry is the directory (including the
%       path from the working directory) that contains the mapped data.
%       Each directory must contain the corresponding workspace in mats_roi
%       as well as the workspaces json_data.mat and json_metadata.mat.
%
%       mats_roi: cell array; each entry is the name of a workspace
%       containing the mapped iEEG data in the corresponding data_dirs
%       array. Workspace name should start with 'map-' and end with '.mat'.
%
%   OUTPUTS
%
%       roi_all: structure containing fields 'roi_data_avg',
%       'roi_data_by_chan', 'roi_prop_resected', and 'x_oid', which are
%       concatenated arrays of those variables from the mats_roi
%       workspaces. Also contains fields 'atlas_tbl', 'meas_name',
%       'meas_settngs', and 'feat_names' to store the corresponding
%       variables from the ROI data workspaces. See map_all_chan2rois for
%       descriptions of the variables.
%
%       meta_all: table of metadata corresponding to the iEEG segments in
%       roi_all; concatenated from the json_metadata.mat files in each
%       data_dirs directory.
%
%       json_all: JSON structure corresponding to the iEEG segments in
%       roi_all; concatenated from the json_data.mat files in each
%       data_dirs directory.
%
%       meas_name, meas_settings, feat_names: measure/feature description 
%       variables from ROI data workspaces
%
%
%   EXAMPLE
%   
%   % map settings - can use to quickly define file names
%   parc = 'scale36';           % parcellation to use
%   psd_name = 'psdBrain22';    % PSD label
%   meas_name = 'rel_bp';       % final measure computed from PSD
%    
%   % directories containing mapped data
%   data_dirs = cell(2,1);
%   data_dirs{1} = 'Documents/pipeline_exports/export_normUK/'; % data for normative map from UK sites
%   data_dirs{2} = 'Documents/pipeline_exports/export_normUS/'; % data for normative maps from US sites
%    
%   % workspaces in the directories containing the iEEG data mapped to the ROI level
%   % note that only the preprocessing labels (preUK vs preUS) differ -
%   % maps are for the same parcellation and measure
%   mats_roi = cell(2,1);
%   mats_roi{1} = ['map-' parc '-preUK-' psd_name '-' meas_name '.mat']; % corresponds to data_dirs{1} - UK data
%   mats_roi{2} = ['map-' parc '-preUS-' psd_name '-' meas_name '.mat']; % corresponds to data_dirs{2} - US data
%
%
%
%   See also MAP_ALL_CHAN2ROIS, DB_EXTRACT_AND_SAVE_JSON_METADATA, 
%   DB_LOAD_JSON_DATA, DB_SAVE_JSON_DATA
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% February 2023


% number of maps to combine
n_maps = length(data_dirs);

% initialise arrays
roi_all = struct();
roi_all.roi_data_avg = [];
roi_all.roi_data_by_chan = [];
roi_all.roi_prop_resected = [];
roi_all.x_oid = [];
meta_all = [];
json_all = [];

% use first map to check that all other files have same atlas and measure
% settings
load([data_dirs{1} '/' mats_roi{1}],'atlas_tbl','meas_settings','meas_name','feat_names')
ref = struct();
ref.atlas_tbl = atlas_tbl; 
ref.meas_settings = meas_settings; 
ref.meas_name = meas_name; 
ref.feat_names = feat_names;
ref_fields = fieldnames(ref);
clearvars atlas_tbl meas_* feat_names

for i=1:n_maps

    % ROI data
    load([data_dirs{i} '/' mats_roi{i}],'roi_data_avg','roi_data_by_chan','roi_prop_resected',...
        'atlas_tbl','meas_name','meas_settings','feat_names','x_oid')

    % checks
    for j=1:length(ref_fields)
        equals_ref = eval(['isequal(ref.(ref_fields{j}),' ref_fields{j} ')']);
        if ~equals_ref
            error([ref_fields{j} ' is not the same across files'])
        end
    end
    
    % concatenate
    roi_all.roi_data_avg = cat(3,roi_all.roi_data_avg,roi_data_avg);
    roi_all.roi_data_by_chan = cat(2,roi_all.roi_data_by_chan,roi_data_by_chan);
    roi_all.roi_prop_resected = cat(2,roi_all.roi_prop_resected,roi_prop_resected);
    roi_all.x_oid = [roi_all.x_oid; x_oid];
    
    % for first map, add atlas and measure info
    if i==1
        roi_all.atlas_tbl = atlas_tbl;
        roi_all.meas_name = meas_name;
        roi_all.meas_settings = meas_settings;
        roi_all.feat_names = feat_names;
    end
    clearvars roi_data* roi_prop_resected atlas_tbl meas_* feat_names x_oid
    

    % Metadata
    load([data_dirs{i} '/json_metadata.mat'],'json_metadata_tbl')
    meta_all = [meta_all; json_metadata_tbl]; clearvars json_metadata_tbl
    
    % JSON data
    load([data_dirs{i} '/json_data.mat'],'json_data')
    json_all = [json_all; json_data]; clearvars json_data
    
end

% check ids of roi data, metadata, and json data match
ids = db_all_json_id(json_all);
if ~isequal(ids,meta_all.x_oid)
    error('Segment IDs of metadata and JSON data do not match')
elseif ~isequal(ids,roi_all.x_oid)
    error('Segment IDs of ROI data and JSON data do not match')
end

