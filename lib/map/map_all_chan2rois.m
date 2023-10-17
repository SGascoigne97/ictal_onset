function [all_roi_data_avg,all_roi_data_by_chan,all_roi_prop_resected,all_roi_data_paths] = ...
    map_all_chan2rois(json_data,parc_names,mat_path,atlas_path,meas_fields,opts)
% MAP_ALL_CHAN2ROIS Map and save channel-level data of all iEEG segments 
% specified by the JSON structure to ROI-level data in the requested 
% parcellations. One workspace is saved per parcellation. 
%
% Applies the "retains" boolean to remove unwanted ROIs.
%
% INPUTS
% 
%   json_data: JSON structure
%
%   parc_names: cell array sof the names of the parcellations; will be used
%   to get atlas info and label saved workspaces. Must be the same order
%   and length as the number of parcellations in the mapping field.
%
%   mat_path: path to folder containing workspaces of channel-level measure
%   to be mapped.
%
%   atlas_path: path to atlas table; will be imported to get atlas info for
%   specified parcellations.
%
%   meas_fields: cell array of field names specifying path to channel-level
%   measure filenames in JSON data; see db_all_json_analysis_fn
%
%   Optional arguments specified as name/value pairs:
%       - 'ROIfield': variable in channel_details to use for mapping info
%       (default: 'ROIids').
%       - 'SaveMaps': boolean of whether to save mapped data in workspaces
%       in the folder mat_path (default: true).
%
% OUTPUTS
%
%   all_roi_data_avg: cell array, one entry per parcellation. Each entry 
%   contains a 3D array of all ROI-level data, ROI x measure feature
%   x iEEG segment (see map_chan2rois for details).
%
%   all_roi_data_by_chan: cell array, one entry per parcellation. Each
%   entry contains a cell array, size ROI x iEEG segment, of the 
%   channel-level data for each ROI (see map_chan2rois for details).
%
%   all_roi_prop_resected: cell array, one entry per parcellation. Each
%   entry contains a numeric array, size ROI x iEEG segment, of the
%   proportion of channels in each ROI that were resected (see 
%   map_chan2rois for details).
%
%   all_roi_data_paths: cell array of complete paths to each saved
%   workspace within the folder mat_path. The corresponding entries of
%   all_roi_data_avg and all_roi_data_by_chan are saved as roi_data_avg and
%   roi_data_by_chan, respectively, along with the parcellation info and
%   the segment database identifiers. 
%
% See also MAP_CHAN2ROIS, DB_ALL_JSON_ANALYSIS_FN
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% January 2023

arguments
    json_data struct
    parc_names cell
    mat_path char
    atlas_path char
    meas_fields cell
    opts.ROIfield char = 'ROIids';
    opts.SaveMaps {mustBeNumericOrLogical} = true;
end


% number of iEEG segments
n_segm = length(json_data);

% number of parcellations
n_atlas = length(parc_names);

% check that number of parcellations matches number of columns in ROIids
% (just check first segment and first channel)
n_col = length(json_data(1).channel_details.(opts.ROIfield){1});
if n_col ~= n_atlas
    error(['Number of atlas names (' num2str(n_atlas) ') must match number of parcellations in ' opts.ROIfield ' (' num2str(n_col) ').'])
end
clearvars n_col

% load atlas info
load(atlas_path,'atlas_tbl')
all_atlas_tbl = atlas_tbl;
clearvars atlas_tbl

% find each atlas in atlas table
all_atlas_idx = zeros(1,n_atlas);
for i=1:n_atlas
    all_atlas_idx(i) = find(strcmp(parc_names{i},all_atlas_tbl.parc_name));
end

% get number of rois for each atlas
all_n_roi = all_atlas_tbl.n_roi(all_atlas_idx);

% get paths to measure workspaces
json_meas_fn = db_all_json_analysis_fn(json_data,meas_fields{:});

% load measure data of first segment to get number of features
load([mat_path '/' json_meas_fn{1}],'meas')
n_feat = size(meas,2);

% initialise cell arrays for storing mapping
% will split if save in matlab workspace
all_roi_data_by_chan = cell(n_atlas,1);     % channel-level data
all_roi_data_avg = cell(n_atlas,1);         % averaged across channels
all_roi_prop_resected = cell(n_atlas,1);    % proportion of each ROI resected
for i=1:n_atlas
    all_roi_data_by_chan{i} = cell(all_n_roi(i),n_segm);
    all_roi_data_avg{i} = nan(all_n_roi(i),n_feat,n_segm);
    all_roi_prop_resected{i} = nan(all_n_roi(i),n_segm);
end

% map each segment to each parcellation
for j=1:n_segm
    load([mat_path '/' json_meas_fn{j}],'meas','eeg_channels')
    opts.ROIfield = 'ROIids';
    for i=1:n_atlas
        [all_roi_data_avg{i}(:,:,j),all_roi_data_by_chan{i}(:,j),all_roi_prop_resected{i}(:,j)] = ...
            map_chan2rois(meas,eeg_channels,json_data(j).channel_details,i,...
            all_n_roi(i),'ROIfield',opts.ROIfield);
    end
end

% apply "retain" boolean
for i=1:n_atlas
    retain = all_atlas_tbl.retain{all_atlas_idx(i)};
    all_roi_data_avg{i} = all_roi_data_avg{i}(retain==1,:,:);
    all_roi_data_by_chan{i} = all_roi_data_by_chan{i}(retain==1,:);
    all_roi_prop_resected{i} = all_roi_prop_resected{i}(retain==1,:);
end
clearvars retain

% split and save 
all_roi_data_paths = cell(n_atlas,1); % will return empty if not saved
if opts.SaveMaps
    % load measure info from one filename to save in workspaces
    load([mat_path '/' json_meas_fn{1}],'meas_settings','meas_name','feat_names');
    
    % get segment ids
    x_oid = db_all_json_id(json_data);
        
    % save each atlas separately so easily to work with in other
    % programming languages
    for i=1:n_atlas
        
        % name for file
        mat_fn = ['map-' parc_names{i}];
        for j=1:length(meas_fields) % add all preprocessing and measure labels
            mat_fn = [mat_fn '-' meas_fields{j}];
        end
        mat_fn = [mat_fn '.mat'];
        
        % read out variables
        roi_data_avg = all_roi_data_avg{i};
        roi_data_by_chan = all_roi_data_by_chan{i};
        roi_prop_resected = all_roi_prop_resected{i};
        
        atlas_tbl = all_atlas_tbl(all_atlas_idx(i),:);
        retain = atlas_tbl.retain{1};
        atlas_tbl.n_roi = sum(retain);
        atlas_tbl.names{1} = atlas_tbl.names{1}(retain==1);
        atlas_tbl.xyz{1} = atlas_tbl.xyz{1}(retain==1,:);
        atlas_tbl.vol{1} = atlas_tbl.vol{1}(retain==1);
        atlas_tbl.dists{1} = atlas_tbl.dists{1}(retain==1,retain==1);
        atlas_tbl.retain{1} = atlas_tbl.retain{1}(retain==1);
        
        atlas_parc_name = atlas_tbl.parc_name{1};
        atlas_n_roi = atlas_tbl.n_roi;
        atlas_names = atlas_tbl.names{1};
        atlas_xyz = atlas_tbl.xyz{1};
        atlas_vol = atlas_tbl.vol{1};
        atlas_dists = atlas_tbl.dists{1};
        
        % save
        % v7 to make python integration easier
        all_roi_data_paths{i} = [mat_path '/' mat_fn];
        all_roi_data_paths{i} = strrep(all_roi_data_paths{i},'//','/'); % removes any double //
        save(all_roi_data_paths{i},'roi_*','atlas_*','meas_settings','meas_name','x_oid','-v7','feat_names');
        
        clearvars roi_data_avg roi_data_by_chan roi_prop_resected atlas_* mat_fn retain
    end
end

