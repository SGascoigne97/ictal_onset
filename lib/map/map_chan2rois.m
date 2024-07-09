function [roi_data_avg,roi_data_by_chan,roi_prop_resected] = ...
    map_chan2rois(chan_data,chan_names,chan_tbl,parc_col,n_roi,opts)
% MAP_CHAN2ROIS Map channel-level data of one iEEG segment to ROI-level 
% data by averaging data across channels.
%
% Only maps channel data to one ROI parcellation - will need to call
% repeatedly to map to multiple parcellations. 
%
% Does not apply any "retain" booleans to restrict ROIs.
%
% Channels without ROI mapping must be removed first during preprocessing
% (see pre_find_chan_with_no_mapping). 
%
% INPUTS
%
%   chan_data: matrix of channel data, channels x features (e.g., band
%   power in different frequency bands).
%
%   chan_names: cell array of channel names that correspond to rows of
%   chan_data.
%
%   chan_tbl: segment's channel_details table (from JSON structure) that
%   contains mapping info (by default, uses table variable "ROIids")
%
%   parc_col: integer specifying which mapping to use in the 
%   channel_details mapping variable (e.g., parc_col = 2 uses 2nd column of
%   ROI IDs)
%
%   n_roi: number of ROIs in the full parcellation (ensures export is the
%   correct size)
%
%   Optional arguments specified as name/value pairs:
%       - 'ROIfield': variable in channel_details to use for mapping info
%       (default: 'ROIids')
%       - 'ResectedField': variable in channel_details to use for resection
%       info (default: 'is_resected5')
%
%
% OUTPUTS
%
%   roi_data_avg: matrix of mapped data, ROIs x features. Will contain NaNs
%   for ROIs that have no corresponding data.
%
%   roi_data_by_chan: cell array of channel data for each ROI. For ROI i
%   containing n channels that have r features, roi_data_by_chan{i} will be
%   size n x r. roi_data_avg(i,:) = mean(roi_data_by_chan{i},1)
%
%   roi_prop_resected: numeric vector of proportion of each ROI resected.
%   Entry i will be NaN if any channels in ROI i are missing resection
%   info or if ROI i has no channels in the segment.
%
% See also MAP_ALL_CHAN2ROIS
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% January 2023
%
%
% Development notes: 
% 
% Would probably be faster to use matrix multiplication to map data, but
% this approach also exports channel-level data organised by ROI.
%
% Could add average type (e.g., mean, median) as optional argument. 

arguments
    chan_data {mustBeNumeric}
    chan_names cell
    chan_tbl table
    parc_col (1,1) {mustBeNumeric}
    n_roi (1,1) {mustBeNumeric}
    opts.ROIfield char = 'ROIids'
    opts.ResectedField char = 'is_resected5';
end
    
% number of features in channel data
n_feat = size(chan_data,2);

% initialise arrays
roi_data_avg = nan(n_roi,n_feat);
roi_data_by_chan = cell(n_roi,1);
roi_prop_resected = nan(n_roi,1);

% rows of channel table that correspond to chan_names
[~,~,idx] = intersect(chan_names,chan_tbl.chan_name,'stable');
if ~isequal(chan_names,chan_tbl.chan_name(idx))
    error('Channel table names do not match measure channel names')
end

% roi ids of channels
roi_ids = chan_tbl.(opts.ROIfield)(idx);    % sort/subset so in same order as chan_names
for i=1:length(roi_ids)                     % ensure row vector
    if size(roi_ids{i},1) > 1 
        roi_ids{i} = roi_ids{i}'; 
    end
end
roi_ids = cell2mat(roi_ids);                % make into array
roi_ids = roi_ids(:,parc_col);              % select desired column 

% resection info of channels
is_resected = chan_tbl.(opts.ResectedField)(idx); % sort/subset so in same order as chan_names

% map
for i=unique(roi_ids)' % only need to loop through IDs present in segment

    % get channels with roi id
    chan_bool = roi_ids == i;
    
    % extract data
    roi_data_by_chan{i} = chan_data(chan_bool,:); 
    
    % average
    roi_data_avg(i,:) = mean(roi_data_by_chan{i},1);
    
    % proportion resected
    roi_prop_resected(i) = mean(is_resected(chan_bool)); 
    
end
    
    
