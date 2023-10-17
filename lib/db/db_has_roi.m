function [has_roi,unique_roi] = db_has_roi(json_data,mapping_col,opts)
% DB_HAS_ROI Determine which iEEG segments have data for which ROIs in a
% parcellation.
%
%   [has_roi,unique_roi] = DB_HAS_ROI(json_data,mapping_col) returns a
%   segment x ROI boolean (has_roi) of which iEEG segments in the JSON 
%   structure json_data have data for which ROIs as well as a cell array of
%   the unique ROI IDs in each segment. The integer mapping_col determines
%   which parcellation (i.e., which element in the field ROIids in the
%   segment's channel details) to use.
%
%   has_roi,unique_roi] = DB_HAS_ROI(...,'ROIfield',roi_field) specifies
%   which field in channel details to use for the mapping information
%   (default: 'ROIids'). 
%
%   has_roi,unique_roi] = DB_HAS_ROI(...,'nROI',n_roi) specifies the number
%   of ROIs in the specified parcellation, which wll determine the number
%   of columns in has_roi. If not specified, the number of ROIs is set to
%   the largest ROI ID number across all iEEG segments for the specified 
%   parcellation (which may be smaller than the full number of ROIs if the
%   latter ROIs do not have coverage in the data).
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% January 2023


arguments
    json_data struct
    mapping_col (1,1) {mustBeNumeric}
    opts.ROIfield char = 'ROIids'
    opts.nROI (1,1) = 0
end


% number of iEEG segments
n_segm = length(json_data);

% initialise array for unique_roi
unique_roi = cell(n_segm,1);

% remove channels with no mapping info
% will allow mapping cell array to be converted to a numeric array
json_data = pre_find_chan_with_no_mapping(json_data,opts.ROIfield);
for i=1:n_segm
    rm_bool = ismember(json_data(i).channel_details.chan_name,json_data(i).bad_chan_no_mapping);
    json_data(i).channel_details = json_data(i).channel_details(~rm_bool,:);
end

% get unique list of ROI IDs for each segment
for i=1:n_segm
    chan_tbl = json_data(i).channel_details;
    if ~isempty(chan_tbl)                                   % check that has 1+ channels with mapping info
        roi_ids = chan_tbl.(opts.ROIfield);
        for j=1:length(roi_ids)                             % ensure row vector
            if size(roi_ids{j},1) > 1
                roi_ids{j} = roi_ids{j}';
            end
        end
        roi_ids = cell2mat(roi_ids);                        % make into array
        unique_roi{i} = unique(roi_ids(:,mapping_col));     % select column and get unique list
    end
end    

% if number of ROIs isn't provided, get max ID number
if opts.nROI == 0
    for i=1:n_segm
        opts.nROI = max(opts.nROI,max(unique_roi{i}));
    end
end

% initalise array for has_roi
has_roi = zeros(n_segm,opts.nROI);

% 1 if has ROI
for i=1:n_segm
    has_roi(i,unique_roi{i}) = 1;
end

% logical
has_roi = logical(has_roi);