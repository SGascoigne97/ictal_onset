function [roi_onset_binary] = ch2roi_thresh(roi_onset_mat, ch2roi_map_mat, chan_to_roi_thresh_type, chan_to_roi_thresh)

% Applying threshold for including a region when converting from channels
% input:
%   - roi_onset_mat: matrix denoting the number of channels per region are
%                    included
%   - ch2roi_map_mat: channel to ROI conversion matrix
%   - chan_to_roi_thresh_type: type of threshold (count of channels or
%                              proportion of channels)
%   - chan_to_roi_thresh: threshold to include the region

% output
%   - roi_onset_binary: full data table with additional patient information

% Sarah Jane Gascoigne
% 25/05/2023

arguments
   roi_onset_mat
   ch2roi_map_mat
   chan_to_roi_thresh_type (1,1) string {mustBeMember(chan_to_roi_thresh_type, ["count", "prop"])}
   chan_to_roi_thresh (1,1) double {mustBeNumeric}  % Maybe add a requirement here that if "count", 
                                                    % it must be an integer and if "prop" it must be between 0.1 and 1
end

roi_onset_binary = zeros(size(roi_onset_mat,1), size(roi_onset_mat,2));

for sz = 1:size(roi_onset_mat,2)
    sz_roi_chan_count = roi_onset_mat(:,sz);

    if chan_to_roi_thresh_type == "count"
        roi_onset_binary(sz_roi_chan_count >= chan_to_roi_thresh,sz) = 1;
        roi_onset_binary(sz_roi_chan_count < chan_to_roi_thresh &...
            sz_roi_chan_count<0,sz) = NaN;
        roi_onset_binary(sz_roi_chan_count == 0,sz) = 0;
    else 
        min_chan = ceil(sum(ch2roi_map_mat,2)*chan_to_roi_thresh);
        roi_onset_binary(sz_roi_chan_count >= min_chan,sz) = 1;
        roi_onset_binary(sz_roi_chan_count < min_chan & ...
            sz_roi_chan_count > 0, sz) = NaN;
        roi_onset_binary(sz_roi_chan_count == 0,sz) = 0;

    end
end
