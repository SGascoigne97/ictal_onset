function [json_data, rm_bool] = db_rm_segm_no_channel_details(json_data)
% DB_RM_SEGM_NO_CHANNEL_DETAILS Checks for segments with no (empty) or very
% limited (one channel) channel details and removes them from the JSON
% structure.
%
% This step allows downstream formatting of each segment's channel details.
%
%   [json_data, rm_bool] = DB_RM_SEGM_NO_CHANNEL_DETAILS(json_data) removes
%   segments with 0-1 channels in the channel_details field from the JSON
%   structure json_data. The boolean rm_bool indicates which segments were
%   removed.
%
%   See also DB_LOAD_JSON_DATA, DB_FORMAT_JSON_DATA
% 
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% February 2023

% number of iEEG segments
n_segm = length(json_data);

% initialise array for storing segments to remove
rm_bool = zeros(n_segm,1);

% find segments with no/very limited channel details
for i=1:n_segm
    if isempty(json_data(i).channel_details)                    % empty channel details
        rm_bool(i) = 1;
    else
        rm_bool(i) = length(json_data(i).channel_details) <= 1; % 1 or fewer channels
    end
end

% as logical
rm_bool = logical(rm_bool);

% subset JSON structure
json_data = json_data(~rm_bool);