function json_data = pre_find_chan_with_no_mapping(json_data,mapping_field)
% PRE_FIND_CHAN_WITH_NO_MAPPING Adds list of channels with no mapping
% information to JSON structure.
%
% Finds any channels with mapping entries that are empty (either [] or ''). 
%
% Note that channels that have no mapping information because they are not
% present in channel_details are not found by this function; use
% pre_find_miss_chan instead.
%
%   json_data = PRE_FIND_CHAN_WITH_NO_MAPPING(json_data) returns a
%   modified JSON structure, json_data, with a new field 
%   'bad_chan_no_mapping' that contains a cell array of channels that do
%   not having mapping information - i.e., the channel's entry in
%   channel_details.ROIids is empty.
%
%   json_data = PRE_FIND_CHAN_WITH_NO_MAPPING(json_data,mapping_field)
%   specifies the field in channel details that contains the mapping 
%   information. Default field is 'ROIids' if mapping_field is empty or 
%   not provided. 
%
% See also PRE_FIND_MISS_CHAN
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% January 2023

arguments
    json_data struct
    mapping_field = 'ROIids'
end

% default option for mapping_field if empty
if isempty(mapping_field) 
    mapping_field = 'ROIids';
end

disp('FINDING CHANNELS WITH NO MAPPING INFORMATION')
disp(['Mapping field: ' mapping_field])

% number of iEEG segments
n_segm = length(json_data);

% find any channels with no mapping in each segment 
for i=1:n_segm
    
    % segment's channel details
    chan_tbl = json_data(i).channel_details; 
    
    % get names of channels with empty mapping info
    chan_mapping = chan_tbl.(mapping_field);
    is_empty = cellfun(@isempty,chan_mapping);      % boolean of empty entries
    chan_no_mapping = chan_tbl.chan_name(is_empty); % channel names
    
    % add to json_data
    json_data(i).bad_chan_no_mapping = chan_no_mapping;
    
    % number of channels with no mapping data
    disp(['segment ' num2str(i) ': ' num2str(length(json_data(i).bad_chan_no_mapping)) '/'...
        num2str(length(chan_tbl.chan_name)) ' channels in channel details have no mapping info'])
end


