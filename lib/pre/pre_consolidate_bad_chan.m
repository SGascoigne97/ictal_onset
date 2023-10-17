function json_data = pre_consolidate_bad_chan(json_data,bad_chan_fields)
% PRE_CONSOLIDATE_BAD_CHAN Consolidates each iEEG segment's lists of "bad"
% channels into a single list in the JSON structure that can then be
% referenced during preprocessing.
%
%   json_data = PRE_CONSOLIDATE_BAD_CHAN(json_data,bad_chan_fields)
%   returns a modified JSON structure, json_data, with a new field 
%   "bad_chan_all" that contains a cell array of all "bad" channels found
%   across multiple fields in the JSON structure (bad_chan_fields).
%   Duplicates are removed from bad_chan_all. 
%   By default, bad_chan_fields are 
%       "bad_chan_filtered' (from pre_filt_chan_by_feature)
%       'bad_chan_detected' (from pre_all_detect_bad_chan)
%       'bad_chan_mat' (from pre_add_bad_chan_to_json)
%       'bad_chan_miss' (from pre_find_miss_chan)
%       'bad_chan_no_mapping (from pre_find_chan_with_no_mapping)
%   A warning is thrown if any bad_chan_fields fields are not present in 
%   the JSON structure.
%
%   See also PRE_FILT_CHAN_BY_FEATURE, PRE_DETECT_BAD_CHAN, 
%   PRE_ALL_DETECT_BAD_CHAN, PRE_ADD_BAD_CHAN_TO_JSON, PRE_FIND_MISS_CHAN,
%   PRE_FIND_CHAN_WITH_NO_MAPPING, VIS_BAD_CHAN
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% December 2022

arguments
    json_data struct
    bad_chan_fields cell = {'bad_chan_filtered',...
        'bad_chan_detected','bad_chan_mat','bad_chan_miss','bad_chan_no_mapping'}
end

% json structure field names
json_fields = fieldnames(json_data);

% any fields not in json structure?
miss_fields = setdiff(bad_chan_fields,json_fields,'stable');

if ~isempty(miss_fields)
    miss_fields = cellfun(@string,miss_fields) + ' ';
    warning(['Following fields are missing and will not be added to final bad channels list: ' miss_fields{:}])
end

% only use fields present in both
bad_chan_fields = intersect(bad_chan_fields,json_fields);
n_fields = length(bad_chan_fields);

% number of iEEG segments
n_segm = length(json_data);

% consolidate bad channels into single list
for i=1:n_segm
    all_bad = {};
    
    % get channel names from each list
    for j=1:n_fields
        all_bad = [all_bad; json_data(i).(bad_chan_fields{j})];
    end
    
    % save 
    json_data(i).bad_chan_all = unique(all_bad);
end


