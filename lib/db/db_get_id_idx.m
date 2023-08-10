function id_idx = db_get_id_idx(id,json_data)
% DB_GET_ID_IDX  Find the index of an iEEG segment in the exported database
% JSON structure using the iEEG segment's database ID.
%
%   id_idx = DB_GET_ID_IDX(id,json_data) returns the index, id_idx, in the
%   database structure json_data based on the iEEG segment's database ID, 
%   id. 
%
%   See also DB_ALL_JSON_ID
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% November 2022

arguments
    id (1,:) char
    json_data struct
end

% get database IDs of all iEEG segments in database structure
all_id = db_all_json_id(json_data);

% find index of specified ID
id_idx = find(strcmp(id,all_id));

% check that exactly one index was found
if isempty(id_idx)
    error('iEEG segment ID was not found in database structure')
elseif length(id_idx)>1
    error('Multiple entries with iEEG segment ID were found in database structure')
end