function json_id = db_all_json_id(json_data)
% DB_ALL_JSON_ID  Create an array of all the iEEG segment database IDs in 
% the iEEG json structure.
%
%   json_id = DB_ALL_JSON_ID(json_data) creates a cell array, json_id, of
%   all the iEEG segemnts in the json structure json_data, where json_id{i}
%   is json_data(i).x_id.x_oid.
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% November 2022

arguments
    json_data struct
end

% number of iEEG segments
n_segm = length(json_data);

% initalise array
json_id = cell(n_segm,1);
for i = 1:n_segm
    json_id{i} = json_data(i).x_id.x_oid;
end