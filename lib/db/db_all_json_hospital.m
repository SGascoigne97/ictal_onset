function json_hospital = db_all_json_hospital(json_data)
% DB_ALL_JSON_HOSPITAL  Create an array of all the iEEG segment database 
% hospital names in the iEEG json structure.
%
%   json_hospital = DB_ALL_JSON_HOSPITAL(json_data) creates a cell array, 
%   json_hospital, of the hospital names of the all iEEG segments in the 
%   json structure json_data, where json_hospital{i} is
%   json_data(i).patient_details.Hospital
%
%   See also DB_EXTRACT_JSON_METADATA
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% January 2023

arguments
    json_data struct
end

% number of iEEG segments
n_segm = length(json_data);

% initalise array
json_hospital = cell(n_segm,1);
for i = 1:n_segm
    json_hospital{i} = json_data(i).patient_details.Hospital;
end