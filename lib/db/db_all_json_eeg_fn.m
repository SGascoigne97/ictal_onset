function json_eeg_fn = db_all_json_eeg_fn(json_data)
% DB_ALL_JSON_EEG_FN  Create an array of all the iEEG segment database 
% file names (i.e., paths to iEEG mat files) in the iEEG json structure.
%
%   json_eeg_fn = DB_ALL_JSON_EEG_FN(json_data) creates a cell array, 
%   json_eeg_fn, of the file names of all iEEG segemnts in the json 
%   structure json_data, where json_eeg_fn{i} is json_data(i).eeg_fn
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
json_eeg_fn = cell(n_segm,1);
for i = 1:n_segm
    json_eeg_fn{i} = json_data(i).eeg_fn;
end