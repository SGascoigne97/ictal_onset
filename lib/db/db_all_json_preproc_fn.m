function json_preproc_fn = db_all_json_preproc_fn(json_data,preproc_label)
% DB_ALL_JSON_PREPROC_FN Get list of the MATLAB workspace filenames for
% previously preprocessed iEEG data for all iEEG segments in the JSON
% structure.
%
% Filenames must be stored in the fn_analysis field in the JSON structure.
%
% Can also use db_all_json_analysis_fn, which this function calls.
%
%   json_preproc_fn = DB_ALL_JSON_PREPROC_FN(json_data,preproc_label)
%   returns a cell array of filenames, json_preproc_fn, that give the paths
%   to the workspaces containing the preprocessed iEEG with label 
%   preproc_label for each iEEG segment in json_data. 
%
%   See also DB_ALL_JSON_ANALYSIS_FN, DB_ALL_JSON_PSD_FN
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% November 2022

arguments
    json_data struct
    preproc_label char
end

% get filenames using generic function
json_preproc_fn = db_all_json_analysis_fn(json_data,preproc_label,'preproc');
