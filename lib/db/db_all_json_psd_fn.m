function json_psd_fn = db_all_json_psd_fn(json_data,preproc_label,psd_label)
% DB_ALL_JSON_PSD_FN Get list of the MATLAB workspace filenames for
% previously computed PSDs for all iEEG segments in the JSON structure.
%
% Filenames must be stored in the fn_analysis field in the JSON structure.
%
% Can also use db_all_json_analysis_fn, which this function calls.
%
%   json_psd_fn = DB_ALL_JSON_PSD_FN(json_data,preproc_label,psd_label)
%   returns a cell array of filenames, json_psd_fn, that give the paths
%   to the workspaces containing the PSDs with label psd_label that were
%   computed from preprocessed iEEG with label preproc_label for each iEEG 
%   segment in json_data. 
%
%   See also DB_ALL_JSON_ANALYSIS_FN, DB_ALL_JSON_PREPROC_FN
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% November 2022

arguments
    json_data struct
    preproc_label char
    psd_label char
end

% get filenames using generic function
json_psd_fn = db_all_json_analysis_fn(json_data,preproc_label,psd_label,'psd');

