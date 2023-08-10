function json_analysis_fn = db_all_json_analysis_fn(json_data,varargin)
% DB_ALL_JSON_ANALYSIS_FN Get list of the MATLAB workspace filenames that
% contain any previously computed measure for all iEEG segments in the JSON
% structure.
%
% Filenames must be stored in the fn_analysis field in the JSON structure.
%
%   json_analysis_fn = DB_ALL_JSON_ANALYSIS_FN(json_data,...) returns a
%   cell array of filenames, json_analysis_fn, that give the paths to the
%   workspaces containing the requested measure for each iEEG segment in
%   json_data. Each additional argument specifies the folder/field name
%   path to the measure's workspace. 
%   
%   EXAMPLES:
%
%       db_all_json_analysis_fn(json_data,'preUK','preproc') returns the
%       preprocessed iEEG filenames stored in
%       json_data.fn_analysis.preUK.preproc and is equivalent to
%       db_all_json_preproc_fn(json_data,'preUK')
%
%       db_all_json_analysis_fn(json_data,'preUK','psd_Brain22','rel_bp')
%       returns the rel_bp measure filenames stored in
%       json_data.fn_analysis.preUK.psd_Brain22.rel_bp (i.e., relative band
%       power computed from the PSD with label "psd_Brain22" computed from
%       the preprocessed iEEG with label "preUK").
%
%   See also DB_ALL_JSON_PREPROC_FN, DB_ALL_JSON_PSD_FN
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% January 2023
    
% number of iEEG segments
n_segm = length(json_data);

% initalise array
json_analysis_fn = cell(n_segm,1);

% get filenames
for i = 1:n_segm
    json_analysis_fn{i} = getfield(json_data(i),'fn_analysis',varargin{:});
end
