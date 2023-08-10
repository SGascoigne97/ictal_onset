function [folder_structure,workspace_names] = db_decompose_filenames(json_data)
% DB_DECOMPOSE_FILENAMES Extracts filename structure from database iEEG
% workspaces. 
%
% Folder and workspace names can then be used to save downstream
% workspaces.
%
%   [folder_structure,workspace_names] = DB_DECOMPOSE_FILENAMES(json_data)
%   returns the folder structure (excluding "eeg/") and the iEEG workspace
%   names of all iEEG segments in the JSON structure json_data. 
%
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% December 2022

% get complete filenames
json_eeg_fn = db_all_json_eeg_fn(json_data);
n = length(json_eeg_fn);

% find first index of workspace name
expression = '[\w-]*.mat';
idx = regexpi(json_eeg_fn,expression);

% extract parts of the file name
workspace_names = cell(n,1);
folder_structure = cell(n,1);
for i=1:n
    workspace_names{i} = json_eeg_fn{i}(idx{i}:end);
    folder_structure{i} = json_eeg_fn{i}(1:idx{i}-5); % also removes eeg/ 
end


