function json_data_path = db_save_json_data(json_data,mat_path,opts)
% DB_SAVE_JSON_DATA Save MATLAB JSON structure in a MATLAB workspace.
% Recommended after editing JSON structure. Will overwrite any file with
% the same name. JSON structure will be saved as variable "json_data"
% unless another filename is specified.
%
% Recommend saving in same folder as corresponding pipeline analysis 
% exports.
%
%   json_data_path = DB_SAVE_JSON_DATA(json_data,mat_path) saves the JSON
%   structure json_data as 'json_data.mat' in the folder mat_path and
%   returns the full path to the new workspace.
%
%   json_data_path = DB_SAVE_JSON_DATA(...,'JSONfilename',json_filename) 
%   saves the JSON structure json_data using the specified filename, 
%   json_filename. json_filename does not need to include the .mat
%   extension.
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% January 2023

arguments
    json_data struct
    mat_path char
    opts.JSONfilename char = 'json_data';
end
[folder_structure, ~] = db_decompose_filenames(json_data);

if ~exist([mat_path '/' folder_structure{1}],'dir')
        mkdir([mat_path '/' folder_structure{1}])
end

json_data_path = [mat_path '/' folder_structure{1} opts.JSONfilename '.mat'];   % full path to file
json_data_path = strrep(json_data_path,'//','/');           % remove any double / 
json_data_path = strrep(json_data_path,'.mat.mat','.mat');  % remove any double .mat 
save(json_data_path,'json_data');                           % save