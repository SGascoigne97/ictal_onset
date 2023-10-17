function json_metadata_path = db_extract_and_save_json_metadata(json_data,mat_path,opts)
% DB_EXTRACT_AND_SAVE_JSON_METADATA Extract and save all specified metadata 
% from the JSON structure.
%
% Metadata is saved as both a table and as individual variables in the
% workspace. 
%
%   json_metadata_path = DB_EXTRACT_AND_SAVE_JSON_METADATA(json_data,mat_path)
%   extracts all default metadata variables for the iEEG segments in the
%   JSON structure json_data and saves the metadata as "json_metadata.mat"
%   in the folder mat_path. The function returns the full path to the saved
%   workspace.
%
%   json_metadata_path = DB_EXTRACT_AND_SAVE_JSON_METADATA(...'MetadataTypes',metadata_types)
%   uses a cell array, metadata_types, to specify which types of metadata 
%   to extract and save. Metadata types must be an option in 
%   db_extract_json_metadata.
%
%   json_metadata_path = DB_EXTRACT_AND_SAVE_JSON_METADATA(...'MetadataFilename',metadata_filename)
%   saves the metadata workspace as metadata_filename. The filename does
%   not need to include the file extension .mat.
%
% See also DB_EXTRACT_JSON_METADATA, DB_EXTRACT_JSON_METADATA_TABLE
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% January 2023

arguments
    json_data struct
    mat_path char
    opts.MetadataTypes = {}; % uses defaults of db_extract_json_metadata_table
    opts.MetadataFilename = 'json_metadata';
end

% get metadata
if isempty(opts.MetadataTypes)
    json_metadata_tbl = db_extract_json_metadata_table(json_data);
else
    json_metadata_tbl = db_extract_json_metadata_table(json_data,opts.MetadataTypes);
end

% save table
% using v7 to facilitate python integration
json_metadata_path = [mat_path '/' opts.MetadataFilename '.mat'];   % full path to file
json_metadata_path = strrep(json_metadata_path,'//','/');           % remove any double / 
json_metadata_path = strrep(json_metadata_path,'.mat.mat','.mat');  % remove any double .mat
save(json_metadata_path,'json_metadata_tbl','-v7');

% also read out all variables and append to workspace
all_var = json_metadata_tbl.Properties.VariableNames;
for i=1:length(all_var)
    eval([all_var{i} '= json_metadata_tbl.(all_var{i});']);
    save(json_metadata_path,all_var{i},'-append');
end
    