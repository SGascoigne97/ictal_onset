function json_metadata_tbl = db_extract_json_metadata_table(json_data,metadata_types)
% DB_EXTRACT_JSON_METADATA_TABLE Extract all specified metadata from the
% JSON structure and store in a table.
%
%   json_metadata_tbl = DB_EXTRACT_JSON_METADATA_TABLE(json_data) extracts
%   all default metadata types and stores them in the table 
%   json_metadata_tbl. The name of each variable in the table is the same
%   as its metadata type name. 
%
%   json_metadata_tbl = DB_EXTRACT_JSON_METADATA_TABLE(json_data,metadata_types)
%   extracts the metadata types listed in the cell array metadata_types.
%   Each type must be an option for db_extract_json_metadata.
%
%   See also DB_EXTRACT_JSON_METADATA, DB_EXTRACT_AND_SAVE_JSON_METADATA
%
% Gabrielle M. Schroeder
% CNNP Lab, Newcastle University
% January 2023

arguments
    json_data struct
    metadata_types cell = {...
        'x_oid','segm_num',...
        'eeg_duration','Hospital','patID','patientDOB',...
        'patientSex','patientAge_onset','patientAge_exam_est',...
        'exam_id','exam_date','outcome_ILAE1','outcome_ILAE',...
        'outcome_Engel1','outcome_Engel','outcome_year','outcome_date',...
        'outcome1_minus_exam_year',...
        'op_type','op_lobe','op_side','op_pathology'};
end

% number of variables to extract
n_var = length(metadata_types);

% initialise table
json_metadata_tbl = table();

% extract each variable and save in table
for i=1:n_var
    json_metadata_tbl.(metadata_types{i}) = ...
        db_extract_json_metadata(json_data,metadata_types{i});
end


    