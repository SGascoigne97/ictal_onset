function [onset_output] = add_meta(onset_output, json_data)
% input:
%   - onset_output: % table with seizure onsets for subject
%   - pat_data: table of patient information (from which we will extract
%               metadata)

% output
%   - onset_output: full data table with additional patient information

% Sarah Jane Gascoigne & Yujiang Wang
% 25/05/2023

arguments
    onset_output 
    json_data
end

% Surgery_outcome
% Surgery year
% Outcome year
% Op type

ncol_exist = size(onset_output,2);

treatment_details = json_data(1,:).treatment_details;
treatment_year = str2double(treatment_details.treatment_date.x_date(1:4));

sz_types = strings(size(json_data,1),1);
for sz = 1:size(json_data,1)
    seizure_details = json_data(sz,:).seizure_details;
    if isfield(seizure_details, "sz_type_ilae")
        sz_types(sz) = seizure_details.sz_type_ilae;
    else
        sz_types(sz) = "unknown";
    end
end

if size(json_data,1) == 1
    seizure_details = json_data.seizure_details;
    if isfield(seizure_details, "sz_type_ilae")
        sz_types = seizure_details.sz_type_ilae;
    else
        sz_types = "unknown";
    end
end

onset_output = [onset_output {treatment_details.outcome_ILAE},...
    {treatment_year}, {treatment_details.outcome_year},...
    {treatment_details.op_type}, {cat(1,json_data.eeg_duration)}, {cellstr(sz_types)}];

onset_output.Properties.VariableNames(ncol_exist+(1:6)) = {'Surgery_outcome', ...
    'Surgery year', 'Outcome year', 'Op type', 'duration', 'sz_types'};