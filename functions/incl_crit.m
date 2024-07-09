function [data_tbl, metadata_tbl, sz_count_pat] = incl_crit(data_tbl, opts)
% Relabel patients to match Seizure Severity paper (Gascoigne et al., 2023)

% At the moment, the code only works for one seizure type at a time, this
% can be adjusted to inmclude multiple seizure types 

% input:
%   - data_tbl: full data table
%   - optional inputs
%       - sz_type: seizure type of interest
%       - min_sz_duration: minimum duration for each seizure 
%       - min_sz_count: minimum number of seizures to ahve been recorded per patient

% output
%   - data_tbl: full data table with additional patient information

% outputs are in order of data_tbl, as we assume that input is in order of
%meta_tbl

    arguments
        data_tbl
        opts.sz_type = ["focal", "sg"] % Chosen seizure type
        opts.min_sz_duration (1,1) double = 7 % Minimum duration of seizure to be included (at least 7 seconds are required for PLHG code)
        opts.min_sz_count (1,1) double {mustBeNumeric} = 8 % Minimum number of chosen seizure type for each patient
    end
    
    %fill in optional arguments
    sz_type = opts.sz_type;
    min_sz_duration = opts.min_sz_duration;
    min_sz_count = opts.min_sz_count;
    
    %% Identify patients meeting inclusion criteria

    pat_sz_type = fillmissing(string(data_tbl.ilae_sz_type), 'constant', "N");

    % Remove seizures that are no in the selected subtypes from data_tbl
    incl_data = data_tbl(ismember(pat_sz_type, sz_type),:);

%     if length(sz_type) == 1
%         incl_data = data_tbl(pat_sz_type == sz_type,:); % Seizure type criteria
%     elseif length(sz_type) == 2
%         incl_data = [data_tbl(pat_sz_type == sz_type(1),:);...
%             data_tbl(pat_sz_type == sz_type(2),:)];
%     elseif length(sz_type) == 3
%         incl_data = [data_tbl(pat_sz_type == sz_type(1),:);...
%             data_tbl(pat_sz_type == sz_type(2),:);...
%             data_tbl(pat_sz_type == sz_type(3),:);];
%     end

    incl_data = incl_data(incl_data.duration >= min_sz_duration,:); % Duration criteria

    % Create a table of counts of seizures for all included patients
    sz_count_pat = tabulate(incl_data.patient_id);
    sz_count_pat = sz_count_pat(cell2mat(sz_count_pat(:,2)) >= min_sz_count,:);

    % Remove patients not meeting inclusion criteria from data_tbl
    data_tbl = incl_data(contains(incl_data.patient_id,string(unique(sz_count_pat(:,1)))),:);
    metadata_tbl = data_tbl(:,1:(end-1));


end
