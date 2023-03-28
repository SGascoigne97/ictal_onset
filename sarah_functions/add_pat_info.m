function [data_tbl, metadata_tbl] = add_pat_info(data_tbl,  metadata_tbl)
% Add patient information (surgical success, diagnosis, etc) to data)
% 
% input:
%   - data_tbl: full data table

% output
%   - data_tbl: full data table with additional patient information

    arguments
        data_tbl
        metadata_tbl
    end
    
    % Get treatment outcome and patient metadata
    treat_out = readtable('../treatment_outcome.xlsx');
    pat_info = readtable('../patient_meta.xlsx');
    patients = string(unique(data_tbl.patient_id));
    
    % Merge treament outcome and TLE/eTLE into data
    for pat = 1
        data_tbl.treatment_outcome(string(data_tbl.patient_id) == patients{pat} ) = ...
            treat_out(string(treat_out.ID_Jane) == string(patients(pat)),:).treatment_outcome;
        data_tbl.is_tle(string(data_tbl.patient_id) == patients{pat} ) = ...
            string(pat_info(string(treat_out.ID_Jane) == string(patients(pat)),:).is_TLE);
    end
    
    % Create binary variable of good (ILAE 1-2) or bad (ILAE >=3) outcome
    data_tbl.treatment_binary = double(data_tbl.treatment_outcome <=2); % 1 denotes a good outcome
    data_tbl.treatment_binary(isnan(data_tbl.treatment_outcome)) = NaN;
end