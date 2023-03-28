function [data_tbl, metadata_tbl] = relab_var(data_tbl,  metadata_tbl)
% Relabel patients to match Seizure Severity paper (Gascoigne et al., 2023)
% 
% input:
%   - data_tbl: full data table

% output
%   - data_tbl: full data table with additional patient information

% outputs are in order of data_tbl, as we assume that input is in order of
%meta_tbl

    arguments
        data_tbl
        metadata_tbl
    end
    
    %% Patient IDs were changed in the paper - see pat_id for conversion
    % Read in CSV with patient IDs 
    patients = unique(data_tbl.patient_id);
%     pat_id = readcell('pat_relab.csv');
%     pat_id = pat_id((2:size(pat_id,1)),(2:3));
%     data_tbl.patient_id = string(data_tbl.patient_id);
%     for pat = 1:length(patients)
%         pat_data = data_tbl(string(data_tbl.patient_id) == string(patients(pat)),:);
%         data_tbl(string(data_tbl.patient_id) == string(patients(pat)),:).patient_id = ...
%             repelem(string(pat_id(string(pat_id)==string(patients(pat)),2)), size(pat_data,1),1);
%     end
    
    % Relabel aura to focal
    sz_label = mergecats(data_tbl.ilae_sz_type,{'focal','aura'});
    sz_label = mergecats(sz_label,{'subclin','subclinical'});
    sz_label = mergecats(sz_label,{'na','nan'});
    data_tbl.ilae_sz_type = sz_label;
    metadata_tbl.ilae_sz_type = sz_label;
    
    % Patient IDs were changed in the paper - see pat_id for conversion
    % Read in CSV with patient IDs 
%     pat_id = readcell('pat_relab.csv');
%     pat_id = pat_id((2:size(pat_id,1)),(2:3));
%     data_tbl.patient_id = string(data_tbl.patient_id);
%     for pat = 1:length(patients)
%         pat_data = data_tbl(string(data_tbl.patient_id) == string(patients(pat)),:);
%         data_tbl(string(data_tbl.patient_id) == string(patients(pat)),:).patient_id = ...
%             repelem(string(pat_id(string(pat_id)==string(patients(pat)),2)), size(pat_data,1),1);
%     end
end