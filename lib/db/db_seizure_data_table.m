function metadata_renamed = db_seizure_data_table(json_data, path_pipeline_exports, preproc_label)
    % DB_EXTRACT_JSON_METADATA_TABLE_IMPRINT Extract all specified metadata from the
    % JSON structure and store in a table.
    %
    %   json_metadata_tbl = DB_EXTRACT_JSON_METADATA_TABLEICT(json_data,...
    %   path_pipeline_exports, patient_ids, preproc_label)) extracts
    %   all default metadata types and stores them in the table 
    %   json_metadata_tbl. The name of each variable in the table is the same
    %   as its metadata type name. 
    %
    %   See also DB_EXTRACT_JSON_METADATA,
    %   DB_EXTRACT_AND_SAVE_JSON_METADATA, 
    %   DB_EXTRACT_METADATA_TABLE
    %
    % Sarah J. Gascoigne (adapted from Gabrielle M. Schroeder)
    % CNNP Lab, Newcastle University
    % March 2023
    
    arguments
        json_data struct
        path_pipeline_exports string
        preproc_label string

    end
    
    metadata_types = {
            'x_oid','Hospital','patID','patientDOB',...
            'patientSex','patientAge_onset','patientAge_exam_est',...
            'exam_date','outcome_ILAE1','outcome_ILAE','outcome_year',...
            'outcome_date','op_lobe','op_side'};
    % number of variables to extract
    n_var = length(metadata_types);
    
    % initialise table
    json_metadata_tbl = table();
    
    % extract each variable and save in table
    for i=1:n_var
        json_metadata_tbl.(metadata_types{i}) = ...
            db_extract_json_metadata(json_data,metadata_types{i});
    end
    
    % Add additional variables so data can be run through imprint pipeline
    % Change column names to match with previous data structure
    metadata_renamed = renamevars(json_metadata_tbl, ["x_oid", "patID"],...
        ["segment_id", "patient_id"]);
    
    % Add additional columns of data required for imprint pipeline
    %metadata_renamed.VI_status = cat(1, json_data.VI_status);
    
    % Iterate through segments to obtain information stored within structs
    json_preproc_fn = db_all_json_preproc_fn(json_data, preproc_label);
    for seg = 1:length(json_data)
        sz_details = json_data(seg).seizure_details;
        metadata_renamed.start(seg) =  extractfield(sz_details.sz_datetime, 'x_date');
        metadata_renamed.duration(seg) =  extractfield(sz_details, 'sz_duration');
        if isfield(sz_details, 'sz_type_ilae')
            metadata_renamed.ilae_sz_type(seg) = string(sz_details.sz_type_ilae);
        else
            metadata_renamed.ilae_sz_type(seg) = NaN;
        end

        if isfield(sz_details, 'sz_loss_of_awareness')
            metadata_renamed.loss_of_awareness(seg) = int8(sz_details.sz_loss_of_awareness);
        else
            metadata_renamed.loss_of_awareness(seg) = NaN;
        end

        if ~isempty(json_data(seg).eeg_baseline_awake)
            metadata_renamed.baseline_awake(seg) = int8(json_data(seg).eeg_baseline_awake);
        else
            metadata_renamed.baseline_awake(seg) = NaN;
        end
        
        metadata_renamed.op_type(seg) = string(json_data(seg).treatment_details.op_type);
        metadata_renamed.ilae(seg) = {json_data(seg).treatment_details.outcome_ILAE};
        metadata_renamed.ilae_year(seg) = {json_data(seg).treatment_details.outcome_year};
        metadata_renamed.onset_channels(seg) = {json_data(seg).channel_details.chan_name(json_data(seg).channel_details.is_soz == 1)};
        metadata_renamed.resected_3mm(seg) = {json_data(seg).channel_details.chan_name(json_data(seg).channel_details.is_resected3 == 1)};
        metadata_renamed.resected_5mm(seg) = {json_data(seg).channel_details.chan_name(json_data(seg).channel_details.is_resected5 == 1)};
        
        % Seizure data
        clearvars 'eeg_data' 'eeg_channels' 'eeg_fs';
        load(strjoin([path_pipeline_exports '/' json_preproc_fn{seg}], ''),'eeg_data');
         metadata_renamed.segment_data(seg) = {eeg_data}; 
         metadata_renamed.segment_channel_labels(seg) = {json_data(seg).pre_eeg_channels};
         metadata_renamed.segment_fs(seg) = json_data(seg).pre_eeg_fs;
         metadata_renamed.segment_pre(seg) = 120;
         metadata_renamed.segment_post(seg) = 120;
    
    end

    % Reorder table
    metadata_renamed  = metadata_renamed(:, {'segment_id', 'patient_id', ...
        'start', 'duration', 'ilae_sz_type', 'loss_of_awareness',...
        'baseline_awake', 'op_type', 'ilae', 'ilae_year'...
        'onset_channels', 'resected_3mm', 'resected_5mm', 'segment_pre',...
        'segment_post', 'segment_fs', 'segment_channel_labels',...
        'segment_data'});
    
end