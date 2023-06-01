clear all
close all
data_location = '/media/b5007876/DATA/Data/UCLH_export/imprintDataExportFinal/'; % replace
path_pipeline_exports = [data_location, 'export_ictal'];

% Add paths to all functions required
addpath(genpath('help_functions'))
addpath('lib_biomarkers')
addpath('lib_dataflow')
addpath('sarah_functions')
addpath(genpath('/home/campus.ncl.ac.uk/b5007876/Desktop/Database Code/ieeg-norm-map-pipeline/lib/'))

%%
% List patients with pre-processed data
patients_dir = dir(path_pipeline_exports);
patients = {patients_dir(3:end).name};

% Choose atlas
atl = 2; %1=scale36,2=60,3=125,4=250
% Choose threshold for allowing propagated activity (seconds)
det = 0;

min_sz = 5;
chan_to_roi_thresh_type = "count";
chan_to_roi_thresh = 1; % One channel in region is sufficient to include region

onset_calc_loc = "onset_calcs_overlap"; % Specify folder to store imprint values in
% Need a new folder if using a different subset of the data/different data
% as it will load previous save if folder is not empty

%
% For each patient, compute onset based on imprint, EI, and PLHG
for pat = 1 %:length(patients)
    patient = patients{pat};

    if exist(sprintf('%s/%s.mat', data_location, patient), 'file')
        load(sprintf('%s/%s.mat', data_location, patient));
    else
        fprintf('Patient %s does not have saved data \n', patient)
        continue 
    end
    pat_data = data_export;
    clear data_export
    
    % load json_data
    filelist = dir(fullfile(strcat(path_pipeline_exports, "/", patient, "/"), '**/*.*'));  % get list of files and folders in any subfolder
    filelist = filelist(~[filelist.isdir]);  % remove folders from list
    folder = filelist(strcmp({filelist(:).name},'json_data.mat')).folder;
    load(strcat(folder, '/json_data.mat')); % load the json_data file in the folder
    
    % Apply inclusion criteria
    [pat_data, ~, ~] = incl_crit(pat_data, 'min_sz_count', min_sz);
    % Append patient data to data table (if inclusion criteria are met)
    if size(pat_data,1) < min_sz
        fprintf('Patient %s does not meet inclusion criteria \n', patient)
        continue
    end
    pat_meta = pat_data(:,1:(end-1));

    % Remove seizures from json data that do not meet inclusion criteria
    json_data = json_data(contains(string(extractfield(cat(2,json_data.x_id),...
        'x_oid')), string(pat_data.segment_id)));
    
    % Compute onsets
    mkdir(onset_calc_loc)

    if size(pat_data,1) > 0
        % Compute imprint for all seizures
        [pat_data, pat_meta, cell_imprint,  sz_count_pat] = calc_imprint(pat_data, pat_meta);
    else 
        fprintf("Patient %s does not meet inclusion criteria \n", patient)
        continue
    end

    sz_count = size(pat_data,1);
    if sz_count < min_sz
        fprintf("Patient %s does not meet inclusion criteria \n", patient)
        continue
    end

     %%
    % Check order of segments in data and cell_imprint are the same
    if any(strcmp(string(pat_data.segment_id), string(cell_imprint.segment_id)) == 0)
        [~, id] = ismember(string(cell_imprint.segment_id), string(pat_data.segment_id));
        % Reorder cell_imprint to match data and metadata
        cell_imprint = cell_imprint(id,:);
    end
    
    % Choose atlas
%     for atl = [60, 125]
        % Choose threshold for allowing propagated activity (seconds)
%         for det = [0,1,2]

    %% convert channel data to ROI with mapping
    if all(cellfun(@isempty, json_data(1).channel_details.ROIname) == 1)
        fprintf("Patient %s does not have ROI information", patient)
        continue
    end
    [ch2roi_map_mat,ch_names,roi_names, ch_det_used]=map_to_roi(json_data,pat_data,atl);
    
    %% detect onset based on channels, then convert to onset ROIs
    [auto_det_onset] = compute_onset(pat_data, cell_imprint, 'wdw_sz', 2, 'det', 0);

    %%
    for method = ["imprint", "EI", "PLHG"]
        onset_mat = auto_det_onset(:,sprintf("%s_chan", method));
        onset_mat = onset_mat{1,1}{:};
    
        roi_onset_mat = ch2roi_map_mat*onset_mat;
        roi_onset_binary = ch2roi_thresh(roi_onset_mat, ch2roi_map_mat,...
            chan_to_roi_thresh_type, chan_to_roi_thresh);
        % Add this to the onset_output_table
        auto_det_onset(:,sprintf("%s_roi", method)) =  mat2cell(roi_onset_binary,length(roi_names),sz_count);
    end
    
    %% Convert labelled onset and resection to ROI
    resect_chan = ch_det_used.is_resected5; % Determine which channels were resected
    clo_chan = ch_det_used.is_soz; % Determine which channels were labelled as onset (by clinicians)
    
    % Convert from channels to ROI using matrix multiplication (then apply
    % threshold criteria for including region)
    resect_roi = ch2roi_thresh(ch2roi_map_mat*resect_chan, ch2roi_map_mat,...
        chan_to_roi_thresh_type, chan_to_roi_thresh);
    clo_roi = ch2roi_thresh(ch2roi_map_mat*clo_chan, ch2roi_map_mat,...
        chan_to_roi_thresh_type, chan_to_roi_thresh);
    
    % Add onsets to output table
    onset_output = [{patient}, {pat_data.segment_id}, {ch_names}, {roi_names}, {ch2roi_map_mat}, auto_det_onset,...
        {clo_chan}, {clo_roi}, {resect_chan}, {resect_roi}];
    onset_output.Properties.VariableNames([1:5,13:16]) = ...
        {'Patient_id', 'Segment_ids','channel_names', 'roi_names', 'chan_2_roi_matrix', ...
        'clo_chan', 'clo_roi', 'resected_chan', 'resected_roi'};
    
    %% attach meta data
    onset_output = add_meta(onset_output, json_data);
    
    % Add in parameters 
    onset_output.atl = atl;
    onset_output.det = sprintf("+%d seconds", det);
    onset_output.thresh_roi = {[chan_to_roi_thresh_type, chan_to_roi_thresh]};
    
    % Add patient onset_output to final output table
    if size(onset_output,1) > 0
        if exist('final_output', 'var')
            final_output = [final_output; onset_output];
        else
            final_output = onset_output;
        end
    end

end
%%
% Remove patients with onsets across all regions in most seizures
rm_pat = zeros(size(final_output, 1),1);
for pat = 1:size(final_output,1)
    ons_all_region = nansum(final_output.imprint_roi{pat,1})/size(final_output.imprint_roi{pat,1},1) == 1;
    if sum(ons_all_region)/length(ons_all_region) >= 0.5
        rm_pat(pat) = 1;
    end
end
final_output = final_output(find(~rm_pat),:);
%%
save(sprintf('tables/final_output_updated_code_%d', atl),"final_output")
