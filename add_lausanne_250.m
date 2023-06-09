% Add in onsets for Lausanne 250 atlas too (table already has 120)

data_location = '/media/b5007876/DATA/Data/UCLH_export/imprintDataExportFinal/'; % replace
path_pipeline_exports = [data_location, 'export_ictal'];

chan_to_roi_thresh_type = "count";
chan_to_roi_thresh = 1;

atl = 3; % Lausanne 250

%%
for pat = 6:size(final_output,1)
    patient = final_output.Patient_id{pat};
    fprintf("Patient %s (%d/%d) \n", patient, pat, size(final_output,1))
    pat_onset = final_output(string(final_output.Patient_id) == patient,:);
    sz_count = size(pat_onset.Segment_ids{:},1);

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

    % Compute channel to ROI atlas
    [ch2roi_map_mat,ch_names,roi_names, ch_det_used]=map_to_roi(json_data,pat_data,atl);
    pat_onset.chan_2_roi_matrix_125 = {ch2roi_map_mat};
    pat_onset.roi_names_125 = {roi_names};

     for method = ["clo", "imprint", "EI", "PLHG"]
        onset_mat = pat_onset(:,sprintf("%s_chan", method));
        onset_mat = onset_mat{1,1}{:};
    
        roi_onset_mat = ch2roi_map_mat*onset_mat;
        roi_onset_binary = ch2roi_thresh(roi_onset_mat, ch2roi_map_mat,...
            chan_to_roi_thresh_type, chan_to_roi_thresh);
        % Add this to the onset_output_table
        if method == "clo"
            pat_onset(:,sprintf("%s_roi_125", method)) =  mat2cell(roi_onset_binary,length(roi_names),1);
        else 
            pat_onset(:,sprintf("%s_roi_125", method)) =  mat2cell(roi_onset_binary,length(roi_names),sz_count);
        end
     end

    % Add patient onset_output to final output table
    if exist('final_output_both_atlas', 'var')
        final_output_both_atlas = [final_output_both_atlas; pat_onset];
    else
        final_output_both_atlas = pat_onset;
    end
end