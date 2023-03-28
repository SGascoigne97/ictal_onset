data_location = '/media/b5007876/DATA/Data/UCLH_export/imprintDataExport/'; % replace
path_pipeline_exports = [data_location, 'export_ictal'];

% Add paths to all functions required
addpath(genpath('help_functions'))
addpath('lib_biomarkers')
addpath('lib_dataflow')
addpath('sarah_functions')

% List patients with pre-processed data
patients_dir = dir(path_pipeline_exports);
patients = {patients_dir(3:end).name};

%% For each patient, compute onset based on imprint, EI, and PLHG
for pat = 6:10%length(patients)
    patient = patients{pat};
    load(sprintf('%s/%s.mat', data_location, patient));
    pat_data = data_export;
    pat_meta = pat_data(:,1:(end-1));

    % load json_data
    filelist = dir(fullfile(strcat(path_pipeline_exports, "/", patient, "/"), '**/*.*'));  % get list of files and folders in any subfolder
    filelist = filelist(~[filelist.isdir]);  % remove folders from list
    folder = filelist(strcmp({filelist(:).name},'json_data.mat')).folder;
    load(strcat(folder, '/json_data.mat')); % load the json_data file in the folder
    
    % Apply inclusion criteria
    [pat_data, ~, ~] = incl_crit(pat_data);
    % Append patient data to data table (if inclusion criteria are met
    if size(pat_data,1) < 1
        fprintf('Patient %s does not meet inclusion criteria \n', patient)
        continue
    end
    metadata = pat_data(:,1:(end-1));
    
    % Compute onsets
    mkdir('onset_calcs')
    % Compute imprint for all seizures
    [pat_data, metadata, cell_imprint,  sz_count_pat] = calc_imprint(pat_data, metadata);
    
    % Check order of segments in data and cell_imprint are the same
    if any(strcmp(string(pat_data.segment_id), string(cell_imprint.segment_id)) == 0)
        [~, id] = ismember(string(cell_imprint.segment_id), string(pat_data.segment_id));
        % Reorder cell_imprint to match data and metadata
        cell_imprint = cell_imprint(id,:);
    end
    
    % Compute onset using imprint, EI, and PLHG
    onset_output = compute_onset(pat_data, json_data, cell_imprint, sz_count_pat);
    
    % Compute ROI-based onset for imprint, EI, and PLHG
    for met = ["imprint", "EI", "PLHG"]
        [onset_output] = onset_chan_to_roi(pat_data, json_data, onset_output, 'method', met);
    end
    
    % Add patient onset_output to final output table
    if size(onset_output,1) > 0
        if exist('final_output', 'var')
            final_output = [final_output; onset_output];
        else
            final_output = onset_output;
        end
    end
end

%% CODE UPDATED TO HERE

%% TO DO
% - rewrite code so it runs on one patient at a time then affs to output
% table (as code now uses patient-specific json data


% Plot onsets for each patient and save
patients = onset_output.Patient_id;
% for pat = 1%:length(patients)
%     figure(pat)
%     subplot(181)
%     imagesc(onset_output.Labelled_onset{pat,1})
%     title("CLO")
%     subplot(1,8,2:3)
%     imagesc(onset_output.imprint_roi{pat,1})
%     title("IO")
%     subplot(1,8,4:5)
%     imagesc(onset_output.EI_roi{pat,1})
%     title("EIO")
%     subplot(1,8,6:7)
%     imagesc(onset_output.PLHG_roi{pat,1})
%     title("PLHGO")
%     subplot(188)
%     imagesc(onset_output.Resected{pat,1})
%     title("Resected")
%     sgtitle(sprintf("Patient %s (onset computed channel-wise)", patients(pat)))
% %     mkdir(sprintf('figures/roi_second/%s', patients(pat)))
% %     saveas(gcf,sprintf('figures/roi_second/%s/onset', patients(pat)), 'png')
% end

% Save output table 
%save('onset_output.mat', 'onset_output');

