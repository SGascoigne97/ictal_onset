% clear all
% close all
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
atlases = [2,3]; %1=scale36,2=60,3=125,4=250
% Choose atlas
% Choose threshold for allowing propagated activity (seconds)

min_sz = 1;
min_sz_dur = 9; % minimum seizure duration for inclusion
chan_to_roi_thresh_type = "count";
chan_to_roi_thresh = 1; % One channel in region is sufficient to include region
window_size = 1; % size of sliding window (in seconds)
window_overlap = 7/8; % Overlap between imprint windows
det = 7; % count of windows beyond first detected onset which will also be considered as onset
rec_type = "sec"; % type of threshold used to validate activity is ictal activity (minimum duration) - can use "sec" or "prop"
rec_thresh = 9; % number of seconds for which activity must persist to be labelled as ictal
mad_thresh = 3;
prop_rec = 0.8;
ict_buffer = 10;

onset_calc_loc = "imprint_ons"; % Specify folder to store imprint values in
% Need a new folder if using a different subset of the data/different data
% as it will load previous save if folder is not empty


%% For each patient, compute onset based on imprint, EI, and PLHG
for pat = 1:length(patients)
    patient = patients{pat};

    if exist(sprintf('%s%s.mat', data_location, patient), 'file')
        load(sprintf('%s%s.mat', data_location, patient));
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
    [pat_data, ~, ~] = incl_crit(pat_data, 'min_sz_count', min_sz, ...
        'sz_type', ["focal", "sg", "subclin", "N"], 'min_sz_duration', min_sz_dur);
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
%%
    if size(pat_data,1) > 0
        % Compute imprint using Mahalanobis distance and MAD to identify outliers 
        [pat_data, cell_imprint,  sz_count_pat] = ...
            calc_imprint_mahal(pat_data, "window_size", window_size,...
            "min_sz_count", min_sz, "folder", onset_calc_loc, "mad_thresh", mad_thresh,...
            "rec_thresh", rec_thresh, "window_overlap", window_overlap);
    else 
        fprintf("Patient %s does not meet inclusion criteria \n", patient)
        continue
    end

    sz_count = size(pat_data,1);
    if sz_count < min_sz
        fprintf("Patient %s does not meet inclusion criteria \n", patient)
        continue
    end

    % Remove seizures from json data that were removed when computing
    % imprint
    json_data = json_data(contains(string(extractfield(cat(2,json_data.x_id),...
        'x_oid')), string(pat_data.segment_id)));

% %%
%     % Save figures (activity detected and preictal noise)
%     for s = 1:size(pat_data,1)
%         % EEGs with imprint and onset highlighted
%         plot_eeg_imprints(pat_data(s,:), cell_imprint(s,:), "save_fig",0,...
%             "fig_visible","on", "imprint_changes_label", "mahal_distance_2", "sz", s,...
%             "save_fig_loc",  sprintf("figures/checking_imprint/%s/", string(patient)))
%         % preictal EEG withpreictal abnormalities highlighted
%         plot_preict_abnormalities(pat_data(s,:), cell_imprint(s,:),... 
%             "save_fig",0,"fig_visible","on", "imprint_changes_label",...
%             "mahal_distance_2","sz",s, "save_fig_loc", sprintf("figures/checking_imprint/%s/", string(patient)))
% % 
%         imprint = cell_imprint(s,:).cell_imprint{:};
%         onset_time = find(sum(imprint),1, 'first');
%         if isempty(onset_time)
%             continue
%         else
%             onset = sum(imprint(:,onset_time+[0:7]),2)>=1;
%         end
%         
%         pre_mahal_mad_mat = cell_imprint.cell_pre_features_mad{s,1};
%         mahal_mad_mat = cell_imprint.cell_madscores.mahal_MAD{s,1};
%         %mahal_mad_movsum = cell_imprint.cell_madscores.mahal_MAD_movsum{s,1};
% 
%         %Onset channel focussed plots
%         onset_channels = find(onset);
%         for ch_ind = 1:length(onset_channels)
%             x_max = size(pre_mahal_mad_mat,2) + size(mahal_mad_mat,2); 
%             chan = onset_channels(ch_ind);
%             eeg_dat = pat_data(s,:).segment_data{:};
%             fig = figure(Visible="off");
%             sgtitle(string(pat_data.segment_channel_labels{1,1}(chan)))
%             subplot(311)
%             eeg_dat_ict = eeg_dat(:,120*512:(120+pat_data(s,:).duration)*512);
%             plot((1:size(eeg_dat,2))/(512/8), eeg_dat(chan,:))
%             xline(960)
%             xlim([480, x_max])
%             set(gca, "XTick", 480:80:x_max, "XTickLabel", -60+(0:10:length(480:80:x_max)*10))
%             title("EEG")
%             subplot(312)
%             plot(1:(size(pre_mahal_mad_mat,2) + size(mahal_mad_mat,2)), [pre_mahal_mad_mat(chan,:), mahal_mad_mat(chan,:)])
%             xlim([480, x_max])
%             set(gca, "XTick", 480:80:x_max, "XTickLabel", -60+(0:10:length(480:80:x_max)*10))
%             yline(mad_thresh)
%             xline(960)
%             title("MAD across time")
%             subplot(313)
%             recruitment_threshold = rec_thresh*8;
%             ms_a=movsum( mahal_mad_mat(chan,:)>=mad_thresh,[0 recruitment_threshold-1],2);%forward looking sum
%             ms_b=movsum( mahal_mad_mat(chan,:)>=mad_thresh,[recruitment_threshold-1 0],2);%backward looking sum
%             ms_c=movsum( mahal_mad_mat(chan,:)>=mad_thresh,[(recruitment_threshold/2)-1 (recruitment_threshold/2)],2); % sum across centre
%             imprint_no_movsum = ms_a >= recruitment_threshold*0.8 | ms_b >= recruitment_threshold*0.8 | ms_c >= recruitment_threshold*0.8 ;
%             plot(1:(size(pre_mahal_mad_mat,2) + size(mahal_mad_mat,2)), [zeros(1,size(pre_mahal_mad_mat,2)), imprint(chan,:)])
%             %plot(1:(size(pre_mahal_mad_mat,2) + size(mahal_mad_mat,2)), [zeros(1, size(pre_mahal_mad_mat,2)), mahal_mad_mat(chan,:)]>=2.5)
%             ylim([-0.5,1.5])
%             xlim([480, x_max])
%             set(gca, "XTick", 480:80:x_max, "XTickLabel", -60+(0:10:length(480:80:x_max)*10))
%             xline(960)
%             title(sprintf("Imprint (MAD thresh = %.1f)", mad_thresh))
% % %             subplot(514)
% % %             plot(1:(size(pre_mahal_mad_mat,2) + size(mahal_mad_movsum,2)),...
% % %                 [pre_mahal_mad_mat(chan,:), mahal_mad_movsum(chan,:)])
% % %             xlim([480, x_max])
% % %             set(gca, "XTick", 480:80:x_max, "XTickLabel", -60+(0:10:length(480:80:x_max)*10))
% % %             xline(960)
% % %             yline(mad_thresh)
% % %             title(sprintf("MAD moving median (%.2f seconds forwards and back)", mov_med_val/8))
% % %             subplot(515)
% % %             plot(1:(size(pre_mahal_mad_mat,2) + size(mahal_mad_mat,2)),...
% % %                 [zeros(1, size(pre_mahal_mad_mat,2)), imprint(chan,:)])
% % %             ylim([-0.5,1.5])
% % %             xlim([480, x_max])
% % %             set(gca, "XTick", 480:80:x_max, "XTickLabel", -60+(0:10:length(480:80:x_max)*10))
% % %             xline(960)
% % %             title(sprintf("Imprint using moving median (MAD thresh = %.1f)", mad_thresh))
% %             saveas(fig, sprintf("figures/checking_imprint/%s/mahal_distance_2/%d_%s.png", ...
% %                     patient, s, string(pat_data.segment_channel_labels{1,1}(chan))))
%         end
%     end

     %%
    % Check order of segments in data and cell_imprint are the same
    if any(strcmp(string(pat_data.segment_id), string(cell_imprint.segment_id)) == 0)
        [~, id] = ismember(string(cell_imprint.segment_id), string(pat_data.segment_id));
        % Reorder cell_imprint to match data and metadata
        cell_imprint = cell_imprint(id,:);
    end

    %% convert channel data to ROI with mapping
    if all(cellfun(@isempty, json_data(1).channel_details.ROIname) == 1)
        fprintf("Patient %s does not have ROI information", patient)
        continue
    end
    
    %% detect onset based on channels, then convert to onset ROIs
    [auto_det_onset] = compute_onset(pat_data, cell_imprint, 'wdw_sz', window_size, 'det', det); 
    % Extract resection and CLO channels
    resect_chan = ismember(pat_data.segment_channel_labels{1},...
        pat_data.resected_5mm{1}); % Determine which channels were resected
    clo_chan = ismember(pat_data.segment_channel_labels{1},...
        pat_data.onset_channels{1}); % Determine which channels were labelled as onset (by clinicians)

    %% Create output table
    onset_output = [{patient}, {pat_data.segment_id}, {pat_data.segment_channel_labels{1}},...
        {clo_chan}, {resect_chan}, auto_det_onset];
    onset_output.Properties.VariableNames = ["Patient_id", "Segment_ids",...
        "channel_names","clo_chan","resected_chan","imprint_chan","EI_chan",...
        "PLHG_chan","when_onset"];
   laus = [72, 120, 250];
    for atl = atlases
        % Remove roi columns from auto_det_onset so they can be recreated using the next atlas
        roi_cols = contains(auto_det_onset.Properties.VariableNames,"_roi");
        auto_det_onset(:, roi_cols) = [];
        % Create channels to regions matrix
        [ch2roi_map_mat,ch_names,roi_names, ch_det_used]=map_to_roi(json_data,pat_data,atl);

        for method = ["imprint", "EI", "PLHG"]
            onset_mat = auto_det_onset(:,sprintf("%s_chan", method));
            onset_mat = onset_mat{1,1}{:};
        
            roi_onset_mat = ch2roi_map_mat*onset_mat;
            roi_onset_binary = ch2roi_thresh(roi_onset_mat, ch2roi_map_mat, ...
                chan_to_roi_thresh_type, chan_to_roi_thresh);
            % Add this to the onset_output_table
            auto_det_onset(:,sprintf("%s_roi", method)) =  mat2cell(roi_onset_binary,length(roi_names),sz_count);
        end
        
        %% Convert labelled onset and resection to ROI
        % Convert from channels to ROI using matrix multiplication (then apply
        % threshold criteria for including region)
        resect_roi = ch2roi_thresh(ch2roi_map_mat*resect_chan, ch2roi_map_mat,...
            "prop", 0.25);
        resect_roi(isnan(resect_roi)) = 0;
        clo_roi = ch2roi_thresh(ch2roi_map_mat*clo_chan, ch2roi_map_mat,...
            chan_to_roi_thresh_type, chan_to_roi_thresh);
        
        % Add onsets to output table
        roi_cols = contains(auto_det_onset.Properties.VariableNames,"_roi");
        roi_atl_tbl = [{roi_names}, {ch2roi_map_mat}, auto_det_onset(:,roi_cols),...
            {clo_roi}, {resect_roi}];
        col_names = ['roi_names', 'chan_2_roi_matrix', ...
            auto_det_onset.Properties.VariableNames(roi_cols), 'clo_roi',...
            'resected_roi'];
        roi_atl_tbl.Properties.VariableNames = append(col_names, sprintf("_%d", laus(atl)));

        onset_output = [onset_output, roi_atl_tbl];
    
    end
    
    %% attach meta data
    onset_output = add_meta(onset_output, json_data);
    
    % Add in parameters 
%     onset_output.det = sprintf("+%d seconds", det);
%     onset_output.thresh_roi = {[chan_to_roi_thresh_type, chan_to_roi_thresh]};

    if ~iscell(onset_output.("Outcome year"))
        onset_output.outcome = onset_output.Surgery_outcome;

    elseif length(onset_output.("Outcome year"){:}) > 1
        outcome_id = onset_output.("Outcome year"){:}-onset_output.("Surgery year") == 1;
        onset_output.outcome = onset_output.Surgery_outcome{:}(outcome_id);
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
%%
save('subquestions/final_output.mat',"final_output")

% save('tables/final_output.mat',"final_output")

% % Remove patients with onsets across all regions in most seizures
% rm_pat = zeros(size(final_output, 1),1);
% for pat = 1:size(final_output,1)
%     ons_all_region = nansum(final_output.imprint_roi_250{pat,1})/size(final_output.imprint_roi_250{pat,1},1) == 1;
%     if sum(ons_all_region)/length(ons_all_region) >= 0.5
%         rm_pat(pat) = 1;
%     end
% end
% final_output_clean = final_output(find(~rm_pat),:);

