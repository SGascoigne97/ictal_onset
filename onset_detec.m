%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Seizure onset detection code associated with 
% "Incomplete resection of the icEEG seizure onset zone is not associated 
% with post-surgical outcomes" (Gascoigne, 2024 submission)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Indicate where icEEG is stored
data_location = 'data/';

% Add paths to all functions required
addpath(genpath('help_functions'))
addpath('lib_biomarkers')
addpath('lib_dataflow')
addpath('functions')
addpath(genpath('lib/'))

%%
% List patients with pre-processed data
subjects_dir = dir(data_location);
subjects = {subjects_dir(3:end).name};
subjects = extractBefore(subjects, "_Sz");

% List all parameters for imprint onset detection
min_sz = 1;
min_sz_dur = 9; % minimum seizure duration for inclusion
window_size = 1; % size of sliding window (in seconds)
window_overlap = 7/8; % Overlap between imprint windows
det = 7; % count of windows beyond first detected onset which will also be 
         % considered as onset
rec_type = "sec"; % type of threshold used to validate activity is ictal
                  % activity (minimum duration) - can use "sec" or "prop"
rec_thresh = 9; % number of seconds for which activity must persist to be 
                % labelled as ictal
mad_thresh = 3; % MAD threshold for detection of activity is 3
prop_rec = 0.8; % proportion of the rec_thresh that is required for
                % activity to be considered in onset 
ict_buffer = 10; % number of seconds before clinically labelled onset that
                 % we begin looking for activity (corrects for potentially 
                 % mislabelling of onset time)

plt_show = true; % change to false if you do not wish to see and save onset 
                 % detection and preictal noise removal visualisations
onset_calc_loc = "imprint_ons/"; % Specify folder to store imprint values 
% Need a new folder if using a different subset of the data or different data
% as it will load previous save if folder is not empty

%% For each patient, compute onset based on imprint, EI, and PLHG
for subj = 1:length(subjects)
    subject = subjects{subj};
    fprintf("\n Subject %s", subject)

    if exist(sprintf('%s%s_Sz_preproc.mat', data_location, subject), 'file')
        load(sprintf('%s%s_Sz_preproc.mat', data_location, subject));
    else
        fprintf('\n Subject %s does not have saved data', subject)
        continue 
    end
    
    % Apply inclusion criteria
    [subj_data, ~, ~] = incl_crit(subj_data, 'min_sz_count', min_sz, ...
        'sz_type', "", 'min_sz_duration', min_sz_dur);
    % Append patient data to data table (if inclusion criteria are met)
    if size(subj_data,1) < min_sz
        fprintf('\n Subject %s does not meet inclusion criteria', subject)
        continue
    end
    pat_meta = subj_data(:,1:(end-1));
    
    % Compute onsets
    mkdir(onset_calc_loc)
%% Use imprint to algorithmically detect seizure activity
    if size(subj_data,1) > 0
        % Compute imprint using Mahalanobis distance and MAD to identify outliers 
        [subj_data, cell_imprint,  sz_count_pat] = ...
            calc_imprint_mahal(subj_data, "window_size", window_size,...
            "min_sz_count", min_sz, "folder", onset_calc_loc, "mad_thresh", mad_thresh,...
            "rec_thresh", rec_thresh, "window_overlap", window_overlap);
    else 
        fprintf("\n Subject %s does not meet inclusion criteria", subject)
        continue
    end

    sz_count = size(subj_data,1);
    if sz_count < min_sz
        fprintf("\n Subject %s does not meet inclusion criteria", subject)
        continue
    end

%% Create and save plots 
    if plt_show 
        for s = 1:size(subj_data,1)
            % Determine if the data has channel labels included
            if any(strcmp("segment_channel_labels",subj_data.Properties.VariableNames))
                channel_labels = subj_data.segment_channel_labels{1,1};
            else
                channel_labels = num2cell(1:size(subj_data.segment_data{s},1));
            end

            % EEGs with imprint and onset highlighted
            plot_eeg_imprints(subj_data(s,:), cell_imprint(s,:), "save_fig",plt_show,...
                "fig_visible","on", "imprint_changes_label", "", "sz", s,...
                "save_fig_loc",  sprintf("figures/imprint_plots/%s/", string(subject)), ...
                "channel_labels", channel_labels)
            % preictal EEG withpreictal abnormalities highlighted
            plot_preict_abnormalities(subj_data(s,:), cell_imprint(s,:),... 
                "save_fig",plt_show,"fig_visible","on", "imprint_changes_label",...
                "","sz",s, "save_fig_loc", sprintf("figures/preict_noise/%s/", ...
                string(subject)), "channel_labels", channel_labels)

            imprint = cell_imprint(s,:).cell_imprint{:};
            onset_time = find(sum(imprint),1, 'first');
            if isempty(onset_time)
                continue
            else
                onset = sum(imprint(:,onset_time+[0:7]),2)>=1;
            end
            
            pre_mahal_mad_mat = cell_imprint.cell_pre_features_mad{s,1};
            mahal_mad_mat = cell_imprint.cell_madscores.mahal_MAD{s,1};
    
            % Onset channel focussed plots
            onset_channels = find(onset);
            for ch_ind = 1:length(onset_channels)
                x_max = size(pre_mahal_mad_mat,2) + size(mahal_mad_mat,2); 
                chan = onset_channels(ch_ind);
                eeg_dat = subj_data(s,:).segment_data{:};
                fig = figure(Visible="off");
                sgtitle(string(channel_labels(chan)))
                subplot(311)
                eeg_dat_ict = eeg_dat(:,120*512:(120+subj_data(s,:).duration)*512);
                plot((1:size(eeg_dat,2))/(512/8), eeg_dat(chan,:))
                xline(960)
                xlim([480, x_max])
                set(gca, "XTick", 480:80:x_max, "XTickLabel", -60+(0:10:length(480:80:x_max)*10))
                title("EEG")
                subplot(312)
                plot(1:(size(pre_mahal_mad_mat,2) + size(mahal_mad_mat,2)), [pre_mahal_mad_mat(chan,:), mahal_mad_mat(chan,:)])
                xlim([480, x_max])
                set(gca, "XTick", 480:80:x_max, "XTickLabel", -60+(0:10:length(480:80:x_max)*10))
                yline(mad_thresh)
                xline(960)
                title("MAD across time")
                subplot(313)
                recruitment_threshold = rec_thresh*8;
                ms_a=movsum( mahal_mad_mat(chan,:)>=mad_thresh,[0 recruitment_threshold-1],2);%forward looking sum
                ms_b=movsum( mahal_mad_mat(chan,:)>=mad_thresh,[recruitment_threshold-1 0],2);%backward looking sum
                ms_c=movsum( mahal_mad_mat(chan,:)>=mad_thresh,[(recruitment_threshold/2)-1 (recruitment_threshold/2)],2); % sum across centre
                imprint_no_movsum = ms_a >= recruitment_threshold*0.8 | ms_b >= recruitment_threshold*0.8 | ms_c >= recruitment_threshold*0.8 ;
                plot(1:(size(pre_mahal_mad_mat,2) + size(mahal_mad_mat,2)), [zeros(1,size(pre_mahal_mad_mat,2)), imprint(chan,:)])
                ylim([-0.5,1.5])
                xlim([480, x_max])
                set(gca, "XTick", 480:80:x_max, "XTickLabel", -60+(0:10:length(480:80:x_max)*10))
                xline(960)
                title(sprintf("Imprint (MAD thresh = %.1f)", mad_thresh))
                saveas(fig, sprintf("figures/imprint_plots/%s/sz%d_chan_%s.png", ...
                    subject, s, string(channel_labels(chan))))
         
            end
        end
    end

     %% Check order of segments in data and cell_imprint are the same
    if any(strcmp(string(subj_data.segment_id), string(cell_imprint.segment_id)) == 0)
        [~, id] = ismember(string(cell_imprint.segment_id), string(subj_data.segment_id));
        % Reorder cell_imprint to match data and metadata
        cell_imprint = cell_imprint(id,:);
    end

    %% detect onset based on channels
    [auto_det_onset] = compute_onset(subj_data, cell_imprint, 'wdw_sz', window_size, 'det', det); 

    %% Create output table
    onset_output = [{subject}, {subj_data.segment_id}, {string(channel_labels)'},...
        auto_det_onset];
    onset_output.Properties.VariableNames = ["Patient_id", "Segment_ids",...
        "channel_names","imprint_chan","EI_chan",...
        "PLHG_chan","when_onset"];
    % Add onset_output to final output table
    if size(onset_output,1) > 0
        if exist('final_output', 'var')
            final_output = [final_output; onset_output];
        else
            final_output = onset_output;
        end
    end
    close all % Closes all figures (this will prevent the computer from crashing due to having too many images open)
end

% Save output table
save('tables/final_output.mat',"final_output")