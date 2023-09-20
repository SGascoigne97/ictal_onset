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
window_size = 1;
window_overlap = 7/8;
mad_thresh = 2.5;

for pat = 1%:length(patients)
    
    patient = patients{pat};
    fprintf(" \n Patient %s (%d)", patient, pat)

    basefolder = sprintf('example_onset_detection/%s', extractAfter(patient, 'UCLH'));
    
    load(sprintf('%s%s.mat', data_location, patient));
    data_tbl = data_export;

    % Remove any seizures with a duration below 9 seconds
    data_tbl = data_tbl(data_tbl.duration >= 9,: );
    
    imprints = cell(size(data_tbl,1),1);
    onset_tab = table(repmat(string(patient),size(data_tbl,1),1), repelem("id", size(data_tbl,1),1), nan(size(data_tbl,1),1),...
        cell(size(data_tbl,1),1), 'VariableNames', ["patient_id", "segment_id", "Onset time", "Onset"]);
    
    % path_pipeline_exports = [data_location, 'export_ictal'];
    
    % Set basefolder to store markers
    
    sampling_rate=data_tbl.segment_fs(1); %just using the first, as all sampling was the same in our data after preproc.
    
    linelength_db = LL_db([basefolder '/LL_db/']); %setup folder for all Line Length measures
    linelength_db.add_paramset('wndw_len',sampling_rate*window_size,'wndw_overlap',sampling_rate*window_overlap);
    linelength_db.paramset_tbl                       % display all currently tracked paramsets
    linelength_db.calc(data_tbl,[]);                       % calculate all parametersets for all segments in data
    
    energy_db = Energy_db([basefolder '/Energy_db']);
    energy_db.add_paramset('wndw_len',sampling_rate*window_size,'wndw_overlap',sampling_rate*window_overlap);
    energy_db.paramset_tbl                       % display all currently tracked paramsets
    energy_db.calc(data_tbl,[]);                       % calculate all parametersets for all segments in data
    
    %bands: [1 4; 4 8; 8 13; 13 30; 30 60, 60 100]; 
    bandpower_db = BP_db([basefolder '/BP_db']);
    bandpower_db.add_paramset('wndw_len',sampling_rate*window_size,'wndw_overlap',sampling_rate*window_overlap,'bandbounds',[1 4]);
    bandpower_db.add_paramset('wndw_len',sampling_rate*window_size,'wndw_overlap',sampling_rate*window_overlap,'bandbounds',[4 8]);
    bandpower_db.add_paramset('wndw_len',sampling_rate*window_size,'wndw_overlap',sampling_rate*window_overlap,'bandbounds',[8 13]);
    bandpower_db.add_paramset('wndw_len',sampling_rate*window_size,'wndw_overlap',sampling_rate*window_overlap,'bandbounds',[13 30]);
    bandpower_db.add_paramset('wndw_len',sampling_rate*window_size,'wndw_overlap',sampling_rate*window_overlap,'bandbounds',[30 60]);
    bandpower_db.add_paramset('wndw_len',sampling_rate*window_size,'wndw_overlap',sampling_rate*window_overlap,'bandbounds',[60 100]);
    bandpower_db.paramset_tbl                       % display all currently tracked paramsets
    bandpower_db.calc(data_tbl,[]);   
    
    
    %% Import markers (line length, energy, bandpowers)
    
    calcs_ll = linelength_db.get(1);    % get all calculation outputs in variable. 
    calcs_energy = energy_db.get(1);
    calcs_bp_delta = bandpower_db.get(1);
    calcs_bp_theta = bandpower_db.get(2);
    calcs_bp_alpha = bandpower_db.get(3);
    calcs_bp_beta = bandpower_db.get(4);
    calcs_bp_gamma = bandpower_db.get(5);
    calcs_bp_hgamma = bandpower_db.get(6);
    
    %% create table of markers from which we will compute baseline and detect seizure activity
    val_tbl=[calcs_ll.LL_ms calcs_energy.energy_ms calcs_bp_delta.bp ...
        calcs_bp_theta.bp calcs_bp_alpha.bp calcs_bp_beta.bp calcs_bp_gamma.bp calcs_bp_hgamma.bp];
    
    %% Create a table of seizure data
    sz_mat_tab = table(data_tbl.segment_id, cell(size(val_tbl,1),1), calcs_ll.t_wndw,...
        'VariableNames',["segment_id", "feat_mat", "tw"]);
    for sz = 1:size(val_tbl,1)
        sz_mat=log(cat(3,val_tbl{sz,:})); % Here we log-transform markers 
        sz_mat_tab.feat_mat{sz} = sz_mat;
    end
    
    %% Testing use of mahal distance
    ict_buffer = 10;
    rec_thresh = 9;
    wl = 1/8;
    mov_med_val = (rec_thresh/2)*8;
    
 for sz = 6 %1:size(data_tbl,1)

        fprintf("\n Seizure %d", sz)
        tw = sz_mat_tab.tw{sz,1};
        pre_ids = find(tw<=-1*ict_buffer);
        pre_vals = sz_mat_tab.feat_mat{sz}(:,pre_ids,:);
        ictal_ids = find(tw>-1*ict_buffer & tw<=data_tbl.duration(sz));
        ict_vals = sz_mat_tab.feat_mat{sz}(:,ictal_ids,:);
        
%         %% Remove preictal outliers
%         for chan = 1:size(pre_vals,1)
%             for mark = 1:8
%                 % Remove preictal outliers (e.g. spikes)
%                 pre_vals(chan,:,mark) = filloutliers(pre_vals(chan,:,mark), NaN);
%             end
%         end
        %% For each channel create a distribution of Mahalanobis distances (remove
        % time point and compare distance against all other time points)
        % Create a distribution of Mahalanobis distances (one value per
        % channel, marker, and time point)
        pre_mahal_mat = nan(size(pre_vals,1), size(pre_vals(1,:,1),2));
        for chan = 1:size(pre_vals,1)
            pre_time_points = 1:size(pre_vals(chan,:,1),2);
            for pre_time_point = pre_time_points
                if all(~isnan(squeeze(pre_vals(chan,pre_time_point,:))))
                    ref_dist = squeeze(pre_vals(chan,pre_time_points(pre_time_points~=pre_time_point),:));
                    ref_dist = ref_dist(~isnan(sum(ref_dist,2)),:);
                    pre_mahal_mat(chan, pre_time_point) = mahal(squeeze(pre_vals(chan,pre_time_point,:))', ref_dist);
                end
            end
        end
%% Compute MAD for preictal Mahalanobis distances ( we will remove any time windows where MAD >= mad_thresh)
        mc = -1/(sqrt(2)*erfcinv(3/2)); % fixed factor for MAD score calculation
        pre_mahal_mad_mat = nan(size(pre_vals,1), size(pre_vals(1,:,1),2));
        
        for chan = 1:size(pre_vals,1)
            pre_dist = pre_mahal_mat(chan,:);
            for ict_time_point = ict_time_points
                pre_med = median(pre_dist, 2, 'omitnan');
                pre_smad=mc*mad_rewrite(pre_dist,1,2);
                pre_mahal_mad_mat(chan,:) = (pre_mahal_mat(chan,:)-pre_med)./pre_smad; %score ictal to median & scaled mad
                
            end
        end
        
        %% Compute Mahalanobis distance between ictal observations and preictal distribution
        ict_mahal_mat = nan(size(ict_vals,1), size(ict_vals(1,:,1),2));
        for chan = 1:size(pre_vals,1)
            ict_time_points = 1:size(ict_vals(chan,:,1),2);
            ref_dist = squeeze(pre_vals(chan,pre_mahal_mad_mat(chan,:)< mad_thresh,:));
            ref_dist = ref_dist(~isnan(sum(ref_dist,2)),:);
            for ict_time_point = ict_time_points
                ref_dist = squeeze(pre_vals(chan,:,:));
                ref_dist = ref_dist(~isnan(sum(ref_dist,2)),:);
                ict_mahal_mat(chan, ict_time_point) = mahal(squeeze(ict_vals(chan,ict_time_point,:))', ref_dist);
            end
        end
        
        %% MAD score ictal Mahalanobis distances against preictal Mahalanobis distances
        mc = -1/(sqrt(2)*erfcinv(3/2)); % fixed factor for MAD score calculation
        mahal_mad_mat = nan(size(ict_vals,1), size(ict_vals(1,:,1),2));
        pre_mahal_mad_mat = nan(size(pre_vals,1), size(pre_vals(1,:,1),2));
        
        for chan = 1:size(pre_vals,1)
            ict_time_points = 1:size(ict_vals(chan,:,1),2);
            pre_dist = pre_mahal_mat(chan,:);
            for ict_time_point = ict_time_points
                pre_med = median(pre_dist, 2, 'omitnan');
                pre_smad=mc*mad_rewrite(pre_dist,1,2);
                pre_mahal_mad_mat(chan,:) = (pre_mahal_mat(chan,:)-pre_med)./pre_smad; %score ictal to median & scaled mad
                mahal_mad_mat(chan,:) = (ict_mahal_mat(chan,:)-pre_med)./pre_smad; %score ictal to median & scaled mad
            end
        end
        
        %%
    %     figure()
    %     imagesc(mahal_mad_mat>=mad_thresh)
        %%
    %     markers = ["line length", "energy", "delta", "theta", "alpha", "beta", "low gamma", "high gamma"];
    %     chan = 63;
    %     figure()
    %     tiledlayout(8,8)
    %     for mark1 = 1:8
    %         for mark2 = 1:8
    %             if mark1 < mark2
    %                nexttile
    %                scatter(pre_vals(chan,:,mark1), pre_vals(chan,:,mark2),10,'.')
    %                hold on
    %                scatter(ict_vals(chan,1:10,mark1), ict_vals(chan,1:10,mark2),10,'.')
    %                xlabel(markers(mark1))
    %                ylabel(markers(mark2))
    %                hold off
    %             elseif mark1 == mark2
    %                 nexttile
    %                 axis off
    %             else
    %                 nexttile
    %                 axis off
    %             end
    %     
    %         end
    %     end
        
        %%
        time_points = (0:9)/8;
    %     figure("Position", [0,0,1050,1500])
    %     tiledlayout(5,1)
    %     for chan = 61:65
    %         nexttile
    %         histogram(pre_mahal_mat(chan,:))
    %         title(string(data_tbl(1,:).segment_channel_labels{:}(chan)))
    %         xlabel("Preictal Mahalanobis distance")
    %         ylabel("Frequency")
    %         xlim([0, 45])
    %         for ict_time_point = 1:10
    %             xline(ict_mahal_mat(chan,ict_time_point), 'r',...
    %                 sprintf("t(%d), MAD = %.3f",...
    %                 ict_time_point, mahal_mad_mat(chan, ict_time_point)))
    %         end
    %     end
    %   

        
        %% Add inclusion criteria (activity persisting beyond threshold)
        recruitment_threshold = rec_thresh/wl;
        movmed_mahal_mad_mat = movmedian(mahal_mad_mat, [mov_med_val,mov_med_val], 2);

        ms_a=movsum(movmed_mahal_mad_mat>=mad_thresh,[0 recruitment_threshold-1],2);%forward looking sum
        ms_b=movsum(movmed_mahal_mad_mat>=mad_thresh,[recruitment_threshold-1 0],2);%backward looking sum
        imprint = ms_a >= recruitment_threshold*0.95 | ms_b >= recruitment_threshold*0.95;%combining the backward and forward looking sum has the advantage of not cutting off early activity, or neglecting late activity
        imprints{sz} = imprint;
        
        onset_time = find(sum(imprint),1, 'first');
        if isempty(onset_time)
            onset = [];
            onset_tab(sz,:) = table(string(patient), data_tbl.segment_id(sz), NaN, {onset}, 'VariableNames', ["patient_id", "segment_id", "Onset time", "Onset"]);
        else
            onset = sum(imprint(:,onset_time+[0:7]),2)>=1;
            onset_tab(sz,:) = table(string(patient), data_tbl.segment_id(sz), (onset_time-1)/8, {onset}, 'VariableNames', ["patient_id", "segment_id", "Onset time", "Onset"]);
        end
%%

    sz_cell_imprint = table({imprint}, data_tbl.segment_id(sz), 'VariableNames', ["cell_imprint", "segment_id"]);
      %%  
    % Plot EEG with imprint overlaid

    if sum(sum(sz_cell_imprint.cell_imprint{:}))>0
        plot_eeg_imprints(data_tbl(sz,:), sz_cell_imprint, "save_fig",1, "save_fig_loc",...
                sprintf("figures/checking_imprint/%s", patient),...
                "imprint_changes_label", "mahal_distance", "sz", sz)
    end

    onset_channels = find(onset);
    for ch_ind = 1:length(onset_channels)
        x_max = size(pre_mahal_mad_mat,2) + size(mahal_mad_mat,2); 
        chan = onset_channels(ch_ind);
        eeg_dat = data_tbl(sz,:).segment_data{:};
        fig = figure();
        sgtitle(string(data_tbl.segment_channel_labels{1,1}(chan)))
        subplot(511)
        eeg_dat_ict = eeg_dat(:,120*512:(120+data_tbl(sz,:).duration)*512);
        plot((1:size(eeg_dat,2))/(512/8), eeg_dat(chan,:))
        xline(960)
        xlim([480, x_max])
        set(gca, "XTick", 480:80:x_max, "XTickLabel", -60+(0:10:length(480:80:x_max)*10))
        title("EEG")
        subplot(512)
        plot(1:(size(pre_mahal_mad_mat,2) + size(mahal_mad_mat,2)), [pre_mahal_mad_mat(chan,:), mahal_mad_mat(chan,:)])
        xlim([480, x_max])
        set(gca, "XTick", 480:80:x_max, "XTickLabel", -60+(0:10:length(480:80:x_max)*10))
        yline(mad_thresh)
        xline(960)
        title("MAD across time")
        subplot(513)
        plot((120*8)+(1:size(imprint,2)), imprint(chan,:))
        ylim([-0.5,1.5])
        xlim([480, x_max])
        set(gca, "XTick", 480:80:x_max, "XTickLabel", -60+(0:10:length(480:80:x_max)*10))
        xline(960)
        title(sprintf("Imprint (MAD thresh = %.1f)", mad_thresh))
        subplot(514)
        plot(1:(size(pre_mahal_mad_mat,2) + size(mahal_mad_mat,2)),...
            movmedian([pre_mahal_mad_mat(chan,:), mahal_mad_mat(chan,:)], [mov_med_val,mov_med_val]))
        xlim([480, x_max])
        set(gca, "XTick", 480:80:x_max, "XTickLabel", -60+(0:10:length(480:80:x_max)*10))
        xline(960)
        yline(mad_thresh)
        title(sprintf("MAD moving median (%.2f seconds forwards and back)", mov_med_val/8))
        subplot(515)
        plot(1:(size(pre_mahal_mad_mat,2) + size(mahal_mad_mat,2)),...
            movmedian([pre_mahal_mad_mat(chan,:), mahal_mad_mat(chan,:)], [mov_med_val,mov_med_val])>=mad_thresh)
        ylim([-0.5,1.5])
        xlim([480, x_max])
        set(gca, "XTick", 480:80:x_max, "XTickLabel", -60+(0:10:length(480:80:x_max)*10))
        xline(960)
        title(sprintf("Imprint using moving median (MAD thresh = %.1f)", mad_thresh))
        saveas(fig, sprintf("figures/checking_imprint/%s/mahal_distance/%d_%s.png", ...
                patient, sz, string(data_tbl.segment_channel_labels{1,1}(chan))))
    end
       
%%
     

 end

 %%
channel = "EEG GB_02-REF";
chan = strcmp(data_tbl.segment_channel_labels{sz}, channel);
  x_max = size(pre_mahal_mad_mat,2) + size(mahal_mad_mat,2); 
        eeg_dat = data_tbl(sz,:).segment_data{:};
        fig = figure();
        tiledlayout(5,1)
        sgtitle(string(data_tbl.segment_channel_labels{1,1}(chan)))
        nexttile
        eeg_dat_ict = eeg_dat(:,120*512:(120+data_tbl(sz,:).duration)*512);
        plot((1:size(eeg_dat,2))/(512/8), eeg_dat(chan,:))
        xline(960)
        xlim([480, x_max])
        set(gca, "XTick", 480:80:x_max, "XTickLabel", -60+(0:10:length(480:80:x_max)*10))
        title("EEG")
        nexttile
        plot(1:(size(pre_mahal_mad_mat,2) + size(mahal_mad_mat,2)), [pre_mahal_mad_mat(chan,:), mahal_mad_mat(chan,:)])
        xlim([480, x_max])
        set(gca, "XTick", 480:80:x_max, "XTickLabel", -60+(0:10:length(480:80:x_max)*10))
        yline(mad_thresh)
        xline(960)
        title("MAD across time")
        nexttile
        plot((120*8)+(1:size(imprint,2)), imprint(chan,:))
        ylim([-0.5,1.5])
        xlim([480, x_max])
        set(gca, "XTick", 480:80:x_max, "XTickLabel", -60+(0:10:length(480:80:x_max)*10))
        xline(960)
        title(sprintf("Imprint (MAD thresh = %.1f)", mad_thresh))
        nexttile
        plot(1:(size(pre_mahal_mad_mat,2) + size(mahal_mad_mat,2)),...
            movmedian([pre_mahal_mad_mat(chan,:), mahal_mad_mat(chan,:)], [mov_med_val,mov_med_val]))
        xlim([480, x_max])
        set(gca, "XTick", 480:80:x_max, "XTickLabel", -60+(0:10:length(480:80:x_max)*10))
        xline(960)
        yline(mad_thresh)
        title(sprintf("MAD moving median (%.2f seconds forwards and back)", mov_med_val/8))
        nexttile
        plot(1:(size(pre_mahal_mad_mat,2) + size(mahal_mad_mat,2)),...
            movmedian([pre_mahal_mad_mat(chan,:), mahal_mad_mat(chan,:)], [mov_med_val,mov_med_val])>=mad_thresh)
        ylim([-0.5,1.5])
        xlim([480, x_max])
        set(gca, "XTick", 480:80:x_max, "XTickLabel", -60+(0:10:length(480:80:x_max)*10))
        xline(960)
        title(sprintf("Imprint using moving median (MAD thresh = %.1f)", mad_thresh))
       
        %%


%      %% If you'd like to look at a channel in more detail
%    channel = "EEGSFG4";
%         chan = strcmp(data_tbl.segment_channel_labels{sz}, channel);
%         eeg_dat = data_tbl(sz,:).segment_data{:};
%         figure()
%         sgtitle(string(data_tbl.segment_channel_labels{1,1}(chan)))
%         subplot(511)
%         eeg_dat_ict = eeg_dat(:,120*512:(120+data_tbl(sz,:).duration)*512);
%         plot((1:size(eeg_dat,2))/(512/8), eeg_dat(chan,:))
%         xline(960)
%         xlim([0, size(pre_mahal_mad_mat,2) + size(mahal_mad_mat,2)])
%         title("EEG")
%         subplot(512)
%         plot(1:(size(pre_mahal_mad_mat,2) + size(mahal_mad_mat,2)), [pre_mahal_mad_mat(chan,:), mahal_mad_mat(chan,:)])
%         xlim([0, size(pre_mahal_mad_mat,2) + size(mahal_mad_mat,2)])
%         yline(mad_thresh)
%         xline(960)
%         title("MAD across time")
%         subplot(513)
%         plot((120*8)+(1:size(imprint,2)), imprint(chan,:))
%         ylim([-0.5,1.5])
%         xlim([0, size(pre_mahal_mad_mat,2) + size(mahal_mad_mat,2)])
%         xline(960)
%         title(sprintf("Imprint (MAD thresh = %.1f)", mad_thresh))
%         subplot(514)
%         plot(1:(size(pre_mahal_mad_mat,2) + size(mahal_mad_mat,2)),movmedian([pre_mahal_mad_mat(chan,:), mahal_mad_mat(chan,:)], [mov_med_val,mov_med_val]))
%         xlim([0, size(pre_mahal_mad_mat,2) + size(mahal_mad_mat,2)])
%         xline(960)
%         yline(mad_thresh)
%         title(sprintf("MAD moving median (%.2f seconds forwards and back)", mov_med_val/8))
%         subplot(515)
%         plot(1:(size(pre_mahal_mad_mat,2) + size(mahal_mad_mat,2)), movmedian([pre_mahal_mad_mat(chan,:), mahal_mad_mat(chan,:)], [mov_med_val,mov_med_val])>=4)
%         ylim([-0.5,1.5])
%         xlim([0, size(pre_mahal_mad_mat,2) + size(mahal_mad_mat,2)])
%         xline(960)
%         title(sprintf("Imprint using moving median (MAD thresh = %.1f)", mad_thresh))


   % pause
   %close all
end

    %%
    
%     cell_imprint = table(imprints, data_tbl.segment_id, 'VariableNames', ["cell_imprint", "segment_id"]);
%     % Remove seizures with no imprint
%     missing_onset = cellfun(@isempty, onset_tab.Onset);
%     onset_tab = onset_tab(~missing_onset,:);
%     cell_imprint = cell_imprint(~missing_onset,:);
%     data_tbl = data_tbl(~missing_onset,:);
%     
%     % Plot EEG with imprint overlaid
%     plot_eeg_imprints(data_tbl, cell_imprint, "save_fig", 0, "save_fig_loc",...
%             sprintf("figures/checking_imprint/%s", patient),...
%             "imprint_changes_label", "mahal_distance")
%     