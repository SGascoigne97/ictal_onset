
atl = 60;
plot_atl = 'lausanne250'; %'lausanne120_aseg';
det = 0;
param_set = final_output(final_output.atl == 60 & ...
    final_output.det == sprintf("+%d seconds", det),:);


%%
% Create colour map for plotting seizure onset (as proportion of seizures)
    cm = [0.9,0.9,0.9; %  Not included in onsets
        1,1,0; % Onset in 1/3 of seizures (yellow)
        0.85,0.325,0.098; % Onset in 2/3 of seizures (orange)
        1,0,0]; % Onset in all seizures (red)
    cm = interp1(cm, 1:0.01:size(cm,1));
%% 
% 
loc_ons = cell(length(final_output.Patient_id), 4);
loc_ons = cell2table(loc_ons, 'VariableNames', {'patient_id', 'onset_prop', 'resec', 'outcome'});
loc_ons.patient_id = final_output.Patient_id;


for pat = 1:length(final_output.Patient_id)
    patient = string(final_output.Patient_id(pat));
    pat_onset = final_output(final_output.Patient_id == patient,:);
%     figure(pat)
%     subplot(1,3,1)
%     imagesc(pat_onset.imprint_roi{:})
%     set(gca, 'YTick', 1:length(pat_onset.ROI_ids{1,1}) , ...
%         'YTickLabels', pat_onset.ROI_ids{1,1})
%     title("Imprint")
%     subplot(1,3,2)
%     imagesc(pat_onset.EI_roi{:})
%     set(gca, 'YTick', [])
%     title("Imprint")
%     title("EI")
%     subplot(1,3,3)
%     imagesc(pat_onset.PLHG_roi{:})
%     set(gca, 'YTick', [])
%     title("PLHG")
%     sgtitle(sprintf("Patient %s", patient))

    % Plot variability in seizure onsets on brain
    unq_roi = pat_onset.roi_names{:};
    unq_roi = strrep(unq_roi, 'r.', 'ctx-rh-');
    unq_roi = strrep(unq_roi, 'l.', 'ctx-lh-');

    col_sch = sum(pat_onset.imprint_roi{:},2)/ size(pat_onset.imprint_roi{:},2);

%     plotBrain(unq_roi, col_sch, cm, 'atlas', plot_atl, 'limits', [0,1], ...
%         'savePath', char(sprintf('./figures/onset_var/%s', patient)))
%     plotBrain(unq_roi, pat_onset.resected_roi{:}, cm, 'atlas', plot_atl, 'limits', [0,1], ...
%         'savePath', char(sprintf('./figures/onset_var/%s_resec', patient)))

    loc_ons(pat,:).onset_prop = {col_sch};
    loc_ons(pat,:).resec = pat_onset.resected_roi;

    pat_surg_out = pat_onset.Surgery_outcome{:};
    year_1 = pat_onset.("Outcome year"){:} - pat_onset.("Surgery year") == 1;
    loc_ons(pat,:).outcome = num2cell(pat_surg_out(year_1));
    loc_ons.sz_count(pat) = num2cell(length(pat_onset.Segment_ids{:}));
end


%%
percentle = 95; 

for pat = 1:size(loc_ons,1)
    ons = loc_ons(pat,:).onset_prop{:};
    freq_ons = zeros(length(ons),1);
    freq_ons(find(ons>= prctile(ons, percentle))) =...
        1;
    resec_freq_ons = freq_ons + loc_ons(pat,:).resec{:} == 2;
    loc_ons.resec_freq_ons(pat) = sum(resec_freq_ons)/sum(freq_ons);
end

%%
loc_ons = loc_ons( cat(1,loc_ons.outcome{:}) ~= 8,:);
figure()
boxplot(loc_ons.resec_freq_ons, loc_ons.outcome)


%% Patient 950 resection before preprocessing
patient = "950";
    % load json_data
    filelist = dir(fullfile(strcat(path_pipeline_exports, "/UCLH", patient, "/"), '**/*.*'));  % get list of files and folders in any subfolder
    filelist = filelist(~[filelist.isdir]);  % remove folders from list
    folder = filelist(strcmp({filelist(:).name},'json_data.mat')).folder;
    load(strcat(folder, '/json_data.mat')); % load the json_data file in the folder

    chan_det = json_data.channel_details;
    chan_missing_roi =  cellfun(@isempty,chan_det.ROIname);
    chan_det = chan_det(~chan_missing_roi,:);
    resec_chan = chan_det.is_resected5;
    roi_names =  cat(2, chan_det.ROIname{:});
    incl_roi = roi_names(2,:);
    unq_roi = unique(incl_roi, 'stable');

    [resec_bin, resec_names] = chan_to_roi_crit(resec_chan, incl_roi, unq_roi, "threshold", NaN);

    unq_roi = strrep(unq_roi, 'r.', 'ctx-rh-');
    unq_roi = strrep(unq_roi, 'l.', 'ctx-lh-');

    plotBrain(unq_roi, resec_bin, cm, 'atlas', plot_atl, 'limits', [0,1])