
atl = 60;
plot_atl = 'lausanne120_aseg';
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
for pat = 3%1:length(param_set.Patient_id)
    patient = param_set.Patient_id(pat);
    pat_onset = param_set(param_set.Patient_id == patient,:);
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
    unq_roi = pat_onset.ROI_ids{:};
    unq_roi = strrep(unq_roi, 'r.', 'ctx-rh-');
    unq_roi = strrep(unq_roi, 'l.', 'ctx-lh-');

    col_sch = sum(pat_onset.imprint_roi{:},2)/ length(pat_onset.Segment_ids{1,1});

    plotBrain(unq_roi, col_sch, cm, 'atlas', plot_atl, 'limits', [0,1], ...
        'savePath', char(sprintf('./figures/onset_var/%s', patient)))
     plotBrain(unq_roi, pat_onset.resected_roi{:}, cm, 'atlas', plot_atl, 'limits', [0,1], ...
        'savePath', char(sprintf('./figures/onset_var/%s_resec', patient)))
end

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