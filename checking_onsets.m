%% Checking imprint onsets for example patients
% Direct to location of saved subject data
data_location = '/media/b5007876/DATA/Data/UCLH_export/imprintDataExportFinal/'; % replace
path_pipeline_exports = [data_location, 'export_ictal'];
basefolder = 'example_onset_detection';

addpath('lib_biomarkers/')
addpath('lib_dataflow/')
addpath(genpath('lib/'))

final_output_struct = load('tables/final_output.mat'); % Here any subjects without follow-up have been removed
final_output = final_output_struct.final_output;
clear final_output_struct
%%
% Select a subject to create figures for
% Choose atlas
% Choose threshold for allowing propagated activity (seconds)
det = 8;
min_sz = 1;
min_sz_dur = 9; % minimum seizure duration for inclusion
chan_to_roi_thresh_type = "count";
chan_to_roi_thresh = 1; % One channel in region is sufficient to include region
wind_overlap = 7/8; % Overlap between imprint windows
rec_type = "sec";
rec_thresh = 3;
mad_thresh = 5;
for pat = 1%:size(final_output,1) %[1, 7, 12, 39, 60, 45, 22, 59]
    pat_id = string(final_output(pat,:).Patient_id);
    atlases = [2,3]; %1=scale36,2=60,3=125,4=250

    % We will show a subsection of EEG channels (25) for clearer visualisations
    % We must include any onset channels so imprint onset across seizures has
    % no empty columns
%     incl_channels = sum(final_output(pat,:).imprint_chan{:}, 2) > 0;
%     remaining_channels = find(~incl_channels);
%     
%     rng(1)
%     random_channels = remaining_channels(randsample(length(remaining_channels), 25 - sum(incl_channels)));
%     incl_channels(random_channels) = 1;
    incl_channels = 1:size(final_output(pat,:).imprint_chan{:},1);
    
    %% Compute imprint
    pat_data = load(sprintf('%s/%s.mat', data_location, pat_id));
    pat_data = pat_data.data_export;
    
    pat_onset = final_output(final_output.Patient_id == pat_id,:);
    pat_data = pat_data(ismember(pat_data.segment_id, pat_onset.Segment_ids{1,1}),:);
    pat_meta = pat_data(:, 1:end);
    
    [~, ~, cell_imprint,  ~] = ...
                           calc_imprint(pat_data, pat_meta, 'window_overlap',...
                           wind_overlap,'folder',basefolder, ...
                           'min_sz_count', min_sz, 'mad_thresh',mad_thresh,...
                           'rec_type', "sec", 'rec_thresh', 3);
    
    cell_imprint = cell_imprint(find(ismember(cell_imprint.segment_id, pat_onset.Segment_ids{1,1})),:);
  
    plot_eeg_imprints(pat_data, cell_imprint, "save_fig",1, "save_fig_loc",...
        sprintf("figures/checking_imprint/%s",string(pat_onset.Patient_id)),...
        "imprint_changes_label", "remove_preict_outliers")
end