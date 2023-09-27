%% Code for creating figure one
% Direct to location of saved subject data
data_location = '/media/b5007876/DATA/Data/UCLH_export/imprintDataExportFinal/'; % replace
path_pipeline_exports = [data_location, 'export_ictal'];
basefolder = 'example_onset_detection';

addpath('lib_biomarkers/')
addpath('lib_dataflow/')
addpath(genpath('lib/'))

%% Load data and Lausanne atlases
final_output_struct = load('tables/final_output.mat'); % Here any subjects without follow-up have been removed
final_output = final_output_struct.final_output;
clear final_output_struct
addpath(genpath('sarah_functions'))

load('roi_info/ATLAS.mat')
%%
% Select a subject to create figure for
pat = 2;
pat_id = string(final_output(pat,:).Patient_id);
atlases = [2,3]; %1=scale36,2=60,3=125,4=250
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

% We will show a subsection of EEG channels (25) for clearer visualisations
% We must include any onset channels so imprint onset across seizures has
% no empty columns
incl_channels = sum(final_output(pat,:).imprint_chan{:}, 2) > 0;
remaining_channels = find(~incl_channels);

rng(1)
random_channels = remaining_channels(randsample(length(remaining_channels), 25 - sum(incl_channels)));
incl_channels(random_channels) = 1;
incl_chan_id = find(incl_channels);

%% Organise channels into order of regions
roi_names = pat_onset.roi_names_120{:}; 
% List channel ids in region order
chan_ind = [];
roi_ind = [];
for roi = 1:length(roi_names)
    chan_ind = [chan_ind, find(final_output.chan_2_roi_matrix_120{1,1}(roi,:))];
    roi_ind = [roi_ind, repmat(roi, 1, length(find(final_output.chan_2_roi_matrix_120{1,1}(roi,:))))];
end
roi_tbl = table(chan_ind', roi_ind', 'VariableNames', ["channel_id", "roi_id"]);

% roi_tbl is my channel order reference table
roi_tbl_included = roi_tbl(ismember(roi_tbl.channel_id, incl_chan_id),:);

incl_chan_ordered = flipud(roi_tbl_included.channel_id);

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
%%
xlims = nan(size(pat_data,1), 2);
xlims([3,7],:) = [0, 1200; 0, 960];
for sz =[3,7]
    f = figure();
    f.Position = [100,100, 800, 500];
    eeg_dat = pat_data.segment_data{sz,1};
    imprint = cell_imprint.cell_imprint{sz,1};
    imprint_bin = [zeros(size(imprint,1),60*8), imprint];
    imprint_bin = repelem(imprint_bin(incl_chan_ordered,:),1000,1);
    for chan = 1:25
        top = 900 + [(chan-1)*1000];
        bottom = chan*1000;
        imprint_bin(top:bottom,:) = 0;
    end
    imprint_bin = [imprint_bin; zeros(500,size(imprint_bin,2))];
    imprint_bin = imprint_bin(500:end,:);

    % Extract segments of EEG which are included in imprint
    % Add preictal and ictal segments to imprint
    imprint_full = [zeros(size(imprint,1),120*8), imprint, zeros(size(imprint,1),120*8)];
    imprint_eeg = repelem(imprint_full, 1, round(size(eeg_dat,2)/size(imprint_full,2)));
    % Add in zeros at each side so the sizes are equivalent 
    dif = size(eeg_dat,2) - size(imprint_eeg,2);
    if dif > 0
        imprint_eeg = [zeros(size(imprint_eeg,1), round(dif/2)), imprint_eeg, zeros(size(imprint_eeg,1), round(dif/2)-1)];
    elseif dif < 0
        dif = abs(dif);
        imprint_eeg = imprint_eeg(:, round(dif/2):end-round(dif/2));
    end

    imprint_eeg(imprint_eeg == 0) = NaN;

    imprint_eeg_dat = eeg_dat.*imprint_eeg;

    hold on
    imagesc(imprint_bin, 'AlphaData', 0.5)
    colormap([1,1,1;0,0.2,1])
    vis_eeg(eeg_dat(incl_chan_ordered,(512*60):(512*(180+pat_data.duration(sz)))),...
        512/8, "ChannelNames", pat_data.segment_channel_labels{1,1}(incl_chan_ordered),...
        "PlotNewFig", false, "Color", [0.2,0.2,0.2], "Offset", 1000);
    vis_eeg(imprint_eeg_dat(incl_chan_ordered,(512*60):(512*(180+pat_data.duration(sz)))),...
        512/8, "ChannelNames", pat_data.segment_channel_labels{1,1}(incl_chan_ordered),...
        "PlotNewFig", false, "Color", [0,0,0.7], "Offset", 1000);
    hold off

    xline(60*8, LineWidth=2, Color="red")
    xline(find(sum(imprint_bin),1, 'first'), LineWidth=2, Color="blue")
    xline((60+pat_data.duration(sz))*8, LineWidth=2, Color="red")
    set(gca, "XTick", 0:240:1680, "XTickLabel", -60:30:150)
    ylim([-1000 26000])
    xlim(xlims(sz,:))
    
    saveas(f,sprintf("figures/paper_figures/Figure 1/%s_sz%d.png", string(pat_onset.Patient_id), sz))
end

%% Panel C: Channel-wise imprint onset 
onsets = final_output(pat,:).imprint_chan{:}(incl_chan_ordered,:);
f = figure();
heatmap(flipud(onsets), 'CellLabelColor','none')
colormap([0.8,0.8,0.8;1,0.5,0.5])
colorbar off
saveas(f,sprintf("figures/paper_figures/Figure 1/%s_all_chan_onsets.svg", string(pat_onset.Patient_id)))

%% Panel D: ROI-120 imprint onset (regions included in channels from Panels A-C)
f = figure();
onsets = final_output(pat,:).imprint_roi_120{:}(unq_rois(1:10),:);
subplot(1,5,1:4)
heatmap(onsets, 'CellLabelColor','none')
colorbar off
subplot(1,5,5)
heatmap(double(mean(final_output(pat,:).imprint_roi_120{:}(unq_rois(1:10),:),2)>=0.5), 'CellLabelColor','none')
colormap([0.8,0.8,0.8;1,0.5,0.5])
colorbar off
saveas(f,sprintf("figures/paper_figures/Figure 1/%s_all_roi_onsets.svg", string(pat_onset.Patient_id)))


%% Panel E: Onset regions for two example seizures highlighted on brain
roi_names = pat_onset.roi_names_120{:};
roi_names = strrep(roi_names, 'r.', 'ctx-rh-');
roi_names = strrep(roi_names, 'l.', 'ctx-lh-');
cm = [0.8,0.8,0.8;1,0.5,0.5];

plotBrain(roi_names, double(mean(final_output(pat,:).imprint_roi_120{:},2)>=0.5),...
    cm, 'atlas', 'lausanne120_aseg', 'limits', [0,1], 'savePath', 'figures/paper_figures/Figure 1/consensus')

% plotBrain(roi_names, pat_onset.imprint_roi_120{:}(:,3),cm, 'atlas', 'lausanne120_aseg', 'limits', [0,1])%, 'savePath', 'figures/paper_figures/sz3')
% plotBrain(roi_names, pat_onset.imprint_roi_120{:}(:,7),cm, 'atlas', 'lausanne120_aseg', 'limits', [0,1], 'savePath', 'figures/paper_figures/sz7')



%% 
test = table(chan_ind', (1:length(chan_ind))')