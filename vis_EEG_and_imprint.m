%%
det = 8;
min_sz = 1;
min_sz_dur = 9; % minimum seizure duration for inclusion
chan_to_roi_thresh_type = "count";
chan_to_roi_thresh = 1; % One channel in region is sufficient to include region
wind_overlap = 7/8; % Overlap between imprint windows
rec_type = "sec";
rec_thresh = 3;
mad_thresh = 5;

onset_calc_loc = "imprint_ons"; 

patient = "UCLH1005";
onset_output = final_output(string(final_output.Patient_id) == patient,:);

load(sprintf('%s/%s.mat', data_location, patient));
pat_data = data_export;
[pat_data, ~, ~] = incl_crit(pat_data, 'min_sz_count', 1, sz_type = ["focal", "sg", "subclin", "N"]);

pat_meta = pat_data(:,1:(end-1));

[pat_data, pat_meta, cell_imprint,  sz_count_pat] = ...
    calc_imprint(pat_data, pat_meta, "window_overlap", wind_overlap,...
    "folder",onset_calc_loc, "min_sz_count", min_sz,...
    "rec_type", rec_type, "rec_thresh",rec_thresh);
%%
% 
% ons_ind = find(onset_output.imprint_chan{:}(:,3)+onset_output.imprint_chan{:}(:,6));
% example_chan_count = 20; %sum(onset_output.imprint_chan{:}(:,sz)); %length(onset_output.channel_names{:}); 
%     % We will pull out 20 example channels (including onset channels)
% chans = 1:size(sz_data,1);
% chans(ons_ind) = [];
% rng(1)
% vis_chans = [randsample(chans,example_chan_count-length(ons_ind)), ons_ind'];
% vis_binary = zeros(size(sz_data,1),1);
% vis_binary(vis_chans) = 1;
% vis_binary = logical(vis_binary);

% We will show a subsection of EEG channels (25) for clearer visualisations
% We must include any onset channels so imprint onset across seizures has
% no empty columns
incl_channels = sum(final_output(pat,:).imprint_chan{:}, 2) > 0;
remaining_channels = find(~incl_channels);

rng(1)
random_channels = remaining_channels(randsample(length(remaining_channels), 25 - sum(incl_channels)));
incl_channels(random_channels) = 1;

for sz = [3,6]%:size(pat_data,1)
    sz_row = pat_data(sz,:);
    sz_data = sz_row.segment_data{:};
    fs = sz_row.segment_fs;
    t = 1:size(sz_data,2);
    
    
    
    clrs = repelem(0, size(sz_data,1),3);
    %clrs(ons_ind,:) = [1,0,0];
    clrs = clrs(incl_channels,:);
    
    % Plot example channels
    opts.plot_labels = true;
    opts.clrs = clrs;
    chan_names = onset_output.channel_names{:}(incl_channels);
    opts.labels = chan_names;
    opts.offset = 1000;
    
%     figure(1)
% %     subplot(1,2,1)
%     onset_time = pat_onset.when_onset{:}(sz,1)/8;
%     vis_seg_t = (512*(100)):(512*(130+onset_time));
%     vis_plot_eeg((vis_seg_t/512)-120,sz_data(vis_binary,vis_seg_t),opts);
%     xline(onset_time, 'b', LineWidth=2)
%     xline(0,'g',  LineWidth=2)
%     
%     %Plot imprint
%     figure(2)
% %     subplot(1,2,2)
%     imagesc(cell_imprint{sz,1}{:}(vis_binary, 1:((onset_time+10)*8)))%1:onset_time*8+20))
%     cm = [1,1,1; 1,0,0];
%     colormap(cm)
%     set(gca, "YTick", [], "XTick",[])
% %     set(gca, "YTick", [], "XTick", 1:((onset_time+20)*8), "XTickLabel", [])%(1:((onset_time+20)*8))/8)
% %     grid on


    figure()
    subplot(1,2,1)
    onset_time = onset_output.when_onset{:}(sz,1)/8;
    vis_seg_t = (512*(1)):(512*(120+pat_data(sz,:).duration));
    vis_plot_eeg((vis_seg_t/512)-120,sz_data(incl_channels,vis_seg_t),opts);
    xline(0,'g',  LineWidth=2)
    xline(onset_time, 'b', LineWidth=2)
    xline(pat_data(sz,:).duration, 'g',  LineWidth=2)
    
    %Plot imprint
    %figure(sz+size(pat_data,1))
    subplot(1,2,2)
    imagesc(cell_imprint{sz,1}{:}(incl_channels,(onset_time*8):end))
%     set(gca, "XTick", (1:size(cell_imprint{sz,1}{:},2)/8)*8, "XTickLabel", ...
%         (1:size(cell_imprint{sz,1}{:},2)/8),"YTick", 1:example_chan_count,...
%         "YTickLabel",onset_output.channel_names{:}(vis_binary), ...
%         'TickLength',[0, 0])
 set(gca, "YTick", [],'TickLength',[0, 0])
    cm = [1,1,1; 0,0,1];
    colormap(cm)
end

%%
figure()
subplot(2,3,1)
imagesc(test_output(string(test_output.Patient_id) == patient,:).clo_roi_120{:})
set(gca, "YTick", 1:length(test_output(string(test_output.Patient_id) == patient,:).roi_names_120{:}), "YTickLabel", test_output(string(test_output.Patient_id) == patient,:).roi_names_120{:})
subplot(2,3,2:3)
imagesc(test_output(string(test_output.Patient_id) == patient,:).imprint_roi_120{:})
set(gca, "YTick", [],  "XTick", [])

subplot(2,3,4)
imagesc(test_output(string(test_output.Patient_id) == patient,:).clo_roi_250{:})
set(gca, "YTick", 1:length(test_output(string(test_output.Patient_id) == patient,:).roi_names_250{:}), "YTickLabel", test_output(string(test_output.Patient_id) == patient,:).roi_names_250{:})
subplot(2,3,5:6)
imagesc(test_output(string(test_output.Patient_id) == patient,:).imprint_roi_250{:})
set(gca, "YTick", [],  "XTick", [])


%%
roi_names = pat_onset.roi_names_120{:};
roi_names = strrep(roi_names, 'r.', 'ctx-rh-');
roi_names = strrep(roi_names, 'l.', 'ctx-lh-');
cm = [0.9,0.9,0.9; %  recorded but no channels remaining (-1: pale grey)
        0,0,1];% onset (blue)
cm = interp1(cm, 1:0.01:size(cm,1));
    
plotBrain(roi_names, pat_onset.imprint_roi_120{:}(:,3),cm, 'atlas', 'lausanne120_aseg', 'limits', [0,1])
plotBrain(roi_names, pat_onset.imprint_roi_120{:}(:,11),cm, 'atlas', 'lausanne120_aseg', 'limits', [0,1])
plotBrain(roi_names, pat_onset.clo_roi_120{:},cm, 'atlas', 'lausanne120_aseg', 'limits', [0,1])
