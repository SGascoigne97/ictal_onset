patient = "UCLH964";
data_location = '/media/b5007876/DATA/Data/UCLH_export/imprintDataExportFinal/'; % replace
path_pipeline_exports = [data_location, 'export_ictal'];

pat_data = load(sprintf('%s/%s.mat', data_location, patient));
pat_data = pat_data.data_export;

pat_onset = final_output(final_output.Patient_id == patient,:);
pat_data = pat_data(ismember(pat_data.segment_id, pat_onset.Segment_ids{1,1}),:);



%%

[~, ~, cell_imprint, ~] = calc_imprint(pat_data, pat_data(:,1:(end-1)));

cell_imprint = cell_imprint(find(ismember(cell_imprint.segment_id, pat_onset.Segment_ids{1,1})),:);

%%
for sz = 1:12%length(final_output.Segment_ids{1,1})  
    figure(sz)
    eeg_dat = pat_data.segment_data{sz,1};
    subplot(2,4,1:3)
    vis_eeg(eeg_dat(:,(512*120):(512*(120+pat_data.duration(sz)))),...
        512, "ChannelNames",...
        pat_data.segment_channel_labels{1,1}, "PlotNewFig", false);
%     set(gca, 'xtick',[], 'xticklabel',[])
%     xlabel("")
    subplot(2,4,5:7)
    imagesc(cell_imprint.cell_imprint{sz,1})
    yticks = 1:length(pat_data.segment_channel_labels{1,1});
    xlabel("Time since marked onset (s)")
    set(gca, 'YTick', yticks, 'YTickLabel', ...
        pat_data.segment_channel_labels{1,1})
    title("Imprint")
    subplot(2,4,8)
    imagesc(pat_onset.imprint_chan{1,1}(:,sz))
    set(gca,'xtick',[], 'ytick', [])
    title("Imprint onset")
    
    sgtitle(sprintf("Patient %s seizure %s channel-level imprint", patient, final_output.Segment_ids{1,1}{sz}))
end



%%
figure()
subplot(2,1,1)
imagesc(final_output.imprint_chan{1,1})
subplot(2,1,2)
imagesc(final_output.imprint_chan{2,1})