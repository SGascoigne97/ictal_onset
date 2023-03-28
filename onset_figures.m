% For each patient, create and save plots displaying onset as binary map
% and onset times labelled on EEGs
for pat = 3%:length(final_output.Patient_id)
    close all
    % Pull out patient-specific onset from table
    onset_output = final_output(pat,:);
    patient = sprintf('UCLH%s', onset_output.Patient_id);
    % Load patient data (this is where we find their EEGs)
    load(sprintf('%s/%s.mat', data_location, patient));
    data = data_export;
    % Some seizures were removed in the onset detction script, remove these
    % seizures here
    data = data(contains(string(data.segment_id), string(onset_output.Segment_ids{:})),:);

    % Create a folder to store figures in
    mkdir(sprintf('figures/%s', patient))
    
    % First figure shows channel-wise seizure onsets
    figure(1)
    subplot(1,7,2:3)
    imagesc(onset_output.Imprint_onset{:})
    set(gca,'ytick',1:1:(size(data.segment_channel_labels{1,1},1)),...
        'yticklabel',data.segment_channel_labels{1,1})
    title('Imprint')
    subplot(1,7,4:5)
    imagesc(onset_output.EI_onset{:})
    set(gca,'ytick',[],'yticklabel',[])
    title('EI')
    subplot(1,7,6:7)
    imagesc(onset_output.PLHG_onset{:})
    set(gca,'ytick',[],'yticklabel',[])
    title('PLHG')
    sgtitle(sprintf('%s (channel-wise)', onset_output.Patient_id))
    saveas(gcf, sprintf('figures/%s/channel_wise_onset', patient), 'png' )
    
    % Second figure shows ROI-wise seizure onsets
    figure(2)
    subplot(1,7,1)
    imagesc(onset_output.Labelled_onset{:})
    set(gca,'ytick',1:1:size(onset_output.ROI_ids{1,1},1),...
        'yticklabel',onset_output.ROI_ids{1,1}, 'xtick', [])
    title('CLO')
    subplot(1,7,2:3)
    imagesc(onset_output.imprint_roi{:})
    set(gca,'ytick',[],'yticklabel',[])
    title('Imprint')
    subplot(1,7,4:5)
    imagesc(onset_output.EI_roi{:})
    set(gca,'ytick',[],'yticklabel',[])
    title('EI')
    subplot(1,7,6:7)
    imagesc(onset_output.PLHG_roi{:})
    set(gca,'ytick',[],'yticklabel',[])
    title('PLHG')
    sgtitle(sprintf('%s (ROI-wise)', onset_output.Patient_id))
    saveas(gcf, sprintf('figures/%s/ROI_wise_onset', patient), 'png' )

    % Finally we plot each seizure's EEG with CLO and automatically
    % detected onset times plotted on top
    opts.plot_labels = true;
    opts.labels= data.segment_channel_labels{1,1};
    opts.offset = 600;
    opts.clrs = repmat([0.5 0.5 0.5], length(opts.labels),1);
    for i = 1:size(data,1)
        figure(i+2)
        vis_plot_eeg((1:size(data.segment_data{i,1},2))/512,data.segment_data{i,1}, opts);
        xline(120, '-k', 'Labelled Onset', 'LabelHorizontalAlignment','left')
        xline(onset_output.when_onset{1,1}(i,:)+120,'--r',{'Imprint','EI','PLHG'})
        title(sprintf('Patient %s seizure %s', patient, string(onset_output.Segment_ids{1}(i))))
        saveas(gcf, sprintf('figures/%s/%s', patient, string(onset_output.Segment_ids{1}(i))), 'png' )
    end
end