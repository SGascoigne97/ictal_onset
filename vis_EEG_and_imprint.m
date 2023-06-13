%%
pat_onset = final_output(string(final_output.Patient_id) == patient,:);

for sz = 8%:size(pat_data,1)
    sz_row = pat_data(sz,:);
    sz_data = sz_row.segment_data{:};
    fs = sz_row.segment_fs;
    t = 1:size(sz_data,2);
    
    example_chan_count = 15;
    % We will pull out 15 example channels (including onset channels)
    ons_ind = find(pat_onset.imprint_chan{:}(:,sz));
    chans = 1:size(sz_data,1);
    chans(ons_ind) = [];
    rng(7)
    vis_chans = [randsample(chans,example_chan_count-length(ons_ind)), ons_ind'];
    vis_binary = zeros(size(sz_data,1),1);
    vis_binary(vis_chans) = 1;
    vis_binary = logical(vis_binary);
    
    clrs = repelem(0, size(sz_data,1),3);
    %clrs(ons_ind,:) = [1,0,0];
    clrs = clrs(vis_binary,:);
    
    % Plot example channels
    opts.plot_labels = true;
    opts.clrs = clrs;
    chan_names = pat_onset.channel_names{:}(vis_binary);
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


    figure(3)
    %subplot(1,2,1)
    onset_time = pat_onset.when_onset{:}(sz,1)/8;
    vis_seg_t = (512*(90)):(512*(120+pat_data(sz,:).duration));
    vis_plot_eeg((vis_seg_t/512)-120,sz_data(vis_binary,vis_seg_t),opts);
    xline(0,'g',  LineWidth=2)
    xline(onset_time, 'b', LineWidth=2)
    xline(pat_data(sz,:).duration, 'g',  LineWidth=2)
    
    %Plot imprint
    figure(4)
    %subplot(1,2,2)
    imagesc(cell_imprint{sz,1}{:}(vis_binary,:))
    set(gca, "YTick", [], "XTick",[])
    cm = [1,1,1; 1,0,0];
    colormap(cm)
end

%%
figure(5)
imagesc(final_output.imprint_roi_250{:})
colormap(cm)
