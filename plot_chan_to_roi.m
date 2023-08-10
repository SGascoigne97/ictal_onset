for pat = 1:size(final_output,1)

    if length(final_output(pat,:).Segment_ids{:}) <3
        continue
    end

    fig = figure("Position", [10,10,900,1800]);
    subplot(9,9,[1,10,19])
    imagesc(final_output(pat,:).imprint_chan{:}(:,1))
    axis off
    clim([0,3])
    
    subplot(9,9,[2:4,11:13,20:22])
    chan_2_roi_vis_1 = final_output(pat,:).imprint_chan{:}(:,1) + 2*final_output(pat,:).chan_2_roi_matrix_120{:}';
    vis_1 = chan_2_roi_vis_1;
    for roi = 1:size(chan_2_roi_vis_1,2)
        vis_1(find(chan_2_roi_vis_1(:,roi)==3, 1,'first')+1:end,roi) = 1;
    end
    vis_1(chan_2_roi_vis_1==3) = 3;
    imagesc(vis_1)
    title('Seizure One')
    axis off
    clim([0,3])
    
    subplot(9,9,29:31)
    imagesc(final_output(pat,:).imprint_roi_120{:}(:,1)')
    axis off
    clim([0,3])
    
    subplot(9,9,[6,15,24])
    imagesc(final_output(pat,:).imprint_chan{:}(:,3))
    axis off
    clim([0,3])
    
    subplot(9,9,[7:9,16:18,25:27])
    chan_2_roi_vis_3 = final_output(pat,:).imprint_chan{:}(:,3) + 2*final_output(pat,:).chan_2_roi_matrix_120{:}';
    vis_3 = chan_2_roi_vis_3;
    for roi = 1:size(chan_2_roi_vis_1,2)
        vis_3(find(chan_2_roi_vis_3(:,roi)==3, 1,'first')+1:end,roi) = 1;
    end
    vis_3(chan_2_roi_vis_3==3) = 3;
    imagesc(vis_3)
    title('Seizure Three')
    axis off
    clim([0,3])
    
    subplot(9,9, 34:36)
    imagesc(final_output(pat,:).imprint_roi_120{:}(:,3)')
    axis off
    clim([0,3])
    
    subplot(9,9, [39:43, 48:51, 57:60, 66:69, 75:78])
    imagesc(final_output(pat,:).imprint_roi_120{:}')
    axis square
    set(gca, "XTick", 1:length(final_output(pat,:).roi_names_120{:}),...
        "XTickLabel", final_output(pat,:).roi_names_120{:},...
        'TickLength',[0, 0])
    clim([0 3])
    
    sgtitle(final_output(pat,:).Patient_id{:})
    
    saveas(fig, sprintf('figures/onset_var/%s.png',final_output(pat,:).Patient_id{:}))
    saveas(fig, sprintf('figures/onset_var/%s.svg',final_output(pat,:).Patient_id{:}))
end
