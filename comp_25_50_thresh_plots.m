for pat = 1:size(final_output,1)
    pat_onset = final_output(pat,:);
    roi_120 = pat_onset.roi_names_120{:};
    
    roi_120 = strrep(roi_120, 'r.', 'ctx-rh-');
    roi_120 = strrep(roi_120, 'l.', 'ctx-lh-');
    
    imprint_120 = pat_onset.imprint_roi_120{:};
    
    prop_sz = sum(imprint_120,2)/size(imprint_120,2);
    
    over_thresh_25 = prop_sz>=0.25;
    over_thresh_50 = prop_sz>=0.5;
    
    comp_thresh = over_thresh_25 + over_thresh_50;
    cm = [.7,.7,.7;0,1,0;1,1,0];
     cm = interp1(cm, 1:0.01:size(cm,1));
    plotBrain(roi_120, comp_thresh, cm, "atlas", 'lausanne120_aseg',...
        'savePath', sprintf('figures/rapid_prototype/comp_onset_across_thresh/vis_brain/%s_imprint',pat_onset.Patient_id{:}))
    plotBrain(roi_120, pat_onset.clo_roi_120{:}, [.7,.7,.7; 1,0,0], "atlas", 'lausanne120_aseg',...
        'savePath', sprintf('figures/rapid_prototype/comp_onset_across_thresh/vis_brain/%s_clo',pat_onset.Patient_id{:}))
end

%%
col_names =["Perc", "Perc_z", "Jacc", "Jacc_z", "Coh", "Coh_z"];
jitt = rand(size(comp_meth_with_clo_25.across_sz.roi_120.imprint,1),1)/10;
%%
for col = col_names
    f = figure("Position",[10,10,900,900]);
    thresh_25 = comp_meth_with_clo_25.across_sz.roi_120.imprint.(sprintf("%s",col));
    thresh_50 = comp_meth_with_clo_50.across_sz.roi_120.imprint.(sprintf("%s",col));
    subplot(3,3,[2,3,5,6])
    scatter(thresh_25+jitt, thresh_50+jitt, "b", "filled")
    l = lsline();
    set(l, 'color', 'b')

    if ~contains(col, "z")
        if col ~= "Coh"
            x_lim = [0,1.1];
        else
            x_lim = [-0.5,1.1];
        end
        bin_wid = 0.1;
    else 
        x_lim = [min(min(thresh_25),min(thresh_50))*1.15, max(max(thresh_25),max(thresh_50))*1.15];
        bin_wid = 0.5;
        hold on
        patch([-2,2,2,-2],[-2,-2,2,2],[0.7,0.7,0.7],'FaceAlpha',.3, 'EdgeColor', "none")
        text(0, 1.5, "Chance", "HorizontalAlignment", "center")
        hold off
        xline(2, '--')
        yline(2, '--')
    end

    xlabel("Region in 25% of seizures")
    ylabel("Region in 50% of seizures")
    xlim(x_lim)
    ylim(x_lim)
    subplot(3,3,[1,4])
    histogram(thresh_50, "BinWidth",bin_wid, "FaceColor","b")
    xlim(x_lim)
    camroll(90)
    axis off
    subplot(3,3,[8,9])
    histogram(thresh_25, "BinWidth",bin_wid, "FaceColor","b")
    xlim(x_lim)
    set(gca, 'YDir','reverse')
    axis off
    sgtitle(sprintf("Comparing onset across (%s) using 25%% and 50%% threshold",strrep(col,'_',' ')))
    saveas(f, sprintf('figures/rapid_prototype/comp_onset_across_thresh/%s.png', col))
    %saveas(f, sprintf('figures/rapid_prototype/comp_onset_across_thresh/%s.svg', col))

end


for col = col_names
    f = figure("Position",[10,10,900,900]);
    thresh_25 = comp_meth_with_clo_25.across_sz.roi_120.imprint.(sprintf("%s",col));
    thresh_50 = comp_meth_with_clo_50.across_sz.roi_120.imprint.(sprintf("%s",col));
    out = comp_meth_with_clo_25.across_sz.roi_120.imprint.Outcome;
    subplot(3,3,[2,3,5,6])
    scatter(thresh_25(out<3)+jitt(out<3), thresh_50(out<3)+jitt(out<3), "b", "filled")
    hold on
    scatter(thresh_25(out>2)+jitt(out>2), thresh_50(out>2)+jitt(out>2), "r", "filled")
    hold off
    l = lsline();
    set(l(1), 'color', 'b')
    set(l(2), 'color', 'r')
    if ~contains(col, "z")
        if col ~= "Coh"
            x_lim = [0,1.1];
        else
            x_lim = [-0.5,1.1];
        end
        bin_wid = 0.1;
    else 
        x_lim = [min(min(thresh_25),min(thresh_50))*1.25, max(max(thresh_25),max(thresh_50))*1.25];
        bin_wid = 0.5;

        hold on
        patch([-2,2,2,-2],[-2,-2,2,2],[0.7,0.7,0.7],'FaceAlpha',.3, 'EdgeColor', "none")
        text(0, 1.5, "Chance", "HorizontalAlignment", "center")
        hold off
        xline(2, '--')
        yline(2, '--')
    end

    xlabel("Region in 25% of seizures")
    ylabel("Region in 50% of seizures")
    xlim(x_lim)
    ylim(x_lim)
    subplot(3,3,[1,4])
    histogram(thresh_50(out<3), "BinWidth",bin_wid, "FaceColor","b")
    camroll(90)
    hold on
    histogram(thresh_50(out>2), "BinWidth",bin_wid, "FaceColor","r")
    hold off
    xlim(x_lim)
    axis off
    subplot(3,3,[8,9])
    histogram(thresh_25(out<3), "BinWidth",bin_wid, "FaceColor","b")
    hold on
    histogram(thresh_25(out>2), "BinWidth",bin_wid, "FaceColor","r")
    hold off
    xlim(x_lim)
    set(gca, 'YDir','reverse')
    axis off
    sgtitle(sprintf("Comparing onset across (%s) using 25%% and 50%% threshold",strrep(col,'_',' ')))
    saveas(f, sprintf('figures/rapid_prototype/comp_onset_across_thresh/%s_outcome.png', col))
    %saveas(f, sprintf('figures/rapid_prototype/comp_onset_across_thresh/%s_outcome.svg', col))

end

