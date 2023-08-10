for pat = 1:size(final_output, 1)
    pat_onset = final_output(pat,:);
    outcome_id = pat_onset.("Outcome year"){:}-pat_onset.("Surgery year") == 1;
    final_output.outcome(pat) = pat_onset.Surgery_outcome{:}(outcome_id);
end

final_output = final_output(final_output.outcome~=8,:);
%%
% low_bound = 0.01;
% up_bound = 0.25;
bound_tbl = [repelem([0.01,0.1,0.2,0.25,1/3],7);...
    repmat([0.25,1/3,0.5,2/3,0.75,0.8,1],1,5)]';
chan_or_roi = "roi_120";

for bound = 1%[1, 3, 16]%:size(bound_tbl,1)
    low_bound = bound_tbl(bound,1);
    up_bound = bound_tbl(bound,2);

    low_bound_print = round(low_bound*100);
    up_bound_print = round(up_bound*100);


    count_outside = nan(1,size(final_output, 1));
    prop_outside = nan(1,size(final_output, 1));
    prop_low_resec = nan(1,size(final_output, 1));
    prop_low_out_resec = nan(1,size(final_output, 1));
    coh_low_resec = nan(1,size(final_output, 1));
    ratio_outside = nan(1,size(final_output, 1));
    none_between_low_high = nan(1,size(final_output, 1));
    some_between_low_high = nan(1,size(final_output, 1));

    for pat = 1:size(final_output, 1)
        pat_onset = final_output(pat,:);
        imprint = pat_onset.(sprintf("imprint_%s",chan_or_roi)){:};
        imprint = imprint(:,sum(imprint)~=0);
        resec = pat_onset.(sprintf("resected_%s",chan_or_roi)){:};
        
        consensus_high = double((sum(imprint,2)/size(imprint,2))>=up_bound);
        consensus_low = double((sum(imprint,2)/size(imprint,2))>=low_bound);
        in_low_not_high = consensus_low-consensus_high;

        if sum(in_low_not_high) == 0
            none_between_low_high(pat) = pat_onset.outcome;
            continue
        end
        some_between_low_high(pat) = pat_onset.outcome;
        
        % figure("Position", [0,0,1800,900])
        % subplot(1,5,1:2)
        % imagesc(imprint)
        % title("Channel-wise imprint onset")
        % subplot(153)
        % imagesc(consensus_high)
        % title("Consensus Onset (50% threshold)")
        % subplot(154)
        % imagesc(consensus_low)
        % title("Consensus Onset (25% threshold)")
        % subplot(155)
        % imagesc(consensus_low-consensus_high)
        % title("Channels in 25%<s<50% of onsets")
        % sgtitle(sprintf("%s (ILAE %d)", pat_onset.Patient_id{:}, pat_onset.outcome))
        
        count_outside(pat) = sum(in_low_not_high); % Count of regions in lower threshold consensus onset but not upper threshold
        prop_outside(pat) = count_outside(pat)/sum(consensus_low); % Proportion of regions in lower threshold consensus onset but not upper threshold
        ratio_outside(pat) = count_outside(pat)/sum(consensus_high); % Ratio of the number of regions in lower threshold consensus onset but not upper threshold against in upper threshold only 
        prop_low_resec(pat) = sum((consensus_low+resec)==2)/sum(consensus_low); % Proportion of regions in lower threshold consensus onset but not upper threshold that were resected
        prop_low_out_resec(pat) = sum((consensus_low+(2*resec))==1)/sum(consensus_low); % Count of the number of onsets in lower threshold consensus onset outside of resection
        coh_low_resec(pat) = cohensKappa(logical(consensus_low), logical(resec)); % Consensus between lower threshold consensus onset and resection
    end
    
    fig1 = figure("Position",[0,0,1800,450]);%, 'Visible','off');
    subplot(151)
    boxchart(double(final_output.outcome>2), count_outside, "MarkerStyle","none")
    hold on
    swarmchart(double(final_output.outcome>2), count_outside, 'filled')
    hold off
    set(gca, "XTick", [0,1], "XTickLabel", ["ILAE 1-2","ILAE 3+"])
    title(sprintf("Count in %d%%<s<%d%% of onsets", low_bound_print, up_bound_print))

    subplot(152)
    boxchart(double(final_output.outcome>2), prop_outside, "MarkerStyle","none")
    hold on
    swarmchart(double(final_output.outcome>2), prop_outside, 'filled')
    hold off
    set(gca, "XTick", [0,1], "XTickLabel", ["ILAE 1-2","ILAE 3+"])
    title(sprintf("Proportion (wrt total recorded) \n in %d%%<s<%d%% of onsets", low_bound_print, up_bound_print))

    subplot(153)
    boxchart(double(final_output.outcome>2), ratio_outside, "MarkerStyle","none")
    hold on
    swarmchart(double(final_output.outcome>2), ratio_outside, 'filled')
    hold off
    set(gca, "XTick", [0,1], "XTickLabel", ["ILAE 1-2","ILAE 3+"])
    title(sprintf("Ratio (wrt size of %d%% consensus) \n in %d%%<s<%d%% of onsets", up_bound_print, low_bound_print, up_bound_print))

%     subplot(154)
%     boxchart(double(final_output.outcome>2), prop_low_resec, "MarkerStyle","none")
%     hold on
%     swarmchart(double(final_output.outcome>2), prop_low_resec, 'filled')
%     hold off
%     set(gca, "XTick", [0,1], "XTickLabel", ["ILAE 1-2","ILAE 3+"])
%     title(sprintf("Percentage of regions \n in %d%%<s<%d%% of onsets resected", low_bound_print, up_bound_print))

    subplot(154)
    boxchart(double(final_output.outcome>2), prop_low_out_resec, "MarkerStyle","none")
    hold on
    swarmchart(double(final_output.outcome>2), prop_low_out_resec, 'filled')
    hold off
    set(gca, "XTick", [0,1], "XTickLabel", ["ILAE 1-2","ILAE 3+"])
    title(sprintf("Proportion of regions \n in >%d%% of \n onsets outside resection", low_bound_print))

    subplot(155)
    boxchart(double(final_output.outcome>2), coh_low_resec, "MarkerStyle","none")
    hold on
    swarmchart(double(final_output.outcome>2), coh_low_resec, 'filled')
    hold off
    set(gca, "XTick", [0,1], "XTickLabel", ["ILAE 1-2","ILAE 3+"])
    title(sprintf("Cohen's kappa between regions \n in %d%%<s<%d%% of onsets and resection",low_bound_print, up_bound_print))
    sgtitle(sprintf("%s (Low %d%%, High %d%%)",chan_or_roi,low_bound_print, up_bound_print))

    %saveas(fig1, sprintf("figures/rapid_prototype/comparing_consensus_thresh/box_comp_%d_%d.png",low_bound_print, up_bound_print ))

    none_tab = crosstab(double(none_between_low_high(~isnan(none_between_low_high))>2))./crosstab(double(final_output.outcome)>2);
    some_tab = crosstab(double(some_between_low_high(~isnan(some_between_low_high))>2))./crosstab(double(final_output.outcome)>2);
    fig2 = figure("Position", [0,0,1800,450])%, 'Visible','off');
    bar([none_tab, some_tab], 'stacked')
    set(gca, "XTick", [1,2], "XTickLabel", ["ILAE 1-2", "ILAE 3+"], "TickLength", [0,0])
    legend([sprintf("None in %d-%d", low_bound_print, up_bound_print), sprintf(">1 in %d-%d", low_bound_print, up_bound_print)])
    sgtitle(sprintf("%s (Low %d%%, High %d%%)",chan_or_roi,low_bound_print, up_bound_print))
    %saveas(fig2, sprintf("figures/rapid_prototype/comparing_consensus_thresh/bar_comp_%d_%d.png",low_bound_print, up_bound_print))

end
