addpath(genpath('sarah_functions/'))
%% Concordance between onsets (both CLO and automatically detected)
final_output = final_output_all;

summary_meths = ["median", "max"];
chan_or_roi = "roi_120";
for summ = 1:2
    summary_meth = summary_meths(summ);
    roi_tab = table(final_output.Patient_id, final_output.Segment_ids, final_output.(sprintf("clo_%s", chan_or_roi)),...
    final_output.(sprintf("imprint_%s", chan_or_roi)),final_output.(sprintf("EI_%s", chan_or_roi)), ...
    final_output.(sprintf("PLHG_%s", chan_or_roi)), final_output.(sprintf("resected_%s", chan_or_roi)), final_output.outcome);
    col_names = ["CLO", "Imprint", "EI", "PLHG", "Resec"];
    roi_tab.Properties.VariableNames = ["Patient_id", "Segment_ids", col_names, "outcome"];
    % If we are using median we need each individual to have 5 recorded
    % seizures
    if summary_meth == "median"
        roi_tab = roi_tab(cellfun(@iscell, roi_tab.Segment_ids),:);
        roi_tab = roi_tab(cellfun(@length, roi_tab.Segment_ids)>=5,:);
    end
    
    concord_mat = zeros(length(col_names),length(col_names), size(roi_tab,1));
    for pat = 1:size(roi_tab,1)
        pat_roi = roi_tab(pat,:);
        for col1 = 1:length(col_names)
            for col2 = 1:length(col_names)
                mat1 = pat_roi.(sprintf("%s", col_names(col1))){:};
                mat2 = pat_roi.(sprintf("%s", col_names(col2))){:};
                concord = zeros(max(size(mat1,2), size(mat2,2)),1);
                if size(mat1,2) == size(mat2,2)
                    for sz = 1:size(mat1,2)
                        concord(sz) = cohensKappa(mat1(:,sz), mat2(:,sz));
                    end
                elseif size(mat1,2) == 1 & size(mat2,2) > 1
                    for sz = 1:max(size(mat1,2), size(mat2,2))
                        concord(sz) = cohensKappa(mat1, mat2(:,sz));
                    end
                elseif size(mat1,2) > 1 & size(mat2,2) == 1
                    for sz = 1:max(size(mat1,2), size(mat2,2))
                        concord(sz) = cohensKappa(mat1(:,sz), mat2);
                    end
                end
                % Store concordance (if between CLO and resec) or median concordance across
                % seizures (when comparison includes automatic detection method)
                if length(concord) == 1
                    concord_mat(col1,col2,pat) = concord;
                else
                    if summary_meth == "median"
                        concord_mat(col1,col2,pat) = median(concord);
                    else
                        concord_mat(col1,col2,pat) = max(concord);
                    end
                end
            end
        end
    end
    
    
    f = figure('Position', [0,0,1800,900]);
    tiledlayout(5,5)
    for col1 = 1:length(col_names)
        for col2 = 1:length(col_names)
            nexttile
            if col1 <= col2
                if col1 ~= col2
                    concord_array = squeeze(concord_mat(col1,col2,:));
                    histogram(concord_array(concord_array>0.6), 'BinWidth',0.1,'FaceColor','g')
                    hold on
                    histogram(concord_array(concord_array<=0.6), 'BinWidth',0.1, 'FaceColor',[0.8500 0.3250 0.0980])
                    histogram(concord_array(concord_array<0.4), 'BinWidth',0.1, 'FaceColor', 'r')
                    hold off
                    xlim([-0.5 1])
                    xline(0,'--', 'no agreement', 'LabelHorizontalAlignment','center')
                    xline(0.2,'--', 'none to slight', 'LabelHorizontalAlignment','center')
                    xline(0.4,'--', 'fair', 'LabelHorizontalAlignment','center')
                    xline(0.6,'--', 'moderate', 'LabelHorizontalAlignment','center')
                    xline(0.8, '--', 'substantial', 'LabelHorizontalAlignment','center')
                    ylim([0 15])
                    xlabel("ILAE Outcome")
                    ylabel("Concordance")
                    title(sprintf("%s and %s", col_names(col1), col_names(col2)))
                else
                    scatter(1,1)
                    text(0,0,col_names(col1))
                    xlim([-.05, .05])
                    ylim([-.05, .05])
                    axis off
                end
            else
                concord_array = squeeze(concord_mat(col1,col2,:));
                scatter(ones(1,3),1:3, [], [0,1,0; 0.85,0.325,0.098;1,0,0], 'filled')
                text(2*ones(1,3),1:3,[string(length(concord_array(concord_array>0.6))),...
                    string(length(concord_array(concord_array<=0.6 & concord_array>0.4))),...
                    string(length(concord_array(concord_array<=0.4)))])
                 xlim([-1 3])
                 ylim([-3 6])
                axis off
            end
        end
    end
    sgtitle(sprintf("Within-subject %s Cohen's Kappa (n=%d)", summary_meth, size(roi_tab,1)))
    saveas(f, sprintf('figures/rapid_prototype/within_patient_%s_cohen_kappa_%s.png', summary_meth, chan_or_roi))
end
%%

%% Concordance between seizure onsets within the same onset detection method
chan_or_roi = "chan";
final_output = final_output_all(cellfun(@iscell, final_output_all.Segment_ids),:);
final_output = final_output(cellfun(@length, final_output.Segment_ids)>4,:);


roi_tab = table(final_output.Patient_id, final_output.Segment_ids,...
    final_output.(sprintf("imprint_%s", chan_or_roi)),...
    final_output.(sprintf("EI_%s", chan_or_roi)), ...
    final_output.(sprintf("PLHG_%s", chan_or_roi)), final_output.outcome);

col_names = ["Imprint", "EI", "PLHG"];
roi_tab.Properties.VariableNames = ["Patient_id", "Segment_ids", col_names, "outcome"];
% We require at least 5 seizures to be included in this analysis
roi_tab = roi_tab(cellfun(@iscell, roi_tab.Segment_ids),:);
roi_tab = roi_tab(cellfun(@length, roi_tab.Segment_ids)>=5,:);

within_method_med_concord_mat = zeros(length(col_names), size(roi_tab,1));
within_method_max_concord_mat = zeros(length(col_names), size(roi_tab,1));

within_method_concord_tab = roi_tab;

%%
for pat = 1:size(roi_tab,1)
    pat_roi = roi_tab(pat,:);
    for col = 1:length(col_names)
        mat = pat_roi.(sprintf("%s", col_names(col))){:};
        concord = nan(max(size(mat,2), size(mat,2)),max(size(mat,2), size(mat,2)));
        for sz1 = 1:size(mat,2)
            for sz2 = 1:size(mat,2)
                if sz1<sz2
                    concord(sz1,sz2) = cohensKappa(mat(:,sz1), mat(:,sz2));
                end
            end
        end

            % Compute median and maximum concordance on patient and onset
            % detection method level
            within_method_med_concord_mat(col,pat) = median(concord, 'all', 'omitnan');
            within_method_max_concord_mat(col,pat) = max(max(concord, [],'omitnan'), [], 'omitnan');
            concord_array = concord(:);
            within_method_all_sz(col,pat) = {concord_array(~isnan(concord_array))};
    end
end
   
    
    
f2 = figure('Position', [0,0,1800,900]);
tiledlayout(2,3)
for col = 1:length(col_names)
    nexttile
    concord_array = squeeze(within_method_med_concord_mat(col,:));
    histogram(concord_array(concord_array>0.6), 'BinWidth',0.1,'FaceColor','g')
    hold on
    histogram(concord_array(concord_array<=0.6), 'BinWidth',0.1, 'FaceColor',[0.8500 0.3250 0.0980])
    histogram(concord_array(concord_array<0.4), 'BinWidth',0.1, 'FaceColor', 'r')
    hold off
    xlim([-0.5 1])
    xline(0,'--', 'no agreement', 'LabelHorizontalAlignment','center')
    xline(0.2,'--', 'none to slight', 'LabelHorizontalAlignment','center')
    xline(0.4,'--', 'fair', 'LabelHorizontalAlignment','center')
    xline(0.6,'--', 'moderate', 'LabelHorizontalAlignment','center')
    xline(0.8, '--', 'substantial', 'LabelHorizontalAlignment','center')
    ylim([0 25])
    xlabel("ILAE Outcome")
    ylabel("Concordance")
    title(sprintf("Median %s", col_names(col)))
end
sgtitle(sprintf("Within-subject median/max of Cohen's Kappa between seizures (n=%d)", size(roi_tab,1)))

for col = 1:length(col_names)
    nexttile
    concord_array = squeeze(within_method_max_concord_mat(col,:));
    histogram(concord_array(concord_array>0.6), 'BinWidth',0.1,'FaceColor','g')
    hold on
    histogram(concord_array(concord_array<=0.6), 'BinWidth',0.1, 'FaceColor',[0.8500 0.3250 0.0980])
    histogram(concord_array(concord_array<0.4), 'BinWidth',0.1, 'FaceColor', 'r')
    hold off
    xlim([-0.5 1])
    xline(0,'--', 'no agreement', 'LabelHorizontalAlignment','center')
    xline(0.2,'--', 'none to slight', 'LabelHorizontalAlignment','center')
    xline(0.4,'--', 'fair', 'LabelHorizontalAlignment','center')
    xline(0.6,'--', 'moderate', 'LabelHorizontalAlignment','center')
    xline(0.8, '--', 'substantial', 'LabelHorizontalAlignment','center')
    ylim([0 25])
    xlabel("ILAE Outcome")
    ylabel("Concordance")
    title(sprintf("Maximim %s", col_names(col)))
end
saveas(f2, sprintf('figures/rapid_prototype/within_patient_%s_seizure_concordance_%s.svg', summary_meth, chan_or_roi))
%%
det_methods = ["Imprint", "EI", "PLHG"];
clr = [1,0,0; 0,1,0; 0,0,1];
sz_counts = cellfun(@length,roi_tab.Segment_ids);
for pat = 1:size(roi_tab,1)
    within_pat_concord = within_method_all_sz(:,pat);

    f = figure('Position', [0,0,1800,900], 'Visible','off');
    tiledlayout(3,1)
    for meth = 1:3
        nexttile
        histogram(within_pat_concord{meth,1}, 'BinWidth', 0.1, 'FaceColor',clr(meth,:))
        title(det_methods(meth))
        xlim([-0.2, 1.2])
        xline([0.6,0.4],'--')
    end
  
    sgtitle(sprintf("%s (ILAE %d) between-seizure concordance (n=%d)",roi_tab.Patient_id{pat}, roi_tab.outcome(pat), sz_counts(pat)))
    saveas(f, sprintf('figures/rapid_prototype/within_sz_concord/%s_between_sz_concord.svg', roi_tab.Patient_id{pat}))
end
% Need to clean up these plots but it is looking interesting
% How can we analyse this data??

%% Plot all comparison measures across channels and ROIs
for chan_or_roi = ["chan", "roi_120", "roi_250"]
    jacc = nan(size(final_output,1),1);
    perc = jacc;
    coh = jacc;
    jacc_z = jacc;
    perc_z = jacc;
    coh_z = jacc;
    
    n_perm = 1000;
    for pat = 1:size(final_output,1)
        pat_onset = final_output(pat,:);
        clo = pat_onset.(sprintf("clo_%s", chan_or_roi)){:};
        if sum(clo) == 0
            continue
        end
        resec = pat_onset.(sprintf("resected_%s", chan_or_roi)){:};
        if sum(resec) == 0
            continue
        end
        jacc(pat) = jaccard(clo,resec);
        perc(pat) = sum(resec + clo == 2)/sum(clo);
        coh(pat) = cohensKappa(clo,resec);
    
        % Perform permutation test 
        jacc_perm = zeros(1,n_perm);
        perc_perm = jacc_perm; 
        for perm = 1:n_perm
            rng(perm)
            perm_clo = clo(randperm(length(clo)));
            rng(perm+1)
            perm_resec = resec(randperm(length(clo)));
            jacc_perm(perm) = jaccard(perm_clo,perm_resec);
            perc_perm(perm) = sum(perm_resec + perm_clo == 2)/sum(perm_clo);
            coh_perm(perm) = cohensKappa(perm_clo,perm_resec);
        end
    
        jacc_z(pat) = (jacc(pat) - mean(jacc_perm))/std(jacc_perm);
        perc_z(pat) = (perc(pat) - mean(perc_perm))/std(perc_perm);
        coh_z(pat) = (coh(pat) - mean(coh_perm))/std(coh_perm);
    
    end
    
    f1 = figure('Position', [0,0,1800,900]);
    subplot(321)
    histogram(perc, 'BinWidth',0.1)
    xlim([0 1])
    title("Distribution of Percentage of CLO resected")
    
    subplot(322)
    histogram(perc_z, 'BinWidth',0.5)
    xlim([-3 6.5])
    title("Distribution of Percentage of CLO resected (z)")
    
    subplot(323)
    histogram(jacc, 'BinWidth',0.1)
    xlim([0 1])
    title("Distribution of Jaccard's index between CLO and resection")
    
    subplot(324)
    histogram(jacc_z, 'BinWidth',0.5)
    title("Distribution of Jaccard's index between CLO and resection (z)")
    xlim([-3 6.5])
    
    subplot(325)
    histogram(coh, 'BinWidth',0.1)
    title("Distribution of Cohen's kappa between CLO and resection")
    xlim([-0.5 1])
    
    subplot(326)
    histogram(coh_z, 'BinWidth',0.5)
    title("Distribution of Cohen's kappa between CLO and resection (z)")
    xlim([-3 6.5])
    sgtitle(sprintf("Comparing CLO and resection (%s)", chan_or_roi))
    saveas(f1, sprintf('figures/clo_vs_resec/%s.png',chan_or_roi))
    
    f2 = figure('Position', [0,0,1800,900]);
    subplot(321)
    histogram(perc(final_output.outcome <=2), 'BinWidth',0.1)
    hold on
    histogram(perc(final_output.outcome >2), 'BinWidth',0.1)
    hold off
    xlim([0 1])
    title("Distribution of Percentage of CLO resected")
    
    subplot(322)
    histogram(perc_z(final_output.outcome <=2), 'BinWidth',0.5)
    hold on
    histogram(perc_z(final_output.outcome >2), 'BinWidth',0.5)
    hold off
    xlim([-3 6.5])
    legend(["ILAE 1-2", "ILAE 3+"])
    title("Distribution of Percentage of CLO resected (z)")
    
    subplot(323)
    histogram(jacc(final_output.outcome <=2), 'BinWidth',0.1)
    hold on
    histogram(jacc(final_output.outcome >2), 'BinWidth',0.1)
    hold off
    xlim([0 1])
    title("Distribution of Jaccard's index between CLO and resection")
    
    subplot(324)
    histogram(jacc_z(final_output.outcome <=2), 'BinWidth',0.5)
    hold on
    histogram(jacc_z(final_output.outcome >2), 'BinWidth',0.5)
    hold off
    title("Distribution of Jaccard's index between CLO and resection (z)")
    xlim([-3 6.5])
    
    subplot(325)
    histogram(coh(final_output.outcome <=2), 'BinWidth',0.1)
    hold on
    histogram(coh(final_output.outcome >2), 'BinWidth',0.1)
    hold off
    title("Distribution of Cohen's kappa between CLO and resection")
    xlim([-0.5 1])
    
    subplot(326)
    histogram(coh_z(final_output.outcome <=2), 'BinWidth',0.5)
    hold on
    histogram(coh_z(final_output.outcome >2), 'BinWidth',0.5)
    hold off
    title("Distribution of Cohen's kappa between CLO and resection (z)")
    xlim([-3 6.5])
    sgtitle(sprintf("Comparing CLO and resection (%s) across outcomes", chan_or_roi))
    saveas(f2, sprintf('figures/clo_vs_resec/%s_outcomes.png',chan_or_roi))

    f3 = figure('Position', [0,0,1800,900]);
    subplot(321)
    boxchart(final_output.outcome, perc, "MarkerStyle", "none")
    hold on
    swarmchart(final_output.outcome, perc, 'filled')
    hold off
    xlabel("Outcome")
    title("Percentage of CLO resected")
    ylim([0 1])
    xlim([0 6])
    
    subplot(322)
    boxchart(final_output.outcome, perc_z, "MarkerStyle", "none")
    hold on
    swarmchart(final_output.outcome, perc_z, 'filled')
    hold off
    xlabel("Outcome")
    title("Percentage of CLO resected (z)")
    ylim([-3 6.5])
    xlim([0 6])

    subplot(323)
    boxchart(final_output.outcome, jacc, "MarkerStyle", "none")
    hold on
    swarmchart(final_output.outcome, jacc, 'filled')
    hold off
    xlabel("Outcome")
    ylim([0 1])
    xlim([0 6])    
    title("Jaccard's index between CLO and resection")
    
    subplot(324)
    boxchart(final_output.outcome, jacc_z, "MarkerStyle", "none")
    hold on
    swarmchart(final_output.outcome, jacc_z, 'filled')
    hold off
    xlabel("Outcome")
    title("Jaccard's index between CLO and resection (z)")
    ylim([-3 6.5])
    xlim([0 6])
    
    subplot(325)
    boxchart(final_output.outcome, coh, "MarkerStyle", "none")
    hold on
    swarmchart(final_output.outcome, coh, 'filled')
    hold off
    xlabel("Outcome")
    title("Cohen's kappa between CLO and resection")
    ylim([-0.5 1])
    xlim([0 6])
    
    subplot(326)
    boxchart(final_output.outcome, coh_z, "MarkerStyle", "none")
    hold on
    swarmchart(final_output.outcome, coh_z, 'filled')
    hold off
    title("Cohen's kappa between CLO and resection (z)")
    ylim([-3 6.5])
    xlim([0 6])

    sgtitle(sprintf("Comparing CLO and resection (%s) across outcomes", chan_or_roi))
    saveas(f3, sprintf('figures/clo_vs_resec/%s_outcome_scatter.svg',chan_or_roi))
end