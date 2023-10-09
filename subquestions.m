

%%


%% Sub question 3: Is resecting a larger proportion of the consensus onset associated with more favourable outcomes?
onset_across_titles = ["subject", "onset across"];
det_meths = ["imprint", "EI"];
comp_measures = ["Perc", "Perc_z"];

onset_across = 1;

save_loc = "figures/subquestions/q3/";
file_type = "svg";

for parc = ["chan", "roi_120"] 
     if parc == "chan"
        comp_measures = "Perc";
    else
        comp_measures = ["Perc", "Perc_vol"];
    end
    half_violin(comp_meth_with_resec, det_meths,comp_measures, "parc", parc,...
        "onset_across", onset_across, "vis_plot", "on",...
        "save_fig", 1, "save_loc", save_loc,  "file_type", file_type)
end

for parc = ["chan", "roi_120"]
    fig = figure("Position",[10,10,900,900]);
    tiledlayout(1,2)
    for det_meth = det_meths
        nexttile
        data_tab = comp_meth_with_resec.across_sz.(sprintf(parc)).(sprintf(det_meth));
        data = data_tab.Perc;
        out = data_tab.Outcome;
       
        out_str = string();
        out_str(out<3) = "Good";
        out_str(out>2) = "Bad";
        resec_str = string();
        resec_str(data>=0.5) = "Yes";
        resec_str(data<0.5) = "No";
        [contingency_tab, ~, p, lab] = crosstab(resec_str, out_str);
        out_tab = table(resec_str', out_str');
        out_tab.Properties.VariableNames = ["resec_50", "Outcome"];
        if any(any(contingency_tab <5))
            [~, p, ~] = fishertest(contingency_tab);
            test_name = "Fishers exact";
        else 
            test_name = "Chi-square";
        end
        heatmap(out_tab, "resec_50", "Outcome")
        colorbar off
        title(sprintf("%s (%s p=%.2f)", det_meth, test_name, p))
    end
     sgtitle(sprintf("%s-wise comparison against resection (%s)", parc, onset_across_titles(onset_across+1)))
     saveas(fig, sprintf("%s%s_contingency.%s", save_loc, parc, file_type))
      
end



%% Additional: Testing association between size of onset and size of resection across outcome groups
final_output = final_output_all(cellfun(@sum, final_output_all.resected_chan)>0,:);
final_output = final_output(cellfun(@sum, final_output.clo_chan)>0,:);

regions_resec = log(cellfun(@sum,final_output.resected_roi_120));
regions_clo = log(cellfun(@sum,final_output.clo_roi_120));

[rho_f, p_f] = corr(regions_clo(final_output.outcome<3), regions_resec(final_output.outcome<3), "type","Pearson");
[rho_u, p_u] = corr(regions_clo(final_output.outcome>2), regions_resec(final_output.outcome>2), "type","Pearson");

figure("Position", [10,10,900,700])
subplot(3,3,[2,3,5,6])
gscatter(regions_resec, regions_clo, final_output.outcome>2)
text(3,1.45,sprintf("r = %.2f, p = %.2f", rho_f, p_f),'rotation',19)
text(3,0.8,sprintf("r = %.2f, p = %.2f", rho_u, p_u) ,'rotation',-2)
ylim([0,2.5])
xlim([0,5])
lsline()
legend({"ILAE 1-2", "ILAE 3+"})
ylabel("Log(Number of regions in CLO)")
xlabel("Log(Number of regions resected)")
subplot(3,3,[1,4])
histogram(regions_clo(final_output.outcome<3), "BinWidth",0.5)
hold on
histogram(regions_clo(final_output.outcome>2), "BinWidth",0.5)
hold off
xlim([0,2.5])
set(gca,'view',[90 -90], 'YDir','reverse')
axis off
subplot(3,3,8:9)
histogram(regions_resec(final_output.outcome<3), "BinWidth",0.5)
hold on
histogram(regions_resec(final_output.outcome>2), "BinWidth",0.5)
hold off
xlim([0,5])
set(gca, 'YDir','reverse')
axis off

%% Subquestion 6:  Are more diffuse onsets associated with less favorable surgical outcomes?
final_output = final_output_all(final_output_all.outcome ~= 8,:);

save_loc = "figures/subquestions/q6/";
file_type = "svg";

clear comp_tab
for subj = 1:size(final_output,1)
    subj_onset = final_output(subj,:);
    imprint = subj_onset.imprint_roi_120{:}; 
    regions = subj_onset.roi_names_120{:};

    % Create subsection of touching matrix for relevant regions
    % (n regions in hemisphere x n regions in hemisphere)
    region_index = find_preserved(names,regions);
    imprint = logical(imprint);
    clear onset_reg_count_sz
    clear max_dist
    clear onset_vol
    clear subj_tab

    max_dist = nan(size(imprint,2), 1);
    onset_vol = max_dist;
    onset_reg_count = max_dist;

    subj_max_dist = max(max(dists(region_index,region_index)));
    subj_max_vol = sum(vols(region_index));

    % Seizure specific process
    for sz = 1:size(imprint,2)% Assess one seizure at a time 
        % Extend imprint to include all regions in relevant hemisphere
        imprint_full = zeros(length(names),1);
        imprint_full(region_index(imprint(:,sz))) = 1;
    
        onset_reg_count(sz) = sum(imprint_full==1);
        max_dist(sz) = max(max(dists(imprint_full==1,imprint_full==1)));
        onset_vol(sz) = sum(vols(imprint_full==1));
    end

    subj_tab = table(repmat(subj_onset.Patient_id, size(imprint,2),1),...
        (1:size(imprint,2))', onset_reg_count, max_dist, onset_vol, ...
        repmat(subj_max_dist, size(imprint,2),1),... 
        repmat(subj_max_vol, size(imprint,2),1),... 
        'VariableNames', ["Subj_id","sz","onset_reg_count",...
        "max_dist","onset_vol", "subj_max_dist", "subj_max_vol"]);

    if exist('comp_tab', 'var')
        comp_tab = [comp_tab; subj_tab];
    else 
        comp_tab = subj_tab;
    end
end

% Create subject level table for distanes and volumes
tab_fill = nan(size(final_output,1),1);
subj_lvl_comp = table(final_output.Patient_id, tab_fill, cell(length(tab_fill),1),...
    cell(length(tab_fill),1), tab_fill, tab_fill,...
        'VariableNames', ["Subj_id","sz_count","max_dist", "onset_vol",...
        "subj_max_dist", "subj_max_vol"]);
for subj = 1:size(final_output,1)
    subj_dist = comp_tab(comp_tab.Subj_id == string(final_output(subj,:).Patient_id),:);
    subj_lvl_comp(subj,:).sz_count = size(subj_dist,1);
    for mark = ["max_dist", "onset_vol","subj_max_dist", "subj_max_vol"]
        if contains(mark, "subj")
            rep_vals = subj_dist.(sprintf(mark)); % From previous table, this values was repeated for each recorded seizure
            subj_lvl_comp(subj,:).(sprintf("%s", mark)) = rep_vals(1);
        else
            % Store raw distance and volume measures
            subj_lvl_comp(subj,:).(sprintf("%s", mark))= {subj_dist.(sprintf(mark))};
        end
    end
end

% Create summary table (with median and max for each subject)
med_dist = nan(size(final_output,1),1);
summary_tab = table(string(final_output.Patient_id), cellfun(@length,subj_lvl_comp.max_dist),...
    final_output.outcome, med_dist, med_dist, med_dist, med_dist,...
    subj_lvl_comp.subj_max_dist, subj_lvl_comp.subj_max_vol,...
    'VariableNames', {'Patient_id', 'sz_count', 'outcome', 'med_max_dist',...
    'med_onset_vol', 'max_max_dist', 'max_onset_vol', 'subj_max_dist', 'subj_max_vol'});

for mark = ["max_dist", "onset_vol"]
    for summ_meas = ["med", "max"]
        if summ_meas == "med"
            summ = @median;
        else 
            summ = @max;
        end
        summary_tab.(sprintf("%s_%s", summ_meas, mark)) = cellfun(summ, subj_lvl_comp.(sprintf(mark)));
        summary_tab.(sprintf("subj_%s_prop_%s",summ_meas, extractAfter(mark,"_"))) = ...
            summary_tab.(sprintf("%s_%s", summ_meas, mark))./summary_tab.(sprintf("subj_max_%s",extractAfter(mark,"_")));
    end
end

% For median, we only want to include subjects with at least 5 seizures
% recorded
for col = find(contains(summary_tab.Properties.VariableNames, "med"))
    vals = table2array(summary_tab(:,col));
    vals(summary_tab.sz_count <5) = NaN;
    summary_tab(:,col) = table(vals);
end

% Create violing plots
half_violin_upd(summary_tab,  ["med_max_dist","max_max_dist",...
    "subj_med_prop_dist","subj_max_prop_dist"],...
    double(summary_tab.outcome>2), "save_fig",1,...
    "file_type","svg", "save_loc",sprintf("%sdist_comp", save_loc), "grp_names", ["ILAE 1-2", "ILAE 3+"])

half_violin_upd(summary_tab,  ["med_onset_vol","max_onset_vol",...
    "subj_med_prop_vol","subj_max_prop_vol"],...
    double(summary_tab.outcome>2), "save_fig",1,...
    "file_type","svg", "save_loc",sprintf("%svol_comp", save_loc), "grp_names", ["ILAE 1-2", "ILAE 3+"])

save('tables/subquestions/dist_vol_tab.mat', 'summary_tab')

%% Subquestion 7: Is variablity in seizure onsets associated with surgical outcomes?
%final_output = final_output_all(final_output_all.outcome ~= 8,:);

save_loc = "../figures/subquestions/onset_var/";
file_type = "svg";

subj_lvl_comp_clean = subj_lvl_comp(subj_lvl_comp.sz_count>4,:);
subj_lvl_comp_clean.Properties.VariableNames{1} = 'Patient_id';

if ~any(contains(subj_lvl_comp_clean.Properties.VariableNames, "outcome"))
    subj_lvl_comp_clean = innerjoin(subj_lvl_comp_clean, final_output(:,[1,35]), 'Keys', 'Patient_id');
end

fig = figure();
tiledlayout(2,2)

for comp = ["MAD", "variance"]
    if comp == "MAD"
        comp_func = @mad;
    else
        comp_func = @var;
    end
    for column = ["max_dist", "onset_vol"]
        nexttile
        col_vals = subj_lvl_comp_clean.(sprintf(column));
        summ_vals = cellfun(comp_func, col_vals);
       
        boxchart(double(subj_lvl_comp_clean.outcome>2), summ_vals, "MarkerStyle", "none")
        hold on
        swarmchart(double(subj_lvl_comp_clean.outcome>2), summ_vals, "filled")
        hold off
        set(gca, "XTick", [0,1], "XTickLabel", ["ILAE 1-2", "ILAE 3+"])
        title(sprintf("%s %s", comp, strrep(column, "_", " ")))
    end
end
sgtitle("Seizure onset variability")
%saveas(fig, sprintf("%s_onset_var.%s", save_loc, file_type))


fig = figure();
tiledlayout(2,2)

for comp = ["MAD", "variance"]
    if comp == "MAD"
        comp_func = @mad;
    else
        comp_func = @var;
    end
    for column = ["max_dist", "onset_vol"]
        nexttile
        col_vals = subj_lvl_comp_clean.(sprintf(column));

        prop_vals = col_vals;
        for subj = 1:length(col_vals)
            prop_vals{subj,:} = col_vals{subj,:}/subj_lvl_comp_clean(subj,:).(sprintf("subj_max_%s", extractAfter(column,"_")));
        end
        
        summ_vals = cellfun(comp_func, prop_vals);
       
        boxchart(double(subj_lvl_comp_clean.outcome>2), summ_vals, "MarkerStyle", "none")
        hold on
        swarmchart(double(subj_lvl_comp_clean.outcome>2), summ_vals, "filled")
        hold off
        set(gca, "XTick", [0,1], "XTickLabel", ["ILAE 1-2", "ILAE 3+"])
        title(sprintf("%s %s", comp, strrep(column, "_", " ")))
    end
end
sgtitle("Seizure onset variability (values from proportions)")
%saveas(fig, sprintf("%s_onset_var(prop).%s", save_loc, file_type))

n_perm = 5000;
resec_thresh = 0.5;
comp_measures = ["max_dist", "onset_vol"]

%         for comp_measure = comp_measures
%             f = figure();
%             f.Position = [10,10,1000,500];
%             tiledlayout(1,3)
%             data_tab = comp_meth_with_resec.across_sz.(sprintf(chan_or_roi)).(sprintf(det_meth));
%             expl = data_tab.(sprintf(comp_measure));
%             out = data_tab.Outcome>2;
%             mod = fitglm(expl, out,'Distribution','binomial','Link','logit');
% 
%             probs = mod.Fitted.Probability;
%             [X,Y,T,AUC] = perfcurve(out,probs,1);
%         
%             perm_auc = nan(n_perm,1);
%             for perm = 1:n_perm
%                 rng(perm)
%                 expl_perm = expl(randperm(length(expl)));
%                 mod_perm = fitglm(expl_perm, out,...
%                     'Distribution','binomial','Link','logit');
%                 probs_perm = mod_perm.Fitted.Probability;
%                 [~,~,~,AUC_perm] = perfcurve(out,probs_perm,1);
%                 perm_auc(perm) = AUC_perm;
%             end
%         
%             nexttile
%             plot(X,Y)
%             hold on
%             plot([0,1], [0,1])
%             hold off
%             xlabel('False positive rate') 
%             ylabel('True positive rate')
%             title(sprintf('AUC = %.3f', AUC))
%         
%             nexttile
%             histogram(perm_auc, 'BinWidth', 0.025)
%             hold on
%             xline(AUC)
%             hold off
%             title(sprintf("p = %.3f", mean(AUC<perm_auc)))
% 
%             nexttile
%             scatter(log(probs./(1-probs)), expl)
%             lsline()
%             title("Linearity assumption")
% 
%             sgtitle(strrep(sprintf("%s %s %s", chan_or_roi, det_meth,comp_measure), "_", " "))
%             saveas(f, sprintf("%s%s_%s_AUC.svg", save_loc))
%   
%         end
% 
% 






%% Alternative subquestion 7: Is concordance between seizure onsets associated with outcomes
% Concordance within seizures
concord_mat_tab = table(final_output.Patient_id, final_output.outcome,...
    nan(size(final_output,1),1), cell(size(final_output,1),1),...
    nan(size(final_output,1),1), nan(size(final_output,1),1),...
    nan(size(final_output,1),1),'VariableNames',...
    ["Subj_id", "Outcome", "Sz_count", "concord_mat",...
    "mean_concord", "med_concord", "max_concord"]);
for subj = 1:size(final_output,1)
    subj_onset = final_output(subj,:);
    for det_meth = ["imprint", "EI"]
        onset = subj_onset.(sprintf("%s_roi_120", det_meth)){:};
        sz_count = size(onset,2);
        concord_mat_tab(subj,:).Sz_count = sz_count;
        concordance = nan(sz_count,sz_count);
        for sz1 = 1:sz_count
            for sz2 = 1:sz_count
                if sz1<sz2
                    concordance(sz1,sz2) = cohensKappa(onset(:,sz1), onset(:,sz2));
                end

            end
        end 
        concord_mat_tab(subj,:).concord_mat = {concordance};
        concord_mat_tab(subj,:).mean_concord = mean(mean(concordance, 'omitnan'), 'omitnan');
        concord_mat_tab(subj,:).med_concord = median(median(concordance, 'omitnan'), 'omitnan');
        concord_mat_tab(subj,:).max_concord = max(max(concordance));
    end
end


concord_mat_tab = concord_mat_tab(concord_mat_tab.Sz_count>=5,:);

%%
half_violin_upd(concord_mat_tab, ["mean_concord", "med_concord", "max_concord"],...
    double(concord_mat_tab.Outcome>2), "y_lim", [-0.55, 1.05], "save_fig",1,...
    "file_type","svg", "save_loc",sprintf("%ssz_concord", save_loc), "grp_names", ["ILAE 1-2", "ILAE 3+"])

save('tables/subquestions/concord_mat_tab.mat', 'concord_mat_tab')
%%
figure(1)
tiledlayout(2,2)
for col_name = ["mean_concord", "med_concord"]
    nexttile
    histogram(concord_mat_tab(concord_mat_tab.Outcome<3,:).(sprintf(col_name)), "BinWidth", 0.2)
    xlim([-0.5 1])
    title(sprintf("%s ILAE 1-2", strrep(col_name, "_", " ")))
    nexttile
    histogram(concord_mat_tab(concord_mat_tab.Outcome>2,:).(sprintf(col_name)), "BinWidth", 0.2)
    xlim([-0.5 1])
    title(sprintf("%s ILAE 3+", strrep(col_name, "_", " ")))
end


%%
median(concord_mat_tab(concord_mat_tab.Outcome<3,:).mean_concord, 'omitnan')
median(concord_mat_tab(concord_mat_tab.Outcome>2,:).mean_concord, 'omitnan')


median(concord_mat_tab(concord_mat_tab.Outcome<3,:).med_concord, 'omitnan')
median(concord_mat_tab(concord_mat_tab.Outcome>2,:).med_concord, 'omitnan')


median(concord_mat_tab(concord_mat_tab.Outcome<3,:).max_concord, 'omitnan')
median(concord_mat_tab(concord_mat_tab.Outcome>2,:).max_concord, 'omitnan')



