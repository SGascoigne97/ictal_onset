final_output_all = load('tables/final_output_all_sz.mat');
final_output_all = final_output_all.final_output;
addpath(genpath('sarah_functions'))

% Remove patients with no labelled CLO or no outcome 
final_output = final_output_all(final_output_all.outcome ~= 8,:);


%% Sub question 1: Does seizure onset tend to be resected?
det_meths = ["clo", "imprint", "EI"];
comp_measures = ["Perc", "Perc_z"];

onset_acr_thresh = 0.5;
% Create an empty table to store output
%fprintf("%s\n", onset_across_title(onset_across))
atl = [72,120,250,500];

% Create structure with all comparisons between onset (CLO and automatic)
% and resection
parcs = ["chan", "roi_" + atl([2,3])];

for parc = parcs
    vols = atlas.vol{atl_ind};
    for det_meth = det_meths
        comp_meth_tab = compute_resec_comp(final_output,atlas(atl_ind,:),vols,...
            "chan_or_roi", parc, "det_meth", det_meth,...
            "onset_acr_thresh",onset_acr_thresh,"onset_across",onset_across);
        comp_meth_tab.outcome_cat = categorical(comp_meth_tab.Outcome>2,[0,1], ["ILAE 1-2", "ILAE 3+"]);
        comp_meth_with_resec.(sprintf("%s", onset_across_title(onset_across))).(sprintf("%s", parc)).(sprintf("%s", det_meth)) = comp_meth_tab;
    end
end

save_loc = "figures/subquestions/q1/";
file_type = "svg";

for parc = ["chan", "roi_120"]
    if parc == "chan"
        comp_measures = "Perc";
    else
        comp_measures = ["Perc", "Perc_vol"];
    end
    half_violin(comp_meth_with_resec, det_meths, comp_measures, "parc", parc,...
        "onset_across", 1, "vis_plot", "on",...
        "save_fig", 1, "save_loc", save_loc,  "file_type", file_type)
end

%%
var_names = {'Parc','Det','All','Total','Fav','Fav_total', 'Unfav', 'Unfav_total'};
onset_resec_tab = array2table(nan(0,8),'VariableNames',var_names);


for parc = ["chan", "roi_120", "roi_250"]
    for det = det_meths
        data_tab = comp_meth_with_resec.across_sz.(sprintf("%s", parc)).(sprintf("%s", det));
        data = data_tab.Perc;
        out = data_tab.Outcome;
        
        acr = sum(data>=0.5);
        fav = sum(data(out<3)>=0.5);
        unfav = sum(data(out>2)>=0.5);

        param_tab = array2table([parc, det], "VariableNames", var_names(1:2));
        val_tab = array2table([acr, length(data), fav, length(data(out<3)), unfav, length(data(out>2))], "VariableNames", var_names(3:end));
        onset_resec_tab = [onset_resec_tab;param_tab, val_tab];
    end
end

onset_resec_tab.All_perc = onset_resec_tab.All./onset_resec_tab.Total;
onset_resec_tab.Fav_perc = onset_resec_tab.Fav./onset_resec_tab.Fav_total;
onset_resec_tab.Unfav_perc = onset_resec_tab.Unfav./onset_resec_tab.Unfav_total;

%% Sub question 2: Is resecting a larger proportion of the clinically 
% labelled onset (based on volume) associated with more favourable outcomes?

save_loc = "figures/subquestions/q2/";

comp_outcome = 0;
onset_across = 1;
onset_across_title = ["per_sz", "across_sz"];
comp_outcome_title = ["", "_comp_outcome"];
n_perm = 1000;

det_meths = ["clo", "imprint", "EI"];

onset_acr_thresh = 0.5;
% Create an empty table to store output
fprintf("%s\n", onset_across_title(onset_across))
chan_or_roi = "roi_120";

final_output = final_output_all;

for det_meth = det_meths
    comp_meth_tab = compute_resec_vol_comp(final_output, atlas,...
        "chan_or_roi", chan_or_roi, "det_meth", det_meth, "n_perm",n_perm,...
        "onset_acr_thresh",onset_acr_thresh,"onset_across",onset_across);
    
        %comp_meth_tab = comp_meth_tab(comp_meth_tab.Sz_count >4,:);
        comp_meth_with_resec_vol.(sprintf("%s", onset_across_title(onset_across+1))).(sprintf("%s", chan_or_roi)).(sprintf("%s", det_meth)) = comp_meth_tab;
end
%%

atl = [72,120,250,500];
atl_ind = 2;
chan_or_roi = sprintf("roi_%d", atl(atl_ind));
vols = atlas.vol{atl_ind};


for det_meth = det_meths
    comp_meth_tab = compute_resec_comp(final_output,atlas(atl_ind,:),vols,...
        "chan_or_roi", chan_or_roi, "det_meth", det_meth, "n_perm",n_perm,...
        "onset_acr_thresh",onset_acr_thresh,"onset_across",onset_across);
        comp_meth_with_resec.(sprintf("%s", onset_across_title(onset_across))).(sprintf("%s", chan_or_roi)).(sprintf("%s", det_meth)) = comp_meth_tab;
    
    comp_meth_tab.outcome_cat = categorical(comp_meth_tab.Outcome>2,[0,1], ["ILAE 1-2", "ILAE 3+"]);
    
   
    figure("Position", [100,1000,800,800])
    sgtitle(det_meth)
    subplot(3,3,[2,3,5,6])
    hold on
    scatter(comp_meth_tab(comp_meth_tab.outcome_cat == "ILAE 1-2",:).Perc, comp_meth_tab(comp_meth_tab.outcome_cat == "ILAE 1-2",:).Perc_vol, 'b',"filled", "XJitter","randn", "YJitter","randn", "XJitterWidth",0.1,"YJitterWidth",0.1)
    scatter(comp_meth_tab(comp_meth_tab.outcome_cat == "ILAE 3+",:).Perc, comp_meth_tab(comp_meth_tab.outcome_cat == "ILAE 3+",:).Perc_vol, 'r',"filled", "XJitter","randn", "YJitter","randn", "XJitterWidth",0.1,"YJitterWidth",0.1)
    plot([0,1], [0,1], 'k--')   
    hold off
    xlabel("Perc (binary)")
    ylabel("Perc (vol)")
    xlim([0,1.1])
    ylim([0,1.1])

    subplot(3,3,[1,4])
    hold on
    histogram(comp_meth_tab(comp_meth_tab.outcome_cat == "ILAE 1-2",:).Perc_vol, "FaceColor", 'b', "BinWidth",0.05)
    histogram(comp_meth_tab(comp_meth_tab.outcome_cat == "ILAE 3+",:).Perc_vol, "FaceColor", 'r', "BinWidth",0.05)
    hold off
    xlim([0,1.1])
    set(gca, 'view',[90, -90],'YDir','reverse')
    axis off
  

    subplot(3,3,[8,9])
    hold on
    histogram(comp_meth_tab(comp_meth_tab.outcome_cat == "ILAE 1-2",:).Perc, "FaceColor", 'b', "BinWidth",0.05)
    histogram(comp_meth_tab(comp_meth_tab.outcome_cat == "ILAE 3+",:).Perc, "FaceColor", 'r', "BinWidth",0.05)
    hold off
    xlim([0,1.1])
    set(gca, 'YDir','reverse')
    axis off
end

%% 
onset_across = 1;
comp_measures = ["Perc", "Perc_vol"];

for parc = ["roi_120"]%, "roi_250"] 
    half_violin(comp_meth_with_resec, det_meths, comp_measures, "parc", parc,...
        "onset_across", onset_across, "vis_plot", "on",...
        "save_fig", 0)
end

%%
resec_thresh = 0.5;
for comp_measure = comp_measures
    fig = figure("Position",[100,1000, 1200, 400]);
    tiledlayout(1,length(det_meths))
    for det_meth = det_meths
        nexttile
        data_tab = comp_meth_with_resec.across_sz.roi_120.(sprintf(det_meth));
        data = data_tab.(sprintf(comp_measure));
        out = data_tab.Outcome;
       
        out_str = string();
        out_str(out<3) = "Good";
        out_str(out>2) = "Bad";
        resec_str = string();
        resec_str(data>=resec_thresh) = "Yes";
        resec_str(data<resec_thresh) = "No";
        [contingency_tab, ~, p, lab] = crosstab(resec_str, out_str);
        out_tab = table(resec_str', out_str');
        out_tab.Properties.VariableNames = ["resec_thresh", "Outcome"];
        if any(any(contingency_tab <5))
            [~, p, ~] = fishertest(contingency_tab);
            test_name = "Fishers exact";
        else 
            test_name = "Chi-square";
        end
        heatmap(out_tab, "resec_thresh", "Outcome")
        colorbar off
        title(sprintf("%s (%s p=%.2f)", det_meth, test_name, p))
    end
    sgtitle(sprintf("%s (thresh = %d%%)", strrep(comp_measure, "_", " "),resec_thresh*100))
    saveas(fig, sprintf("%scontingency_%s_%d.%s", save_loc,comp_measure,resec_thresh*100,file_type))
end

% OLD CODE
% onset_across_titles = ["subject", "onset across"];
% det_meths = ["clo"];
% 
% onset_across = 1;
% 
% comp_measures = ["Perc", "Perc_z"];
% 
% for parc = ["chan", "roi_120"] 
%     half_violin(comp_meth_with_resec, det_meths,comp_measures, "parc", parc,...
%         "onset_across", onset_across, "vis_plot", "on",...
%         "save_fig", 0, "save_loc", save_loc,  "file_type", file_type)
% end
% 
% fig = figure();
% tiledlayout(1,2)
% for parc = ["chan", "roi_120"]
%     for det_meth = det_meths
%         nexttile
%         data_tab = comp_meth_with_resec.across_sz.(sprintf(parc)).(sprintf(det_meth));
%         data = data_tab.Perc;
%         out = data_tab.Outcome;
%        
%         out_str = string();
%         out_str(out<3) = "Good";
%         out_str(out>2) = "Bad";
%         resec_str = string();
%         resec_str(data>=0.5) = "Yes";
%         resec_str(data<0.5) = "No";
%         [contingency_tab, ~, p, lab] = crosstab(resec_str, out_str);
%         out_tab = table(resec_str', out_str');
%         out_tab.Properties.VariableNames = ["resec_50", "Outcome"];
%         if any(any(contingency_tab <5))
%             [~, p, ~] = fishertest(contingency_tab);
%             test_name = "Fishers exact";
%         else 
%             test_name = "Chi-square";
%         end
%         heatmap(out_tab, "resec_50", "Outcome")
%         colorbar off
%         title(sprintf("%s (%s p=%.2f)", parc, test_name, p))
%     end
%      sgtitle(sprintf("Comparison against resection (%s)", onset_across_titles(onset_across+1)))
%      %saveas(fig, sprintf("figures/rapid_prototype/onset_resec/%s_contingency_%s.png", parc, onset_across_titles(onset_across+1)))
%       
% end

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

%% Subquestion 4: Does a larger resection tend to result in more favourable outcomes?
% Is a larger resection (i.e., more regions resected) associated with better
% post-surgical outcomes?
% Remove subjects with no recorded regions resected
final_output = final_output_all(final_output_all.outcome ~= 8,:);
final_output = final_output(cellfun(@sum, final_output.resected_chan)>0,:);

save_loc = "figures/subquestions/q4/";
file_type = "svg";

vols = atlas.vol{2};
names = atlas.name{2};

% We will consider the count of regions and the total volume
resec_count = cellfun(@sum,final_output.resected_roi_120);
resec_vol = nan(size(final_output,1),1);

for pat = 1:size(final_output)
    pat_onset = final_output(pat,:);
    regions = pat_onset.roi_names_120{:};
    resec = pat_onset.resected_roi_120{:};
    resec_regions = regions(logical(resec));
    resec_reg_atlas = contains(names, resec_regions);
    resec_vol(pat) = sum(vols(resec_reg_atlas));
end

fig = figure("Position",[10,10,500,900]);
subplot(3,2,1:2)
boxchart(double(final_output.outcome>2), resec_count, 'MarkerStyle','none')
hold on
swarmchart(double(final_output.outcome>2), resec_count, 'filled')
hold off
set(gca, "XTick", [0,1], "XTickLabel", ["ILAE 1-2", "ILAE 3+"])
title("Comparing resection region count across outcome groups")
subplot(3,2,3:4)
boxchart(double(final_output.outcome>2), log(resec_count), 'MarkerStyle','none')
hold on
swarmchart(double(final_output.outcome>2), log(resec_count), 'filled')
hold off
set(gca, "XTick", [0,1], "XTickLabel", ["ILAE 1-2", "ILAE 3+"])
title("Comparing log(resection size) across outcome groups")
subplot(3,2,5)
histogram(log(resec_count(final_output.outcome<3)))
title("Log(resection size) (ILAE 1-2)")
xlim([0 3])
ylim([0 14])
subplot(3,2,6)
histogram(log(resec_count(final_output.outcome>2)))
title("Log(resection size) (ILAE 3+)")
xlim([0 3])
ylim([0 14])
[~,p, ~,st] = ttest(log(resec_count(final_output.outcome<3)), log(resec_count(final_output.outcome>2)));
sgtitle(sprintf("t(%d)=%.3f, p=%.3f",st.df ,st.tstat,p))
saveas(fig, sprintf("%s%s_comp_resec_size_count.%s", save_loc, parc, file_type))

fig2 = figure("Position",[10,10,500,900]);
subplot(2,2,1:2)
boxchart(double(final_output.outcome>2), resec_vol, 'MarkerStyle','none')
hold on
swarmchart(double(final_output.outcome>2), resec_vol, 'filled')
hold off
set(gca, "XTick", [0,1], "XTickLabel", ["ILAE 1-2", "ILAE 3+"])
title("Comparing resection volume across outcome groups")
subplot(2,2,3)
histogram(resec_vol(final_output.outcome<3), "BinWidth",5*10^4)
title("Resection size (ILAE 1-2)")
ylim([0 14])
xlim([0,3*10^5])
subplot(2,2,4)
histogram(resec_vol(final_output.outcome>2), "BinWidth",5*10^4)
title("Resection size (ILAE 3+)")
ylim([0 14])
xlim([0,3*10^5])
[~,p, ~,st] = ttest(resec_vol(final_output.outcome<3), resec_vol(final_output.outcome>2));
sgtitle(sprintf("t(%d)=%.3f, p=%.3f",st.df ,st.tstat,p))
saveas(fig2, sprintf("%s%s_comp_resec_vol.%s", save_loc, parc, file_type))



%% Subquestion 5: Does a larger onset tend to result in less favourable outcomes?
final_output = final_output_all(final_output_all.outcome ~= 8,:);
final_output = final_output(cellfun(@sum, final_output.clo_chan)>0,:);

save_loc = "figures/subquestions/q5/";
file_type = "svg";

vols = atlas.vol{2};
names = atlas.name{2};

% We will consider the count of regions and the total volume
clo_count = cellfun(@sum,final_output.clo_roi_120);
clo_vol = nan(size(final_output,1),1);

for pat = 1:size(final_output)
    pat_onset = final_output(pat,:);
    regions = pat_onset.roi_names_120{:};
    clo = pat_onset.clo_roi_120{:};
    resec_regions = regions(logical(clo));
    resec_reg_atlas = contains(names, resec_regions);
    clo_vol(pat) = sum(vols(resec_reg_atlas));
end

fig = figure("Position",[10,10,500,900]);
subplot(3,2,1:2)
boxchart(double(final_output.outcome>2), clo_count, 'MarkerStyle','none')
hold on
swarmchart(double(final_output.outcome>2), clo_count, 'filled')
hold off
set(gca, "XTick", [0,1], "XTickLabel", ["ILAE 1-2", "ILAE 3+"])
title("Comparing resection region count across outcome groups")
subplot(3,2,3:4)
boxchart(double(final_output.outcome>2), log(clo_count), 'MarkerStyle','none')
hold on
swarmchart(double(final_output.outcome>2), log(clo_count), 'filled')
hold off
set(gca, "XTick", [0,1], "XTickLabel", ["ILAE 1-2", "ILAE 3+"])
title("Comparing log(CLO size) across outcome groups")
subplot(3,2,5)
histogram(log(clo_count(final_output.outcome<3)))
title("Log(CLO size) (ILAE 1-2)")
xlim([0 3])
ylim([0 14])
subplot(3,2,6)
histogram(log(clo_count(final_output.outcome>2)))
title("Log(CLO size) (ILAE 3+)")
xlim([0 3])
ylim([0 14])
[~,p, ~,st] = ttest(log(clo_count(final_output.outcome<3)), log(clo_count(final_output.outcome>2)));
sgtitle(sprintf("t(%d)=%.3f, p=%.3f",st.df ,st.tstat,p))
saveas(fig, sprintf("%s%s_comp_clo_count.%s", save_loc, parc, file_type))

fig2 = figure("Position",[10,10,500,900]);
subplot(3,2,1:2)
boxchart(double(final_output.outcome>2), clo_vol, 'MarkerStyle','none')
hold on
swarmchart(double(final_output.outcome>2), clo_vol, 'filled')
hold off
set(gca, "XTick", [0,1], "XTickLabel", ["ILAE 1-2", "ILAE 3+"])
title("Comparing resection volume across outcome groups")
subplot(3,2,3:4)
boxchart(double(final_output.outcome>2), log(clo_vol), 'MarkerStyle','none')
hold on
swarmchart(double(final_output.outcome>2), log(clo_vol), 'filled')
hold off
set(gca, "XTick", [0,1], "XTickLabel", ["ILAE 1-2", "ILAE 3+"])
title("Comparing log(CLO volume) across outcome groups")
subplot(3,2,5)
histogram(log(clo_vol(final_output.outcome<3)))
title("Log(CLO volume) (ILAE 1-2)")
ylim([0 14])
subplot(3,2,6)
histogram(log(clo_vol(final_output.outcome>2)))
title("Log(CLO volume) (ILAE 3+)")
ylim([0 14])
[~,p, ~,st] = ttest(log(clo_vol(final_output.outcome<3)), log(clo_vol(final_output.outcome>2)));
sgtitle(sprintf("t(%d)=%.3f, p=%.3f",st.df ,st.tstat,p))
saveas(fig2, sprintf("%s%s_comp_clo_vol.%s", save_loc, parc, file_type))

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
final_output = final_output_all(final_output_all.outcome ~= 8,:);

save_loc = "figures/subquestions/q7/";
file_type = "svg";

subj_lvl_comp_clean = subj_lvl_comp(subj_lvl_comp.sz_count>4,:);


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
       
        boxchart(double(final_output.outcome(incl_subj)>2), summ_vals, "MarkerStyle", "none")
        hold on
        swarmchart(double(final_output.outcome(incl_subj)>2), summ_vals, "filled")
        hold off
        set(gca, "XTick", [0,1], "XTickLabel", ["ILAE 1-2", "ILAE 3+"])
        title(sprintf("%s %s", comp, strrep(column, "_", " ")))
    end
end
sgtitle("Seizure onset variability")
saveas(fig, sprintf("%s_onset_var.%s", save_loc, file_type))


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
       
        boxchart(double(final_output.outcome(incl_subj)>2), summ_vals, "MarkerStyle", "none")
        hold on
        swarmchart(double(final_output.outcome(incl_subj)>2), summ_vals, "filled")
        hold off
        set(gca, "XTick", [0,1], "XTickLabel", ["ILAE 1-2", "ILAE 3+"])
        title(sprintf("%s %s", comp, strrep(column, "_", " ")))
    end
end
sgtitle("Seizure onset variability (values from proportions)")
saveas(fig, sprintf("%s_onset_var(prop).%s", save_loc, file_type))


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
    title(sprintf("%s ILAE 1-2", strrep(col_name, "_", " ")))
    nexttile
    histogram(concord_mat_tab(concord_mat_tab.Outcome>2,:).(sprintf(col_name)), "BinWidth", 0.2)
    title(sprintf("%s ILAE 3+", strrep(col_name, "_", " ")))
end


