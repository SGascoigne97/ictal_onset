%% Sub question 1: Does seizure onset tend to be resected?
det_meths = ["clo", "imprint", "EI"];
comp_measures = ["Perc", "Perc_z"];

for parc = ["chan", "roi_120"]
    half_violin(comp_meth_with_resec, det_meths,comp_measures, "parc", parc,...
        "onset_across", 1, "vis_plot", "on",...
        "save_fig", 0, "save_loc", save_loc,  "file_type", file_type)
end

var_names = {'Parc','Det','All','Total','Fav','Fav_total', 'Unfav', 'Unfav_total'};
onset_resec_tab = array2table(zeros(0,8),'VariableNames',var_names);


for parc = ["chan", "roi_120"]
    for det = det_meths
        data_tab = comp_meth_with_resec.across_sz.(sprintf("%s", parc)).(sprintf("%s", det));
        data = data_tab.Perc;
        out = data_tab.Outcome;
        
        acr = sum(data>0.5);
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

%% Sub question 2: Is resecting a larger proportion of the clinically labelled onset associated with more favourable outcomes?
final_output = final_output_all;

onset_across = 1;
onset_across_titles = ["across_sz", "per_sz"];
onset_across_title = onset_across_titles(onset_across);
n_perm = 1000;

det_meths = ["clo", "imprint", "EI"];

onset_acr_thresh = 0.5;
% Create an empty table to store output
fprintf("%s\n", onset_across_title(onset_across))
atl = [72,120,250,500];

for atl_ind = [2,3]
    chan_or_roi = sprintf("roi_%d", atl(atl_ind));
    vols = atlas.vol{atl_ind};
    
    for det_meth = det_meths
        comp_meth_tab = compute_resec_comp(final_output,atlas(atl_ind,:),vols,...
            "chan_or_roi", chan_or_roi, "det_meth", det_meth, "n_perm",n_perm,...
            "onset_acr_thresh",onset_acr_thresh,"onset_across",onset_across);
        comp_meth_tab.outcome_cat = categorical(comp_meth_tab.Outcome>2,[0,1], ["ILAE 1-2", "ILAE 3+"]);
        comp_meth_with_resec.(sprintf("%s", onset_across_title(onset_across))).(sprintf("%s", chan_or_roi)).(sprintf("%s", det_meth)) = comp_meth_tab;
    end
end

% Create half-violin plots
comp_measures = ["Perc", "Perc_vol"];

for parc = ["roi_120"]%, "roi_250"] 
    half_violin(comp_meth_with_resec, det_meths, comp_measures, "parc", parc,...
        "onset_across", onset_across, "vis_plot", "on",...
        "save_fig", 0)
end


for comp_measure = comp_measures
    fig = figure("Position",[100,1000, 1200, 400]);
    tiledlayout(1,length(det_meths))
    for det_meth = det_meths
        nexttile
        data_tab = comp_meth_with_resec.consensus.roi_120.(sprintf(det_meth));
        data = data_tab.(sprintf(comp_measure));
        out = data_tab.outcome_cat;

        resec_str = string();
        resec_str(data>=0.5) = "Yes";
        resec_str(data<0.5) = "No";
        [contingency_tab, ~, p, lab] = crosstab(resec_str, out);
        out_tab = table(resec_str', out);
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
    sgtitle(strrep(comp_measure, "_", " "))
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

for parc = ["chan", "roi_120"] 
    half_violin(comp_meth_with_resec, det_meths,comp_measures, "parc", parc,...
        "onset_across", onset_across, "vis_plot", "on",...
        "save_fig", 0, "save_loc", save_loc,  "file_type", file_type)
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
     %saveas(fig, sprintf("figures/rapid_prototype/onset_resec/%s_contingency_%s.png", parc, onset_across_titles(onset_across+1)))
      
end

%% Subquestion 4: Does a larger resection tend to result in more favourable outcomes?
% Is a larger resection (i.e., more regions resected) associated with better
% post-surgical outcomes?
final_output_all = load('tables/final_output_all_sz.mat');
final_output_all = final_output_all.final_output;

% Add surgical outcome to table
for pat = 1:size(final_output_all, 1)
    pat_onset = final_output_all(pat,:);
    outcome_id = pat_onset.("Outcome year"){:}-pat_onset.("Surgery year") == 1;
    final_output_all.outcome(pat) = pat_onset.Surgery_outcome{:}(outcome_id);
end

% Remove subjects with no surgical outcome listed
final_output_all = final_output_all(final_output_all.outcome ~= 8,:);
% Remove subjects with no recorded regions resected
final_output = final_output_all(cellfun(@sum, final_output_all.resected_chan)>0,:);

regions_resec = cellfun(@sum,final_output.resected_roi_120);

figure("Position",[10,10,500,900])
subplot(3,2,1:2)
boxchart(double(final_output.outcome>2), regions_resec, 'MarkerStyle','none')
hold on
swarmchart(double(final_output.outcome>2), regions_resec, 'filled')
hold off
set(gca, "XTick", [0,1], "XTickLabel", ["ILAE 1-2", "ILAE 3+"])
title("Comparing resection size across outcome groups")
subplot(3,2,3:4)
boxchart(double(final_output.outcome>2), log(regions_resec), 'MarkerStyle','none')
hold on
swarmchart(double(final_output.outcome>2), log(regions_resec), 'filled')
hold off
set(gca, "XTick", [0,1], "XTickLabel", ["ILAE 1-2", "ILAE 3+"])
title("Comparing log(resection size) across outcome groups")
subplot(3,2,5)
histogram(log(regions_resec(final_output.outcome<3)))
title("Log(resection size) (ILAE 1-2)")
xlim([0 3])
ylim([0 14])
subplot(3,2,6)
histogram(log(regions_resec(final_output.outcome>2)))
title("Log(resection size) (ILAE 3+)")
xlim([0 3])
ylim([0 14])
[h,p, ~,st] = ttest(log(regions_resec(final_output.outcome<3)), log(regions_resec(final_output.outcome>2)))

%% Subquestion 5: Does a larger onset tend to result in less favourable outcomes?
final_output_all = load('tables/final_output_all_sz.mat');
final_output_all = final_output_all.final_output;

% Add surgical outcome to table
for pat = 1:size(final_output_all, 1)
    pat_onset = final_output_all(pat,:);
    outcome_id = pat_onset.("Outcome year"){:}-pat_onset.("Surgery year") == 1;
    final_output_all.outcome(pat) = pat_onset.Surgery_outcome{:}(outcome_id);
end

% Remove subjects with no surgical outcome listed
final_output_all = final_output_all(final_output_all.outcome ~= 8,:);
% Remove subjects with no recorded regions resected
final_output = final_output_all(cellfun(@sum, final_output_all.clo_chan)>0,:);

regions_clo = cellfun(@sum,final_output.clo_roi_120);

figure("Position",[10,10,500,900])
subplot(3,2,1:2)
boxchart(double(final_output.outcome>2), regions_clo, 'MarkerStyle','none')
hold on
swarmchart(double(final_output.outcome>2), regions_clo, 'filled')
hold off
set(gca, "XTick", [0,1], "XTickLabel", ["ILAE 1-2", "ILAE 3+"])
title("Comparing CLO size across outcome groups")
subplot(3,2,3:4)
boxchart(double(final_output.outcome>2), log(regions_clo), 'MarkerStyle','none')
hold on
swarmchart(double(final_output.outcome>2), log(regions_clo), 'filled')
hold off
set(gca, "XTick", [0,1], "XTickLabel", ["ILAE 1-2", "ILAE 3+"])
title("Comparing log(CLO size) across outcome groups")
subplot(3,2,5)
histogram(log(regions_clo(final_output.outcome<3)))
title("Log(clo size) (ILAE 1-2)")
xlim([0 3])
ylim([0 14])
subplot(3,2,6)
histogram(log(regions_clo(final_output.outcome>2)))
title("Log(clo size) (ILAE 3+)")
xlim([0 3])
ylim([0 14])
[h,p, ~,st] = ttest(log(regions_clo(final_output.outcome<3)), log(regions_clo(final_output.outcome>2)))

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

%% Subquestion 6: Is contiguity of seizure onset associated with surgical outcomes?
figure()
subplot(2,2,1)
imagesc(touching.touching)
title("Touching")
subplot(2,2,2)
imagesc(touching.wm_cortical_to_cortical)
title("WM cortical to cortical")
subplot(2,2,3)
imagesc(touching.wm_subcortical_to_cortical)
title("WM subcortical to cortical")
subplot(2,2,4)
imagesc(touching.wm_subcortical_to_subcortical)
title("WM subcortical to subcortical")


%% Working with Nathan's code
% DO NOT USE THIS NOW - DOESN'T WORK!!!!
addpath(genpath('Brewermap'))
addpath(genpath('Nathan Code'))
% load atlas
load('/home/campus.ncl.ac.uk/b5007876/Desktop/Ictal_onset/roi_info/ATLAS.mat', 'atlas');
load('/home/campus.ncl.ac.uk/b5007876/Desktop/Ictal_onset/roi_info/lhSurf.mat', 'lf', 'lv');
load('/home/campus.ncl.ac.uk/b5007876/Desktop/Ictal_onset/roi_info/rhSurf.mat', 'rf', 'rv');

% load connectivity matricies
load('sFrances.mat', 's');
s_table = struct2table(s);

% load touching matrix 36 60 125 250
% scale 125 touching matrix left hand side not symetric???????
scale = 60;
[names, touching, roiIndex] = load_touching_matrix(scale);

groups = nan(1,size(final_output,1));
connection_types = groups;

% DOESN'T WORK YET - IT IS SEPARATING HIPPOCAMPUS AND AMYGDALA :(

for pat = 39%1:size(final_output,1)
    patient = string(final_output(pat,:).Patient_id);
   
    conData = s_table.scale60(strcmp(s_table.subjID,extractAfter(patient, 'UCLH')));
    if isempty(conData) | all(cellfun(@isempty,conData))
        fprintf("%s not included in s_table \n", patient)
        continue
    end
    
    pat_onset = final_output_all(string(final_output_all.Patient_id) == patient,:);
    
    % create white matter touching matrix
    wm_names = conData{1}.names;
    wm_names = regexp(wm_names,"[^']*",'match','once');
    wm_touching = logical(conData{1}.adjCount);
    
    % match touching matrix to wm touching matrix
    idx = find_preserved(names,wm_names);
    idx2 = find_preserved(wm_names,names);
    touching_wm_size = zeros(size(wm_touching));
    touching_wm_size(idx2,idx2) = touching(idx,idx);
    
    region_names = pat_onset.(sprintf("roi_names_%d", scale*2)){:}';
    region_indexes = find_preserved(wm_names,region_names);
    
    % Create consensus imprint onset
    imprint_onsets = pat_onset.(sprintf("imprint_roi_%d", scale*2)){:};
    imprint_consensus = mean(imprint_onsets,2)>=0.5;
    
    seizureRoute = SeizureRoute(size(wm_names,1),1, region_indexes, touching_wm_size, wm_touching);
    seizureRoute.newRegions(region_indexes(imprint_consensus),1)

    seizureRoutes{pat} = seizureRoute;

    % See plotting limits to show relevant hemispheres
    if any(contains(region_names,"l.")) & any(contains(region_names,"r."))
        % Bilateral placement
        y_lim = [0 450];
        x_lim = [-500 1300];

    elseif any(contains(region_names,"l.")) & ~any(contains(region_names,"r."))
        % Left hemisphere placement
        y_lim = [230 450]; % If only left hemisphere
        x_lim = [-250 900];

    elseif ~any(contains(region_names,"l.")) & any(contains(region_names,"r."))
        % Right hemisphere placement
        y_lim = [0 210];
        x_lim = [-250 900];
    end

    groups(pat) = seizureRoutes{pat}.n_touching_groups;
    connection_types(pat) = max(seizureRoutes{pat}.region_connections);

    fig = figure();
    seizureRoutes{pat}.simpleBrainPlot(atlas, roiIndex, 1, 1);
    ylim(y_lim)
    title(sprintf("%s (ILAE %d) %d touching group(s)", patient, pat_onset.outcome, seizureRoute.n_touching_groups))
    pause(2)
    %saveas(fig, sprintf("figures/rapid_prototype/contiguity/%s.png", patient))

end

figure()
subplot(121)
histogram(groups(final_output.outcome <3))
hold on
histogram(groups(final_output.outcome >2))
hold off
ylim([0,12])
title("Group count")
subplot(122)
histogram(connection_types(final_output.outcome <3))
hold on
histogram(connection_types(final_output.outcome >2))
xlim([0,5])
ylim([0,12])
title("Connection type (max)")
legend({"ILAE 1-2", "ILAE 3+"})
%% Subquestion 7: Is the presence/prevalence of diffuse onsets associated with with post-surgical outcomes?
figure()
tiledlayout(7,7)
for pat = 1:size(final_output,1)
nexttile
imagesc(final_output(pat,:).imprint_roi_120{:})
end

%%
pat = 5;
pat_onset = final_output(pat,:);
imprint = pat_onset.imprint_roi_120{:}; 
figure()
for sz = 1:size(imprint,2)
    subplot(size(imprint,2),1,sz)
    plotBrain_NE(pat_onset.roi_names_120{:}, imprint(:,sz), 'cm',[0.7,0.7,0.7;0,0.7,0.7]);
    colorbar off
    ylim([0,210])
end
load('Nathan Code/SimpleRouteExample/touching_struct.mat', 'touching')
gm_touching = touching.touching;

%% Subquestion 8: Is a wider epileptogenic network associated with less favourable outcomes? (regions in any/few onsets)



