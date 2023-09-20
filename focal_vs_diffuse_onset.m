final_output_all = load('tables/final_output_all_sz.mat');
final_output_all = final_output_all.final_output;
load('roi_info/ATLAS.mat')

addpath(genpath('sarah_functions'))
addpath(genpath('Nathan Code'))

save_loc = "figures/subquestions/q6/";
file_type = "svg";

% Remove patients with no labelled CLO or no outcome 
final_output = final_output_all(final_output_all.outcome ~= 8,:);

% Extract region information for all regions in Lausanne 120 atlas
names = atlas.name{2};
xyz = atlas.xyz{2};
dists = atlas.dists{2};
vols = atlas.vol{2};

% Load touching matrices
load('Nathan Code/SimpleRouteExample/touching_struct_threshold_100.mat', 'touching')
gm_touching = touching.touching;
wm_sc_sc = touching.wm_subcortical_to_subcortical;
wm_sc_c = touching.wm_subcortical_to_cortical;

% Create matrix with grey matter cortical-cortical connections and white
% matter subcortical-subcortical connections
connec = gm_touching+wm_sc_sc+wm_sc_c;
connec = connec - diag(diag(connec));

% Set colour maps
cm = [0,0,0;.7,.7,.7;1,0,0];
cm = interp1(cm, 1:0.01:size(cm,1)); % low black, mid grey, high red

cm1 = [0.7,0.7,0.7;0,0.9,0.9;0,0.7,0.7];
cm1 = interp1(cm1, 1:0.01:size(cm1,1));

%% Choose subject to assess

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

%%
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

%% Need to determine the best way to incorporate variation in onset distance and volume into comparisons
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
%%
half_violin_upd(summary_tab,  ["med_max_dist","max_max_dist",...
    "subj_med_prop_dist","subj_max_prop_dist"],...
    double(summary_tab.outcome>2), "save_fig",1,...
    "file_type","svg", "save_loc",sprintf("%svol_comp", save_loc), "grp_names", ["ILAE 1-2", "ILAE 3+"])

half_violin_upd(summary_tab,  ["med_onset_vol","max_onset_vol",...
    "subj_med_prop_vol","subj_max_prop_vol"],...
    double(summary_tab.outcome>2), "save_fig",1,...
    "file_type","svg", "save_loc",sprintf("%svol_comp", save_loc), "grp_names", ["ILAE 1-2", "ILAE 3+"])