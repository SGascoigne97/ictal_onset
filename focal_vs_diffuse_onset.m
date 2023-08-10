final_output_all = load('tables/final_output_all_sz.mat');
final_output_all = final_output_all.final_output;
load('roi_info/ATLAS.mat')

addpath(genpath('sarah_functions'))
addpath(genpath('Nathan Code'))

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
incl_thresh = 0.75;

plot_grph = 1;
for subj = 65%1:size(final_output,1)
    subj_onset = final_output(subj,:);
    imprint = subj_onset.imprint_roi_120{:}; 
    regions = subj_onset.roi_names_120{:};
    
    % Isolate only the recorded hemisphere to make computations simpler 
    % and more efficient
    if any(contains(regions,"l.")) & any(contains(regions,"r."))
        % Bilateral placement
        incl_reg = repelem(1,length(names),1);
        y_lim = [0,440];
        xyz_view = [0 90 0];
    elseif any(contains(regions,"l.")) & ~any(contains(regions,"r."))
        % Left hemisphere placement
        incl_reg = contains(names,"Left"|"l.");
        y_lim = [230,440];
        xyz_view = [-90 0 45];
    elseif ~any(contains(regions,"l.")) & any(contains(regions,"r."))
        % Right hemisphere placement
        incl_reg = contains(names,"Right"|"r.");  
        y_lim = [0,210];
        xyz_view = [90 0 25];
    end
    
    % Keep only the names of the regions in the relavant hemisphere
    names_hem = names(find(incl_reg));
    xyz_hem = xyz(find(incl_reg),:);
    dists_hem = dists(find(incl_reg),find(incl_reg));
    vols_hem = vols(find(incl_reg));
    
    % Create subsection of touching matrix for relevant regions
    % (n regions in hemisphere x n regions in hemisphere)
    region_index = find_preserved(names_hem,regions);
    imprint = logical(imprint);
    
    % For graph theory 'learning' (ask Iain for correct nomenclature), we will
    % exclude recorded regions that we know are non-onset
    clear onset_reg_count_sz
    clear max_dist
    clear onset_vol
    clear subj_tab

    max_dist_orig = nan(size(imprint,2), 1);
    max_dist_upd = max_dist_orig;
    onset_vol_orig = max_dist_orig;
    onset_vol_upd = max_dist_orig;

    onset_reg_count_orig = max_dist_orig;
    onset_reg_count_upd = max_dist_orig;

    % Seizure specific process
    for sz = 1:size(imprint,2)% Assess one seizure at a time 
        clear g
        if plot_grph == 1
            f1 = figure();
            f1.Position = [2000,1000,500,500];
            tiledlayout(2,1)
            sgtitle(sprintf("%s sz %d", string(subj_onset.Patient_id),sz))
            f2 = figure();
            f2.Position = [2500,1000,500,500];
            tiledlayout(2,1)
            sgtitle(sprintf("%s sz %d", string(subj_onset.Patient_id),sz))
        end
    
        % Extend imprint to include all regions in relevant hemisphere
        imprint_full = zeros(length(names_hem),1);
        imprint_full(region_index(imprint(:,sz))) = 1;
        % We give recorded regions that are known to be non-onset a value of -1 so
        % they can be removed prior to creating graphs
        rec_non_onset =  zeros(length(names_hem),1);
        rec_non_onset(region_index(imprint(:,sz)==0)) = 1;
    
        % We will give each known onset regions value of 1 after each update
        % We will give recorded non-onset regions value of 0 after each update
        
        % Create subsection of touching matrix to only include relevant hemisphere
        % and exclude recorded non-onset regions
        connec_hem = connec(find(incl_reg),find(incl_reg));
    
        g = graph(connec_hem);

        for edge = 1:length(g.Edges.Weight)
            g.Edges.Weight(edge) = mean([imprint_full(g.Edges.EndNodes(edge,1)), imprint_full(g.Edges.EndNodes(edge,2))]);
        end
    
        imprint_full_plt = imprint_full;
        imprint_full_plt(find(rec_non_onset)) = -1;

        if plot_grph == 1
            set(0, 'CurrentFigure', f1)
            nexttile
            h = plot(g,'LineWidth',(g.Edges.Weight+0.01)*3,'EdgeLabel',{},...
                'NodeLabel', {},'XData',xyz_hem(:,1),'YData',xyz_hem(:,2),...
                'ZData',xyz_hem(:,3));
            colormap(cm)
            h.EdgeCData = g.Edges.Weight;
            h.NodeCData = imprint_full_plt;
            title("Original")
            axis off
            clim([-1 1])
            view(xyz_view);
            xlabel("x")
            ylabel("y")
          
            set(0, 'CurrentFigure', f2)
            if min(imprint_full_plt)==-1 & max(imprint_full_plt)==1
                cm1 = [0,0,0;1,1,1;0,0.7,0.7];
                cm1 = interp1(cm1, 1:0.01:size(cm1,1));
            elseif min(imprint_full_plt)==0 & max(imprint_full_plt)==1
                cm1 = [1,1,1;0,0.7,0.7];
                cm1 = interp1(cm1, 1:0.01:size(cm1,1));
            end
            nexttile;
            plotBrain_NE(names_hem, imprint_full_plt ,'cm',cm1);
            colorbar off
            title("Original")
            ylim(y_lim)
        end

        upd_ons = nan(height(g.Nodes),1);
        for node = 1:height(g.Nodes)
            if imprint_full_plt(node) == 0
                node_vals = g.Edges.Weight(g.Edges.EndNodes(:,1) == node | g.Edges.EndNodes(:,2) == node);
                upd_ons(node) = sum(node_vals >= 0.5) >= length(node_vals)*incl_thresh;
            end
        end
            
        upd_ons(logical(imprint_full))= 1;
        upd_ons(logical(rec_non_onset)) = -1;

        for edge = 1:length(g.Edges.Weight)
            g.Edges.Weight(edge) = all([upd_ons(g.Edges.EndNodes(edge,1)), upd_ons(g.Edges.EndNodes(edge,2))]==1);
        end
        
        if plot_grph == 1
            set(0, 'CurrentFigure', f1)
            nexttile
            h = plot(g,'LineWidth',(g.Edges.Weight+0.01)*3,'EdgeLabel',{},...
                'NodeLabel', {},'XData',xyz_hem(:,1),'YData',xyz_hem(:,2),...
                'ZData',xyz_hem(:,3));
            
            colormap(cm)
            h.EdgeCData = g.Edges.Weight;
            h.NodeCData = upd_ons;
            title("Updated")
            axis off
            clim([-1 1])
            view(xyz_view);

            upd_ons_plot = upd_ons;
            upd_ons_plot(upd_ons == 1 & imprint_full ==0) = 0.5;

            if min(upd_ons)==-1 & max(upd_ons)==1
                cm2 = [0,0,0;0,0,0;1,1,1;0,0.9,0.9;0,0.7,0.7];
                cm2 = interp1(cm2, 1:0.01:size(cm2,1));
            elseif min(upd_ons)==0 & max(upd_ons)==1
                cm2 = [1,1,1;0,0.9,0.9;0,0.7,0.7];
                cm2 = interp1(cm2, 1:0.01:size(cm2,1));
            end
         
            set(0, 'CurrentFigure', f2)
            nexttile;
            plotBrain_NE(names_hem, upd_ons_plot ,'cm', cm2);
            colorbar off
            title("Updated")
            ylim(y_lim)
        end
    
        onset_reg_count_orig(sz) = sum(imprint_full==1);
        onset_reg_count_upd(sz) = sum(upd_ons==1);

        max_dist_orig(sz) = max(max(dists_hem(imprint_full==1,imprint_full==1)));
        max_dist_upd(sz) = max(max(dists_hem(upd_ons==1,upd_ons==1)));

        onset_vol_orig(sz) = sum(vols_hem(imprint_full==1));
        onset_vol_upd(sz) = sum(vols_hem(upd_ons==1));
 
    end

    subj_tab = table(repmat(subj_onset.Patient_id, size(imprint,2),1),...
        (1:size(imprint,2))', onset_reg_count_orig, onset_reg_count_upd, ...
        max_dist_orig, max_dist_upd, onset_vol_orig, onset_vol_upd,...
        'VariableNames', ["Subj_id","sz","onset_reg_count_orig",...
        "onset_reg_count_upd", "max_dist_orig", "max_dist_upd",...
        "onset_vol_orig","onset_vol_upd"]);

    if exist('comp_tab', 'var')
        comp_tab = [comp_tab; subj_tab];
    else 
        comp_tab = subj_tab;
    end

end

%%

diff = comp_tab.onset_reg_count_upd-comp_tab.onset_reg_count_orig;
diff_vol = comp_tab.onset_vol_upd-comp_tab.onset_vol_orig;
diff_dist = comp_tab.max_dist_upd-comp_tab.max_dist_orig;

figure()
subplot(311)
histogram(diff(diff~=0))
title("Change in count of regions")
subplot(312)
histogram(diff_dist(diff~=0), "BinWidth",1)
title("Change in maximum distance between regions")
subplot(313)
histogram(diff_vol(diff~=0), "BinWidth", 10000)
title("Change in volume")

%%

diff_thresh = 50;
diff_thresh_vol = 50000;
% Setting threshold at 50mm for now

for subj = 1:size(final_output,1)
    subj_dist = comp_tab(comp_tab.Subj_id == string(final_output(subj,:).Patient_id),:);
    for typ = ["orig", "upd"]
        for mark = ["max_dist", "onset_vol"]
            if mark == "max_dist"
                thresh = 50;
            else
                thresh = 50000;
            end
            % Store raw distance and volume measures
            subj_lvl_comp.(sprintf(typ)).(sprintf("%s", mark))(subj,1) = {subj_dist.(sprintf("%s_%s", mark, typ))};

            subj_lvl_comp.(sprintf(typ)).(sprintf("present_diffuse_%s", mark))(subj,1) = sum(subj_dist.(sprintf("%s_%s", mark, typ))> thresh)>0;
            subj_lvl_comp.(sprintf(typ)).(sprintf("sz_diffuse_%s", mark))(subj,1) = {subj_dist.(sprintf("%s_%s", mark, typ))>thresh};
        end
    end
end
% %%
% for typ = ["orig", "upd"]
%         for mark = ["max_dist", "onset_vol"]
% 
%             [tab,chi,p] = crosstab(final_output.outcome>2, subj_lvl_comp.(sprintf(typ)).(sprintf("present_diffuse_%s", mark)));
%             figure()
%             heatmap(tab,...
%                 "YData",["ILAE 1-2", "ILAE 3+"], "XData", ["No diffuse onsets", ">=1 diffuse onset"])
%             title(sprintf("Count of subjects with at >=1 diffuse onset (%s, %s) \n chi^2 = %.2f, p = %.3f", typ, mark,chi,p))
% 
%         end
% end

%% Need to determine the best way to incorporate variation in onset distance and volume into comparisons
med_dist_o = nan(length(subj_lvl_comp.upd.max_dist),1);
med_vol_o = med_dist;
med_dist_u = med_dist;
med_vol_u = med_dist;

for subj = 1:length(med_dist)
    if length(subj_lvl_comp.orig.max_dist{subj}) >2
        med_dist_o(subj) = median(subj_lvl_comp.orig.max_dist{subj});
        med_vol_o(subj) = median(subj_lvl_comp.orig.onset_vol{subj});
    
        med_dist_u(subj) = median(subj_lvl_comp.upd.max_dist{subj});
        med_vol_u(subj) = median(subj_lvl_comp.upd.onset_vol{subj});
    else
        med_dist_o(subj) = nan;
        med_vol_o(subj) = nan;
    
        med_dist_u(subj) = nan;
        med_vol_u(subj) = nan;
    end
end

summary_tab = table(string(final_output.Patient_id), final_output.outcome, med_dist_o, ...
    med_dist_u, med_vol_o, med_vol_u, 'VariableNames', {'Patient_id', ...
    'outcome', 'med_dist_orig', 'med_dist_upd', 'med_vol_orig', 'med_vol_upd'});

%%
summary_tab = summary_tab(~isnan(summary_tab.med_dist_orig),:);
comp_var = ["med_dist", "med_vol"];
onset_type = ["orig", "upd"];
subjs = size(summary_tab.(sprintf("med_dist_%s", ons_typ)),1);

offset = (subjs+10)/2;
max_x = (offset*((length(onset_type))-1)) +10;
x_lim = [-10, max_x];

f = figure(1);
f.Position = [100,1000,900,450];
tiledlayout(1,2)


for var = comp_var
    nexttile
    hold on
    for ons_ind = 1:length(onset_type)
        ons_typ = onset_type(ons_ind);
        data = summary_tab.(sprintf("%s_%s", var, ons_typ));
        out = summary_tab.outcome;
        out = out(~isnan(data));
        data = data(~isnan(data));

        bin_wd = 5*(max(data)-min(data))/length(data);
    
%         if var == "med_dist"
%             bin_wd = 8;
%         else
%             bin_wd = 10000;
%         end
                   
        g = data(out<3);
        b = data(out>2);
                  
        min_val_g = round(min(g)*(1/bin_wd))*bin_wd;
        min_val_b = round(min(b)*(1/bin_wd))*bin_wd;
        max_val_g = round(max(g)*(1/bin_wd))*bin_wd;
        max_val_b = round(max(b)*(1/bin_wd))*bin_wd;
        
        if min_val_g > min(g)
            min_val_g = min_val_g - bin_wd;
        end 
        if min_val_b > min(b)
            min_val_b = min_val_b - bin_wd;
        end
        if max_val_g < max(g)
            max_val_g = max_val_g + bin_wd;
        end
        if max_val_b < max(b)
            max_val_b = max_val_b + bin_wd;
        end
        
        if min_val_g == max_val_g
            max_val_g = max_val_g+bin_wd;
        end
        if min_val_b == max_val_b
            max_val_b = max_val_b+bin_wd;
        end
        
        hist_vals_g = histcounts(g, min_val_g:bin_wd:max_val_g);
        % Add scatter points
        non_zero_grps = hist_vals_g(hist_vals_g~=0);
        jitt = [];
        for grp_sz = non_zero_grps
            jitt = [jitt; rand(grp_sz,1)*grp_sz*0.8];
        end
        scatter((ons_ind-1)*offset+jitt, sort(g), [], [0.3010 0.7450 0.9330], 'filled');
        
        hist_vals_b = histcounts(b, min_val_b:bin_wd:max_val_b);
        % Add scatter points
        non_zero_grps = hist_vals_b(hist_vals_b~=0);
        jitt = [];
        for grp_sz = non_zero_grps
            jitt = [jitt; rand(grp_sz,1)*grp_sz*0.8];
        end
        scatter((ons_ind-1)*offset-jitt, sort(b), [],[0.8500 0.3250 0.0980], 'filled');
        if length(hist_vals_g) >1
                smooth_data_g = 1.1*(smoothdata(interp1(hist_vals_g, 1:0.01:length(hist_vals_g))));
            fill((ons_ind-1)*offset+[0 smooth_data_g 0], ...
                rescale((1:length([0 smooth_data_g 0]))/length([0 smooth_data_g 0]),...
                min(g) - (bin_wd/10), max(g) + (bin_wd/10)), [0.3010 0.7450 0.9330], 'FaceAlpha',0.3,...
                'LineStyle','none')
        end 
        
        if length(hist_vals_b) > 1
            smooth_data_b = 1.1*(smoothdata(interp1(hist_vals_b, 1:0.01:length(hist_vals_b))));
            fill((ons_ind-1)*offset -[0 smooth_data_b 0],...
                rescale((1:length([0 smooth_data_b 0]))/length([0 smooth_data_b 0]),...
                min(b) - (bin_wd/10), max(b) + (bin_wd/10)), [0.8500 0.3250 0.0980], 'FaceAlpha',0.3,...
                'LineStyle','none')

        end
    end

    if var == "med_dist"
        yline(25, '--')
    else 
        yline(25^3, '--')
    end

    set(gca, "XTick", (0:(length(onset_type))-1)*offset, "XTickLabel", onset_type)
    title(strrep(var, "_", " "))
end
