final_output_all = load('tables/final_output_all_sz.mat');
final_output_all = final_output_all.final_output;
load('roi_info/ATLAS.mat')

addpath(genpath('sarah_functions'))
addpath(genpath('Nathan Code'))

% Remove patients with no labelled CLO or no outcome 
final_output = final_output_all(final_output_all.outcome ~= 8,:);


xyz = atlas.xyz{2};
dists = atlas.dists{2};
vols = atlas.vol{2};
incl_thresh = 1/3;

% Load touching matrices
load('Nathan Code/SimpleRouteExample/touching_struct.mat', 'touching')
gm_touching = touching.touching;
subcort_touching = touching.wm_subcortical_to_subcortical;

% Create matrix with grey matter cortical-cortical connections and white
% matter subcortical-subcortical connections
touch_gm_subcort = gm_touching+subcort_touching;
touch_gm_subcort = touch_gm_subcort - diag(diag(touch_gm_subcort));

% Extract region names for all regions in Lausanne 120 atlas
names = atlas.name{2};

% Set colour maps
cm = [0,0,0;.7,.7,.7;1,0,0];
cm = interp1(cm, 1:0.01:size(cm,1)); % low black, mid grey, high red

cm1 = [0.7,0.7,0.7;0,0.9,0.9;0,0.7,0.7];
cm1 = interp1(cm1, 1:0.01:size(cm1,1));
    
cm2 = [.1,.1,.1;.1,.1,.1;1,1,1;0,0.9,0.9;0,0.7,0.7];
cm2 = interp1(cm2, 1:0.01:size(cm2,1));

%% Choose subject to assess

plot_grph = 0;
for subj = 2:size(final_output,1)
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
    
    % Seizure specific process
    for sz = 1:size(imprint,2)% Assess one seizure at a time 
        clear g

        if plot_grph == 1
            f1 = figure();
            f1.Position = [2000,1000,500,1800];
            tiledlayout(5,1)
            sgtitle(sprintf("%s sz %d", string(subj_onset.Patient_id),sz))
            f2 = figure();
            f2.Position = [2500,1000,500,1000];
            tiledlayout(5,1)
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
        touch_gm_subcort_hem = touch_gm_subcort(find(incl_reg),find(incl_reg));
    
        g = graph(touch_gm_subcort_hem);
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
            nexttile;
            plotBrain_NE(names_hem, imprint_full_plt ,'cm',cm2);
            colorbar off
            title("Original")
            ylim(y_lim)
        end
%     
        onset_reg_count_sz(sz,1) = sum(imprint_full);
        for step = 1:4
            avg_weight = nan(1,height(g.Nodes));
            for node = 1:height(g.Nodes)
                node_vals = g.Edges.Weight(g.Edges.EndNodes(:,1) == node | g.Edges.EndNodes(:,2) == node);
                if length(node_vals)>1 % if the region is connected to at least two edges, we compute the mean
                    avg_weight(node) = mean(node_vals);
                else % if the region is connected to at fewer than two edges, connected onset node (only connection) will mean that this regions is automatically extended into in the second update step
                    avg_weight(node) = node_val(node); 
                end
            end
            
            avg_weight(logical(imprint_full))= 1;
            avg_weight(logical(rec_non_onset)) = 0;
            
            for edge = 1:length(g.Edges.Weight)
                g.Edges.Weight(edge) = mean([avg_weight(g.Edges.EndNodes(edge,1)), avg_weight(g.Edges.EndNodes(edge,2))]);
            end

            upd_node_val_plt = double(avg_weight>0.5);
            upd_node_val_plt(logical(rec_non_onset)) = -1;

            upd_node_brain_plt = upd_node_val_plt;
            upd_node_brain_plt(logical(avg_weight>0.5 & avg_weight<1)) = 0.5;

            if plot_grph == 1
            
                set(0, 'CurrentFigure', f1)
                nexttile
                h = plot(g,'LineWidth',(g.Edges.Weight+0.01)*3,'EdgeLabel',{},...
                    'NodeLabel', {},'XData',xyz_hem(:,1),'YData',xyz_hem(:,2),...
                    'ZData',xyz_hem(:,3));
                colormap(cm)
                h.EdgeCData = g.Edges.Weight;
                
                h.NodeCData = upd_node_val_plt;
                title(sprintf("Update step %d",step))
                axis off
                clim([-1 1])
                view(xyz_view);
             
                set(0, 'CurrentFigure', f2)
                nexttile;
                plotBrain_NE(names_hem, upd_node_brain_plt ,'cm', cm2);
                colorbar off
                title(sprintf("Update step %d", step))
                ylim(y_lim)
            end
    
            onset_reg_count(step) = sum(upd_node_val_plt==1);
            onset_reg_count_sz(sz,step+1) = onset_reg_count(step);
            onset_reg_vol_sz(sz,step+1) = sum(vols_hem(upd_node_val_plt==1));
            
        end

        upd_onset = upd_node_val_plt;
        upd_onset(upd_onset == -1) = 0;
    
        max_dist(sz) = max(max(dists_hem(logical(upd_onset),logical(upd_onset))));
        onset_vol(sz) = sum(vols_hem(logical(upd_onset)));
 
    end

    subj_tab = table(repmat(subj_onset.Patient_id, size(imprint,2),1),...
        (1:size(imprint,2))', max_dist', onset_vol', 'VariableNames', ["Subj_id","sz","max_dist","onset_vol"]);

    if exist('dist_tab', 'var')
        dist_tab = [dist_tab; subj_tab];
    else 
        dist_tab = subj_tab;
    end

    subj_tab = table(repmat(subj_onset.Patient_id, size(imprint,2),1),...
        (1:size(imprint,2))', onset_reg_count_sz, onset_reg_vol_sz, 'VariableNames', ["Subj_id","sz","n_regions_upd", "region_vol_upd"]);

    if exist('robust_tab', 'var')
        robust_tab = [robust_tab; subj_tab];
    else 
        robust_tab = subj_tab;
    end
end

%%
diff_thresh = 50;
diff_thresh_vol = 30000;
% Setting threshold at 20mm for now
for subj = 1:size(final_output,1)
    subj_dist = dist_tab(dist_tab.Subj_id == string(final_output(subj,:).Patient_id),:);
    subj_prop_diffuse(subj) = sum(subj_dist.max_dist>diff_thresh)/size(subj_dist,1);
    subj_present_diffuse(subj) = sum(subj_dist.max_dist>diff_thresh)>0;
    subj_prop_diffuse_vol(subj) = sum(subj_dist.onset_vol>diff_thresh_vol)/size(subj_dist,1);
    subj_present_diffuse_vol(subj) = sum(subj_dist.onset_vol>diff_thresh_vol)>0;


end

[tab,chi,p] = crosstab(final_output.outcome>2, double(subj_prop_diffuse>0.5));
    
figure()
heatmap(tab,...
    "YData",["ILAE 1-2", "ILAE 3+"], "XData", ["<=50% diffuse", ">50% diffuse"])

title(sprintf("Count of subjects with at least 50%% diffuse onsets (>%dmm) \n chi^2 = %.2f, p = %.3f", diff_thresh,chi,p))

[tab,chi,p] = crosstab(final_output.outcome>2, double(subj_present_diffuse));
figure()
heatmap(tab,...
    "YData",["ILAE 1-2", "ILAE 3+"], "XData", ["No diffuse", ">=1 diffuse onset"])

title(sprintf("Count of subjects with diffuse onsets (>%dmm) \n chi^2 = %.2f, p = %.3f", diff_thresh, chi,p))

[tab,chi,p] = crosstab(final_output.outcome>2, double(subj_prop_diffuse_vol>0.5));
figure()
heatmap(tab,...
    "YData",["ILAE 1-2", "ILAE 3+"], "XData", ["<=50% diffuse", ">50% diffuse"])

title(sprintf("Count of subjects with at least 50%% diffuse onsets (vol>%d) \n chi^2 = %.2f, p = %.3f", diff_thresh_vol,chi,p))

[tab,chi,p] = crosstab(final_output.outcome>2, double(subj_present_diffuse_vol));
figure()
heatmap(tab,...
    "YData",["ILAE 1-2", "ILAE 3+"], "XData", ["No diffuse", ">=1 diffuse onset"])

title(sprintf("Count of subjects with diffuse onsets (vol>%d) \n chi^2 = %.2f, p = %.3f", diff_thresh_vol, chi,p))




%%
diff_reg_count = zeros(size(robust_tab,1),9);

for step = 1:9
    diff_reg_count(:,step) = (robust_tab.n_regions_upd(:,step+1) - robust_tab.n_regions_upd(:,1))./robust_tab.n_regions_upd(:,1);
end

figure()
subplot(211)
hold on
for subj = 1:size(final_output,1)
    plot(robust_tab(string(robust_tab.Subj_id) == string(final_output.Patient_id{subj}),:).n_regions_upd')
end
hold off
xlabel("Number of update steps")
ylabel("Count of regions in onset")
set(gca, "XTick", 1:11, "XTickLabel", 0:10)
title("Comparing number of onset regions at each update step") 

subplot(212)
boxchart(diff_reg_count, 'MarkerStyle', 'none')
hold on
swarmchart(1:9,diff_reg_count, 'k','filled')
hold off

xlabel("Number of update steps")
ylabel("Ratio between count of regions in onset between steps i and original onset") 
title("Comparing number of onset regions between update steps and original (ratio)")

          
   

















% %% Load onset table
% final_output_all = load('tables/final_output_all_sz.mat');
% final_output_all = final_output_all.final_output;
% 
% % % Add in patient outcome column
% % for subj = 1:size(final_output_all, 1)
% %     subj_onset = final_output_all(subj,:);
% %     outcome_id = subj_onset.("Outcome year"){:}-subj_onset.("Surgery year") == 1;
% %     final_output_all.outcome(subj) = subj_onset.Surgery_outcome{:}(outcome_id);
% % end
% 
% % Remove patients with no labelled CLO or no outcome 
% final_output = final_output_all(final_output_all.outcome ~= 8,:);
% 
% % Load atlas
% load('roi_info/ATLAS.mat')
% 
% addpath(genpath('Nathan Code'))
% 
% %% Choose subject to assess
% subj = 1;
% subj_onset = final_output(subj,:);
% imprint = subj_onset.imprint_roi_120{:}; 
% regions = subj_onset.roi_names_120{:};
% 
% % Load touching matrices
% load('Nathan Code/SimpleRouteExample/touching_struct.mat', 'touching')
% gm_touching = touching.touching;
% subcort_touching = touching.wm_subcortical_to_subcortical;
% 
% % Create matrix with grey matter cortical-cortical connections and white
% % matter subcortical-subcortical connections
% touch_gm_subcort = gm_touching+subcort_touching;
% touch_gm_subcort = touch_gm_subcort - diag(diag(touch_gm_subcort));
% 
% % Extract region names for all regions in Lausanne 120 atlas
% names = atlas.name{2};
% 
% % Isolate only the recorded hemisphere to make computations simpler 
% % and more efficient
% if any(contains(regions,"l.")) & any(contains(regions,"r."))
%     % Bilateral placement
%     incl_reg = repelem(1,length(names),1);
%     y_lim = [0,440];
% elseif any(contains(regions,"l.")) & ~any(contains(regions,"r."))
%     % Left hemisphere placement
%     incl_reg = contains(names,"Left"|"l.");
%     y_lim = [230,440];
% elseif ~any(contains(regions,"l.")) & any(contains(regions,"r."))
%     % Right hemisphere placement
%     incl_reg = contains(names,"Right"|"r.");  
%     y_lim = [0,210];
% end
% 
% % Keep only the names of the regions in the relavant hemisphere
% names_hem = names(find(incl_reg));
% 
% % Create subsection of touching matrix for relevant regions
% % (n regions in hemisphere x n regions in hemisphere)
% region_index = find_preserved(names_hem,regions);
% imprint = logical(imprint);
% 
% % For graph theory 'learning' (ask Iain for correct nomenclature), we will
% % exclude recorded regions that we know are non-onset
% 
% %% Seizure specific process
% 
% % There does not appear to be any changes between update steps - need to
% % work out why this is the case!!!! SARAH 
% 
% 
% %for sz = 1:5%size(imprint,2)% Assess one seizure at a time 
% sz = 1;
% 
%     % Extend imprint to include all regions in relevant hemisphere
%     imprint_full = zeros(length(names_hem),1);
%     imprint_full(region_index(imprint(:,sz))) = 1;
%     % We give recorded regions that are known to be non-onset a value of -1 so
%     % they can be removed prior to creating graphs
%     imprint_full(region_index(imprint(:,sz)==0)) = -1;
%     % Remove recorded non-onset regions
%     imprint_full_clean = imprint_full(imprint_full~=-1);
%     
%     % Create subsection of touching matrix to only include relevant hemisphere
%     % and exclude recorded non-onset regions
%     touch_gm_subcort_hem = touch_gm_subcort(find(incl_reg),find(incl_reg));
%     touch_gm_subcort_clean = touch_gm_subcort_hem(imprint_full~=-1,imprint_full~=-1);
%     
%     names_clean = names_hem(imprint_full~=-1);
% 
%      xyz = atlas.xyz{2};
%      xyz_hem = xyz(find(incl_reg),:);
%      xyz_clean = xyz_hem(imprint_full~=-1,:);
%     
%     % Any variable with "_clean" has had recorded non-onset regions removed as
%     % we know that they were notincluded in the onset.
%     % We are only interested in unrecorded regions which MAY have been involved
%     % in the onset. 
%     
%     incl_thresh = 1/3;
%     
%     cm1 = [0.7,0.7,0.7;0,0.9,0.9;0,0.7,0.7];
%     cm1 = interp1(cm1, 1:0.01:size(cm1,1));
%     
%     cm2 = [.1,.1,.1;.1,.1,.1;1,1,1;0,0.9,0.9;0,0.7,0.7];
%     cm2 = interp1(cm2, 1:0.01:size(cm2,1));
%     
%     figure(1)
%     tiledlayout(5,2)
%     % Create graph (undirected) based on touching matrix
%     % Create a graph (undirected) between onset and unknown regions
%     g = graph(touch_gm_subcort_clean);
%     upd_onset = nan(height(g.Nodes),5);
%     
%     f1 = figure();
%     tiledlayout(5,1)
%     sgtitle(string(subj_onset.Patient_id))
%     f2 = figure();
%     tiledlayout(5,1)
%     sgtitle(string(subj_onset.Patient_id))
%     
%     
%     for upd_step = 1:5
%         % Set baseline graph as the graph from the previous step
%         if upd_step>1
% %             base_g = upd_g;
%         else % First step, we set baseline graph as the original graph
%             upd_g = g;
%             upd_onset = imprint_full_clean;
%         end
%         
%     
%         for edge = 1:length(g.Edges.Weight)
%             upd_g.Edges.Weight(edge) = mean([upd_onset(upd_g.Edges.EndNodes(edge,1)), upd_onset(upd_g.Edges.EndNodes(edge,2))]);
%         end
%         % Compute weight of edges as the average of weights of the connected edges
%         avg_weight = nan(height(upd_g.Nodes),1);
%         for node = 1:height(upd_g.Nodes)
%             node_vals = upd_g.Edges.Weight(upd_g.Edges.EndNodes(:,1) == node | upd_g.Edges.EndNodes(:,2) == node);
%             if length(node_vals)>1 % if the region is connected to at least two edges, we compute the mean
%                 avg_weight(node,upd_step) = mean(node_vals);
%             else % if the region is connected to at fewer than two edges, connected onset node (only connection) will mean that this regions is automatically extended into in the second update step
%                 avg_weight(node,upd_step) = 0; 
%             end
%         end
% 
%         % I think this is where the issue is which is preventing updates -
%         % ask Iain to check
%     
%         upd_g_output{upd_step} = upd_g;
%         
%         onset_upd = avg_weight(:,upd_step)>incl_thresh;
%         addit_regions = upd_onset;
%         addit_regions(onset_upd==1 & upd_onset==0) = 0.5;
%     
%         %upd_onset_mat(:,upd_step) = addit_regions;
%     
%         set(0, 'CurrentFigure', f1)
%         nexttile
%         h = plot(upd_g,'LineWidth',(upd_g.Edges.Weight+0.01)*3, "NodeLabel",{},'XData',xyz_clean(:,1),'YData',xyz_clean(:,2),'ZData',xyz_clean(:,3));
%         colormap(cm1)
%         h.EdgeCData = upd_g.Edges.Weight;
%         h.NodeCData = addit_regions;
%         title(sprintf("Update step %d", upd_step))
%         axis off
%         view([200 110 0]);
%         xlabel("x")
%         ylabel("y")
%         
%     
%         set(0, 'CurrentFigure', f2)
%         nexttile;
%         onset_upd_hem = [addit_regions; repmat(-1,sum(imprint_full==-1),1)];
%         names_hem_upd = cellstr([string(names_clean); string(names_hem(imprint_full==-1))]);
%         plotBrain_NE(names_hem_upd, onset_upd_hem ,'cm',cm2);
%         colorbar off
%         title(sprintf("Update step %d", upd_step))
%         ylim(y_lim)
%     
%     end
%     
%     %saveas(f1, sprintf('figures/rapid_prototype/onset_graph/%s_graph_%d.png',string(subj_onset.Patient_id), sz))
%     %saveas(f2, sprintf('figures/rapid_prototype/onset_graph/%s_brain_%d.png',string(subj_onset.Patient_id), sz))
% %end
% %%
% 
% 
% % %% Create graph (undirected) based on touching matrix
% % % Create a graph (undirected) between onset and unknown regions
% % g = graph(touch_gm_subcort_clean);
% % %g.Nodes = names_clean;
% % 
% % % First update step
% % % Set weights as the average across end nodes
% % for edge = 1:length(g.Edges.Weight)
% %      g.Edges.Weight(edge) = mean([imprint_full_clean(g.Edges.EndNodes(edge,1)), imprint_full_clean(g.Edges.EndNodes(edge,2))]);
% % end
% % % Compute weight of edges as the average of weights of the connected edges
% % avg_weight = nan(height(g.Nodes),1);
% % for node = 1:height(g.Nodes)
% %     node_vals = g.Edges.Weight(g.Edges.EndNodes(:,1) == node | g.Edges.EndNodes(:,2) == node);
% %     if length(node_vals)>1 % if the region is connected to at least two edges, we compute the mean
% %         avg_weight(node) = mean(node_vals);
% %     else % if the region is connected to at fewer than two edges, connected onset node (only connection) will mean that this regions is automatically extended into in the second update step
% %         avg_weight(node) = 0; 
% %     end
% % end
% % 
% % % Second update step
% % % Set weights as the average across end nodes
% % g_upd = g;
% % for edge = 1:length(g.Edges.Weight)
% %     g_upd.Edges.Weight(edge) = mean([avg_weight(g_upd.Edges.EndNodes(edge,1)), avg_weight(g_upd.Edges.EndNodes(edge,2))]);
% % end
% % 
% % avg_weight_upd = nan(height(g.Nodes),1);
% % for node = 1:height(g.Nodes)
% %     node_vals = g_upd.Edges.Weight(g_upd.Edges.EndNodes(:,1) == node | g_upd.Edges.EndNodes(:,2) == node);
% %     if length(node_vals)>1
% %         avg_weight_upd(node) = mean(node_vals);
% %     else
% %         avg_weight_upd(node) = 0;
% %     end
% % end
% % 
% % %% Create the updated onset after two update steps (Threshold of 0.3333)
% % incl_thresh = 1/3;
% % onset_upd = avg_weight>incl_thresh;
% % addit_regions_upd = imprint_full_clean;
% % addit_regions_upd(onset_upd==1 & imprint_full_clean==0) = 0.5;
% % 
% % onset_upd2 = avg_weight>incl_thresh;
% % addit_regions_upd2 = imprint_full_clean;
% % addit_regions_upd2(onset_upd2==1 & imprint_full_clean==0) = 0.5;
% % 
% % onset_upd2_hem = [addit_regions_upd2; repmat(-1,sum(imprint_full==-1),1)];
% % names_hem_upd = cellstr([string(names_clean); string(names_hem(imprint_full==-1))]);
% % 
% % figure(3)
% % subplot(2,1,1)
% % % plotBrain_NE([names_clean; names(imprint_full_atlas==-1)], [imprint_full_atlas_clean; repmat(-1,sum(imprint_full_atlas==-1),1)], 'cm',cm);
% % plotBrain_NE(regions,imprint(:,sz), 'cm',[0.2,0.2,0.2;0,0.7,0.7] );
% % colorbar off
% % ylim(y_lim)
% % title("Original")
% % 
% % cm = [.1,.1,.1;.1,.1,.1;1,1,1;0,0.9,0.9;0,0.7,0.7];
% % cm = interp1(cm, 1:0.01:size(cm,1));
% % 
% % subplot(2,1,2)
% % plotBrain_NE(names_hem_upd, onset_upd2_hem , 'cm',cm);
% % colorbar off
% % ylim(y_lim)
% % title("Onset after two onset steps")
% % sgtitle(string(subj_onset.Patient_id))
% % 
% % xyz = atlas.xyz{2};
% % xyz_hem = xyz(find(incl_reg),:);
% % xyz_clean = xyz_hem(imprint_full~=-1,:);
% % 
% % %% Plot graphs of each update step (including xyz coordinates of included regions
% % cm = [0.7,0.7,0.7;0.85,0.325,0.098;1,0,0];
% % cm = interp1(cm, 1:0.01:size(cm,1));
% % colormap(cm)
% % 
% % figure(4)
% % subplot(311)
% % h1 = plot(g,'XData',xyz_clean(:,1),'YData',xyz_clean(:,2),'ZData',xyz_clean(:,3),'LineWidth',((g.Edges.Weight==1)+0.01)*3,"NodeLabel",{});
% % colormap(cm)
% % h1.EdgeCData = g.Edges.Weight==1;
% % h1.NodeCData = imprint_full_clean;
% % title("Original Graph")
% % axis off
% % 
% % subplot(312)
% % h2 = plot(g,'LineWidth',(g.Edges.Weight+0.01)*3, "NodeLabel",{}, 'XData',xyz_clean(:,1),'YData',xyz_clean(:,2),'ZData',xyz_clean(:,3));
% % colormap(cm)
% % h2.EdgeCData = g.Edges.Weight;
% % h2.NodeCData = onset_upd;
% % title("Update step one")
% % axis off
% % 
% % subplot(313)
% % h3 = plot(g_upd,'LineWidth',(g_upd.Edges.Weight+0.01)*3, "NodeLabel",{},'XData',xyz_clean(:,1),'YData',xyz_clean(:,2),'ZData',xyz_clean(:,3));
% % colormap(cm)
% % h3.EdgeCData = g_upd.Edges.Weight;
% % h3.NodeCData = onset_upd2;
% % title("Update step two")
% % axis off