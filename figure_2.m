dists = atlas.dists{2};
vols = atlas.vol{2};
names = atlas.name{2};

subj = 40;
subj_ons = final_output(subj,:);
imprint = logical(subj_ons.imprint_roi_120{:});
regions = subj_ons.roi_names_120{:};
regions_relab = strrep(regions, 'r.', 'ctx-rh-');
regions_relab = strrep(regions_relab, 'l.', 'ctx-lh-');
sz_count = size(imprint,2);

cm = [0.7,0.7,0.7;0,0.2,1];
cm = interp1(cm, 1:0.01:size(cm,1));

f = figure("Position",[10,10,2500,400]);
for sz = 1:size(imprint,2)
    if min(imprint(:,sz)) ~= max(imprint(:,sz))
        sz_region_atlas = contains(names, regions(imprint(:,sz)));
        sz_dists = dists(sz_region_atlas,sz_region_atlas);
        sz_max_dist = max(max(sz_dists));
        sz_vol = sum(vols(sz_region_atlas));
        plotBrain(regions_relab, double(imprint(:,sz)), cm, 'atlas', 'lausanne120_aseg');
%         colorbar off
%         ylim([230 450])
%         title(sprintf("Seizure %d \n (max dist %3.2f, volume %5.f)", sz, sz_max_dist,sz_vol))
    else
        axis off
    end
end
sgtitle(sprintf("%s (ILAE %d)", string(subj_ons.Patient_id), subj_ons.outcome ))
saveas(f, sprintf('figures/paper_figures/Figure 2/%s_onsets.svg', string(subj_ons.Patient_id)))

dist_mat = nan(length(regions),length(regions));
for reg1 = 1:length(regions)
    for reg2 = 1:length(regions)
        if reg1<reg2
            dist_mat(reg1,reg2) = dists(contains(names,regions(reg1)), contains(names,regions(reg2)));
        end
    end
end

% f = figure("Position", [10,10,1000,800]);
% heatmap(dist_mat(1:(end-1),2:end))
% set(gca, "XData",regions(2:end), "YData", regions(1:(end-1)))
% saveas(f, sprintf('figures/paper_figures/Figure 2/%s_distances.svg', string(subj_ons.Patient_id)))

% vol_arr = nan(length(regions),1);
% for reg = 1:length(regions)
%     vol_arr(reg) = vols(contains(names,regions(reg)));
% end
% 
% f = figure("Position", [10,10,400,1000]);
% heatmap(vol_arr, "CellLabelFormat", '%5.f')
% set(gca, "YData", regions)
%saveas(f, sprintf('figures/paper_figures/Figure 2/%s_volumes.svg', string(subj_ons.Patient_id)))

% 
% f = figure("Position",[0,0,2000,2000]);
% tiledlayout(5,5)
% for reg = 1:length(regions)
%     nexttile
%     plot_vals = zeros(length(regions),1);
%     plot_vals(reg) = 1;
%     plotBrain_NE(regions,plot_vals , "cm", [0.7,0.7,0.7;1,0,0]);
%     title(regions(reg))
%     colorbar off
%     ylim([230 450])
% end
%saveas(f, sprintf('figures/paper_figures/Figure 2/%s_highlighted_regions.svg', string(subj_ons.Patient_id)))


resec = logical(subj_ons.resected_roi_120{:});
f = figure("Position",[10,10,2500,400]);
resec_region_atlas = contains(names, regions(resec));
resec_vol = sum(vols(resec_region_atlas));

plotBrain_NE(regions, resec, "cm", [0.7,0.7,0.7;0,0,0]);
colorbar off
ylim([230 450])
sgtitle(sprintf("%s resection", string(subj_ons.Patient_id)))
saveas(f, sprintf('figures/paper_figures/Figure 2/%s_resection.svg', string(subj_ons.Patient_id)))

clo = logical(subj_ons.clo_roi_120{:});
f = figure("Position",[10,10,2500,400]);
clo_region_atlas = contains(names, regions(clo));
clo_vol = sum(vols(clo_region_atlas));

clo_dist = dists(clo_region_atlas,clo_region_atlas);
clo_max_dist = max(max(clo_dist));
plotBrain_NE(regions, clo, "cm", [0.7,0.7,0.7;0,0.2,1]);
colorbar off
ylim([230 450])
sgtitle(sprintf("%s CLO \n (max dist %3.2f, volume %5.f)", string(subj_ons.Patient_id), clo_max_dist,clo_vol))
saveas(f, sprintf('figures/paper_figures/Figure 2/%s_clo.svg', string(subj_ons.Patient_id)))

% Plotting onset and resection
f = figure("Position",[10,10,2500,400]);
tiledlayout(2,3)
for sz = 1:min([6, size(imprint,2)])
    nexttile
    if min(imprint(:,sz)) ~= max(imprint(:,sz))
        sz_region_atlas = contains(names, regions(imprint(:,sz)));
       
        onset_resec =  2*imprint(:,sz)+resec; 
        onset_resec_region_atlas = contains(names, regions(onset_resec==3));
        onset_resec_vol = sum(vols(onset_resec_region_atlas));
        
        plotBrain_NE(regions, onset_resec, "cm", [.7,.7,.7;0,0,0;0,0.2,1;0,0,0]);
        colorbar off
        ylim([230 450])
        title(sprintf("Seizure %d \n (volume resected %5.f)", sz,onset_resec_vol))
    else
        axis off
    end
end
sgtitle(sprintf("%s (ILAE %d)", string(subj_ons.Patient_id), subj_ons.outcome ))
saveas(f, sprintf('figures/paper_figures/Figure 2/%s_onset_resec.svg', string(subj_ons.Patient_id)))

% Plotting CLO and resection
f = figure("Position",[10,10,2500,400]);
clo_region_atlas = contains(names, regions(clo));
clo_vol = sum(vols(clo_region_atlas));
onset_resec =  2*clo+resec; 

onset_resec_region_atlas = contains(names, regions(onset_resec==3));
onset_resec_vol = sum(vols(onset_resec_region_atlas));
        
plotBrain_NE(regions, onset_resec, "cm", [.7,.7,.7;0,0,0;0,0.2,1;0,0,0]);
colorbar off
ylim([230 450])
title(sprintf("(volume resected %5.f)",onset_resec_vol))
saveas(f, sprintf('figures/paper_figures/Figure 2/%s_clo_resec.svg', string(subj_ons.Patient_id)))

%%
sub = 40;
clo = final_output(sub,:).clo_roi_120{:};
imprint = mean(final_output(sub,:).imprint_roi_120{:},2)>0.5;
resect = final_output(sub,:).resected_roi_120{:};
regions = final_output(sub,:).roi_names_120{:};

%% Set brain plot y limits to only show relevant hemisphere

if any(contains(regions, "r.")) & ~any(contains(regions, "l."))
    % Right hemisphere placement
    y_lim = [0, 220];

elseif any(contains(regions, "l.")) & ~any(contains(regions, "r."))
    % Left hemisphere placement
    y_lim = [230, 450];
elseif any(contains(regions, "r.")) & any(contains(regions, "l."))
    % Bilateral placement
    y_lim = [0, 450];
end


% Give each region a colour
n_regions = length(regions);
rng(2)
region_index = randperm(n_regions);

% Create colour map
cm = [1,0,0; 1,1,0; 0,1,0; 0,1,1; 0,0,1];
cm = interp1(cm,  1+(1:n_regions)/(n_regions/4));
cm = [0.8,0.8,0.8;cm];

f = figure("Position",[100,100,1000,1000]);
subplot(311)
plotBrain_NE(regions, clo.*region_index', "cm", cm);
colorbar off
ylim(y_lim)
clim([0, max(clo.*region_index')])
title("CLO")

subplot(312)
plotBrain_NE(regions, imprint.*region_index', "cm", cm);
clim([0, max(imprint.*region_index')])
colorbar off
ylim(y_lim)

title("Consensus imprint")

subplot(313)
plotBrain_NE(regions, resect.*region_index', "cm", cm);
colorbar off
ylim(y_lim)
clim([0, max(resect.*region_index')])
title("Resection")

sgtitle(string(final_output(sub,:).Patient_id))
saveas(f, sprintf("figures/subquestions/%s_regions.svg", string(final_output(sub,:).Patient_id)))


