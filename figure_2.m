dists = atlas.dists{2};
vols = atlas.vol{2};
names = atlas.name{2};

subj = 40;
subj_ons = final_output(subj,:);
imprint = logical(subj_ons.imprint_roi_120{:});
regions = subj_ons.roi_names_120{:};
sz_count = size(imprint,2);

f = figure("Position",[10,10,2500,400]);
tiledlayout(2,3)
for sz = 1:min([6, size(imprint,2)])
    nexttile
    if min(imprint(:,sz)) ~= max(imprint(:,sz))
        sz_region_atlas = contains(names, regions(imprint(:,sz)));
        sz_dists = dists(sz_region_atlas,sz_region_atlas);
        sz_max_dist = max(max(sz_dists));
        sz_vol = sum(vols(sz_region_atlas));

        plotBrain_NE(regions, imprint(:,sz), "cm", [0.7,0.7,0.7;1,0,0]);
        colorbar off
        ylim([230 450])
        title(sprintf("Seizure %d \n (max dist %3.2f, volume %5.f)", sz, sz_max_dist,sz_vol))
    else
        axis off
    end
end
sgtitle(sprintf("%s (ILAE %d)", string(subj_ons.Patient_id), subj_ons.outcome ))
saveas(f, sprintf('figures/subquestions/Figure 2/%s_onsets.svg', string(subj_ons.Patient_id)))

dist_mat = nan(length(regions),length(regions));
for reg1 = 1:length(regions)
    for reg2 = 1:length(regions)
        if reg1<reg2
            dist_mat(reg1,reg2) = dists(contains(names,regions(reg1)), contains(names,regions(reg2)));
        end
    end
end

f = figure("Position", [10,10,1000,800]);
heatmap(dist_mat(1:(end-1),2:end))
set(gca, "XData",regions(2:end), "YData", regions(1:(end-1)))
saveas(f, sprintf('figures/subquestions/Figure 2/%s_distances.svg', string(subj_ons.Patient_id)))

vol_arr = nan(length(regions),1);
for reg = 1:length(regions)
    vol_arr(reg) = vols(contains(names,regions(reg)));
end

f = figure("Position", [10,10,400,1000]);
heatmap(vol_arr, "CellLabelFormat", '%5.f')
set(gca, "YData", regions)
saveas(f, sprintf('figures/subquestions/Figure 2/%s_volumes.svg', string(subj_ons.Patient_id)))


f = figure("Position",[0,0,2000,2000]);
tiledlayout(5,5)
for reg = 1:length(regions)
    nexttile
    plot_vals = zeros(length(regions),1);
    plot_vals(reg) = 1;
    plotBrain_NE(regions,plot_vals , "cm", [0.7,0.7,0.7;1,0,0]);
    title(regions(reg))
    colorbar off
    ylim([230 450])
end
saveas(f, sprintf('figures/subquestions/Figure 2/%s_highlighted_regions.svg', string(subj_ons.Patient_id)))


resec = logical(subj_ons.resected_roi_120{:});
f = figure("Position",[10,10,2500,400]);
resec_region_atlas = contains(names, regions(resec));
resec_vol = sum(vols(resec_region_atlas));

plotBrain_NE(regions, resec, "cm", [0.7,0.7,0.7;0,0,1]);
colorbar off
ylim([230 450])
sgtitle(sprintf("%s resection", string(subj_ons.Patient_id)))
saveas(f, sprintf('figures/subquestions/Figure 2/%s_resection.svg', string(subj_ons.Patient_id)))

clo = logical(subj_ons.clo_roi_120{:});
f = figure("Position",[10,10,2500,400]);
clo_region_atlas = contains(names, regions(clo));
clo_vol = sum(vols(clo_region_atlas));

clo_dist = dists(clo_region_atlas,clo_region_atlas);
clo_max_dist = max(max(clo_dist));
plotBrain_NE(regions, clo, "cm", [0.7,0.7,0.7;1,0.5,0]);
colorbar off
ylim([230 450])
sgtitle(sprintf("%s CLO \n (max dist %3.2f, volume %5.f)", string(subj_ons.Patient_id), clo_max_dist,clo_vol))
saveas(f, sprintf('figures/subquestions/Figure 2/%s_clo.svg', string(subj_ons.Patient_id)))

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
        
        plotBrain_NE(regions, onset_resec, "cm", [.7,.7,.7;.3,.3,.3;1,0,0;.5,0,.5]);
        colorbar off
        ylim([230 450])
        title(sprintf("Seizure %d \n (volume resected %5.f)", sz,onset_resec_vol))
    else
        axis off
    end
end
sgtitle(sprintf("%s (ILAE %d)", string(subj_ons.Patient_id), subj_ons.outcome ))
saveas(f, sprintf('figures/subquestions/Figure 2/%s_onset_resec.svg', string(subj_ons.Patient_id)))

% Plotting CLO and resection
f = figure("Position",[10,10,2500,400]);
clo_region_atlas = contains(names, regions(clo));
clo_vol = sum(vols(clo_region_atlas));
onset_resec =  2*clo+resec; 

onset_resec_region_atlas = contains(names, regions(onset_resec==3));
onset_resec_vol = sum(vols(onset_resec_region_atlas));
        
plotBrain_NE(regions, onset_resec, "cm", [.7,.7,.7;.3,.3,.3;1,0,0;.5,0,.5]);
colorbar off
ylim([230 450])
title(sprintf("(volume resected %5.f)",onset_resec_vol))

saveas(f, sprintf('figures/subquestions/Figure 2/%s_clo_resec.svg', string(subj_ons.Patient_id)))


