figure()
subplot(131)
imagesc(saved_output.imprint_roi{1,1})
set(gca, 'YTick', 1:length(saved_output.ROI_ids{1,1}), 'YTickLabel', saved_output.ROI_ids{1,1})
title("Atlas 125")

subplot(132)
imagesc(final_output.imprint_roi{1,1})
set(gca, 'YTick', 1:length(final_output.ROI_ids{1,1}), 'YTickLabel', final_output.ROI_ids{1,1})
title("Atlas 125 (3sec)")

subplot(133)
imagesc(atl_60_3sec_output.imprint_roi{1,1})
set(gca, 'YTick', 1:length(atl_60_3sec_output.ROI_ids{1,1}), 'YTickLabel', atl_60_3sec_output.ROI_ids{1,1})
title("Atlas 60 (3sec)")

%%
patients = unique(final_output.Patient_id);

for pat = 1:length(patients)
    pat_output = final_output(final_output.Patient_id == patients(pat),:);
    figure('Position', [10 10 1500 1200])
    if size(pat_output,1) == 6
        for i = 1:6
            subplot(2,3,i)
            imagesc(pat_output.imprint_roi{i,1})
            if i == 1 || i == 4
                set(gca, 'YTick', 1:length(final_output.ROI_ids{i,1}), 'YTickLabel', pat_output.ROI_ids{i,1})
            else
                set(gca, 'YTick', [], 'YTickLabel', pat_output.ROI_ids{i,1})
            end
            title(sprintf("Atlas %d, PropThresh %s", pat_output.atl(i), pat_output.det(i)))
        end
    else
        

    end
    sgtitle(pat_output.Patient_id(1))
end

%%
patients = unique(final_output.Patient_id);
atl_group = [60; 60; 60; 125; 125; 125];
thres_group = [0; 1; 2; 0; 1; 2];
sub_pl = [1; 2; 3; 4; 5; 6];
plot_orientation_tab = table(atl_group, thres_group, sub_pl);

for pat = 1:length(patients)
    pat_output = final_output(final_output.Patient_id == patients(pat),:);
    f = figure(pat);
    f.Position = [10 10 1500 1200];
    for i = 1:6
        subpl_output = pat_output(pat_output.atl == atl_group(i) & pat_output.det == sprintf("+%d seconds",thres_group(i) ),:);
        if ~isempty(subpl_output)
            subplot(2,3,i)
            imagesc(subpl_output.imprint_roi{1,1})
            if i == 1 || i == 4
                set(gca, 'YTick', 1:length(subpl_output.ROI_ids{1,1}), 'YTickLabel', subpl_output.ROI_ids{1,1})
            else
                set(gca, 'YTick', [])
            end
            title(sprintf("Atlas %d, PropThresh %s", subpl_output.atl(1), subpl_output.det(1)))
        end
    end

    sgtitle(pat_output.Patient_id(1))
end






%% Comparisons code


% This script creates a table storing all comparison measures for comparing
% seizures with resected regions (comparison = "resection") or between
% pairs of seizures (comparison = "pairwise)

% Onsets are automatically detected using one of three methods:
%   "imprint": Imprint onset (based on Gascoigne et al., 2023)
%   "EI": Epileptogenicity Index (Bartolomei et al., 2008)
%   "PLHG": Phase-Locked High-Gamma (Weiss et al., 2013)
% comparison = "pairwise";
% det_method = "imprint";

load('roi_info/ATLAS_.mat')
load('roi_info/retains.mat')
load('roi_info/atlasinfo.mat')
load('roi_info/ATLAS.mat')
%load('roi_info/lhSurf.mat') 
%load('roi_info/rhSurf.mat')

addpath('sarah_functions/')
addpath(genpath('help_functions/'))


%%

load("final_output_min_5.mat");
%%
% 
% % Ensure there are no duplicate patients
% patients = final_output.Patient_id;
% [~, ind] = unique(patients);
% final_output = final_output(ind,:);

save_table = 0;
save_plot = 0;
chan_or_roi = "roi"; % I need to add Hausdorff to channel-level comparisons
atl_id =  [36; 60; 125; 250];
col_id = [1;2;3;4];
atl_tab = table(atl_id, col_id);
comparison = "resection";
det_method = "imprint";


%%
patients = unique(final_output.Patient_id);
for atl = [60, 125]
    for det = [0, 1, 2]
        if atl == 60 && det == 0
            continue
        end

        close all 
        clear final_comp
        fprintf('%d atlas with threshold %d seconds \n', atl, det)

        comp_output = final_output(final_output.atl == atl & final_output.det == sprintf("+%d seconds", det ),:);
     
        % Iterate through patients
        for pat = 1:length(patients)
            patient = patients(pat);
            pat_onset = comp_output(comp_output.Patient_id == patient,:);

            if isempty(pat_onset)
                continue
            end
                            
            % Compute normalised Jacccard's index and Sorensen-Dice coefficient
            jac_sor_table = calc_jacc_sorr(pat_onset, "det_method",...
                det_method, "comparison", comparison, 'n_perm',1000,...
                "chan_or_roi", chan_or_roi);
            if height(jac_sor_table) == 0
                continue
            end
            pat_comp = jac_sor_table;
            if chan_or_roi == "roi"
                haus_table = calc_haus(pat_onset, ...
                    atlas(table2array(atl_tab(atl_tab.atl_id == atl,"col_id")),:), "det_method", det_method, "comparison", comparison);
                pat_comp = join(jac_sor_table, haus_table);
            end
            
%             if comparison == "resection"  & det_method == "CLO"
%                 onset_resec = pat_onset.labelled_onset_chan{:} + pat_onset.resected_chan{:} == 2;
%                 perc_chan_resec = sum(onset_resec)/sum(pat_onset.labelled_onset_chan{:});
%                 if chan_or_roi == "roi"
%                     pat_comp = join([jac_sor_table, table(perc_chan_resec)], haus_table);
%                 else 
%                     pat_comp = [jac_sor_table, table(perc_chan_resec)];
%                 end
%                     
%             end

            pat_surg_out = pat_onset.Surgery_outcome{:};
            pat_out_year = pat_onset.Outcome_year{:};
            year_1 = pat_out_year - str2double(pat_onset.Surgery_year) == 1;
            pat_comp.Y1_outcome = repmat(pat_surg_out(year_1),size(pat_comp,1),1);
            
            if pat_comp.Y1_outcome == 8
                pat_comp.Y1_outcome = repmat(NaN,size(pat_comp,1),1);
            end
                           
            if exist('final_comp', 'var')
                final_comp = [final_comp; pat_comp];
            else
                final_comp = pat_comp;
            end
        end
          
        writetable(final_comp, sprintf('tables/atl_%d_thresh_%d.csv', atl, det))
        % Save table for future analyses
        if save_table == 1
            save(sprintf('tables/across_pat_%d_%d.mat', atl, det), 'final_comp')
            % Also save as a CSV to open in R
            writetable(final_comp, sprintf('tables/atl_%d_thresh_%d.csv', atl, det))
        end

%         if comparison == "resection"
%             comp_measures = {final_comp.Properties.VariableNames{3:end-1}};
%         else 
%             comp_measures = {final_comp.Properties.VariableNames{4:end-1}};
%         end
% 
%         
%         for comp = 1:length(comp_measures)
%             if det_method ~= "CLO"
%                 beeswarm_across_pat(final_comp, final_output, "det_method",det_method,...
%                     "comparison",comparison, "comp_measure",comp_measures{comp},...
%                     "save_plot",save_plot,"save_location","figures/across_patients/",...
%                     "chan_or_roi", chan_or_roi)
%             end
%         end
    end  
end

