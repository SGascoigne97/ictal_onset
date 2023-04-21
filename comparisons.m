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

load("final_output.mat");

% Ensure there are no duplicate patients
patients = final_output.Patient_id;
[~, ind] = unique(patients);
final_output = final_output(ind,:);

save_table = 1;
save_plot = 1;
chan_or_roi = "chan"; % I need to add Hausdorff to channel-level comparisons

for comparison = "pairwise" %["resection", "pairwise"]
    for det_method = ["EI", "PLHG"] %["CLO", "imprint", "EI", "PLHG"] %
        close all 
        clear final_comp
        fprintf('%s and %s \n', comparison, det_method)
        if comparison == "pairwise" & det_method == "CLO"
            fprintf("Comparison not possible \n")
        else 
                % Iterate through patients
            for pat = 1:length(patients)
                patient = patients(pat);
                pat_onset = final_output(final_output.Patient_id == patient,:);
                                
                % Compute normalised Jacccard's index and Sorensen-Dice coefficient
                jac_sor_table = calc_jacc_sorr(pat_onset, "det_method",...
                    det_method, "comparison", comparison, 'n_perm',1000,...
                    "chan_or_roi", chan_or_roi);
                if height(jac_sor_table) == 0
                    continue
                end
                if chan_or_roi == "roi"
                    haus_table = calc_haus(pat_onset, atlas(3,:), "det_method", det_method, "comparison", comparison);
                    pat_comp = join(jac_sor_table, haus_table);
                end
                pat_comp = jac_sor_table;

                if comparison == "resection"  & det_method == "CLO"
                    onset_resec = pat_onset.labelled_onset_chan{:} + pat_onset.resected_chan{:} == 2;
                    perc_chan_resec = sum(onset_resec)/sum(pat_onset.labelled_onset_chan{:});
                    if chan_or_roi == "roi"
                        pat_comp = join([jac_sor_table, table(perc_chan_resec)], haus_table);
                    else 
                        pat_comp = [jac_sor_table, table(perc_chan_resec)];
                    end
                        
                end

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
              
            % Save table for future analyses
            if save_table == 1
                save(sprintf('tables/across_pat_%s_%s.mat', det_method, comparison), 'final_comp')
                % Also save as a CSV to open in R
                writetable(final_comp, sprintf('tables/across_pat_%s_%s_%s.csv', det_method, comparison, chan_or_roi))
            end
           
            % Beeswarm plot across patients
%             if comparison == "resection"
%                 comp_measures = ["Jaccard_norm", "Sorensen_norm", "Percentage_resec", "perc_chan_resec", "Hausdorff_norm"];
%             elseif comparison == "pairwise"
%                 comp_measures = ["Jaccard_norm", "Sorensen_norm", "Hausdorff_norm"];
%             end

            if comparison == "resection"
                comp_measures = {final_comp.Properties.VariableNames{3:end-1}};
            else 
                comp_measures = {final_comp.Properties.VariableNames{4:end-1}};
            end

            
            for comp = 1:length(comp_measures)
                if det_method ~= "CLO"
                    beeswarm_across_pat(final_comp, final_output, "det_method",det_method,...
                        "comparison",comparison, "comp_measure",comp_measures{comp},...
                        "save_plot",save_plot,"save_location","figures/across_patients/",...
                        "chan_or_roi", chan_or_roi)
                end
            end
        end
    
        
    end
end

%% Compare values using channels and regions (percentage resected)
fun = @(x) x(1);
final_output.Outcome_1 = cellfun(fun, final_output.Surgery_outcome);
final_output.Outcome_1(final_output.Outcome_1 == 8) = NaN;

figure()
scatter3(final_comp.Percentage_resec*100, final_comp.perc_chan_resec*100, final_output.Outcome_1>2)
lsline
xlabel('Region-wise')
ylabel('Channel-wise')
title('Plot of percentage resected using channels and regions')
%legend({'good','bad','good','bad'})