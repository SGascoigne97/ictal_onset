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
load('roi_info/lhSurf.mat') 
load('roi_info/rhSurf.mat')

%%

%load("final_output.mat");

% Ensure there are no duplicate patients
patients = final_output.Patient_id;
[~, ind] = unique(patients);
final_output = final_output(ind,:);


save_table = 1;
save_plot = 1;

for comparison = "pairwise" %["resection", "pairwise"]
    for det_method = ["CLO", "imprint", "EI", "PLHG"] %
        close all 
        clear final_comp
        sprintf('%s and %s', comparison, det_method)
        if comparison == "pairwise" & det_method == "CLO"
            sprintf("Comparison not possible")
        else 
                % Iterate through patients
            for pat = 1:length(patients)
                patient = patients(pat);
                pat_onset = final_output(final_output.Patient_id == patient,:);
                
                % Compute normalised Jacccard's index and Sorensen-Dice coefficient
                jac_sor_table = calc_jacc_sorr(pat_onset, "det_method",...
                    det_method, "comparison", comparison, 'n_perm',1000);
                % Need to write script to compute other measures 
                haus_table = calc_haus(pat_onset, atlas(3,:), "det_method", det_method, "comparison", comparison);
                
                pat_comp = join(jac_sor_table, haus_table);
                
                if exist('final_comp', 'var')
                    final_comp = [final_comp; pat_comp];
                else
                    final_comp = pat_comp;
                end
            end
            % Add in surgical outcome
            final_comp = join(final_comp, final_output(:,[1 9]));
          
            % Save table for future analyses
            if save_table == 1
                save(sprintf('tables/across_pat_%s_%s.mat', det_method, comparison), 'final_comp')
                % Also save as a CSV to open in R
                writetable(final_comp, sprintf('tables/across_pat_%s_%s.csv', det_method, comparison))
            end
           
            % Beeswarm plot across patients
            if comparison == "resection"
                comp_measures = ["Jaccard_norm", "Sorensen_norm", "Percentage_resec", "Hausdorff_norm"];
            elseif comparison == "pairwise"
                comp_measures = ["Jaccard_norm", "Sorensen_norm", "Hausdorff_norm"];
            end
            
            for comp = 1:length(comp_measures)
                if det_method ~= "CLO"
                    beeswarm_across_pat(final_comp, final_output, "det_method",det_method,...
                        "comparison",comparison, "comp_measure",comp_measures(comp),...
                        "save_plot",save_plot,"save_location","figures/across_patients/")
                end
            end
        end
    
        
    end
end

% Need to add code to perform statistical analysis on output tables
% Add code comparing CLO with resection 