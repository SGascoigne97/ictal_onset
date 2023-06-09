% Comparing detected onset with clinically labelled onset 
% load("tables/all_pat_with_focal_2.mat");

%%
det_comp = struct("imprint", [], "EI", [], "PLHG",[]);
sz_prop_thresh = 0.01; 
tab_type = {'all_sz'; 'per_sz'};
for det_method = ["imprint", "EI", "PLHG"]
    clear all_pat_comp
%for i = 1:2
    per_sz_or_all_sz = 'all_sz';
    for pat = 1:size(final_output_both_atlas,1)
    
        patient = final_output_both_atlas.Patient_id{pat};
        pat_onset = final_output_both_atlas(string(final_output.Patient_id) == patient,: );
    %     if all(isnan(pat_onset.clo_chan{:}))
    %         fprintf("Patient %s does not have clinically labelled onset \n", patient)
    %         continue
    %     end
        %%
        pat_imprint_roi = pat_onset.(sprintf("%s_roi", det_method)){:};
        % Remove any seizures with onset in >50% of regions or no onset
        % detected
       
        rm_sz =  sum(pat_imprint_roi,1) == 0 |...
            sum(pat_imprint_roi,1) >= size(pat_imprint_roi,1)/2;
        pat_imprint_roi = pat_imprint_roi(:,~rm_sz);
    
        % Compare each imprint onset against CLO (with permutation test)
        pat_comp_tab = comp_auto_clo(pat_onset.clo_roi{:},...
            pat_imprint_roi, "per_sz_or_all_sz", per_sz_or_all_sz, "sz_prop_thresh", sz_prop_thresh);
        if isempty(pat_comp_tab)
            fprintf("Patient %s has no CLO \n", patient)
            continue
        end

        pat_surg_out = pat_onset.Surgery_outcome{:};
        pat_surg_year = pat_onset.("Surgery year");
        pat_out_year = pat_onset.("Outcome year"){:};
        year_1 = pat_out_year - pat_surg_year == 1;
        pat_comp_tab = [repmat({patient}, size(pat_comp_tab,1),1), pat_comp_tab,...
            repmat({pat_surg_out(year_1)}, size(pat_comp_tab,1),1), ...
            repmat(pat_onset.("Op type")(1), size(pat_comp_tab,1),1)];
        pat_comp_tab.Properties.VariableNames(1) = "Patient_id";
        pat_comp_tab.Properties.VariableNames(end-1) = "Outcome";
        pat_comp_tab.Properties.VariableNames(end) = "Op_type";

        if exist('all_pat_comp', 'var')
            all_pat_comp = [all_pat_comp; pat_comp_tab];
        else
            all_pat_comp = pat_comp_tab;
        end

        roi_names = pat_onset.roi_names{:};
        roi_names = strrep(roi_names, 'r.', 'ctx-rh-');
        roi_names = strrep(roi_names, 'l.', 'ctx-lh-');
        
        
        min_sz = ceil(sz_prop_thresh*size(pat_imprint_roi,2));
        impr_across = summarise_across_szs(pat_imprint_roi, min_sz);

%         col_sch = 2*impr_across+pat_onset.clo_roi{:};
% 
%         if max(col_sch) == 2   
%             cm = [0.9,0.9,0.9; % neither (0:grey)
%                 0,0,1; % CLO not captured (1: blue)
%                 1,1,0]; % in imprint, not clo (2:yellow)
%         else
%             cm = [0.9,0.9,0.9; % neither (0:grey)
%                 0,0,1; % CLO not captured (1: blue)
%                 1,1,0; % in imprint, not clo (2:yellow)
%                  0,1,0];% CLO captured (3:green)
%         end
%           cm = interp1(cm, 1:0.01:size(cm,1));
% 
%         plotBrain(roi_names, col_sch,...
%        cm, 'atlas', 'lausanne120_aseg', 'savepath',...
%         char(sprintf('./figures/clo_vs_imprint/%s', patient)))
    end
%end
%%

det_comp.(sprintf("%s", det_method)) = all_pat_comp;
end
%%
% 
all_pat_comp_per_sz = all_pat_comp(all_pat_comp.All_or_per_sz =='per_sz',:);
all_pat_comp_across_sz = all_pat_comp(all_pat_comp.All_or_per_sz =='all_sz',:);
%%

figure(1)
subplot(2,2,1)
histogram(all_pat_comp_per_sz.("Percentage of CLO captured"), BinWidth=0.1)
title("Percentage of CLO captured")
subplot(2,2,2)
histogram(all_pat_comp_per_sz.Perc_z)
title("Percentage of CLO captured (z-scored)")
xline(0)
subplot(2,2,3)
histogram(all_pat_comp_per_sz.Jaccard)
title("Jaccard's index (z-scored)")
subplot(2,2,4)
histogram(all_pat_comp_per_sz.Jaccard_z)
xline(0)
title("Jaccard's index (z-scored)")
sgtitle("Per seizure per patient")

figure(2)
subplot(2,2,1)
histogram(all_pat_comp_across_sz.("Percentage of CLO captured"), BinWidth=0.1)
title("Percentage of CLO captured")
subplot(2,2,2)
histogram(all_pat_comp_across_sz.Perc_z)
xline(0)
title("Percentage of CLO captured (z-scored)")
subplot(2,2,3)
histogram(all_pat_comp_across_sz.Jaccard, BinWidth=0.1)
title("Jaccard's index (z-scored)")
subplot(2,2,4)
histogram(all_pat_comp_across_sz.Jaccard_z)
xline(0)
title("Jaccard's index (z-scored)")
sgtitle("Region included in at least one seizure")

figure(3)
subplot(2,2,1)
boxchart(categorical(all_pat_comp_per_sz.Patient_id), all_pat_comp_per_sz.("Percentage of CLO captured"))
hold on 
swarmchart(categorical(all_pat_comp_per_sz.Patient_id), all_pat_comp_per_sz.("Percentage of CLO captured"), 'filled')
hold off
title("Percentage of CLO captured")
subplot(2,2,2)
boxchart(categorical(all_pat_comp_per_sz.Patient_id), all_pat_comp_per_sz.Perc_z)
hold on 
swarmchart(categorical(all_pat_comp_per_sz.Patient_id), all_pat_comp_per_sz.Perc_z, 'filled')
yline(0)
hold off
title("Percentage of CLO captured (z-scored)")
subplot(2,2,3)
boxchart(categorical(all_pat_comp_per_sz.Patient_id), all_pat_comp_per_sz.Jaccard)
hold on 
swarmchart(categorical(all_pat_comp_per_sz.Patient_id), all_pat_comp_per_sz.Jaccard, 'filled')

hold off
title("Jaccard's index ")
subplot(2,2,4)
boxchart(categorical(all_pat_comp_per_sz.Patient_id), all_pat_comp_per_sz.Jaccard_z)
hold on 
swarmchart(categorical(all_pat_comp_per_sz.Patient_id), all_pat_comp_per_sz.Jaccard_z, 'filled')
yline(0)
hold off
title("Jaccard's index (z-scored)")
sgtitle("Per seizure per patient")

figure(4)
subplot(221)
boxchart(categorical(all_pat_comp_across_sz.Op_type), all_pat_comp_across_sz.Perc_z)
hold on
swarmchart(categorical(all_pat_comp_across_sz.Op_type), all_pat_comp_across_sz.Perc_z, "filled")
yline(0)
hold off
title("Perc clo (z) across op types")
subplot(222)
boxchart(categorical(all_pat_comp_across_sz.Outcome), all_pat_comp_across_sz.Perc_z)
hold on
swarmchart(categorical(all_pat_comp_across_sz.Outcome), all_pat_comp_across_sz.Perc_z, "filled")
yline(0)
hold off
title("Perc clo (z) across outcomes")
subplot(223)
boxchart(categorical(all_pat_comp_across_sz.Op_type), all_pat_comp_across_sz.Jaccard_z)
hold on
swarmchart(categorical(all_pat_comp_across_sz.Op_type), all_pat_comp_across_sz.Jaccard_z, "filled")
yline(0)
hold off
title("Jaccard (z) across op types")
subplot(224)
boxchart(categorical(all_pat_comp_across_sz.Outcome), all_pat_comp_across_sz.Jaccard_z)
hold on
swarmchart(categorical(all_pat_comp_across_sz.Outcome), all_pat_comp_across_sz.Jaccard_z, "filled")
yline(0)
hold off
title("Jaccard (z) across outcomes")
sgtitle("Comparisons across groups")