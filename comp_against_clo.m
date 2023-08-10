%% Compare CLO against Imprint, EI, and PLHG
% Include comparison of outcome groups

final_output_all = load('tables/final_output_all_sz.mat');
final_output_all = final_output_all.final_output;

% Add in patient outcome column
for pat = 1:size(final_output_all, 1)
    pat_onset = final_output_all(pat,:);
    outcome_id = pat_onset.("Outcome year"){:}-pat_onset.("Surgery year") == 1;
    final_output_all.outcome(pat) = pat_onset.Surgery_outcome{:}(outcome_id);
end

% Remove patients with no labelled CLO or no outcome 
final_output_all = final_output_all(final_output_all.outcome ~= 8,:);
final_output_all = final_output_all(cellfun(@sum, final_output_all.clo_chan)>0,:);
final_output = final_output_all(cellfun(@iscell, final_output_all.Segment_ids),:);

comp_outcome = 0;
onset_across_title = ["per_sz", "across_sz"];
comp_outcome_title = ["", "_comp_outcome"];
n_perm = 1000;
onset_acr_thresh = 0.5;

for onset_across = [0, 1]
    for chan_or_roi = ["chan", "roi_120", "roi_250"]
        fprintf("%s\n", chan_or_roi)
        for det_meth = ["imprint", "EI", "PLHG"]
            fprintf("%s\n", det_meth)
            clear comp_with_clo_tab
            comp_with_clo_tab = compute_clo_comp_measures(final_output,...
                "det_meth",det_meth, "chan_or_roi", chan_or_roi,...
                "onset_across", onset_across, "onset_acr_thresh",...
                onset_acr_thresh,"n_perm",n_perm);

            if onset_across == 1
                %comp_with_clo_tab = comp_with_clo_tab(comp_with_clo_tab.Sz_count >4,:);
                comp_meth_with_clo.(sprintf("%s", onset_across_title(onset_across+1))).(sprintf("%s", chan_or_roi)).(sprintf("%s", det_meth)) = comp_with_clo_tab;
            else
                for summ_meth = ["med", "max"]
                    if summ_meth == "med"
                        summ_method = @median;
                    else 
                        summ_method = @max;
                    end
                    summ_tab = [comp_with_clo_tab(:,1:3), array2table(cellfun(summ_method, table2cell(comp_with_clo_tab(:,4:9))))];
                    summ_tab.Properties.VariableNames = comp_with_clo_tab.Properties.VariableNames;
                    if summ_meth == "med"
                        summ_tab = summ_tab(summ_tab.Sz_count >4,:);
                    end
                    comp_meth_with_clo.(sprintf("%s", onset_across_title(onset_across+1))).(sprintf("%s", chan_or_roi)).(sprintf("%s", det_meth)).(sprintf("%s", summ_meth))=...
                        summ_tab;
                end
                comp_meth_with_clo.(sprintf("%s", onset_across_title(onset_across+1))).(sprintf("%s", chan_or_roi)).(sprintf("%s", det_meth)).all = comp_with_clo_tab;
            end
        end
    end
end
%%

% Check plots for threshold of 0.25 against 0.5

% %% Plot comparison measures comparing CLO against Imprint, EI, and PLHG
% % Include comparison of outcome groups
% 
% final_output_all = load('tables/final_output_all_sz.mat');
% final_output_all = final_output_all.final_output;
% 
% % Add in patient outcome column
% for pat = 1:size(final_output_all, 1)
%     pat_onset = final_output_all(pat,:);
%     outcome_id = pat_onset.("Outcome year"){:}-pat_onset.("Surgery year") == 1;
%     final_output_all.outcome(pat) = pat_onset.Surgery_outcome{:}(outcome_id);
% end
% % Remove patients with no labelled CLO or no outcome 
% final_output_all = final_output_all(final_output_all.outcome ~= 8,:);
% final_output_all = final_output_all(cellfun(@sum, final_output_all.clo_chan)>0,:);
% 
% %%
% comp_outcome = 0;
% onset_across_title = ["per_sz", "across_sz"];
% comp_outcome_title = ["", "_comp_outcome"];
% %chan_or_roi = "roi_120"; % ["chan", "roi_120", "roi_250"];
% n_perm = 1000;
% onset_acr_thresh = 0.25;
% summ_meth = "median";
% % Create an empty table to store output
% col_names = {'Patient_id', 'Sz_count', 'Outcome','Jacc', 'Jacc_z', 'Perc', 'Perc_z', ...
%     'Coh', 'Coh_z'};
% 
% final_output = final_output_all(cellfun(@iscell, final_output_all.Segment_ids),:);
% 
% for onset_across = [0, 1]
%     for chan_or_roi = ["chan", "roi_120", "roi_250"]
%         fprintf("%s\n", chan_or_roi)
%         for det_meth = ["imprint", "EI", "PLHG"]
%             fprintf("%s\n", det_meth)
%             clear comp_with_clo_tab
% 
%             if onset_across == 1
%                  comp_with_clo_tab = zeros(size(final_output,1), length(col_names));
%             else
%                 
%                 comp_with_clo_tab = cell(size(final_output,1), length(col_names));
%             end
%             comp_with_clo_tab = array2table(comp_with_clo_tab, 'VariableNames',...
%             col_names);
%         
%             comp_with_clo_tab.Patient_id = final_output.Patient_id;
%             comp_with_clo_tab.Sz_count = cellfun(@length, final_output.Segment_ids);
%             comp_with_clo_tab.Sz_count(~cellfun(@iscell, final_output.Segment_ids)) = 1;
%             comp_with_clo_tab.Outcome = final_output.outcome;
%         
%             for pat = 1:size(final_output,1)
%                 pat_onset = final_output(pat,:);
%                 % Pull out onset and CLO regions
%                 onset = pat_onset.(sprintf("%s_%s",det_meth, chan_or_roi)){:};
%                 clo = pat_onset.(sprintf("clo_%s", chan_or_roi)){:};
%                          
%                 if sum(clo) == 0
%                     fprintf("%s No CLO included \n", pat_onset.Patient_id{:})
%                     continue
%                 end
%     
%                 if onset_across == 1
%                     onset_acr = sum(onset,2)/size(onset,2) >= onset_acr_thresh;
%     
%                     if sum(onset_across) == 0
%                         fprintf("%s No onset across found \n", pat_onset.Patient_id{:})
%                         continue
%                     end
%                     jacc = jaccard(logical(clo), logical(onset_acr));
%                     perc = sum(clo + onset_acr == 2)/sum(clo);
%                     coh = cohensKappa(logical(clo),logical(onset_acr));
%     
%                     % Compute comparison measures
%                     jacc_perm = zeros(1,n_perm);
%                     perc_perm = jacc_perm; 
%                     coh_perm = jacc_perm;
%                     for perm = 1:n_perm
%                         rng(perm)
%                         perm_onset = onset_acr(randperm(length(onset_acr)));
%                         rng(perm+1)
%                         perm_clo = clo(randperm(length(clo)));
%                         jacc_perm(perm) = jaccard(double(perm_onset),double(perm_clo));
%                         perc_perm(perm) = sum(perm_clo + perm_onset == 2)/sum(perm_clo);
%                         coh_perm(perm) = cohensKappa(logical(perm_onset),logical(perm_clo));
%                     end
%     
%                     jacc_z = (jacc - mean(jacc_perm))/std(jacc_perm);
%                     perc_z = (perc - mean(perc_perm))/std(perc_perm);
%                     coh_z = (coh- mean(coh_perm))/std(coh_perm);
% 
%                     comp_with_clo_tab(pat, 4:9) = table(jacc, jacc_z, perc,...
%                     perc_z, coh, coh_z);
%                 else
%                     jacc = nan(size(onset,2),1);
%                     perc = jacc;
%                     coh = jacc; 
%                     jacc_z = jacc;
%                     perc_z = jacc;
%                     coh_z = jacc; 
%                     for sz = 1:size(onset,2)
%                         sz_onset = onset(:,sz);
%                         if sum(sz_onset) == 0
%                             continue
%                         end
%                         jacc(sz) = jaccard(logical(clo), logical(sz_onset));
%                         perc(sz) = sum(clo + sz_onset == 2)/sum(clo);
%                         coh(sz) = cohensKappa(logical(clo),logical(sz_onset));
%         
%                         % Compute comparison measures
%                         jacc_perm = zeros(1,n_perm);
%                         perc_perm = jacc_perm; 
%                         coh_perm = jacc_perm;
%                         for perm = 1:n_perm
%                             rng(perm)
%                             perm_onset = sz_onset(randperm(length(sz_onset)));
%                             rng(perm+1)
%                             perm_clo = clo(randperm(length(clo)));
%                             jacc_perm(perm) = jaccard(double(perm_onset),double(perm_clo));
%                             perc_perm(perm) = sum(perm_clo + perm_onset == 2)/sum(perm_clo);
%                             coh_perm(perm) = cohensKappa(logical(perm_onset),logical(perm_clo));
%                         end
%         
%                         jacc_z(sz) = (jacc(sz) - mean(jacc_perm))/std(jacc_perm);
%                         perc_z(sz) = (perc(sz) - mean(perc_perm))/std(perc_perm);
%                         coh_z(sz) = (coh(sz)- mean(coh_perm))/std(coh_perm);
%                     end
%                     comp_with_clo_tab(pat, 4:9) = [{jacc}, {jacc_z}, {perc},...
%                     {perc_z}, {coh}, {coh_z}];
%                
%                 end
%                    
%             end
% 
%             if onset_across == 1
%                 comp_meth_with_clo.(sprintf("%s", onset_across_title(onset_across+1))).(sprintf("%s", chan_or_roi)).(sprintf("%s", det_meth)) = comp_with_clo_tab;
%            
%             else
%                 med_tab = [comp_with_clo_tab(:,1:3), array2table(cellfun(@median, table2cell(comp_with_clo_tab(:,4:9))))];
%                 med_tab.Properties.VariableNames = comp_with_clo_tab.Properties.VariableNames;
%                 med_tab = med_tab(med_tab.Sz_count >4,:);
%             
%                 max_tab = [comp_with_clo_tab(:,1:3), array2table(cellfun(@max, table2cell(comp_with_clo_tab(:,4:9))))];
%                 max_tab.Properties.VariableNames = comp_with_clo_tab.Properties.VariableNames;
%              
%                 comp_meth_with_clo.(sprintf("%s", onset_across_title(onset_across+1))).(sprintf("%s", chan_or_roi)).(sprintf("%s", det_meth)).med = med_tab;
%                 comp_meth_with_clo.(sprintf("%s", onset_across_title(onset_across+1))).(sprintf("%s", chan_or_roi)).(sprintf("%s", det_meth)).max = max_tab;
%                 comp_meth_with_clo.(sprintf("%s", onset_across_title(onset_across+1))).(sprintf("%s", chan_or_roi)).(sprintf("%s", det_meth)).all = comp_with_clo_tab;
%             end
%         end
%     end
% end
