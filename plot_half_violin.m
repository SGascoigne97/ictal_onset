%% Compating automatic onsets to CLO
onset_across_titles = [sprintf("Across subject %s",summ_meas), "Onset across seizures"];

save_loc = 'figures/rapid_prototype/onset_clo/';
file_type = "png";


det_meths = ["imprint", "EI", "PLHG"];
comp_measures = ["Perc", "Perc_z", "Jacc", "Jacc_z", "Coh", "Coh_z"];

for parc = ["chan", "roi_120", "roi_250"]
    for onset_across = [0,1]
        if onset_across == 0
            for summ_meas = ["med", "max"]
                half_violin(comp_meth_with_clo, det_meths,comp_measures, "parc", parc,...
                    "onset_across", onset_across, "summ_meas", summ_meas,...
                    "vis_plot", "on", "save_fig", 0, "save_loc", save_loc,...
                    "file_type", file_type)
            end
        else
            half_violin(comp_meth_with_clo, det_meths,comp_measures, "parc", parc,...
                    "onset_across", onset_across, "vis_plot", "on",...
                    "save_fig", 0, "save_loc", save_loc,  "file_type", file_type)
        end
    end
end

%% Comparing onsets to resection
save_loc = 'figures/rapid_prototype/onset_resec/';
file_type = "png";

det_meths = ["clo", "imprint", "EI"];
comp_measures = ["Perc", "Perc_z"];%, "Coh", "Coh_z"]; %"Jacc", "Jacc_z",

for parc = ["chan", "roi_120", "roi_250"]
    for onset_across = 1 %[0,1]
        if onset_across == 0
            for summ_meas = ["med", "max"]
                half_violin(comp_meth_with_resec, det_meths(2:end), comp_measures, "parc", parc,...
                    "onset_across", onset_across, "summ_meas", summ_meas,...
                    "vis_plot", "on", "save_fig",0, "save_loc", save_loc,...
                    "file_type", file_type)
            end
        else
            half_violin(comp_meth_with_resec, det_meths,comp_measures, "parc", parc,...
                    "onset_across", onset_across, "vis_plot", "on",...
                    "save_fig", 0, "save_loc", save_loc,  "file_type", file_type)
        end
    end
end

%%


%%



% parc = "roi_250";
% onset_across = 1;
% 
% summ_meas = "med";
% onset_across_titles = [sprintf("Across subject %s",summ_meas), "Onset across seizures"];
% if onset_across == 1
%     det_meths = ["clo", "imprint", "EI", "PLHG"];
% else
%     det_meths = ["imprint", "EI", "PLHG"];
% end
% 
% fig = figure("Position", [10, 10, 1800, 900], 'Visible','off');
% tiledlayout(3,2)
% rng(6)
% for comp_measure = ["Perc", "Perc_z", "Jacc", "Jacc_z", "Coh", "Coh_z"]
%     if comp_measure == "Perc" | comp_measure == "Jacc"
%         y_lim = [-0.1,1.1];
%         bin_wd = 0.1;
%     elseif comp_measure == "Coh"
%         y_lim = [-0.6,1.1];
%         bin_wd = 0.1;
%     elseif contains(comp_measure, "z")
%         y_lim = [-3, 12]; % May need to tweak this later
%         bin_wd = 1;
%     end
% 
%     % Compute offset between violin plots as the maximum number of subjects
%     % across methods
%     subjs = [];
%     for det = det_meths
%         subjs = [subjs, size(comp_meth.across_sz.chan.(sprintf("%s", det)),1)];
%     end
%     offset = (max(subjs)+10)/2;
%     max_x = (offset*2) +10;
%     x_lim = [-10, max_x];
%     
%     nexttile
%     hold on
%     for det_ind = 1:3 
%         if onset_across == 1
%             data_tab =  comp_meth.across_sz.(sprintf("%s", parc)).(sprintf("%s", det_meths(det_ind)));
%             data = data_tab.(sprintf("%s", comp_measure));
%             out = data_tab.Outcome;
%         else
%             data_tab = comp_meth.per_sz.(sprintf("%s", parc)).(sprintf("%s", det_meths(det_ind))).(sprintf("%s", summ_meas));
%             data = data_tab.(sprintf("%s", comp_measure));
%             out = data_tab.Outcome;
%         end
%         
%         out = out(~isnan(data));
%         data = data(~isnan(data));
%        
%         g = data(out<3);
%         b = data(out>2);
%         min_val_g = round(min(g)*(1/bin_wd))*bin_wd;
%        
%         min_val_b = round(min(b)*(1/bin_wd))*bin_wd;
%         max_val_g = round(max(g)*(1/bin_wd))*bin_wd;
%         max_val_b = round(max(b)*(1/bin_wd))*bin_wd;
%         if min_val_g > min(g)
%             min_val_g = min_val_g - bin_wd;
%         end 
%         if min_val_b > min(b)
%             min_val_b = min_val_b - bin_wd;
%         end
%         if max_val_g < max(g)
%             max_val_g = max_val_g + bin_wd;
%         end
%         if max_val_b < max(b)
%             max_val_b = max_val_b + bin_wd;
%         end
%        
%         hist_vals_g = histcounts(g, min_val_g:bin_wd:max_val_g);
%         % Add scatter points
%         non_zero_grps = hist_vals_g(hist_vals_g~=0);
%         jitt = [];
%         for grp_sz = non_zero_grps
%             jitt = [jitt; rand(grp_sz,1)*grp_sz*0.8];
%         end
%         scatter((det_ind-1)*offset+jitt, sort(g), [], [0.3010 0.7450 0.9330], 'filled');
%     
%         hist_vals_b = histcounts(b, min_val_b:bin_wd:max_val_b);
%         
%         % Add scatter points
%         non_zero_grps = hist_vals_b(hist_vals_b~=0);
%         jitt = [];
%         for grp_sz = non_zero_grps
%             jitt = [jitt; rand(grp_sz,1)*grp_sz*0.8];
%         end
%         scatter((det_ind-1)*offset-jitt, sort(b), [],[0.8500 0.3250 0.0980], 'filled');
%         if length(hist_vals_g) >1
%             smooth_data_g = smoothdata(interp1(hist_vals_g, 1:0.1:length(hist_vals_g)));
%             fill((det_ind-1)*offset+[0 smooth_data_g 0], ...
%                 rescale((1:length([0 smooth_data_g 0]))/length([0 smooth_data_g 0]),...
%                 min_val_g, max_val_g), [0.3010 0.7450 0.9330], 'FaceAlpha',0.3,...
%                 'LineStyle','none')
%         end 
%         if length(hist_vals_b) >1
%             smooth_data_b = smoothdata(interp1(hist_vals_b, 1:0.1:length(hist_vals_b)));
%             fill((det_ind-1)*offset -[0 smooth_data_b 0],...
%                 rescale((1:length([0 smooth_data_b 0]))/length([0 smooth_data_b 0]),...
%                 min_val_b, max_val_b), [0.8500 0.3250 0.0980], 'FaceAlpha',0.3,...
%                 'LineStyle','none')
%         end
%     
%     end
%     set(gca, "XTick", (0:2)*offset, "XTickLabel", det_meths)
%     title(comp_measure)
%     ylim(y_lim)
%     xlim(x_lim)
% 
%     if contains(comp_measure, "z")
%         patch([x_lim(1), x_lim(2), x_lim(2), x_lim(1)], [-2, -2, 2, 2], [0.7,0.7,0.7], "FaceAlpha", 0.3, "LineStyle", "none")
%     elseif comp_measure == "Coh"
%         yline(0.2)
%     end
% end
% sgtitle(sprintf("%s-wise comparison against resection (%s) ", parc, onset_across_titles(onset_across+1)))
% saveas(fig,sprintf('figures/rapid_prototype/onset_resec/violin_%s_%s.png', parc, strrep(onset_across_titles(onset_across+1), ' ', '_')))
% %saveas(fig,sprintf('figures/rapid_prototype/onset_resec/violin_%s_%s.svg', parc, strrep(onset_across_titles(onset_across+1), ' ', '_')))
% 
% 
% %%
% % onset_across_titles = [sprintf("Across subject %s",summ_meas), "Onset across seizures"];
% % for parc = ["chan", "roi_120", "roi_250"]
% %     for onset_across = [0,1]
% %         for summ_meas = ["med", "max"]
% %             if onset_across == 1
% %                 det_meths = ["clo", "imprint", "EI", "PLHG"];
% %             else
% %                 det_meths = ["imprint", "EI", "PLHG"];
% %             end
% %             
% %             fig = figure("Position", [10, 10, 1800, 900], 'Visible','off');
% %             tiledlayout(3,2)
% %             rng(6)
% %             for comp_measure = ["Perc", "Perc_z", "Jacc", "Jacc_z", "Coh", "Coh_z"]
% %                 if comp_measure == "Perc" | comp_measure == "Jacc"
% %                     y_lim = [-0.1,1.1];
% %                     bin_wd = 0.1;
% %                 elseif comp_measure == "Coh"
% %                     y_lim = [-0.6,1.1];
% %                     bin_wd = 0.1;
% %                 elseif contains(comp_measure, "z")
% %                     y_lim = [-3, 12]; % May need to tweak this later
% %                     bin_wd = 1;
% %                 end
% %             
% %                 % Compute offset between violin plots as the maximum number of subjects
% %                 % across methods
% %                 subjs = [];
% %                 for det = det_meths
% %                     subjs = [subjs, size(comp_meth.across_sz.chan.(sprintf("%s", det)),1)];
% %                 end
% %                 offset = (max(subjs)+10)/2;
% %                 max_x = (offset*2) +10;
% %                 x_lim = [-10, max_x];
% %                 
% %                 nexttile
% %                 hold on
% %                 for det_ind = 1:3 
% %                     if onset_across == 1
% %                         data_tab =  comp_meth.across_sz.(sprintf("%s", parc)).(sprintf("%s", det_meths(det_ind)));
% %                         data = data_tab.(sprintf("%s", comp_measure));
% %                         out = data_tab.Outcome;
% %                     else
% %                         data_tab = comp_meth.per_sz.(sprintf("%s", parc)).(sprintf("%s", det_meths(det_ind))).(sprintf("%s", summ_meas));
% %                         data = data_tab.(sprintf("%s", comp_measure));
% %                         out = data_tab.Outcome;
% %                     end
% %                     
% %                     out = out(~isnan(data));
% %                     data = data(~isnan(data));
% %                    
% %                     g = data(out<3);
% %                     b = data(out>2);
% %                     min_val_g = round(min(g)*(1/bin_wd))*bin_wd;
% %                    
% %                     min_val_b = round(min(b)*(1/bin_wd))*bin_wd;
% %                     max_val_g = round(max(g)*(1/bin_wd))*bin_wd;
% %                     max_val_b = round(max(b)*(1/bin_wd))*bin_wd;
% %                     if min_val_g > min(g)
% %                         min_val_g = min_val_g - bin_wd;
% %                     end 
% %                     if min_val_b > min(b)
% %                         min_val_b = min_val_b - bin_wd;
% %                     end
% %                     if max_val_g < max(g)
% %                         max_val_g = max_val_g + bin_wd;
% %                     end
% %                     if max_val_b < max(b)
% %                         max_val_b = max_val_b + bin_wd;
% %                     end
% %                    
% %                     hist_vals_g = histcounts(g, min_val_g:bin_wd:max_val_g);
% %                     % Add scatter points
% %                     non_zero_grps = hist_vals_g(hist_vals_g~=0);
% %                     jitt = [];
% %                     for grp_sz = non_zero_grps
% %                         jitt = [jitt; rand(grp_sz,1)*grp_sz*0.8];
% %                     end
% %                     scatter((det_ind-1)*offset+jitt, sort(g), [], [0.3010 0.7450 0.9330], 'filled');
% %                 
% %                     hist_vals_b = histcounts(b, min_val_b:bin_wd:max_val_b);
% %                     
% %                     % Add scatter points
% %                     non_zero_grps = hist_vals_b(hist_vals_b~=0);
% %                     jitt = [];
% %                     for grp_sz = non_zero_grps
% %                         jitt = [jitt; rand(grp_sz,1)*grp_sz*0.8];
% %                     end
% %                     scatter((det_ind-1)*offset-jitt, sort(b), [],[0.8500 0.3250 0.0980], 'filled');
% %                     if length(hist_vals_g) >1
% %                         smooth_data_g = smoothdata(interp1(hist_vals_g, 1:0.1:length(hist_vals_g)));
% %                         fill((det_ind-1)*offset+[0 smooth_data_g 0], ...
% %                             rescale((1:length([0 smooth_data_g 0]))/length([0 smooth_data_g 0]),...
% %                             min_val_g, max_val_g), [0.3010 0.7450 0.9330], 'FaceAlpha',0.3,...
% %                             'LineStyle','none')
% %                     end 
% %                     if length(hist_vals_b) >1
% %                         smooth_data_b = smoothdata(interp1(hist_vals_b, 1:0.1:length(hist_vals_b)));
% %                         fill((det_ind-1)*offset -[0 smooth_data_b 0],...
% %                             rescale((1:length([0 smooth_data_b 0]))/length([0 smooth_data_b 0]),...
% %                             min_val_b, max_val_b), [0.8500 0.3250 0.0980], 'FaceAlpha',0.3,...
% %                             'LineStyle','none')
% %                     end
% %                 
% %                 end
% %                 set(gca, "XTick", (0:2)*offset, "XTickLabel", det_meths)
% %                 title(comp_measure)
% %                 ylim(y_lim)
% %                 xlim(x_lim)
% %             
% %                 if contains(comp_measure, "z")
% %                     patch([x_lim(1), x_lim(2), x_lim(2), x_lim(1)], [-2, -2, 2, 2], [0.7,0.7,0.7], "FaceAlpha", 0.3, "LineStyle", "none")
% %                 elseif comp_measure == "Coh"
% %                     yline(0.2)
% %                 end
% %             end
% %             sgtitle(sprintf("%s-wise comparison against resection (%s) ", parc, onset_across_titles(onset_across+1)))
% %             saveas(fig,sprintf('figures/rapid_prototype/onset_resec/violin_%s_%s.png', parc, strrep(onset_across_titles(onset_across+1), ' ', '_')))
% %             %saveas(fig,sprintf('figures/rapid_prototype/onset_resec/violin_%s_%s.svg', parc, strrep(onset_across_titles(onset_across+1), ' ', '_')))
% %         end
% %     end
% % end
% % 
