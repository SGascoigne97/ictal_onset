
%% Plot comparison measures between CLO, Imprint, EI, PLHG and resection
% Include comparison of outcome groups

final_output_all = load('tables/final_output_all_sz.mat');
final_output_all = final_output_all.final_output;

% Remove patients with fewer than 5 seizures and those with no surgical
% outcome listed
for pat = 1:size(final_output_all, 1)
    pat_onset = final_output_all(pat,:);
    outcome_id = pat_onset.("Outcome year"){:}-pat_onset.("Surgery year") == 1;
    final_output_all.outcome(pat) = pat_onset.Surgery_outcome{:}(outcome_id);
end
%final_output_all = final_output_all(final_output_all.outcome ~= 8,:);
final_output_all = final_output_all(cellfun(@sum, final_output_all.resected_chan)>0,:);

%
comp_outcome = 0;
onset_across_title = ["per_sz", "across_sz"];
comp_outcome_title = ["", "_comp_outcome"];
n_perm = 1000;
onset_acr_thresh = 0.5;
% Create an empty table to store output
%col_names = {'Patient_id', 'Sz_count', 'Outcome','Jacc', 'Jacc_z', 'Perc', 'Perc_z', ...
    %'Coh', 'Coh_z'};
for onset_across = [0, 1]
    fprintf("%s\n", onset_across_title(onset_across+1))
    for chan_or_roi = ["chan", "roi_120", "roi_250"]
        fprintf("%s\n", chan_or_roi)
        for det_meth = ["clo", "imprint", "EI", "PLHG"]
            fprintf("%s\n", det_meth)
            clear comp_meth_tab

            if onset_across == 0
                final_output = final_output_all(cellfun(@iscell, final_output_all.Segment_ids),:);
                final_output = final_output(cellfun(@length, final_output.Segment_ids)>4,:);
            else
                final_output = final_output_all;
            end
%             if det_meth ~= "clo" 
%                 final_output = final_output_all(cellfun(@iscell, final_output_all.Segment_ids),:);
%                 final_output = final_output(cellfun(@length, final_output.Segment_ids)>4,:);
%             else
%                 final_output = final_output_all;
%             end
        
            comp_meth_tab = compute_resec_comp_measures(final_output,...
                "chan_or_roi", chan_or_roi, "det_meth", det_meth, "n_perm",n_perm,...
                "onset_acr_thresh",onset_acr_thresh,"onset_across",onset_across);

            if onset_across == 1
                %comp_meth_tab = comp_meth_tab(comp_meth_tab.Sz_count >4,:);
                comp_meth_with_resec.(sprintf("%s", onset_across_title(onset_across+1))).(sprintf("%s", chan_or_roi)).(sprintf("%s", det_meth)) = comp_meth_tab;
            else
                if det_meth ~= "clo"
                    for summ_meth = ["med", "max"]
                        if summ_meth == "med"
                            summ_method = @median;
                        else 
                            summ_method = @max;
                        end
                        summ_tab = [comp_meth_tab(:,1:3), array2table(cellfun(summ_method, table2cell(comp_meth_tab(:,4:end))))];
                        summ_tab.Properties.VariableNames = comp_meth_tab.Properties.VariableNames;
                        if summ_meth == "med"
                            summ_tab = summ_tab(summ_tab.Sz_count >4,:);
                        end
                        comp_meth_with_resec.(sprintf("%s", onset_across_title(onset_across+1))).(sprintf("%s", chan_or_roi)).(sprintf("%s", det_meth)).(sprintf("%s", summ_meth))=...
                            summ_tab;
                    end
                    comp_meth_with_resec.(sprintf("%s", onset_across_title(onset_across+1))).(sprintf("%s", chan_or_roi)).(sprintf("%s", det_meth)).all = comp_meth_tab;
                end
            end
        end
    end
end




%%
        
%         %% Create plots 
%         fig = figure('Position',[10,10,1800,900]);
%         tiledlayout(3,2)
%         plot_cols = ["Perc", "Perc_z", "Jacc", "Jacc_z", ...
%             "Coh", "Coh_z"];
%         for col = plot_cols
%             nexttile
%             if contains(col,"z")
%                 bin_width = 0.5;
%                 x_lim = [-3,9];
%             elseif col == "Coh"
%                  bin_width = 0.1;
%                  x_lim = [-.5, 1.1];
%             else
%                 bin_width = 0.1;
%                 x_lim = [-0.1,1.1];
%             end
%             plot_data = comp_meth.(sprintf("%s", chan_or_roi)).(sprintf("%s", det_meth)).(sprintf("%s", col));
%             if comp_outcome == 1
%                 h_good = histogram(plot_data(comp_meth.(sprintf("%s", chan_or_roi)).(sprintf("%s", det_meth)).Outcome <=2),'FaceColor', 'b', 'BinWidth',bin_width);
%                 max_good = max(h_good.Values);
%                 hold on
%                 h_bad = histogram(plot_data(comp_meth.(sprintf("%s", chan_or_roi)).(sprintf("%s", det_meth)).Outcome >2),'FaceColor', 'r', 'BinWidth',bin_width, 'FaceAlpha',0.5);
%                 max_bad = max(h_bad.Values);
%                 max_y = max(max_good, max_bad) + 2;
%                 legend({"ILAE 1-2", "ILAE 3+"}, 'Location', 'northwest')
%                 title(col)
%             else
%                 skw = skewness(plot_data);
%                 hist = histogram(plot_data,'FaceColor', 'b', 'BinWidth',bin_width);
%                 max_y = max(hist.Values);
%                 title(sprintf("%s (skewness: %.2f)", col, skw))
%             end
%             if contains(col,"z")
%                 patch([-2,-2,2,2],[0, max_y, max_y, 0],[.7,.7,.7], 'FaceAlpha', 0.3, 'LineStyle', 'none')
%             end
%             hold off
%             xlim(x_lim)
%             ylim([0 max_y])
%             
%             
%         end
%         sgtitle(sprintf("%s %s (n=%d)", det_meth, chan_or_roi, length(plot_data)))
%         %saveas(fig, sprintf('figures/rapid_prototype/onset_resec/%s_hist_%s.svg',chan_or_roi, det_meth))
%           
% %         fig_2 = figure('Position',[10,10,1800,900]);
% %         tiledlayout(3,2)
% %         for col = plot_cols
% %             nexttile
% %             if contains(col,"z")
% %                 bin_width = 0.5;
% %                 y_lim = [-3,9];
% %             elseif col == "Coh"
% %                  bin_width = 0.1;
% %                  y_lim = [-.5, 1.1];
% %             else
% %                 bin_width = 0.1;
% %                 y_lim = [-0.1,1.1];
% %             end
% %             plot_data = comp_meth.(sprintf("%s", chan_or_roi)).(sprintf("%s", det_meth)).(sprintf("%s", col));
% %             boxchart(double(comp_meth.(sprintf("%s", chan_or_roi)).(sprintf("%s", det_meth)).Outcome >2), plot_data, "MarkerStyle","none")
% %             hold on
% %             swarmchart(double(comp_meth.(sprintf("%s", chan_or_roi)).(sprintf("%s", det_meth)).Outcome >2), plot_data, "filled")
% %             set(gca,'XTick', [0,1], 'XTickLabel', {"ILAE 1-2","ILAE 3+"})
% %             if contains(col,"z")
% %                 patch([-0.5,-0.5,1.5,1.5],[2,-2,-2,2],[0.7,0.7,0.7], 'FaceAlpha', 0.3, 'LineStyle', 'none')
% %             end
% %             hold off
% %             ylim(y_lim)
% %             title(col)
% %         end
% %         sgtitle(sprintf("%s %s (n=%d)", det_meth, chan_or_roi, length(plot_data)))
%         %saveas(fig_2, sprintf('figures/rapid_prototype/onset_resec/%s_box_%s.svg',chan_or_roi, det_meth))


%% Plots outside of main code
for chan_or_roi = ["chan", "roi_120", "roi_250"]
    for det_meth = ["clo", "imprint"]
        fig = figure('Position',[10,10,1800,900]);
        tiledlayout(3,2)
        plot_cols = ["Perc", "Perc_z", "Jacc", "Jacc_z", ...
            "Coh", "Coh_z"];
        for col = plot_cols
            nexttile
            if contains(col,"z")
                bin_width = 1;
                x_lim = [-3,9];
            elseif col == "Coh"
                 bin_width = 0.1;
                 x_lim = [-.5, 1.1];
            else
                bin_width = 0.1;
                x_lim = [-0.1,1.1];
            end
            plot_data = comp_meth_with_resec.(sprintf("%s", chan_or_roi)).(sprintf("%s", det_meth)).(sprintf("%s", col));
            if comp_outcome == 1
                h_good = histogram(plot_data(comp_meth_with_resec.(sprintf("%s", chan_or_roi)).(sprintf("%s", det_meth)).Outcome <=2),'FaceColor', 'b', 'BinWidth',bin_width);
                max_good = max(h_good.Values);
                hold on
                h_bad = histogram(plot_data(comp_meth_with_resec.(sprintf("%s", chan_or_roi)).(sprintf("%s", det_meth)).Outcome >2),'FaceColor', 'r', 'BinWidth',bin_width, 'FaceAlpha',0.5);
                max_bad = max(h_bad.Values);
                max_y = max(max_good, max_bad) + 2;
                legend({"ILAE 1-2", "ILAE 3+"}, 'Location', 'northwest')
                title(col)
            else
                skw = skewness(plot_data);
                hist = histogram(plot_data,'FaceColor', 'b', 'BinWidth',bin_width);
                max_y = max(hist.Values) + 2;
                title(sprintf("%s (skewness: %.2f)", col, skw))
            end
            if contains(col,"z")
                patch([-2,-2,2,2],[0, max_y, max_y, 0],[.7,.7,.7], 'FaceAlpha', 0.3, 'LineStyle', 'none')
                text(0, max_y-1, 'chance', 'Color',[0.2,0.2,0.2], 'HorizontalAlignment','center')
            end
            hold off
            xlim(x_lim)
            ylim([0 max_y])
            
            
        end
        sgtitle(sprintf("%s %s (n=%d)", det_meth, chan_or_roi, length(plot_data)))
        saveas(fig, sprintf('figures/rapid_prototype/onset_resec/%s_hist_%s.png',chan_or_roi, det_meth))
    end
end


% %% Plot comparison measures between Imprint, EI, PLHG and CLO
% % Include comparison of outcome groups
% 
% addpath(genpath('sarah_functions'))
% final_output_all = load('tables/final_output_all_sz.mat');
% final_output_all = final_output_all.final_output;
% 
% % Remove patients with fewer than 5 seizures and those with no surgical
% % outcome listed
% for pat = 1:size(final_output_all, 1)
%     pat_onset = final_output_all(pat,:);
%     outcome_id = pat_onset.("Outcome year"){:}-pat_onset.("Surgery year") == 1;
%     final_output_all.outcome(pat) = pat_onset.Surgery_outcome{:}(outcome_id);
% end
% final_output_all = final_output_all(final_output_all.outcome ~= 8,:);
% final_output_all = final_output_all(cellfun(@sum, final_output_all.resected_chan)>0,:);
% 
% %%
% 
% %chan_or_roi = "roi_120"; % ["chan", "roi_120", "roi_250"];
% n_perm = 1000;
% 
% % Create an empty table to store output
% col_names = {'Patient_id', 'Sz_count', 'Outcome','Jacc', 'Jacc_z', 'Perc', 'Perc_z', ...
%     'Coh', 'Coh_z'};
% 
% for chan_or_roi = ["chan", "roi_120", "roi_250"]
%     fprintf("%s\n", chan_or_roi)
%     for det_meth = ["clo"]%, "imprint", "EI", "PLHG"]
%         fprintf("%s\n", det_meth)
%         clear comp_meth_tab
%         if det_meth ~= "clo" 
%             final_output = final_output_all(cellfun(@iscell, final_output_all.Segment_ids),:);
%             final_output = final_output(cellfun(@length, final_output.Segment_ids)>4,:);
%         else
%             final_output = final_output_all;
%         end
%     
%         comp_meth_tab = zeros(size(final_output,1), length(col_names));
%         comp_meth_tab = array2table(comp_meth_tab, 'VariableNames',...
%         col_names);
%     
%         comp_meth_tab.Patient_id = final_output.Patient_id;
%         comp_meth_tab.Sz_count = cellfun(@length, final_output.Segment_ids);
%         comp_meth_tab.Outcome = final_output.outcome;
%     
%         for pat = 1:size(final_output,1)
%         
%             pat_onset = final_output(pat,:);
%             % Pull out onset and resected regions (Lausanne 120)
%             onset = pat_onset.(sprintf("%s_%s",det_meth, chan_or_roi)){:};
%             resec = pat_onset.(sprintf("resected_%s", chan_or_roi)){:};
%             
%             % Create one onset per individual (regions included in at least half of
%             % their seizures)
%             if det_meth ~= "clo"
%                 onset_across = (sum(onset,2)/size(onset,2)) >= 0.5; 
%             else
%                 onset_across = onset;
%             end
%             
%             if sum(resec) == 0
%                 fprintf("%s No resection included \n", pat_onset.Patient_id{:})
%                 continue
%             elseif sum(onset_across) == 0
%                 fprintf("%s No onset found \n", pat_onset.Patient_id{:})
%                 continue
%             end
%         
%             comp_meth_tab.Jacc(pat) = jaccard(double(onset_across),double(resec));
%             comp_meth_tab.Perc(pat) = sum(resec + onset_across == 2)/sum(onset_across);
%             comp_meth_tab.Coh(pat) = cohensKappa(logical(onset_across),logical(resec));
%             
%             % Compute comparison measures
%             jacc_perm = zeros(1,n_perm);
%             perc_perm = jacc_perm; 
%             coh_perm = jacc_perm;
%             for perm = 1:n_perm
%                 rng(perm)
%                 perm_onset = onset_across(randperm(length(onset_across)));
%                 rng(perm+1)
%                 perm_resec = resec(randperm(length(onset_across)));
%                 jacc_perm(perm) = jaccard(double(perm_onset),double(perm_resec));
%                 perc_perm(perm) = sum(perm_resec + perm_onset == 2)/sum(perm_onset);
%                 coh_perm(perm) = cohensKappa(logical(perm_onset),logical(perm_resec));
%             end
%         
%               comp_meth_tab.Jacc_z(pat) = (comp_meth_tab.Jacc(pat) - mean(jacc_perm))/std(jacc_perm);
%               comp_meth_tab.Perc_z(pat) = (comp_meth_tab.Perc(pat) - mean(perc_perm))/std(perc_perm);
%               comp_meth_tab.Coh_z(pat) = (comp_meth_tab.Coh(pat) - mean(coh_perm))/std(coh_perm);
%         end
%     
%         comp_meth_tab = comp_meth_tab(comp_meth_tab.Jacc + comp_meth_tab.Perc_z ~=0,:);
%         
%         comp_meth.(sprintf("%s", chan_or_roi)).(sprintf("%s", det_meth)) = comp_meth_tab;
%         
%         %% Create plots 
%         fig = figure('Position',[10,10,1800,900]);
%         tiledlayout(3,2)
%         plot_cols = ["Perc", "Perc_z", "Jacc", "Jacc_z", ...
%             "Coh", "Coh_z"];
%         for col = plot_cols
%             nexttile
%             if contains(col,"z")
%                 bin_width = 0.5;
%                 x_lim = [-3,9];
%             elseif col == "Coh"
%                  bin_width = 0.1;
%                  x_lim = [-.5, 1.1];
%             else
%                 bin_width = 0.1;
%                 x_lim = [-0.1,1.1];
%             end
%             plot_data = comp_meth.(sprintf("%s", chan_or_roi)).(sprintf("%s", det_meth)).(sprintf("%s", col));
%             if comp_outcome == 1
%                 h_good = histogram(plot_data(comp_meth.(sprintf("%s", chan_or_roi)).(sprintf("%s", det_meth)).Outcome <=2),'FaceColor', 'b', 'BinWidth',bin_width);
%                 max_good = max(h_good.Values);
%                 hold on
%                 h_bad = histogram(plot_data(comp_meth.(sprintf("%s", chan_or_roi)).(sprintf("%s", det_meth)).Outcome >2),'FaceColor', 'r', 'BinWidth',bin_width, 'FaceAlpha',0.5);
%                 max_bad = max(h_bad.Values);
%                 max_y = max(max_good, max_bad) + 2;
%                 legend({"ILAE 1-2", "ILAE 3+"}, 'Location', 'northwest')
%                 title(col)
%             else
%                 skw = skewness(plot_data);
%                 hist = histogram(plot_data,'FaceColor', 'b', 'BinWidth',bin_width);
%                 max_y = max(hist.Values);
%                 title(sprintf("%s (skewness: %.2f)", col, skw))
%             end
%             if contains(col,"z")
%                 patch([-2,-2,2,2],[0, max_y, max_y, 0],[.7,.7,.7], 'FaceAlpha', 0.3, 'LineStyle', 'none')
%             end
%             hold off
%             xlim(x_lim)
%             ylim([0 max_y])
%             
%             
%         end
%         sgtitle(sprintf("%s %s (n=%d)", det_meth, chan_or_roi, length(plot_data)))
%         %saveas(fig, sprintf('figures/rapid_prototype/onset_resec/%s_hist_%s.svg',chan_or_roi, det_meth))
%           
% %         fig_2 = figure('Position',[10,10,1800,900]);
% %         tiledlayout(3,2)
% %         for col = plot_cols
% %             nexttile
% %             if contains(col,"z")
% %                 bin_width = 0.5;
% %                 y_lim = [-3,9];
% %             elseif col == "Coh"
% %                  bin_width = 0.1;
% %                  y_lim = [-.5, 1.1];
% %             else
% %                 bin_width = 0.1;
% %                 y_lim = [-0.1,1.1];
% %             end
% %             plot_data = comp_meth.(sprintf("%s", chan_or_roi)).(sprintf("%s", det_meth)).(sprintf("%s", col));
% %             boxchart(double(comp_meth.(sprintf("%s", chan_or_roi)).(sprintf("%s", det_meth)).Outcome >2), plot_data, "MarkerStyle","none")
% %             hold on
% %             swarmchart(double(comp_meth.(sprintf("%s", chan_or_roi)).(sprintf("%s", det_meth)).Outcome >2), plot_data, "filled")
% %             set(gca,'XTick', [0,1], 'XTickLabel', {"ILAE 1-2","ILAE 3+"})
% %             if contains(col,"z")
% %                 patch([-0.5,-0.5,1.5,1.5],[2,-2,-2,2],[0.7,0.7,0.7], 'FaceAlpha', 0.3, 'LineStyle', 'none')
% %             end
% %             hold off
% %             ylim(y_lim)
% %             title(col)
% %         end
% %         sgtitle(sprintf("%s %s (n=%d)", det_meth, chan_or_roi, length(plot_data)))
%         %saveas(fig_2, sprintf('figures/rapid_prototype/onset_resec/%s_box_%s.svg',chan_or_roi, det_meth))
%     end
% end
% 
% 
% %% Plots outside of main code
% for chan_or_roi = ["chan", "roi_120", "roi_250"]
%     for det_meth = ["clo", "imprint"]
%         fig = figure('Position',[10,10,1800,900]);
%         tiledlayout(3,2)
%         plot_cols = ["Perc", "Perc_z", "Jacc", "Jacc_z", ...
%             "Coh", "Coh_z"];
%         for col = plot_cols
%             nexttile
%             if contains(col,"z")
%                 bin_width = 1;
%                 x_lim = [-3,9];
%             elseif col == "Coh"
%                  bin_width = 0.1;
%                  x_lim = [-.5, 1.1];
%             else
%                 bin_width = 0.1;
%                 x_lim = [-0.1,1.1];
%             end
%             plot_data = comp_meth.(sprintf("%s", chan_or_roi)).(sprintf("%s", det_meth)).(sprintf("%s", col));
%             if comp_outcome == 1
%                 h_good = histogram(plot_data(comp_meth.(sprintf("%s", chan_or_roi)).(sprintf("%s", det_meth)).Outcome <=2),'FaceColor', 'b', 'BinWidth',bin_width);
%                 max_good = max(h_good.Values);
%                 hold on
%                 h_bad = histogram(plot_data(comp_meth.(sprintf("%s", chan_or_roi)).(sprintf("%s", det_meth)).Outcome >2),'FaceColor', 'r', 'BinWidth',bin_width, 'FaceAlpha',0.5);
%                 max_bad = max(h_bad.Values);
%                 max_y = max(max_good, max_bad) + 2;
%                 legend({"ILAE 1-2", "ILAE 3+"}, 'Location', 'northwest')
%                 title(col)
%             else
%                 skw = skewness(plot_data);
%                 hist = histogram(plot_data,'FaceColor', 'b', 'BinWidth',bin_width);
%                 max_y = max(hist.Values) + 2;
%                 title(sprintf("%s (skewness: %.2f)", col, skw))
%             end
%             if contains(col,"z")
%                 patch([-2,-2,2,2],[0, max_y, max_y, 0],[.7,.7,.7], 'FaceAlpha', 0.3, 'LineStyle', 'none')
%                 text(0, max_y-1, 'chance', 'Color',[0.2,0.2,0.2], 'HorizontalAlignment','center')
%             end
%             hold off
%             xlim(x_lim)
%             ylim([0 max_y])
%             
%             
%         end
%         sgtitle(sprintf("%s %s (n=%d)", det_meth, chan_or_roi, length(plot_data)))
%         saveas(fig, sprintf('figures/rapid_prototype/onset_resec/%s_hist_%s.png',chan_or_roi, det_meth))
%     end
% end
% 
% 
