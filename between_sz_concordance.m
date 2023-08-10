addpath(genpath('sarah_functions/'))
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


%chan_or_roi = "roi_120"; % ["chan", "roi_120", "roi_250"];
n_perm = 1000;
% Create an empty table to store output
col_names = {'Patient_id', 'Sz_count', 'Outcome','Jacc', 'Jacc_z', ...
    'Coh', 'Coh_z'};

final_output = final_output_all(cellfun(@iscell, final_output_all.Segment_ids),:);
final_output = final_output(cellfun(@length, final_output.Segment_ids)>4,:);

for chan_or_roi =["chan","roi_120", "roi_250"]
    fprintf("%s\n", chan_or_roi)
    for det_meth = ["imprint", "EI", "PLHG"]
        fprintf("%s\n", det_meth)
        clear within_method_concord_tab
    
        within_method_concord_tab = cell(size(final_output,1), length(col_names));
        within_method_concord_tab = array2table(within_method_concord_tab, 'VariableNames',...
        col_names);
    
        within_method_concord_tab.Patient_id = final_output.Patient_id;
        within_method_concord_tab.Sz_count = cellfun(@length, final_output.Segment_ids);
        within_method_concord_tab.Outcome = final_output.outcome;
    
        for pat = 1:size(final_output,1)
            pat_onset = final_output(pat,:);
            % Pull out onset and CLO regions
            onset = pat_onset.(sprintf("%s_%s",det_meth, chan_or_roi)){:};

            jacc = nan(size(onset,2),size(onset,2));
            coh = jacc; 
            jacc_z = jacc;
            coh_z = jacc; 

            for sz1 = 1:size(onset,2)
                sz_onset_1 = onset(:,sz1);
                if sum(sz_onset_1) == 0
                    continue
                end

                for sz2 = 1:size(onset,2)
                    sz_onset_2 = onset(:,sz2);
                    if sum(sz_onset_2) == 0
                        continue
                    end

                    if sz1 < sz2
                        jacc(sz1,sz2) = jaccard(sz_onset_1, sz_onset_2);
                        coh(sz1,sz2) = cohensKappa(logical(sz_onset_1), logical(sz_onset_2));
            
                        % Compute comparison measures
                        jacc_perm = zeros(1,n_perm);
                        coh_perm = jacc_perm;
                        for perm = 1:n_perm
                            rng(perm)
                            perm_onset_1 = sz_onset_1(randperm(length(sz_onset_1)));
                            rng(perm+1)
                            perm_onset_2 = sz_onset_2(randperm(length(sz_onset_2)));
                            jacc_perm(perm) = jaccard(double(perm_onset_1),double(perm_onset_2));
                            coh_perm(perm) = cohensKappa(logical(perm_onset_1),logical(perm_onset_2));
                        end
        
                        jacc_z(sz1,sz2) = (jacc(sz1,sz2) - mean(jacc_perm))/std(jacc_perm);
                        coh_z(sz1,sz2) = (coh(sz1,sz2)- mean(coh_perm))/std(coh_perm);
                    end
                end
            end
            within_method_concord_tab(pat, 4:7) = [{jacc(:)}, {jacc_z(:)}, ...
                {coh(:)}, {coh_z(:)}];
           
        end
%%
        med_tab = [within_method_concord_tab(:,1:3), array2table(cellfun(@nanmedian, table2cell(within_method_concord_tab(:,4:7))))];
        med_tab.Properties.VariableNames = within_method_concord_tab.Properties.VariableNames;
    
        max_tab = [within_method_concord_tab(:,1:3), array2table(cellfun(@max, table2cell(within_method_concord_tab(:,4:7)),'UniformOutput', false))];
        max_tab.Properties.VariableNames = within_method_concord_tab.Properties.VariableNames;
     
        comp_meth_conc.(sprintf("%s", chan_or_roi)).(sprintf("%s", det_meth)).med = med_tab;
        comp_meth_conc.(sprintf("%s", chan_or_roi)).(sprintf("%s", det_meth)).max = max_tab;
        comp_meth_conc.(sprintf("%s", chan_or_roi)).(sprintf("%s", det_meth)).all = within_method_concord_tab;

        %% Create plots 
        for summ_method = ["med", "max"]
            fig = figure('Position',[10,10,1800,900]);
            tiledlayout(2,2)
            plot_cols = [ "Jacc", "Jacc_z", ...
                "Coh", "Coh_z"];
            for col = plot_cols
                nexttile
                if contains(col,"z")
                    bin_width = 0.5;
                    x_lim = [-3,20];
                elseif col == "Coh"
                     bin_width = 0.1;
                     x_lim = [-.5, 1.1];
                else
                    bin_width = 0.1;
                    x_lim = [-0.1,1.1];
                end

                plot_tab = comp_meth_conc.(sprintf("%s", chan_or_roi)).(sprintf("%s", det_meth)).(sprintf("%s", summ_method));

                if summ_method == "med"
                    plot_tab = plot_tab(~isnan(plot_tab.Jacc),:);
                    plot_data = plot_tab.(sprintf("%s", col));
                elseif summ_method == "max"
                    plot_tab = plot_tab(~cellfun(@isempty, plot_tab.Jacc),:);
                    plot_data = cell2mat(plot_tab.(sprintf("%s", col)));
                end

                h_good = histogram(plot_data(plot_tab.Outcome <=2),'FaceColor', 'b', 'BinWidth',bin_width);
                max_good = max(h_good.Values);
                hold on
                h_bad = histogram(plot_data(plot_tab.Outcome >2),'FaceColor', 'r', 'BinWidth',bin_width, 'FaceAlpha',0.5);
                max_bad = max(h_bad.Values);
                max_y = max(max_good, max_bad) + 2;
                if contains(col,"z")
                    patch([-2,-2,2,2],[0, max_y, max_y, 0],[.7,.7,.7], 'FaceAlpha', 0.3, 'LineStyle', 'none')
                end
                hold off
                xlim(x_lim)
                ylim([0 max_y])
                title(col)
                legend({"ILAE 1-2", "ILAE 3+"}, 'Location', 'northwest')
            end
            sgtitle(sprintf("within-patient %s concordance %s %s (n=%d)", summ_method, det_meth, chan_or_roi, length(plot_data)))
            saveas(fig, sprintf('figures/rapid_prototype/onset_concord/%s_hist_%s_%s.svg',chan_or_roi, det_meth, summ_method))
%           
            fig_2 = figure('Position',[10,10,1800,900]);
            tiledlayout(2,2)
            for col = plot_cols
                nexttile
                if contains(col,"z")
                    bin_width = 0.5;
                    y_lim = [-3,20];
                elseif col == "Coh"
                     bin_width = 0.1;
                     y_lim = [-.5, 1.1];
                else
                    bin_width = 0.1;
                    y_lim = [-0.1,1.1];
                end
        
                plot_tab = comp_meth_conc.(sprintf("%s", chan_or_roi)).(sprintf("%s", det_meth)).(sprintf("%s", summ_method));
        
                if summ_method == "med"
                    plot_tab = plot_tab(~isnan(plot_tab.Jacc),:);
                    plot_data = plot_tab.(sprintf("%s", col));
                elseif summ_method == "max"
                    plot_tab = plot_tab(~cellfun(@isempty, plot_tab.Jacc),:);
                    plot_data = cell2mat(plot_tab.(sprintf("%s", col)));
                end
    
                boxchart(double(plot_tab.Outcome >2), plot_data, "MarkerStyle","none")
                hold on
                swarmchart(double(plot_tab.Outcome >2), plot_data, "filled")
                set(gca,'XTick', [0,1], 'XTickLabel', {"ILAE 1-2","ILAE 3+"})
                if contains(col,"z")
                    patch([-0.5,-0.5,1.5,1.5],[2,-2,-2,2],[0.7,0.7,0.7], 'FaceAlpha', 0.3, 'LineStyle', 'none')
                end
                hold off
                ylim(y_lim)
                title(col)
            end
            sgtitle(sprintf("within-patient %s concordance %s %s (n=%d)", summ_method, det_meth, chan_or_roi, length(plot_data)))
            saveas(fig_2, sprintf('figures/rapid_prototype/onset_concord/%s_box_%s_%s.svg',chan_or_roi, det_meth, summ_method))
        end
        
    end
end




%%
% Violin plots of between-seizure concordance for imprint onsets
for chan_or_roi = ["chan", "roi_120"]
    all_concord_tab = comp_meth_conc.(sprintf("%s",chan_or_roi)).imprint.all;
    all_concord_tab = all_concord_tab(~cellfun(@isempty, all_concord_tab.Jacc),:);
    
    coh = cat(1, all_concord_tab.Coh{:});
    coh_z = cat(1, all_concord_tab.Coh_z{:});
    
    full_list_pat = [];% Need to create string array with each patient ID repeated to match number of seizure combinations
    
    pat_count = cellfun(@length,all_concord_tab.Coh);
    for pat = 1:size(all_concord_tab,1)
        full_list_pat = [full_list_pat; repmat(all_concord_tab.Patient_id(pat),...
            pat_count(pat),1)];
    end 
    
    full_list_pat = full_list_pat(~isnan(coh));
    coh_z = coh_z(~isnan(coh));
    coh = coh(~isnan(coh));
    
    
    patients = all_concord_tab.Patient_id;
    ordered_patients = [patients(all_concord_tab.Outcome <= 2); patients(all_concord_tab.Outcome > 2)];
    
    f = figure('Position', [0,0,1800,900]);
    subplot(211)
    violinplot(coh,full_list_pat, 'GroupOrder', ordered_patients);
    set(gca, 'box', 'off')
    yline(0.4)
    xline(11.5, 'r')
    ylim([-0.5, 1])
    title("Cohen's Kappa")
    subplot(212)
    patch([0, length(patients)+1, length(patients)+1, 0], [-2,-2,2,2],[.7,.7,.7], 'FaceAlpha', 0.3, 'LineStyle', 'none')
    violinplot(coh_z,full_list_pat, 'GroupOrder', ordered_patients);
    set(gca, 'box', 'off')
    ylim([-4, 12])
    xline(11.5, 'r')
    title("Z(Cohen's Kappa)")
    sgtitle(sprintf("Across-seizure concordance (%s)", chan_or_roi))
    saveas(f, sprintf('figures/rapid_prototype/onset_concord/all_pat_violin_%s.svg',chan_or_roi))
end