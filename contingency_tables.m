%% See from line 124 for maximum recurring region (rather than setting threshold) comparison and use of Nathan's plotBrain function

onset_across_titles = ["subject", "onset across"];
det_meths = ["clo", "imprint", "EI"];

for parc = "roi_120" %["chan", "roi_120", "roi_250"]
    for onset_across = 1 %[0,1]
        if onset_across == 0
            det_meths = det_meths(det_meths~="clo");
            for summ_meas = ["med", "max"]
                fig = figure("Position",[10,10,900,900]);
                tiledlayout(2,2)
                for det_meth = det_meths
                    nexttile
                        data_tab = comp_meth_with_resec.per_sz.(sprintf(parc)).(sprintf(det_meth)).(sprintf(summ_meas));
                        data = data_tab.Perc;
                        out = data_tab.Outcome;
                        out_str = string();
                        out_str(out<3) = "Good";
                        out_str(out>2) = "Bad";
                        resec_str = string();
                        resec_str(data>=0.75) = "Yes";
                        resec_str(data<0.75) = "No";
                        [contingency_tab, ~, ~, lab] = crosstab(resec_str, out_str);
                        out_tab = table(resec_str', out_str');
                        out_tab.Properties.VariableNames = ["resec_75", "Outcome"];
                        [~, p, ~] = fishertest(contingency_tab);    
                        heatmap(out_tab, "resec_75", "Outcome")
                        title(sprintf("%s (p=%.2f)", det_meth,p))
                end
                sgtitle(sprintf("%s-wise comparison against resection (%s %s)", parc, onset_across_titles(onset_across+1), summ_meas))
                %saveas(fig, sprintf("figures/rapid_prototype/onset_resec/%s_contingency_%s_%s.png", parc, onset_across_titles(onset_across+1), summ_meas))
            end
        else 
            fig = figure("Position",[10,10,900,900]);
            tiledlayout(2,2)
            for det_meth = det_meths
                nexttile
                data_tab = comp_meth_with_resec.across_sz.(sprintf(parc)).(sprintf(det_meth));
                data = data_tab.Perc;
                out = data_tab.Outcome;
               
                out_str = string();
                out_str(out<3) = "Good";
                out_str(out>2) = "Bad";
                resec_str = string();
                resec_str(data>=0.75) = "Yes";
                resec_str(data<0.75) = "No";
                [contingency_tab, ~, ~, lab] = crosstab(resec_str, out_str);
                out_tab = table(resec_str', out_str');
                out_tab.Properties.VariableNames = ["resec_75", "Outcome"];
                [~, p, ~] = fishertest(contingency_tab);

                heatmap(out_tab, "resec_75", "Outcome")
                title(sprintf("%s (p=%.2f)", det_meth,p))
            end
             sgtitle(sprintf("%s-wise comparison against resection (%s)", parc, onset_across_titles(onset_across+1)))
             %saveas(fig, sprintf("figures/rapid_prototype/onset_resec/%s_contingency_%s.png", parc, onset_across_titles(onset_across+1)))
         end
    end
end

%%
onset_across_titles = ["subject", "onset across"];
det_meths = ["imprint", "EI"];

for parc = ["chan", "roi_120", "roi_250"]
    for onset_across = [0,1]
        if onset_across == 0
            det_meths = det_meths(det_meths~="clo");
            for summ_meas = ["med", "max"]
                fig = figure("Position",[10,10,900,900]);
                tiledlayout(2,2)
                for det_meth = det_meths
                    nexttile
                        data_tab = comp_meth_with_clo.per_sz.(sprintf(parc)).(sprintf(det_meth)).(sprintf(summ_meas));
                        data = data_tab.Perc;
                        out = data_tab.Outcome;
                        out_str = string();
                        out_str(out<3) = "Good";
                        out_str(out>2) = "Bad";
                        resec_str = string();
                        resec_str(data>=0.75) = "Yes";
                        resec_str(data<0.75) = "No";
                        [contingency_tab, ~, ~, lab] = crosstab(resec_str, out_str);
                        out_tab = table(resec_str', out_str');
                        out_tab.Properties.VariableNames = ["resec_75", "Outcome"];
                        [~, p, ~] = fishertest(contingency_tab);    
                        heatmap(out_tab, "resec_75", "Outcome")
                        title(sprintf("%s (p=%.2f)", det_meth,p))
                end
                sgtitle(sprintf("%s-wise comparison against clo (%s %s)", parc, onset_across_titles(onset_across+1), summ_meas))
                saveas(fig, sprintf("figures/rapid_prototype/onset_clo/%s_contingency_%s_%s.png", parc, onset_across_titles(onset_across+1), summ_meas))
            end
        else 
            fig = figure("Position",[10,10,900,900]);
            tiledlayout(2,2)
            for det_meth = det_meths
                nexttile
                data_tab = comp_meth_with_clo.across_sz.(sprintf(parc)).(sprintf(det_meth));
                data = data_tab.Perc;
                out = data_tab.Outcome;
               
                out_str = string();
                out_str(out<3) = "Good";
                out_str(out>2) = "Bad";
                resec_str = string();
                resec_str(data>=0.75) = "Yes";
                resec_str(data<0.75) = "No";
                [contingency_tab, ~, ~, lab] = crosstab(resec_str, out_str);
                out_tab = table(resec_str', out_str');
                out_tab.Properties.VariableNames = ["resec_75", "Outcome"];
                [~, p, ~] = fishertest(contingency_tab);

                heatmap(out_tab, "resec_75", "Outcome")
                title(sprintf("%s (p=%.2f)", det_meth,p))
            end
             sgtitle(sprintf("%s-wise comparison against clo (%s)", parc, onset_across_titles(onset_across+1)))
             saveas(fig, sprintf("figures/rapid_prototype/onset_clo/%s_contingency_%s.png", parc, onset_across_titles(onset_across+1)))
         end
    end
end

%%
% Max overlap concensus
final_output = final_output(cellfun(@sum, final_output.resected_roi_120)~=0,:);
chan_or_roi = "roi_120";
max_prop = nan(1,size(final_output,1));
out = max_prop;
n_reg_at_max = max_prop;
prop_reg_at_max = max_prop;

for pat = 11:20%1:size(final_output,1)
    pat_onset = final_output(pat,:);
    imprint = pat_onset.imprint_roi_120{:};
    resec = double(pat_onset.resected_roi_120{:});
    if size(imprint,2) <2
        continue
    end
    
    prop_sz = sum(imprint,2)/size(imprint,2);

    high = max(prop_sz);
    max_prop(pat) = high;
    n_reg_at_max(pat) = sum(prop_sz == high);
    prop_reg_at_max(pat) = n_reg_at_max(pat)/length(prop_sz);
    out(pat) = pat_onset.outcome;
    high_print = round(high*100);

    max_reg = double(prop_sz>=high);
    pat_roi = pat_onset.roi_names_120{:};

    % See plotting limits to show relevant hemispheres
    if any(contains(pat_roi,"l.")) & any(contains(pat_roi,"r."))
        % Bilateral placement
        y_lim = [0 450];
        x_lim = [-500 1300];

    elseif any(contains(pat_roi,"l.")) & ~any(contains(pat_roi,"r."))
        % Left hemisphere placement
        y_lim = [230 450]; % If only left hemisphere
        x_lim = [-250 900];

    elseif ~any(contains(pat_roi,"l.")) & any(contains(pat_roi,"r."))
        % Right hemisphere placement
        y_lim = [0 210];
        x_lim = [-250 900];
    end
    
    fig = figure("Position", [0,0,1200,900]);

    subplot(3,1,1)
    plotBrain_NE(pat_roi, max_reg ,"cm", [.7,.7,.7; 1,0,0]);
    title(sprintf("In >=%d%% of onsets", high_print))
    colorbar off
    ylim(y_lim)
    
    subplot(3,1,2)
    plotBrain_NE(pat_roi, resec ,"cm", [.7,.7,.7;1,1,0]);
    clim([0,1])
    title("Resected")
    colorbar off
    ylim(y_lim)
    xlim(x_lim)
    
    subplot(3,1,3)
    col_sch = (max_reg*2)+resec; 
    % Resected only: yellow
    % In low not high, not resected: red
    % In low not high, resected: orange
    cm =  [.7,.7,.7;1,1,0;1,0,0;1,0.5,0];
    if max(col_sch) < size(cm,1)
        cm = cm(1:(max(col_sch)+1),:);    
    end
    cm = interp1(cm, 1:0.01:size(cm,1));
    plotBrain_NE(pat_roi, col_sch ,"cm", cm);
    colorbar off
    hold on
    cm_scatter = [1,1,0;1,0.5,0;1,0,0];
    for i= 1:3
        scatter(-100,-100, [], cm_scatter(i,:), 'filled')
    end
    hold off
    ylim(y_lim)
    xlim(x_lim)
    
    legend(["Resected only", "In max onsets, resec",...
        "In max onsets, not resec"], 'Location','northeast')
    
    legend boxoff  
    sgtitle(sprintf("%s, (ILAE %d)", pat_onset.Patient_id{:}, pat_onset.outcome))

end


%%
figure()
subplot(311)
boxchart(double(out>2), max_prop, "MarkerColor","none")
set(gca, "XTick", [0,1], "XTickLabel", ["ILAE 1-2", "ILAE 3+"])
hold on
swarmchart(double(out>2), max_prop, 'filled')
hold off
title("Maximum proportion of seizures with consistent onset (in at least one region)")
subplot(312)
boxchart(double(out>2), n_reg_at_max, "MarkerColor","none")
set(gca, "XTick", [0,1], "XTickLabel", ["ILAE 1-2", "ILAE 3+"])
hold on
swarmchart(double(out>2), n_reg_at_max, 'filled')
hold off
title("Number of regions recurring in most seizures")
subplot(313)
boxchart(double(out>2), prop_reg_at_max, "MarkerColor","none")
set(gca, "XTick", [0,1], "XTickLabel", ["ILAE 1-2", "ILAE 3+"])
hold on
swarmchart(double(out>2), prop_reg_at_max, 'filled')
title("Proportion of regions recurring in most seizures")
hold off
%%
figure()
gscatter(max_prop, prop_reg_at_max,out>2,'filled')
xlabel("Maximum proportion of seizures with >=1 recurring onset region")
ylabel("Proportion of subject recorded regions recurring in maximum count of seizures")
legend(["ILAE 1-2", "ILAE 3+"])

%%
parc = "roi_120";
fig = figure("Position",[10,10,900,900]);
tiledlayout(2,2)
for det_meth = det_meths
    nexttile
    data_tab = comp_meth_with_resec.across_sz.(sprintf(parc)).(sprintf(det_meth));
    data = data_tab.Perc;
    out = data_tab.Outcome;
   
    out_str = string();
    out_str(out<3) = "Good";
    out_str(out>2) = "Bad";
    resec_str = string();
    resec_str(data>=0.75) = "Yes";
    resec_str(data<0.75) = "No";
    [contingency_tab, ~, ~, lab] = crosstab(resec_str, out_str);
    out_tab = table(resec_str', out_str');
    out_tab.Properties.VariableNames = ["resec_75", "Outcome"];
    [~, p, ~] = fishertest(contingency_tab);

    heatmap(out_tab, "resec_75", "Outcome")
    title(sprintf("%s (p=%.2f)", det_meth,p))
end
 sgtitle(sprintf("%s-wise comparison against resection (%s)", parc, onset_across_titles(onset_across+1)))

 %%

 for parc = ["chan", "roi_120", "roi_250"]
     clo_tab = comp_meth_with_resec.across_sz.(sprintf(parc)).clo;
     imprint_tab = comp_meth_with_resec.across_sz.(sprintf(parc)).imprint;
     joined_tab = innerjoin(clo_tab,imprint_tab , 'Keys',["Patient_id", "Sz_count", "Outcome"]);
    
     resec_thresh = .2;
    
    fig = figure("Position",[10,10,900,450]);
    tiledlayout(1,2)
    for det_meth = ["clo", "imprint"]
        nexttile
        data = joined_tab.(sprintf("Coh_%s_tab", det_meth));
        out = joined_tab.Outcome;
       
        out_str = string();
        out_str(out<3) = "Good";
        out_str(out>2) = "Bad";
        resec_str = string();
        resec_str(data>=resec_thresh) = "Yes";
        resec_str(data<resec_thresh) = "No";
        [contingency_tab, ~, p, lab] = crosstab(resec_str, out_str);
        out_tab = table(resec_str', out_str');
        out_tab.Properties.VariableNames = [sprintf("resec_%d",resec_thresh*100), "Outcome"];
        %[~, p, ~] = fishertest(contingency_tab);
    
        heatmap(out_tab, sprintf("resec_%d",resec_thresh*100), "Outcome")
        title(sprintf("%s (p=%.2f)", det_meth,p))
    end
     sgtitle(sprintf("%s-wise comparison against resection (%s)", parc, onset_across_titles(onset_across+1)))
 end

 %%
figure()
subplot(221)
scatter(comp_meth_with_resec.across_sz.roi_120.imprint.Coh, comp_meth_with_resec.across_sz.roi_120.imprint.Adj_Coh,'filled', 'XJitter','density')
lsline()
xlabel("Cohen's \kappa")
ylabel("Adjusted Cohen's \kappa")
xlim([-0.5,1.1])
ylim([-0.5,1.1])
subplot(222)
scatter(comp_meth_with_resec.across_sz.roi_120.imprint.Perc, comp_meth_with_resec.across_sz.roi_120.imprint.Adj_Coh,'filled', 'XJitter','density')
lsline()
xlabel("Percentage resected")
ylabel("Adjusted Cohen's \kappa")
xlim([-0.1,1.1])
ylim([-0.5,1.1])

subplot(223)
scatter(comp_meth_with_resec.across_sz.roi_120.imprint.Coh_z, comp_meth_with_resec.across_sz.roi_120.imprint.Adj_Coh_z,'filled', 'XJitter','density')
lsline()
xlabel("Z(Cohen's \kappa)")
ylabel("Z(Adjusted Cohen's \kappa)")
xlim([-2.5,4])
ylim([-2.5,4])
patch([-2 2 2 -2], [-2 -2 2 2], [.7,.7,.7], "FaceAlpha", 0.3, "LineStyle", "none")
subplot(224)
scatter(comp_meth_with_resec.across_sz.roi_120.imprint.Perc_z, comp_meth_with_resec.across_sz.roi_120.imprint.Adj_Coh_z,'filled', 'XJitter','density')
lsline()
xlabel("Z(Percentage resected)")
ylabel("Z(Adjusted Cohen's \kappa)")
xlim([-2.5,4])
ylim([-2.5,4])
patch([-2 2 2 -2], [-2 -2 2 2], [.7,.7,.7], "FaceAlpha", 0.3, "LineStyle", "none")






% 
% 
% for parc = ["chan", "roi_120", "roi_250"]
%     for onset_across = [0,1]
%         if onset_across == 0
%             det_meths = det_meths(det_meths~="clo");
%             for summ_meas = ["med", "max"]
%                 fig = figure("Position",[10,10,900,900]);
%                 tiledlayout(2,2)
%                 for det_meth = det_meths
%                     nexttile
%                         data_tab = comp_meth_with_clo.per_sz.(sprintf(parc)).(sprintf(det_meth)).(sprintf(summ_meas));
%                         data = data_tab.Perc;
%                         out = data_tab.Outcome;
%                         [contingency_tab, ~, ~, lab] = crosstab(data>=0.75, (out >2)+1);
%                         out_tab = table(discretize(data==1,[0,0.5,1], 'categorical', ["No", "Yes"]), discretize(out>2,[0,0.5,1], 'categorical', ["Good", "Bad"]));
%                         out_tab.Properties.VariableNames = ["Over 75% captured", "Outcome"];
%                         [~, p, ~] = fishertest(contingency_tab);    
%                         heatmap(out_tab, "Over 75% captured", "Outcome")
%                         title(sprintf("%s (p=%.2f)", det_meth,p))
%                 end
%                 sgtitle(sprintf("%s-wise comparison against CLO (%s %s)", parc, onset_across_titles(onset_across+1), summ_meas))
%                 saveas(fig, sprintf("figures/rapid_prototype/onset_clo/%s_contingency_%s_%s.png", parc, onset_across_titles(onset_across+1), summ_meas))
%             end
%         else 
%             fig = figure("Position",[10,10,900,900]);
%             tiledlayout(2,2)
%             for det_meth = det_meths
%                 nexttile
%                 data_tab = comp_meth_with_clo.across_sz.(sprintf(parc)).(sprintf(det_meth));
%                 data = data_tab.Perc;
%                 out = data_tab.Outcome;
%                 [contingency_tab, ~, ~, lab] = crosstab(data>=0.75, (out >2)+1);
%                 out_tab = table(discretize(data==1,[0,0.5,1], 'categorical', ["No", "Yes"]), discretize(out>2,[0,0.5,1], 'categorical', ["Good", "Bad"]));
%                 out_tab.Properties.VariableNames = ["Over 75% captured", "Outcome"];
%                 [~, p, ~] = fishertest(contingency_tab);
% 
%                 heatmap(out_tab, "Over 75% captured", "Outcome")
%                 title(sprintf("%s (p=%.2f)", det_meth,p))
%             end
%              sgtitle(sprintf("%s-wise comparison against CLO (%s)", parc, onset_across_titles(onset_across+1)))
%              saveas(fig, sprintf("figures/rapid_prototype/onset_clo/%s_contingency_%s.png", parc, onset_across_titles(onset_across+1)))
%         end
%     end
% end

% onset_across_titles = ["subject", "onset across"];
% 
% for parc = ["chan", "roi_120", "roi_250"]
%     for onset_across = [0,1]
%         if onset_across == 0
%             det_meths = ["imprint", "EI", "PLHG"];
%             for summ_meas = ["med", "max"]
%                 fig = figure("Position",[10,10,900,900]);
%                 tiledlayout(2,2)
%                 for det_meth = det_meths
%                     nexttile
%                         data_tab = comp_meth.per_sz.(sprintf(parc)).(sprintf(det_meth)).(sprintf(summ_meas));
%                         data = data_tab.Perc;
%                         out = data_tab.Outcome;
%                         [contingency_tab, ~, ~, lab] = crosstab(data==1, (out >2)+1);
%                         out_tab = table(discretize(data==1,[0,0.5,1], 'categorical', ["No", "Yes"]), discretize(out>2,[0,0.5,1], 'categorical', ["Bad", "Good"]));
%                         out_tab.Properties.VariableNames = ["All_resec", "Outcome"];
%     
%                         heatmap(out_tab, "All_resec", "Outcome")
%                         title(det_meth)
%                 end
%                 sgtitle(sprintf("%s-wise comparison against resection (%s %s)", parc, onset_across_titles(onset_across+1), summ_meas))
%                 saveas(fig, sprintf("figures/rapid_prototype/onset_resec/%s_contingency_%s_%s.png", parc, onset_across_titles(onset_across+1), summ_meas))
%             end
%         else 
%             det_meths = ["clo", "imprint", "EI", "PLHG"];
%             fig = figure("Position",[10,10,900,900]);
%             tiledlayout(2,2)
%             for det_meth = det_meths
%                 nexttile
%                 data_tab = comp_meth.across_sz.(sprintf(parc)).(sprintf(det_meth));
%                 data = data_tab.Perc;
%                 out = data_tab.Outcome;
%                 [contingency_tab, ~, ~, lab] = crosstab(data==1, (out >2)+1);
%                 out_tab = table(discretize(data==1,[0,0.5,1], 'categorical', ["No", "Yes"]), discretize(out>2,[0,0.5,1], 'categorical', ["Bad", "Good"]));
%                 out_tab.Properties.VariableNames = ["All_resec", "Outcome"];
% 
%                 heatmap(out_tab, "All_resec", "Outcome")
%                 title(det_meth)
%             end
%              sgtitle(sprintf("%s-wise comparison against resection (%s)", parc, onset_across_titles(onset_across+1)))
%              saveas(fig, sprintf("figures/rapid_prototype/onset_resec/%s_contingency_%s.png", parc, onset_across_titles(onset_across+1)))
% 
%         end
%     end
% end
% 
% %%
% 
% 
% 
% 
% 
% figure()
% tiledlayout(2,2)
% 
% for det_meth = ["clo", "imprint", "EI", "PLHG"]
%     nexttile
% 
%     data_tab = comp_meth.across_sz.roi_120.(sprintf(det_meth));
%     data = data_tab.Perc;
%     out = data_tab.Outcome;
%     
%     %[contingency_tab, chisq,p] = crosstab(data==1, out >2)
%     
%     [contingency_tab, ~, ~, lab] = crosstab(data==1, (out >2)+1)
%     
%     
%     %%
%     out_tab = table(discretize(data==1,[0,0.5,1], 'categorical', ["No", "Yes"]), discretize(out>2,[0,0.5,1], 'categorical', ["Bad", "Good"]));
%     out_tab.Properties.VariableNames = ["All_resec", "Outcome"];
%     %[contingency_tab, ~, ~, lab] = crosstab(out_tab)
%     
%     heatmap(out_tab, "All_resec", "Outcome")
%     title(det_meth)
% end
% 
% %%
% figure()
% tiledlayout(2,2)
% 
% for det_meth = ["imprint", "EI", "PLHG"]
%     nexttile
% 
%     data_tab = comp_meth_with_clo.across_sz.roi_120.(sprintf(det_meth));
%     data = data_tab.Perc;
%     out = data_tab.Outcome;
%     
%     %[contingency_tab, chisq,p] = crosstab(data==1, out >2)
%     
%     [contingency_tab, ~, ~, lab] = crosstab(data==1, (out >2)+1);
%     
%     
%     %%
%     out_tab = table(discretize(data==1,[0,0.5,1], 'categorical', ["No", "Yes"]), discretize(out>2,[0,0.5,1], 'categorical', ["Bad", "Good"]));
%     out_tab.Properties.VariableNames = ["All_resec", "Outcome"];
%     
%     heatmap(out_tab, "All_resec", "Outcome")
%     title(det_meth)
% end