%% Analysis
save_plot = 1;

%% Does the proportion of resected CLO regions differ across patients with 
% favourable and unfavourable outcomes (ILAE 1-2 and 3+, respectively)?

% Using proportion of regions resected (CLO)
comparison = "resection";
det_method = "CLO";
load(sprintf('tables/across_pat_%s_%s.mat', det_method, comparison));
clo_comp = final_comp;

clo_comp.outcome_group = categorical(clo_comp.Y1_outcome > 2);
clo_comp.outcome_group(isnan(clo_comp.Y1_outcome)) = "Unknown";

clo_comp.outcome_group = renamecats(clo_comp.outcome_group, "false", "Good");
clo_comp.outcome_group = renamecats(clo_comp.outcome_group, "true", "Bad");

%% Check distributions for both outcome groups (should we use a T-test or 
% Wilcoxon-rank?)
% Perform relevant test and obtain test statistic and p-value
    % Data is not normally distributed - therefore we will use
    % Wilcoxon-Rank

figure()
tiledlayout(3,1)
nexttile
histogram(clo_comp.Percentage_resec, 'BinWidth', 0.1)
xlim([-.05,1.05])
ylim([0 8])
[p1,~, stats1] = ranksum(clo_comp.Percentage_resec(clo_comp.Y1_outcome <3),...
    clo_comp.Percentage_resec(clo_comp.Y1_outcome >2));
title(sprintf('Comp 1-2 vs. 3+ (r = %.3f, p = %.3f)', stats1.ranksum, p1))
nexttile
histogram(clo_comp.Percentage_resec(clo_comp.Y1_outcome <3), 'BinWidth', 0.1)
xlim([-.05,1.05])
ylim([0 8])
[p2,~, stats2] = ranksum(clo_comp.Percentage_resec(clo_comp.Y1_outcome <3),...
    0.5);
title(sprintf('ILAE 1-2 (r = %.3f, p = %.3f)', stats2.ranksum, p2))
nexttile
histogram(clo_comp.Percentage_resec(clo_comp.Y1_outcome >2), 'BinWidth', 0.1)
xlim([-.05,1.05])
ylim([0 8])
[p2,~, stats2] = ranksum(clo_comp.Percentage_resec(clo_comp.Y1_outcome >2),...
    0.5);
title(sprintf('ILAE 3+ (r = %.3f, p = %.3f)', stats2.ranksum, p2))
sgtitle('Histograms of the proportion of CLO regions resected')
if save_plot == 1
    saveas(gcf,'figures/across_patients/resection_roi/CLO_perc_resec_hist', 'png')
end
% Visual checks
figure()
subplot(121)
boxchart(ones(length(clo_comp.Percentage_resec),1),...
    clo_comp.Percentage_resec)
hold on
swarmchart(ones(length(clo_comp.Percentage_resec),1),...
    clo_comp.Percentage_resec, 'filled')
hold off
ylim([-0.05 1.05])
xlim([0.5 1.5])
yline(0.5,'--')
title('All Patients')
subplot(122)
boxchart(clo_comp.outcome_group, clo_comp.Percentage_resec)
hold on
swarmchart(clo_comp.outcome_group, clo_comp.Percentage_resec, 'filled')
hold off
ylim([-0.05 1.05])
yline(0.5,'--')
title('Comparison against outcome')
sgtitle("Do surgeons tend to remove CLO regions?")
if save_plot ==1
    saveas(gcf,'figures/across_patients/resection_roi/CLO_perc_resec', 'png')
end


%% Can we distinguish between good and bad outcome patients based on 
% comparisons between automatically detected onsets and resected regions?

% Is the overlap with resected larger for ILAE 1-2 and/or smaller 
% for ILAE 3+? 
% Is the distance between onset and resection larger for ILAE 3+ and/or 
% smaller for ILAE 1-2? 


%% Can we distinguish between good and bad outcome patients based on 
% pairwise comparisons between automatically detected onsets?

% For each onset detection method and comparison measure, check 
% distributions for both outcome groups (should we use a T-test or 
% Wilcoxon-rank?)
% Perform relevant test and obtain test statistic and p-value

% Is the overlap between onsets in ILAE 1-2 larger than ILAE 3+?
% Is the distance between onsets in ILAE 1-2 less than ILAE 3+?


%% Is the frequency of few (<25%) or many (>75%) automatically detected
% onset regions resected different across good and bad outcome
% patients?
saveplot = 0;

% Plot histograms of frequencies for both ILAE 1-2 and ILAE 3+
% Visual differences?
few_lim = 0.25;
most_lim = 0.75; 
comparison = "resection";
for det_method = ["imprint"] %, "EI", "PLHG"]
    %det_method = "imprint";
    chan_or_roi = "roi";
    resec_roi_comp = readtable(sprintf('tables/across_pat_%s_%s_%s.csv', det_method, comparison, chan_or_roi));
    % For each onset detection method and comparison measure, check 
    % distributions for both outcome groups (should we use a T-test or 
    % Wilcoxon-rank?)
    % T-tests to compare patients with <100% of onset regions resected -must
    % check histograms first!
    
    patients = unique(resec_roi_comp.Patient_id, 'stable');
    most_few_tab = array2table(zeros(length(unique(resec_roi_comp.Patient_id)),3));
    most_few_tab.Properties.VariableNames = {'Patient_ID', 'Prop_few', 'Prop_most'};
    most_few_tab.Patient_ID = patients;
    
    for pat = 1:length(patients)
        pat_comp = resec_roi_comp(resec_roi_comp.Patient_id == patients(pat),:);
        most_few_tab.Prop_few(pat) = sum(pat_comp.Percentage_resec <=few_lim)/size(pat_comp,1);
        most_few_tab.Prop_most(pat) = sum(pat_comp.Percentage_resec >=most_lim)/size(pat_comp,1);
        most_few_tab.Outcome(pat) = discretize(pat_comp.Y1_outcome(1),[0 2.1 5.1], 'categorical',["Good", "Bad"]);
    end
    
    figure()
    subplot(211)
    swarmchart( most_few_tab.Outcome, most_few_tab.Prop_few, 'filled')
    hold on
    boxchart( most_few_tab.Outcome, most_few_tab.Prop_few)
    hold off
    title('Proportion of seizures with few (<=25%) onset regions resected')
    subplot(212)
    swarmchart( most_few_tab.Outcome, most_few_tab.Prop_most, 'filled')
    hold on
    boxchart( most_few_tab.Outcome, most_few_tab.Prop_most)
    hold off
    title('Proportion of seizures with most (>=75%) onset regions resected')
    sgtitle(sprintf("%s", det_method))
    if save_plot ==1
        saveas(gcf,sprintf('figures/across_patients/resection_roi/most_few_%s', det_method), 'png')
    end


end

% For each onset detection method and comparison measure, perform logistic 
% regression to predict outcome based on the frequency of few and most
% onset regions resected 

% TO DO
%% Do patients with focal onsets tend to have favourable outcomes? 
% Similarly, do patients with diffuse onsets tend to have unfavourable 
% outcomes? 

% Using CLO for now (would need to summarise across seizures for each
% automatic onset detection method)
% Categorise onsets based on raw channel xyz coordinates 

% Need to determine a threshold for distinction between focal and diffuse

% Create contingency table:
%         | Favourable  | Unfavourable
%_________|_____________|_____________
% Focal   |             |
%_________|_____________|_____________
% Diffuse |             |
%         |             |

% Perform Chi^2 test 
% Perform post-hoc test (where are the differences?)