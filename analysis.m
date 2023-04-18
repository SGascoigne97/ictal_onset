%% Analysis

%% Does the proportion of resected CLO regions differ across patients with 
% favourable and unfavourable outcomes (ILAE 1-2 and 3+, respectively)?

% Using proportion of regions resected (CLO)
comparison = "resection";
det_method = "CLO";
load(sprintf('tables/across_pat_%s_%s.mat', det_method, comparison));
comp_table = final_comp; % check this 
comp_table.Surgery_outcome = cell2mat(comp_table.Surgery_outcome); % need to change this line

% Check distributions for both outcome groups (should we use a T-test or 
% Wilcoxon-rank?)
% Perform relevant test and obtain test statistic and p-value

% Is the proportion of regions resected larger for ILAE 1-2 and/or smaller 
% for ILAE 3+?

% Visual checks
figure(1)
subplot(121)
beeswarm(ones(length(comp_table.Percentage_resec),1),...
    comp_table.Percentage_resec, 3,'sort_style','up','overlay_style','sd',...
        'dot_size', 2, 'corral_style', 'gutter');
ylim([-0.05 1.05])
xlim([0.5 1.5])
yline(0.5,'--')
title('All Patients')
subplot(122)
beeswarm(comp_table.Surgery_outcome, comp_table.Percentage_resec, 3,...
    'sort_style','up','overlay_style','sd',...
    'dot_size', 2, 'corral_style', 'gutter');
xlim([-0.5 1.5])
ylim([-0.05 1.05])
yline(0.5,'--')
title('Comparison against outcome')
sgtitle("Do surgeons tend to remove CLO regions?")
if save_plot ==1
    saveas(gcf,'figures/across_patients/CLO/perc_resec', 'png')
end


%% Can we distinguish between good and bad outcome patients based on 
% comparisons between automatically detected onsets and resected regions?

% For each onset detection method and comparison measure, check 
% distributions for both outcome groups (should we use a T-test or 
% Wilcoxon-rank?)
% Perform relevant test and obtain test statistic and p-value

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
% onset regions resected different across good and bad outcome patients?

% Plot histograms of frequencies for both ILAE 1-2 and ILAE 3+
% Visual differences?

% For each onset detection method and comparison measure, perform logistic 
% regression to predict outcome based on the frequency of few and most
% onset regions resected 


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