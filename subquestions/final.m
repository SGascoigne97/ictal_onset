%% Code for results in "More complete resection of the seizure onset zone
% based on iEEG is not directly associated to post-surgical outcomes"
% (S J Gascoigne et al., 2023)

%% Load data and Lausanne atlases
%final_output_struct = load('../tables/final_output.mat'); % Here any subjects without follow-up have been removed
%final_output = final_output_struct.final_output;
% clear final_output_struct
addpath(genpath('../sarah_functions'))

load("final_output.mat")

load('../roi_info/ATLAS.mat')
final_output = final_output(final_output.outcome ~= 8,:);

%% Set parameters for analyses
%chan_or_roi = "roi_120"; % All analyses are based on Lausanne 120 atlas
    n_perm = 1000; % All AUC p-values are computed based on permutation test 
                   % with 5000 permutations
    out_thresh = 2;%2; % ILAE 1-2 is considered as a 'favourable' outcome whilst 
                    % ILAE 3+ is considered as an 'unfavourable' outcome
    consensus_thresh = 0.5; % Set threshold for inclusion of regions in 
                            % consensus onset, here we will include regions 
                            % present in at least half of the subject's seizures
    most_resec_thresh = 0.5; % Set threshold above which subjects are considered 
                             % as having 'most' of their onset resected
    diffuse_thresh = 0.75;
    det_meths = ["clo", "imprint", "EI"]; % List onset detection methods you 
                                          % are interested in analysing
                                          
    
    % Set parameters for saving figures
    save_fig = 1; % Each figure that is created will be saved 
    save_folder = '../figures/paper_figures/Figure 3/'; % Location where subfolders will
                                              % be created to store figures and tables
    
    file_type = "svg"; % Figures will be saved as svgs so they can be pulled 
                       % into illustrator for paper figures
    
for parc = ["chan", "roi_250", "roi_120"]
    [chan_or_roi, n_perm, out_thresh, consensus_thresh, most_resec_thresh,...
        det_meths, save_fig, save_folder, file_type, out_grps] =...
        argument_validation(final_output, parc, n_perm, out_thresh, consensus_thresh,...
        most_resec_thresh, det_meths, save_fig, save_folder, file_type);
    
    %% Additional parameters used throughout analyses (not to be adjused)
    atl_inds = 2*(atlas.scale');
    if contains(chan_or_roi, "roi")
        atl = atlas(atl_inds == str2double(extractAfter(chan_or_roi, "_")),:);
        vols = atl.vol{:};
        names = atl.name{:};
        dists = atl.dists{:};
    else
        atl = NaN;
        vols = NaN;
        names = NaN;
        dists = NaN;
    end
    
    %% Add outcome category column to final_output table
    % Outcome category (binary based on selected outcome threshold)
    final_output.outcome_cat = categorical(final_output.outcome>out_thresh,[0,1], out_grps);
    
    %% 3.1 IOLA onsets have moderate agreement with CLO in most subjects
    
    % We will consider concordance between the consensus onset and the maximum
    % consensus between the IOLA onset for any one seizure and the clinically
    % labelled onset (CLO)
    
    concordance_with_clo
    % concord_tab: subject level concordance between consensus onset and CLO
    % and maximum concordance between CLO and any one IOLA onset
    
    % auc_clo_comp: table of AUCs and associated p-values for distinguishing
    % surgical outcome groups based on concordance 
    
    %% 3.2 Seizure onset tends to be resected, however resecting a larger
    % proportion of the SOZ is not associated with more favourable surgical
    % outcomes
    
    % We will compute the proportion of the seizure onset zone (CLO and IOLA)
    % that was subsequently resected then compare across outcome groups
    
    is_onset_resected
    
    %% 3.3 Larger onsets are not associated with surgical outcomes
    
    % Here we will compute the size of seizure onsets (both CLO and
    % automatically captured) using both count of regions and volume (based on
    % controls) and compare across surgical outcome groups
    
    larger_onset
    
    %% 3.4 Larger resection is not associated with surgical outcomes
    
    % Here we will compute the size of resections using both count of regions
    % and volume (if using ROIS, based on controls) and compare across surgical
    % outcome groups
    
    larger_resection
    
    %% 3.5 More diffuse onsets are not associated with surgical outcomes
    
    % Here we will estimate the diffusivity of seizure onsets based on the
    % maximum distance between the centers of regions in onset
    
    if contains(chan_or_roi, "roi") % This analysis is only possible when using regions
        more_diffuse
    end
end



%% SUPPLEMENTARY
%% S2 Subject Metadata
% Here we will look at various metadata variables to report in the 
% supplementary results. Further, we will investigate if any metadata
% groups have significant differences in surgical outcomes
% have confounding effect on surgical outcomes
metadata

%% S3.2 Optimising MAD threshold τ
% Here we will display how the MAD threshold for identifying abnormal iEEG
% activity was tuned for use in IOLA
%../scan_mad_thresh

%% S4.X  Concordance between seizure onsets is not associated with surgical outcomes

