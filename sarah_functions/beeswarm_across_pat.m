function [] = beeswarm_across_pat(final_comp, onset_output, opts)
% Compute average EEG across regions of interest for each patient

% Sarah J Gascoigne 13/02/2023

% input:
%   - final_comp
%   - onset_output

%   - optional inputs
%       - det_method: state the method of onset detection you'd like to use
%           ("imprint", "EI", "PLHG")
%       - comparison: state whether the comparisons are based on
%           comparisons between onset and resected regions ("resction") or
%           between pairs of seizure onsets ("pairwise")
%       - save_plot: 1 to save, 0 to not save
%       - save_location: string indicating where to save the plot 

% output
%   - comp_table: table containing raw and normalised Jaccard's index and
%       Sorensen-Dice coefficient

% Could add an argument to decide if we should plot Jaccard/Sorensen
% against count of regions to show the difference

    arguments
        final_comp
        onset_output
        opts.det_method (1,1) string {mustBeMember(opts.det_method, ["imprint", "EI", "PLHG"])} = "imprint" 
        opts.comparison (1,1) string {mustBeMember(opts.comparison, ["resection", "pairwise"])} = "resection"
        opts.chan_or_roi (1,1) string {mustBeMember(opts.chan_or_roi, ["chan", "roi"])} = "roi"
        opts.comp_measure (1,1) string {mustBeMember(opts.comp_measure,...
            ["Jaccard","Jaccard_norm", "Sorensen", "Sorensen_norm", "Percentage_resec", "Hausdorff", "Hausdorff_norm"])} = "Jaccard_norm" 
        opts.save_plot (1,1) double = 0
        opts.save_location = 'figures/across_patients/';
    end
    
    %fill in optional arguments
    det_method = opts.det_method;
    comparison = opts.comparison; 
    chan_or_roi = opts.chan_or_roi;
    comp_measure = opts.comp_measure;
    save_plot = opts.save_plot;
    save_location = opts.save_location;
    
    % Beeswarm plot across patients
    % We will be organising our beeswarm plots by outcome (good/bad) and by
    % median within these groups
    val_tbl = final_comp(:,ismember(final_comp.Properties.VariableNames,...
        {'Patient_id', 'Jaccard', 'Jaccard_norm', 'Sorensen',...
        'Sorensen_norm', 'Percentage_resec','Hausdorff', 'Hausdorff_norm', 'Y1_outcome'}));
    pat_median = grpstats(val_tbl, "Patient_id", "median");

%     if comparison == "resection"
%         val_tbl = final_comp; %(:,{'Patient_id', 'Jaccard', 'Jaccard_norm', 'Sorensen',...
%             %'Sorensen_norm', 'Percentage_resec','Hausdorff', 'Hausdorff_norm'});
%         pat_median = grpstats(val_tbl, "Patient_id", "median");
%         pat_median = join(pat_median, onset_output(:,[1 9]));
%     elseif comparison == "pairwise"
%         val_tbl = final_comp;%(:,{'Patient_id', 'Jaccard', 'Jaccard_norm', 'Sorensen',...
%             %'Sorensen_norm', 'Hausdorff', 'Hausdorff_norm'});
%         pat_median = grpstats(val_tbl, "Patient_id", "median");
%         pat_median = join(pat_median, onset_output(:,[1 9]));
%     end
    
    % Pull out column of comparison measure of interest
    column_index = find(string(final_comp.Properties.VariableNames)...
        == comp_measure);
    val_col = table2array(final_comp(:,column_index));

    outcome = zeros(1,size(onset_output,1));
    
    % Organise patients into good and bad outcomes
    for pat = 1:size(onset_output,1)
        outcome(pat) = pat_median.median_Y1_outcome(pat);
    end
    good_mean = pat_median(outcome <= 2,:);
    bad_mean = pat_median(outcome > 2,:);
    
    % Within groups, organise by mean 
    [~, Ig] = sort(good_mean.(sprintf("median_%s", comp_measure)), 'descend');
    [~, Ib] = sort(bad_mean.(sprintf("median_%s", comp_measure)), 'descend');
      
    % Find indeces to organise by mean value
    ind = [Ig; (size(good_mean,1)+Ib)];
    
    pat_id_double = str2double(final_comp.Patient_id);
    patient_ids = [str2double(good_mean.Patient_id); str2double(bad_mean.Patient_id)];
    reorder_pat_label = nan(length(pat_id_double), 1);
    for x = 1:length(patient_ids)
        i = patient_ids(ind(x));
        reorder_pat_label(pat_id_double == i) = x;
    end
    
    figure()
    boxchart(reorder_pat_label, val_col, 'MarkerStyle','none')
    hold on
    swarmchart(reorder_pat_label, val_col,'filled', 'XJitterWidth',1)
    hold off
    xlim([-1 (length(outcome)+1)])
    xticks(1:length(outcome))
    xticklabels(patient_ids(ind))
    xtickangle(45)
    xlabel("Patient ID")

    if comp_measure == "Hausdorff_norm"
        ylim([0 1])
    elseif comp_measure == "Percentage_resec"
        ylim([0 1])
    else
        yline(0)
    end
    
    xline(size(good_mean,1)+0.5)
    title(sprintf('%s %s across patients (based on %s)', comparison, strrep(comp_measure,'_',' '), det_method))
    
    if save_plot == 1
        mkdir(sprintf('%s%s_%s', save_location, comparison, chan_or_roi))
        saveas(gcf,sprintf('%s%s_%s/%s_all_pat_%s', save_location, comparison, chan_or_roi, comp_measure, det_method), 'png')
    end
end