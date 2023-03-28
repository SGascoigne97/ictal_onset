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
        opts.comp_measure (1,1) string {mustBeMember(opts.comp_measure,...
            ["Jaccard","Jaccard_norm", "Sorensen", "Sorensen_norm", "Percentage_resec", "Hausdorff", "Hausdorff_norm"])} = "Jaccard_norm" 
        opts.save_plot (1,1) double = 0
        opts.save_location = 'figures/across_patients/';
    end
    
    %fill in optional arguments
    det_method = opts.det_method;
    comparison = opts.comparison; 
    comp_measure = opts.comp_measure;
    save_plot = opts.save_plot;
    save_location = opts.save_location;
    
    % Beeswarm plot across patients
    % We will be organising our beeswarm plots by outcome (good/bad) and by
    % mean within these groups
    if comparison == "resection"
        val_tbl = final_comp(:,{'Patient_id', 'Jaccard', 'Jaccard_norm', 'Sorensen',...
            'Sorensen_norm', 'Percentage_resec','Hausdorff', 'Hausdorff_norm'});
        pat_mean = grpstats(val_tbl, "Patient_id", "mean");
        pat_mean = join(pat_mean, onset_output(:,[1 9]));
    elseif comparison == "pairwise"
        val_tbl = final_comp(:,{'Patient_id', 'Jaccard', 'Jaccard_norm', 'Sorensen',...
            'Sorensen_norm', 'Hausdorff', 'Hausdorff_norm'});
        pat_mean = grpstats(val_tbl, "Patient_id", "mean");
        pat_mean = join(pat_mean, onset_output(:,[1 9]));
    end
    
    % Pull out column of comparison measure of interest
    column_index = find(string(final_comp.Properties.VariableNames)...
        == comp_measure);
    val_col = table2array(final_comp(:,column_index));
    
    % Organise patients into good and bad outcomes
    good_mean = pat_mean(cell2mat(pat_mean.Surgery_outcome) == 1,:);
    bad_mean = pat_mean(cell2mat(pat_mean.Surgery_outcome) == 0,:);
    
    % Within groups, organise by mean 
    [~, Ig] = sort(good_mean.(sprintf("mean_%s", comp_measure)), 'descend');
    [~, Ib] = sort(bad_mean.(sprintf("mean_%s", comp_measure)), 'descend');
      
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
    beeswarm(reorder_pat_label, val_col, 3,'sort_style','up','overlay_style','sd',...
        'dot_size', 2, 'corral_style', 'gutter');
    xlim([-1 16])
    if comp_measure == "Hausdorff_norm"
        ylim([0 1])
    elseif comp_measure == "Percentage_resec"
        ylim([0 1])
    end
    
    xline(size(good_mean,1)+0.5)
    title(sprintf('%s %s across patients (based on %s)', comparison, strrep(comp_measure,'_',' '), det_method))
    
    if save_plot == 1
        saveas(gcf,sprintf('%s%s/%s_all_pat_%s', save_location, comparison, comp_measure, det_method), 'png')
    end
end