load('roi_info/ATLAS_.mat')
load('roi_info/retains.mat')
load('roi_info/atlasinfo.mat')
load('roi_info/ATLAS.mat')
load('roi_info/lhSurf.mat') 
load('roi_info/rhSurf.mat')
addpath('figures')
channel_folder = '../../channels_ROI/';
%%



% Organise patients by good and bad surgical outcome
patients = onset_output.Patient_id;
good_out = cell2mat(onset_output.Surgery_outcome) == 1;
bad_out = ~good_out;
good_pat = patients(good_out);
bad_pat = patients(bad_out);
patients = [good_pat; bad_pat];
onset_measures = ["imprint", "EI", "PLHG"];
comp_measures = ["jac", "sor","haus", "perc"];
av_method = ["median", "mean", "max"];

vals_table_resec = array2table(zeros(length(onset_measures)*...
    length(comp_measures)*length(av_method),3), ...
    'VariableNames',["Onset Measure", "Comparison Measure", ...
    "Averaging Method"]); %, "Values"
vals_table_resec.("Onset Measure") = repelem(onset_measures', ...
    length(comp_measures)*length(av_method),1);
vals_table_resec.("Comparison Measure") = repmat(repelem(comp_measures',...
    length(av_method),1),length(onset_measures),1);
vals_table_resec.("Averaging Method") = repmat(av_method',...
    length(onset_measures)*length(comp_measures), 1);

freq_perc = array2table(zeros(length(onset_measures),1),...
    'VariableNames',"Onset Measure");
freq_perc.("Onset Measure") = onset_measures';
for ons = 1:3
    onset_measure = onset_measures(ons);
    jac = [];
    sor = [];
    perc = [];
    haus = [];
    % Load onset table (if computed and saved previously)
    % onset = load("onset_output");
    % onset_output = onset.onset_output;
    all_pat_index = {};
    for pat = 1:length(onset_output.Patient_id)
        % Choose patient to assess
        patient = patients(pat); 

        % Obtain channel information for this patient (stored under old
        % identification)
        channel_info = load(sprintf('%s%s/channels.mat',channel_folder,patient));
        % Separate out seizures for this patient
        pat_data = data(string(data.patient_id) == patient,:);
        % Capture the names of ROIs included for this patient
        unq_roi = onset_output(onset_output.Patient_id == patient,:).ROI_ids{1,1};
        % Obtain the indicies for each ROI from the atlas 
        rois = contains(string(atlas.name{3,1}), unq_roi);
        rois = find(rois);
        
        % Extract onset information
        pat_onset = onset_output(onset_output.Patient_id == patient,:);
        
        % Each patient has labelled onset from clinicians
        labelled_onset =  rois.*cell2mat(pat_onset.Labelled_onset);
        labelled_onset = labelled_onset(labelled_onset~=0);
        
        % Extract onset based on selected onset metric
        if onset_measure == "imprint"
            pat_onset = onset_output(onset_output.Patient_id ==...
                patient,:).imprint_roi{1,1}; 
        elseif onset_measure == "EI"
            pat_onset = onset_output(onset_output.Patient_id ==...
                patient,:).EI_roi{1,1}; 
        elseif onset_measure == "PLHG"
            pat_onset = onset_output(onset_output.Patient_id ==...
                patient,:).PLHG_roi{1,1}; 
        end
        
        % remove seizures with no onset detected
        pat_onset = pat_onset(:,sum(pat_onset)~=0);
        sz_count = size(pat_onset,2);
        pat_data = pat_data(sum(pat_onset)~=0,:);
        this_pat_index = repmat(pat, sz_count,1);
        all_pat_index(end+1) = {this_pat_index};
       
        pat_table = onset_output(onset_output.Patient_id ==...
                patient,:);
       
        %%
        % Extract matrix of Euclidean distances between ROIs
        roi_distances = cell2mat(atlas.dists(3,:));
        
        % Identify rois within the atlas (indecies)
        unq_binary = contains(string(atlas.name{3,1}), unq_roi);
        % Extract xyz location of regions included in the atlas
        xyz = atlas.xyz{3,1};

        % Identify resected regions
        resec_roi_binary = contains(string(atlas.name{3,1}), unq_roi(find(pat_table.Resected{1,1})));
        resec_xyz = xyz(find(resec_roi_binary),:);
    
        %% 
        % Compute Jaccard's index and Sorenson's-Dice coefficient to quantify 
        % overlap between onset and resection
        % In this calculation we only consider the regions recorded
        jac_ind = nan(1, sz_count);
        sor_dice = nan(1, sz_count);
        perc_rec = nan(1, sz_count);
        hausdorff = nan(1, sz_count);
        
        for i=1:sz_count
            onset_roi_binary = contains(string(atlas.name{3,1}), unq_roi(find(pat_onset(:,i))));
            jac_ind(i) = jaccard(resec_roi_binary, onset_roi_binary);
            sor_dice(i) = (2*jac_ind(i))/(1+jac_ind(i));
            % Compute the percentage of onset regions that were resected
            perc_rec(i) = sum((onset_roi_binary + resec_roi_binary)==2)/sum(onset_roi_binary);
        end

        % Jaccard and Sorensen are biased by the size of the onset region
        % and resection, we will scale these values by the size of the
        % resection for now
        scaled_jac_ind = jac_ind; %/sum(resec_roi_binary);
        scaled_sor_dice = sor_dice; %/sum(resec_roi_binary);
    
        % Store overlap metrics between onset and resection 
        jac = [jac; scaled_jac_ind'];
        sor = [sor; scaled_sor_dice'];
        perc = [perc; perc_rec'];
       
        % Compute Hausdorff distance between onset regions and resection
        resec_xyz = xyz(find(resec_roi_binary),:);
        for sz = 1:sz_count
            onset_roi_binary = contains(string(atlas.name{3,1}), unq_roi(find(pat_onset(:,sz))));
            onset_xyz = xyz(find(onset_roi_binary),:);
            [hd D] = HausdorffDist(resec_xyz, onset_xyz); %, [],'visualize');
            hausdorff(sz) = hd/max(max(atlas.dists{3,1}));
        end
        % Store Hausdorff's distance between onset and resection 
        haus = [haus; hausdorff'];
    
    end
    pat_index = cat(1,all_pat_index{:});
    %% Create a table storing Jaccard, Sorensen, percentage resected, and Hausdorff
    comp_resec_table = table(pat_index, jac, sor, perc, haus);
    writetable(comp_resec_table, sprintf('comp_all_pat_%s.csv', onset_measure))
    for com = 1:length(comp_measures)
        % organise patients by median value within good and bad outcome
        comp_measure = comp_measures(com);
        mean_good = zeros(1,sum(good_out));
        mean_bad = zeros(1,sum(bad_out));
        
        % Pull out relevant value from table
        val_ind = find(string(comp_resec_table.Properties.VariableNames)...
            == comp_measure);
        val_col = comp_resec_table(:,val_ind);
        val_col = val_col{:,1};
        % Compute mean for each patient within outcome groups, this will be used
        % for visualisation
        
        for i = 1:sum(good_out)
            mean_good(i) = mean(val_col(pat_index == i));
        end
        for i = 1:sum(bad_out)
            mean_bad(i) = mean(val_col(pat_index == i+sum(good_out)));
        end
        
        [Vg, Ig] = sort(mean_good, 'descend');
        [Vb, Ib] = sort(mean_bad, 'descend');
        
        % Find indeces to organise by median value
        ind = [Ig (sum(good_out)+Ib)];
        arr = [];
        
        test = pat_index;
        
        for x = 1:size(onset_output,1)
            i = ind(x);
            test(pat_index == i) = x;
        end
        
        figure()
        beeswarm(test, val_col, 3,'sort_style','up','overlay_style','sd',...
            'dot_size', 2, 'corral_style', 'gutter');
        xlim([-1 16])
        ylim([0,1])
        xline(sum(good_out)+0.5)
        title(sprintf('%s across patients (based on %s)', comp_measure, onset_measure))
        %saveas(gcf,sprintf('figures/roi_second/%s_all_pat_%s', comp_measure, onset_measure), 'png')

        %% Compare summary statisticics across good and bad outcome patients (max, median, mean)
        for i = 1:length(patients)
            median_val(i) = median(val_col(pat_index == i));
            mean_val(i) = mean(val_col(pat_index == i));
            max_val(i) = max(val_col(pat_index == i));
        end

        vals = [median_val; mean_val; max_val];
        for av = 1:length(av_method)
            % Add array of values summarised across patients
            row_index = find(vals_table_resec.("Onset Measure") == onset_measures(ons) & ...
            vals_table_resec.("Comparison Measure") == comp_measures(com) & ...
            vals_table_resec.("Averaging Method") == av_method(av)) ;
            vals_table_resec.Values(row_index) = mat2cell(vals(av,:),1,length(patients));
        end
        
        if comp_measures(com) == "perc"
        % Compute the number of instances of 0% and 100% resected in each
        % patient
            for pat = 1:length(patients)
                percentages = val_col(pat_index == pat);
                freq_zero(pat) = sum(percentages == 0)/length(percentages);
                freq_hundred(pat) = sum(percentages == 1)/length(percentages);
            end
            perc_index = find(freq_perc.("Onset Measure") == onset_measures(ons));
            freq_perc.("Frequency 0")(perc_index) = mat2cell(freq_zero, 1, length(patients));
            freq_perc.("Frequency 100")(perc_index) = mat2cell(freq_hundred, 1, length(patients));
               
        end
        





    end
end
 

% Compare values between good and bad outcome patients and store in
% vals_table_resec
for i = 1:size(vals_table_resec,1)
    comp_vals = vals_table_resec.Values{i,1};
    [p, ~, stats] = ranksum(comp_vals(1:sum(good_out)), ...
        comp_vals(sum(good_out)+1:end),"alpha",0.05, "method", "approximate");
    vals_table_resec.z_val(i) = abs(stats.zval)/sqrt(length(patients));
    vals_table_resec.p_val(i) = p;
end

%% Compare frequencies of 0% and 100% resected across good and bad outcomes

figure()

for ons = 1:length(onset_measures)
    zero_freq = cell2mat(freq_perc.("Frequency 0")(ons));
    hundred_freq = cell2mat(freq_perc.("Frequency 100")(ons));
    [p, ~, stats] = ranksum(zero_freq(1:sum(good_out)), zero_freq(sum(good_out)+1:end),...
        "method", "approximate", "tail", "left");
    freq_perc.("Ranksum 0")(ons) = abs(stats.zval)/sqrt(length(patients));
    freq_perc.("p 0")(ons) = p;
    [p, ~, stats] = ranksum(hundred_freq(1:sum(good_out)), hundred_freq(sum(good_out)+1:end),...
        "method", "approximate", "tail", "right");
    freq_perc.("Ranksum 100")(ons) = abs(stats.zval)/sqrt(length(patients));
    freq_perc.("p 100")(ons) = p;

    subplot(2,2,ons)
    title(onset_measures(ons))
    boxplot([hundred_freq(1:sum(good_out)), hundred_freq(sum(good_out)+1:end)],[ones(length(hundred_freq(1:sum(good_out))),1); zeros(length(hundred_freq(sum(good_out)+1:end)),1)])
   ylim([-0.05 1.05])
end
sgtitle('Comparing the proportion of seizures with 100% of onset regions resected beteen good and bad outcome patients')
%saveas(gcf,'Comp 100 0', 'png')

% 
% figure()
% subplot(3,1,ons)
%     histogram(hundred_freq(1:sum(good_out)), 'FaceColor','green', 'BinWidth',0.05)
%     hold on
%     histogram(hundred_freq(sum(good_out)+1:end), 'FaceColor', 'red', 'BinWidth',0.05)
%     hold off
%     xlim([0,1])
%     ylim([0,8])
%     xlabel('Proportion of seizures')
%     ylabel("Frequency")
%     title(sprintf("%s", onset_measures(ons)))