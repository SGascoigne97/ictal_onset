function [comp_table] = calc_jacc_sorr(pat_onset, opts)
% Compute average EEG across regions of interest for each patient

% Sarah J Gascoigne 10/02/2023

% input:
%   - pat_onset: onset_output table for one specific patient 
%   - optional inputs
%       - det_method: state the method of onset detection you'd like to use
%           clinical labels ("CLO") or automatic detection ("imprint", 
%           "EI", "PLHG")
%       - comparison: state whether the comparisons are based on
%           comparisons between onset and resected regions ("resection") or
%           between pairs of seizure onsets ("pairwise")
%       - n_perm: number of permutations to complete for 'chance' distribution
%       - tau: constant term to prevent division by zero

% Note that CLO can only be used in resection comparisons

% output
%   - comp_table: table containing raw and normalised Jaccard's index,
%       Sorensen-Dice coefficient, and percentage of onset regions resected

% Could add an argument to decide if we should plot Jaccard/Sorensen
% against count of regions to show the difference

    arguments
        pat_onset % onset_output table for one specific patient 
        opts.det_method (1,1) string {mustBeMember(opts.det_method, ["CLO", "imprint", "EI", "PLHG"])} = "imprint" 
        opts.chan_or_roi (1,1) string {mustBeMember(opts.chan_or_roi, ["chan", "roi"])} = "roi" 
        opts.comparison (1,1) string {mustBeMember(opts.comparison, ["resection", "pairwise"])} = "resection" 
        opts.n_perm (1,1) double = 100 
        opts.tau (1,1) double = 10^-16 
    end
    
    %fill in optional arguments
    det_method = opts.det_method;
    chan_or_roi = opts.chan_or_roi;
    comparison = opts.comparison;
    n_perm = opts.n_perm;
    tau = opts.tau;

    % Normalising Jaccard's Index and Sorensen Dice
    % Both JI and SD are highly influenced by the size of the data (i.e., the
    % count of onset and resected regions) - we will use a permutation test to
    % adjust for this


    if det_method == "CLO"
        % Extract the clinically labelled onset array
        onset_binary = cell2mat(pat_onset.(sprintf('labelled_onset_%s',chan_or_roi)));
        % skip patient if missing CLO data
        if min(isnan(onset_binary)) == 1
            fprintf("Patient %s does not have CLO data \n", pat_onset.Patient_id)
            comp_table = array2table(nan(0,3), 'VariableNames',...
            {'N_onset_regions', 'Jaccard', 'Sorensen'});
            return
        end
        sz_count = size(onset_binary,2);
    else
        % Extract the automatically detected onsets matrix
        onset_binary = cell2mat(pat_onset.(sprintf('%s_%s',det_method,chan_or_roi)));
        % Extract seizure IDs
        seizure_ids = string(pat_onset.Segment_ids{1,1});
        
        % Remove seizures with no onset detected
        % This should have been done earlier, double check code!
%         onset_binary = onset_binary(:,sum(onset_binary) > 0);
%         seizure_ids = seizure_ids(sum(onset_binary) > 0);
         sz_count = size(onset_binary,2);
    end

    if comparison == "resection"
        % Extract resected region for this patient
        resected_binary = cell2mat(pat_onset.(sprintf('resected_%s',chan_or_roi)));
        nan_resec = isnan(resected_binary); % Compute outside loop as regions will be removed leading to error
        comp_table = array2table(nan(sz_count,3), 'VariableNames',...
            {'N_onset_regions', 'Jaccard', 'Sorensen'});
        % Iterate through seizures and compute Jaccard's index, then
        % normalise based on 'chance' distribution
        for sz = 1:sz_count
            sz_onset = onset_binary(:,sz);
            if det_method == "CLO"
                comp_table.Seizure_id(sz) = "CLO";
            else
                comp_table.Seizure_id(sz) = seizure_ids(sz);
            end
            
            % If any regions are ignored, remove them from analysis now 
            loc_nan = isnan(sz_onset)+nan_resec >0;
            sz_onset = sz_onset(~loc_nan);
            resected_binary_rm_nan = resected_binary(~loc_nan);

            % Add the number of onset regions to the output table 
            comp_table.N_onset_regions(sz) = sum(sz_onset);

              % Add in patient id
            comp_table.Patient_id(sz) = pat_onset.Patient_id(1);

             %If a seizure does not have a detected onset, skip seizure
            if sum(sz_onset, 'omitnan') == 0
                comp_table.Jaccard(sz) = NaN;
                comp_table.Jaccard_norm(sz) = NaN;
                comp_table.Sorensen(sz) = NaN;
                comp_table.Sorensen_norm(sz) = NaN;
                comp_table.Percentage_resec(sz) = NaN;
                continue
            end

            comp_table.Jaccard(sz) = jaccard(sz_onset,resected_binary_rm_nan);
            comp_table.Sorensen(sz) = (2*comp_table.Jaccard(sz))/(1+comp_table.Jaccard(sz));
            perm_jac = zeros(1, n_perm);
            perm_sor = zeros(1, n_perm);
            
            for perm = 1:n_perm
                rng(perm)
                seizure_table = array2table([resected_binary_rm_nan sz_onset], 'VariableNames', {'Resected', 'Onset'});
                seizure_shuffled = seizure_table;
                % Shuffle both onset and resected matrices 
                seizure_shuffled.Resected = seizure_shuffled.Resected(randperm(size(seizure_shuffled,1)));
                seizure_shuffled.Onset = seizure_shuffled.Onset(randperm(size(seizure_shuffled,1)));
                
                % Compute Jaccard's index on the shuffled data
                perm_jac(perm) = jaccard(seizure_shuffled.Onset, seizure_shuffled.Resected);
                perm_sor(perm) = (2*perm_jac(perm))/(1+perm_jac(perm));
            end
            
            % Compute normalised Jaccard's index 
            comp_table.Jaccard_norm(sz) = (comp_table.Jaccard(sz) - mean(perm_jac))/(std(perm_jac)+tau);
            comp_table.Sorensen_norm(sz) = (comp_table.Sorensen(sz) - mean(perm_sor))/(std(perm_sor)+tau);

            % Compute the percentage of onset regions resected
            comp_table.Percentage_resec(sz) = sum(sz_onset+resected_binary_rm_nan==2)/sum(sz_onset);
        
        end

        if max(comp_table.Jaccard_norm) == 10^(16)
            comp_table = comp_table(comp_table.Jaccard_norm < 10^(16),:);
%             % Compute maximum Jaccard (that is not 10^16)
%             max_jacc = max(comp_table(comp_table.Jaccard_norm <10^(16),:).Jaccard_norm);
%             corrected_jac = max_jacc/comp_table.Jaccard((comp_table.Jaccard_norm == max_jacc));
%             comp_table.Jaccard_norm(comp_table.Jaccard_norm == 10^(16)) = corrected_jac;
%             comp_table.Sorensen_norm(comp_table.Jaccard_norm == 10^(16)) = ...
%                 (2*comp_table.Jaccard(comp_table.Jaccard_norm == 10^(16)))...
%                 /(1+comp_table.Jaccard(comp_table.Jaccard_norm == 10^(16)));
        end
        % Reorder table to help with readability
        comp_table = comp_table(:,{'Patient_id', 'Seizure_id',... % 'N_onset_regions'
            'Jaccard', 'Jaccard_norm', 'Sorensen', 'Sorensen_norm',...
            'Percentage_resec'});

    elseif comparison == "pairwise"
        if det_method == "CLO"
            print("Pairwise comparison is not possible for clinically labelled onsets")
            comp_table = [];
        else
            % Normalisation for pairwise comparisons
            % Use a permutation test to correct for varied size of seizure onsets
            % Look at the relationship between the size of the detected onset and
            % Jaccard's index
            comp_table = array2table(nan(sum(1:(sz_count-1)),4), 'VariableNames', ...
                {'N_onset_regions_sz1', 'N_onset_regions_sz2', 'Jaccard', 'Sorensen'});
            count = 0;
            for sz1 = 1:sz_count
                for sz2 = 1:sz_count
                    if sz1<sz2
                        count = count+1;
                        perm_jac = zeros(1, n_perm);
                        perm_sor = zeros(1, n_perm);
                        sz_onset1 = onset_binary(:,sz1);
                        sz_onset2 = onset_binary(:,sz2);
                        % Remove any ignored regions
                        na_rm = isnan(sz_onset1) + isnan(sz_onset2);

                        sz_onset1 = sz_onset1(~na_rm);
                        sz_onset2 = sz_onset2(~na_rm);
        
                        if sum(sz_onset1) ~= 0 && sum(sz_onset2) ~= 0
                            comp_table.Seizure1_id(count) = seizure_ids(sz1);
                            comp_table.Seizure2_id(count) = seizure_ids(sz2);
                            comp_table.N_onset_regions_sz1(count) = sum(sz_onset1);
                            comp_table.N_onset_regions_sz2(count) = sum(sz_onset2);
                            comp_table.Jaccard(count) = jaccard(sz_onset1,sz_onset2);
                            comp_table.Sorensen(count) = (2*comp_table.Jaccard(count))/(1+comp_table.Jaccard(count));
            
                            for perm = 1:n_perm
                                rng(perm)
                                seizure_table = array2table([sz_onset1 sz_onset2], 'VariableNames', {'Onset_1', 'Onset_2'});
                                seizure_shuffled = seizure_table;
                                seizure_shuffled.Resected = seizure_shuffled.Onset_1(randperm(size(seizure_shuffled,1)));
                                seizure_shuffled.Onset = seizure_shuffled.Onset_2(randperm(size(seizure_shuffled,1)));
                                
                                % Compute Jaccard's index
                                perm_jac(perm) = jaccard(seizure_shuffled.Onset_1, seizure_shuffled.Onset_2);
                                perm_sor(perm) = (2*perm_jac(perm))/(1+perm_jac(perm));
                            end
                            
                            % Compute normalised Jaccard's index 
                            comp_table.Jaccard_norm(count) = (comp_table.Jaccard(count) - mean(perm_jac))/(std(perm_jac)+tau);
                            comp_table.Sorensen_norm(count) = (comp_table.Sorensen(count) - mean(perm_sor))/(std(perm_sor)+tau);
    
                            % Add in patient ID
                            comp_table.Patient_id(count) = pat_onset.Patient_id;
                        end

                    end
                end
            end
        end
        % Reorder table to help with readability
        comp_table = comp_table(:,{'Patient_id', 'Seizure1_id', 'Seizure2_id',...%'N_onset_regions_sz1', 'N_onset_regions_sz2',
            'Jaccard', 'Jaccard_norm', 'Sorensen', 'Sorensen_norm'});
        
    end
end
